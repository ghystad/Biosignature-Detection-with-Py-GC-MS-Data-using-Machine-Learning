
##"Creating graphs for pyrolysis gas chromatography-mass spectrometry data analyzed
## by random forest classification method"

## Data consist of 150 samples of organic molecular mixtures analyzed by 
## pyrolysis gas chromatography-mass spectrometry (py-GC-MS).
## Each data set has the scan numbers in the rows and the mass to charge ratio (m/z)
## values in the columns.
## The data values are the raw intensities as measured by py-GC-MS.
## Each sample has either an abiotic (A) or biotic (B) origin. 
## Purpose is to determine classification rules to predict whether the samples are
## biotic or abiotic.
## Purpose is also to determine the features, m/z and scan numbers, that discriminate
## between the biotic and the abiotic groups.
## The data are preprocessed before we apply random forest classification model.
## In this file we create different graphs regarding the importance variables
## from the random forest model and Monte Carlo simulations.


library(dplyr)         # for dataframe computation
library(MALDIquant)    # for chemometrics processing
library(caret)         # for machine learning
library(mlr3)          # for machine learning
library("mlr3verse")   # for machine learning
library("mlr3learners")# for machine learning
library("mlr3tuning")  # for machine learning
library("data.table")  # for rbindlist
library("ggplot2")     # for plots
library(rgl)           # for 3D graphs
library(plotly)        # for interactive graphs


#Data preparation
setwd("C:/Desktop/Templeton_grant_research/Second_paper_material/Data2022")
species=c(rep("A",5),"B", "A","A","B","B","A","C", rep("B",5),rep("C",4),rep("B",4),"A",
          rep("B",11),"A","A","B","A","A","B","B","A",rep("B",3),"A","B",rep("A",6),
          rep("B",3),"A","A","B","B",rep("A",4),"B",rep("A",5),"B",rep("A",10),
          rep("B",4),rep("A",5),rep("B",3),rep("A",3),rep("B",4),rep("A",8),"B",
          "A","A","B","C","C","A",rep("B",3),"A","B", "B","A","B",rep("A",3),"B","A",
          rep("A",7),"B","A","B")
species2=as.numeric(as.factor(species))
ind=which(species2=="3")
species_n=species2[-ind]


#Reading in the data is based on the code found in:
## "How to import multiple .csv files simultaneously in R and create a data frame"
## by Dahna, A., datascience+, March 09, 2019
#https://datascienceplus.com/how-to-import-multiple-csv-files-simultaneously-in-r-and-create-a-data-frame/
species_names <- list.files()
species_names=species_names[-ind]
z=lapply(species_names, read.delim)  #read in all of 134 datasets

setwd("C:/Desktop/Templeton_grant_research/Second_paper_material/Data2023")
species_names2 <- list.files()
z2=lapply(species_names2, read.delim)  #read in all of 16 datasets

NN=700 #number of m/z values
mass=seq(50,NN,1) #m/z 
MM=3240 #number of scans
ndim=length(species_names)
ndim2=length(species_names2)

#Scan number:m/z:intensity values for 134 samples
M=list()
for(i in 1:ndim){
	colnames(z[[i]])="mass"
	#remove commas
	z[[i]]=data.frame(do.call("rbind", strsplit(as.character(z[[i]]$mass), ",", 
	                                            fixed = TRUE)))
	z[[i]]=data.frame(lapply(z[[i]],as.numeric))
	colnames(z[[i]])=c("scan",as.character(seq(50,NN,1)))
	z[[i]]=z[[i]] %>% slice(1:MM)       #selects the first MM rows  
	M[[i]]=z[[i]]
}

#Scan number:m/z:intensity values for 16 samples
M2=list()
for(i in 1:ndim2){
	colnames(z2[[i]])="mass"
	#remove commas
	z2[[i]]=data.frame(do.call("rbind", strsplit(as.character(z2[[i]]$mass), ",",
	                                             fixed = TRUE)))
	z2[[i]]=data.frame(lapply(z2[[i]],as.numeric))
	colnames(z2[[i]])=c("scan",as.character(seq(50,NN,1)))
	z2[[i]]=z2[[i]] %>% slice(1:MM)       #selects the first MM rows
	M2[[i]]=z2[[i]]
}

M_new=c(M,M2)
N=length(M_new)
species_n2=rep(2,16)
species_names_new=c(species_names,species_names2)

y=c(species_n,species_n2)
y=factor(y,labels=c("A","B"))

#########################################################################################
#Preprocessing
#Detect the significant peaks as local max above four times the signal to noise ratio

MZ=151  #number of m/z values to use
#Create Chromatograms for each sample and each m/z value inside each sample
sample_list=list()
sample_location_list=list()
suppressWarnings({
for (i in 1:N){
    S=list()
	  for(j in 1:MZ){
	     S[[j]] = createMassSpectrum(mass=seq(1,MM,1), intensity=unlist(M_new[[i]][,j+1]),
	     metaData=list(name="Chrom"))  
	  }
   chrom = transformIntensity(S,method="sqrt") #stabilize the variance
   chrom = smoothIntensity(chrom, method="MovingAverage",halfWindowSize=5)#smooth the data
   chrom = removeBaseline(chrom, method="SNIP") #remove the baseline

   # Put the processed chromatograms back into a dataframe
   processed_chrom_list=list()
	  for (k in 1:MZ){
	      processed_chrom_list[[k]] = as.numeric(intensity(chrom[[k]]))
	   }

   processed_mass_dataframe = as.data.frame(do.call(rbind, processed_chrom_list))
   Ma=max(processed_mass_dataframe)
   Mi=min(processed_mass_dataframe)
   #Normalize across sample
   processed_mass_dataframe = t((processed_mass_dataframe - Mi)/(Ma - Mi))
   processed_mass_dataframe = as.data.frame(processed_mass_dataframe)
   S2=list()
	  for(t in 1:MZ){
	     S2[[t]] = createMassSpectrum(mass=seq(1,MM,1),
	                                  intensity=processed_mass_dataframe[, t],
	     metaData=list(name="Chrom_normalized"))  
	  }

    peaks = detectPeaks(S2, method="MAD", halfWindowSize=20, SNR=4)
    peak_list=list()
  	for (tt in 1:MZ){
         v=numeric(MM)
         scan_number=mass(peaks[[tt]])
         v[scan_number] = intensity(peaks[[tt]])
	       peak_list[[tt]] = v
	 }
 	 processed_peaks = t(as.data.frame(do.call(rbind, peak_list)))
   row.names(processed_peaks)=c(paste0("R", 1:MM))
   colnames(processed_peaks)=c(paste0("M", 50:(MZ+50-1)))
	 processed_peaks2=as.data.frame(processed_peaks) %>% 
   mutate(bin = cut(seq(1,MM,1),breaks=180,dig.lab = 6)) %>%  # bin the peaks by scan #
   group_by(bin) %>%
   summarise_all(max)
	 sample_list[[i]] = processed_peaks2
   sample_location_list[[i]] = processed_peaks
}
})

#Scan # and mass spectrum for each sample
mass_scan_list_loc=list()
for(i in 1:N){
    sm_loc=as.numeric(unlist(sample_location_list[[i]])) 
    mass_scan_list_loc[[i]]=sm_loc
}

#Sample vs mass/scan numbers
#Put the samples into a dataframe
data_mass_scan_loc = do.call(rbind, mass_scan_list_loc) 

bin2=as.character(seq(1,3240,1))
MS=as.character(seq(50,(MZ+50-1),1))

colnames(data_mass_scan_loc) = paste(outer(bin2, MS, paste, sep = ';'))

data_mass_scan_new=as.data.frame(ifelse(data_mass_scan_loc > 0,1,0))
MZ_name=c(paste0(";",50:(MZ+50-1)))
MZ_name2=c(paste0(50:(MZ+50-1)))
MZ_name3=c(paste0(".",50:(MZ+50-1)))
bin3 = unique(cut(seq(1,MM,1),breaks=180,dig.lab = 6))

#Finding the endpoints in bin4 is based on code found in 
#"Obtain endpoints from interval that is factor variable"

#https://stackoverflow.com/questions/40665240/obtain-endpoints-from-interval-that-is-factor-variable
#Stack Overflow:
bin4=unique(as.numeric(unlist( strsplit( gsub( "[][(]" , "", levels(bin3)) , ","))))

#Determine the mean scan number for each bin and each 
#m/z value across the training samples
mean_scan_name = list()
for (i in 1:MZ){
	  data_mass_scan_new2 = data_mass_scan_new %>% select(ends_with(MZ_name[i])) 
	  scan_name=sub(basename(colnames(data_mass_scan_new2)),pattern = MZ_name[i], 
	                replacement = "" , fixed = TRUE)
	  colnames(data_mass_scan_new2) = scan_name #name the columns with the scan numbers
	  count_elements=function(x){
	    y=sum(x)
    }
   #Count the number of elements for each scan number across the training samples
   n_elements=as.numeric(apply(as.matrix(data_mass_scan_new2),2,count_elements)) 
   #Repeat the nonzero scan numbers with its frequency
   vec=as.numeric(rep(colnames(data_mass_scan_new2), times=n_elements)) 
   scan_dataframe=list() #scan numbers for each bin for the ith m/z value selected 
   indd1=1:floor(bin4[2])
   ind_L1=which(vec %in% indd1)
   Lt1=vec[ind_L1]
   #Mean number of scan number in the first bin
   scan_mean1=if(length(Lt1)>0){round(mean(Lt1)) 
	         }else {-1}
   DD1=data.frame(matrix(ncol = 1, nrow =  1)) 
   colnames(DD1)= paste(outer(as.character(scan_mean1),MZ_name2[i], paste, sep = ';'))
   scan_dataframe[[1]]=DD1
	 for (j in 2:(length(bin4)-1)){
			 indd=(floor(bin4[j])+1):floor(bin4[j+1])
        ind_L=which(vec %in% indd)
        Lt=vec[ind_L]
			  #Mean number of scan number in each bin
        scan_mean=if(length(Lt)>0){round(mean(Lt)) 
				           }else {-j} #to give empty columns a unique name
        DD=data.frame(matrix(ncol = 1, nrow =  1)) 
			  colnames(DD)= paste(outer(as.character(scan_mean),MZ_name2[i], 
			                            paste, sep = ';'))
			  scan_dataframe[[j]]=DD
		}
   #The mean scan number for each bin of scans for the ith m/z value and each sample	
   mean_scan_name[[i]]=do.call(cbind, scan_dataframe) 
}

# To get the column names of the interleaved m/z values and scan numbers as columns
data_mass_scan_interl=do.call(cbind, mean_scan_name) 
names_new=colnames(data_mass_scan_interl)

#Scan # and mass spectrum for each sample
mass_scan_list=list()
for(i in 1:N){
   sm=as.numeric(unlist(sample_list[[i]]))[181:(180*(MZ+1))]     
   mass_scan_list[[i]]=sm
}

#Sample vs mass/scan numbers
data_preprocessed = do.call(rbind, mass_scan_list) #put the samples into a dataframe
colnames(data_preprocessed ) = names_new
data_preprocessed  = as.data.frame(data_preprocessed )

count_nonzero=function(x){(length(which(x > 0)))/(dim(data_preprocessed)[2])}
#Ratio of nonzero feature values for each observation
Perc_Nnonzero=apply(data_preprocessed ,1,count_nonzero) 

data_preprocessed  = data_preprocessed  %>% mutate(Perc_Nnonzero)
colnames(data_preprocessed) = make.names(colnames(data_preprocessed))

#Removing variables with near zero variance, nearZeroVar, is based on the book
#"Applied Predictive Modeling" by Kuhn, M. and Johnson, K.,
#Springer, New York, 2013 

#Detect features with near zero variance
near_zero_variance = nearZeroVar(data_preprocessed) 
#Remove features with near zero variance
data_preprocessed  = data_preprocessed[, -near_zero_variance] 
dim(data_preprocessed) #new dimension 
#Contains all samples with the selected features and corresponding y-values (A or B)
training_transformed=data.frame(data_preprocessed,y=y)



#########################################################################################
#Machine learning
#The machine learning R-code is based on the mlr3 library and the book
## https://mlr3book.mlr-org.com/, "Flexible and Robust Machine Learning Using mlr3 in R"
#by Lang, M. et al

#Create a task 
task=as_task_classif(training_transformed, target = "y", positive="B")

#Using stratified random sampling for each resampling
task$col_roles$stratum = "y"


#Determine importance variables
#Use the results from the random forest model in previous file with
#mtry=1688 and num.trees=2658

learner2 = lrn("classif.ranger", id = "rf",predict_type="prob",importance = "impurity")
learner2$param_set$values$mtry =        1688
learner2$param_set$values$num.trees  =     2658  

set.seed(99)
learner2$train(task)

#Importance variables in ranked order
filter = flt("importance", learner = learner2)
filter$calculate(task)
dff=as.data.table(filter)

#Monte Carlo (MC) simulations using random forest from previous file
df.list=readRDS("C:/Desktop/Templeton_grant_research/Second_paper_material/Materials_for_submission/New_excel_4.3.3/df.list_Monte_Carlo.RDS")
df.test = readRDS("C:/Desktop/Templeton_grant_research/Second_paper_material/Materials_for_submission/New_excel_4.3.3/df.test_Monte_Carlo.RDS")



#####Determine the importance variables that have the lowest average rank across the
#Monte Carlo samples and that are among the original importance variables from the
#random forest model

K=(dim(training_transformed)[2]-1)
variables=colnames(training_transformed)[1:K]
ID=seq(1,K,1)
variables2=data.frame(variables,ID)

feature_rank = list() #Arrange by features and include their ranks
for (i in 1:50){
  ind2=2*i-1
  fn = as.character(unlist(as.data.frame(df.list)[,ind2]))
  feature_rank2 = list()
	for (j in 1:K){
	  ID2=which(variables2[,1] %in% fn[j])
	  feature_rank2[[j]] = data.frame(ID2,j)
	  colnames(feature_rank2[[j]]) = c("Feature","Rank")
	}
  feature_rank3 = do.call(rbind, feature_rank2) #put each feature in ranked order
  feature_rank[[i]] = (feature_rank3 %>% arrange(Feature))[,2] 
}

#Rank for each feature across the MC samples
feature_rank4 = data.frame(ID,do.call(cbind, feature_rank))
#Calculate the mean rank for each feature across the MC samples
feature_rank4 =feature_rank4 %>% mutate(mean_rank = apply(feature_rank4[,2:51],1,mean))
#Order the features according to the lowest mean rank across the MC samples
feature_rank5= feature_rank4 %>% arrange(mean_rank) 

DIMM=60
vari = feature_rank5[1:DIMM,1]
#Ranked importance variables average across the MC samples
average_imp=variables2[vari,] 
# The first DIMM importance variables in the original sample
base=as.character(unlist(dff[1:DIMM,1]))
#Lowest ranked importance variables in the original sample that are also among the lowest
#ranked importance variables across the MC samples
Av=intersect(base,average_imp[,1]) 

y2=c(1,1,1,1,1,3,1,1,3,3,1,rep(2,9),1,3,3,rep(2,9),1,1,3,1,1,2,2,1,3,2,3,1,2,rep(1,6),
     2,3,2,1,1,2,3,rep(1,4),2,rep(1,5),2,rep(1,10),3,2,2,2,rep(1,5),2,2,2,1,1,1,3,3,3,3,
     rep(1,8),3,1,1,3,1,2,3,2,1,2,2,1,3,1,1,1,2,rep(1,8),3,1,3,rep(2,13),3,2,3)
y2=factor(y2,labels=c("A","B","C"))
#A: Abiotic
#B: Contemporary biotic
#C: Altered biotic

indexA=which(y2=="A")
indexB=which(y2=="B")
indexC=which(y2=="C")
lA=length(indexA) #number of abiotic samples
lB=length(indexB) #number of contemporary biotic samples
lC=length(indexC) #number of altered biotic samples
abiotic=species_names_new[indexA] #abiotic sample names
bioticB=species_names_new[indexB] #contemporary biotic sample names
bioticC=species_names_new[indexC] #altered biotic sample names

####################################################################################
#######Calculate the proportion of samples that have the importance variables sorted
#######by abiotic, biotic contemporary, and altered biotic 
calc_proportion=function(Av) {
#Intensity values for the importance variables in the vector Av for each species
dA=list()
dB=list()
dC=list()
for(i in 1:length(Av)){
	dA2=training_transformed %>% select(Av[i]| "y")   %>% filter(y2=="A")
	dA[[i]] = dA2[,1]
	dB2=training_transformed %>% select(Av[i]| "y")   %>% filter(y2=="B")
	dB[[i]] = dB2[,1]
	dC2=training_transformed %>% select(Av[i]| "y")   %>% filter(y2=="C")
	dC[[i]] = dC2[,1]
}

data_A = do.call(cbind, dA)
colnames(data_A) = Av
data_B = do.call(cbind, dB)
colnames(data_B) = Av
data_C = do.call(cbind, dC)
colnames(data_C) = Av

#Put non-zero values = 1 
data_AA=ifelse(data_A[,1:(length(Av))]>0,1,0)
data_BB=ifelse(data_B[,1:(length(Av))]>0,1,0)
data_CC=ifelse(data_C[,1:(length(Av))]>0,1,0)

data_ABC = rbind(data_AA,data_BB,data_CC)

data_stat=data.frame(apply(data_AA,2,mean),apply(data_BB,2,mean),apply(data_CC,2,mean),
                     apply(data_ABC,2,mean))
colnames(data_stat)=c("Pr A", "Pr B(contemporary) ","Pr B(altered)", "Prop.total")
return(data_stat)
}

data_stat=calc_proportion(Av)


#########################################################################################
######Find the distribution of proportion of abiotic samples, biotic samples,
######contemporary biotic and altered biotic samples
allA=training_transformed %>% filter(y=="A")
allB=training_transformed %>% filter(y=="B")
allB_Cont=training_transformed[,1:(K-1)]%>% slice(indexB) #omitting Perc_Nnonzero feature 
allB_Alt=training_transformed[,1:(K-1)]%>% slice(indexC)

allA_bin=ifelse(allA[,1:(K-1)]>0,1,0)
allB_bin=ifelse(allB[,1:(K-1)]>0,1,0)
allB_Cont_bin=ifelse(allB_Cont>0,1,0)
allB_Alt_bin=ifelse(allB_Alt>0,1,0)

# Proportion of abiotic samples that contain each feature 
prop_A= apply(allA_bin,2,mean) 
# Proportion of biotic samples that contain each feature 
prop_B= apply(allB_bin,2,mean) 
# Proportion of contemporary biotic samples that contain each feature 
prop_B_Cont= apply(allB_Cont_bin,2,mean) 
# Proportion of altered biotic samples that contain each feature 
prop_B_Alt = apply(allB_Alt_bin,2,mean)  

#Median number of features for abiotic, biotic,
#contemporary  biotic and altered biotic samples, respectively
median(prop_A)
median(prop_B)
median(prop_B_Cont)
median(prop_B_Alt)


##########################################################################################
#Proportion of samples that contain the lowest ranked importance variables in ranked order 

D_A = data_stat[,1]
D_B_Cont = data_stat[,2]
D_B_Alt = data_stat[,3]
D_B = (data_stat[,2] * lB+ data_stat[,3] * lC)/(lB+lC)

DIM=length(Av)
D_A_array=array(c(D_A), dim=c(DIM,1))
D_B_array=array(c(D_B), dim=c(DIM,1))
D_B_Cont_array=array(c(D_B_Cont), dim=c(DIM,1))
D_B_Alt_array=array(c(D_B_Alt), dim=c(DIM,1))

#These importance variables are on average found in the following proportion of samples:
mean(D_A)
mean(D_B)
mean(D_B_Cont)
mean(D_B_Alt)

#######################################################################################
#####Distribution of the lowest ranked importance variables in ranked order 

allA2=allA %>% select(all_of(Av))
allB2=allB  %>%  select(all_of(Av))
allB_Cont2 =  allB_Cont %>% select(all_of(Av))
allB_Alt2 = allB_Alt  %>% select(all_of(Av))

Nrow = lA+(lB+lC)*2
size_Dbar = Nrow*length(Av)
Dbar=matrix(numeric(size_Dbar),nrow=Nrow)
Dbar=as.data.frame(Dbar)

#DIM=length(Av)
for(i in 1:DIM){
	DbarA=data.frame(allA2[,i])
	colnames(DbarA)=c(Av[i])

	DbarB=data.frame(allB2[,i])
	colnames(DbarB)=c(Av[i])

	DbarBCont=data.frame(allB_Cont2[,i])
	colnames(DbarBCont)=c(Av[i])

	DbarBAlt=data.frame(allB_Alt2[,i])
	colnames(DbarBAlt)=c(Av[i])

	Dbar[,i]=rbind(DbarA,DbarB,DbarBCont, DbarBAlt)
}

Type=c(rep("A",lA),rep("B",(lB+lC)),rep("B(Cont.)",lB),rep("B(Alt.)",lC))
Dbar=data.frame(Dbar,Type)
colnames(Dbar)=c(Av,"Type")
########################################################################################
######Split the vector Xscan#.m/z into the scan# component and the m/z value component

#K is the number of features 
#Omit Perc_nonzero feature
datadf=as.character(unlist(dff[c(seq(1,580,1),seq(582,K,1)),1]))

#Splitting a vector string based on code found in "Splitting Strings in
## R programming – strsplit() method":
## https://www.geeksforgeeks.org/splitting-strings-in-r-programming-strsplit-method/
# 11 November, 2021
datadf_st = strsplit(datadf, split = "[.]+")
datadf_st2=list()
for (i in 1:length(datadf_st)){
datadf_st2[[i]] = strsplit(datadf_st[[i]],split = '""')
}

#Get the first element of a list based on code found in "R list get first item of each
#element":
#https://stackoverflow.com/questions/44176908/r-list-get-first-item-of-each-element ,
#stack overflow
datadf_st21= unlist(sapply(datadf_st2, function(x) x[1]))

#gsub based on code found in "Remove Character From String in R":
#https://sparkbyexamples.com/r-programming/remove-character-from-string-in-r/ 
#by Nelamali,N., March 27, 2024
xx = as.numeric(gsub('[X]','',datadf_st21)) #Scan numbers
yy= as.numeric(unlist(sapply(datadf_st2, function(x) x[2]))) #m/z values
######################################################################################
#Split the scan# and m/z values for the ranked importance variables into two separate
## components for different number of features.

data_split=function(Av){
datadf3=Av

datadf_st3 = strsplit(datadf3, split = "[.]+")
datadf_st23=list()
for (i in 1:length(datadf_st3)){
datadf_st23[[i]] = strsplit(datadf_st3[[i]],split = '""')
}


datadf_st213= unlist(sapply(datadf_st23, function(x) x[1]))

xx2 = as.numeric(gsub('[X]','',datadf_st213))
yy2= as.numeric(unlist(sapply(datadf_st23, function(x) x[2])))
return(list(xx2=xx2,yy2=yy2))
}
xx2=as.numeric(unlist(data_split(Av)[1]))
yy2=as.numeric(unlist(data_split(Av)[2]))

####################################################################################
split_name = strsplit(species_names_new, split = "[.]+")
split_name2=list()
for (i in 1:length(split_name)){
split_name2[[i]] = strsplit(split_name[[i]],split = '""')
}

#Name of samples
split_name_new= unlist(sapply(split_name2, function(x) x[1]))
split_name_new =  unlist( strsplit( gsub( "[][(]" , "", split_name_new) , "3d"))

#Select the first 30 ranked importance variables
data_new=training_transformed %>%  select(all_of(Av[1:30]))

y3=factor(y2,labels=c("Abiotic","Biotic (contemporary)","Biotic (altered)"))
#PCA
pr.out3=prcomp(data_new, scale=T)
summary(pr.out3)

#########################################################################################
######Create a dataframe with scan#, m/z, intensity coordinates, sample type, and name for
## the 48 ranked importance variables
Av3=Av

data_Av3=training_transformed %>% select(all_of(Av3))
sample_import=data.frame(species_names_new,y2,data_Av3)

#Coordinates for scan# and m/z for the first 48 ranked importance variables
cord=data.frame(xx2,yy2)
colnames(cord)=c("scan #", "m/z")

#Intensity values for the first 48 ranked importance variables for each sample
ss=list()
for(i in 1:150){
	ss[[i]] = sample_import[i,]
}

sample_frame = as.data.frame(t(do.call(rbind,ss)))
colnames(sample_frame)=species_names_new
#The 48 ranked importance variables with their scan# and m/z coordinates and intensity 
#value for each sample
cord2=cbind(cord,sample_frame[3:(length(Av3)+2),]) 
#Code for converting dataframe to numeric based on code found in
#"Code for converting entire data frame to numeric":
#https://stackoverflow.com/questions/60288057/code-for-converting-entire-data-frame-to-numeric , 
#stack overflow
cord2=mutate_all(cord2, function(x) as.numeric(as.character(x)))

#Create a dataframe for abiotic samples with scan#, m/z value coordinates, intensity,
## sample type
data_mA=list()
for (i in 1:lA){
      data_mA[[i]]= data.frame(xx2,yy2,as.numeric(data_Av3[indexA[i],]))
      colnames(data_mA[[i]])=c("Scan","Mass_to_charge_ratio","Intensity")
}
typeA=factor(rep("Abiotic",lA*length(Av3)))
data3DA=do.call(rbind,data_mA)
data3DA=data.frame(data3DA,typeA)
colnames(data3DA)=c("Scan","Mass_to_charge_ratio","Intensity","Type")
NameA=factor(rep(split_name_new[indexA],each=length(Av3)))#Abiotic sample names

#Create a dataframe for contemporary biotic samples with scan#, m/z value coordinates,
## intensity, sample type
data_mB=list()
for (i in 1:lB){
     data_mB[[i]]= data.frame(xx2,yy2,as.numeric(data_Av3[indexB[i],]))
     colnames(data_mB[[i]])=c("Scan","Mass_to_charge_ratio","Intensity")
} 
typeB=factor(rep("Biotic (cont.)",lB*length(Av3)))
data3DB=do.call(rbind,data_mB)
data3DB=data.frame(data3DB,typeB)
colnames(data3DB)=c("Scan","Mass_to_charge_ratio","Intensity","Type")
NameB=factor(rep(split_name_new[indexB],each=length(Av3)))

#Create a dataframe for altered biotic samples with scan#, m/z value coordinates,
## intensity, sample type
data_mC=list()
for (i in 1:lC){
     data_mC[[i]]= data.frame(xx2,yy2,as.numeric(data_Av3[indexC[i],]))
    colnames(data_mC[[i]])=c("Scan","Mass_to_charge_ratio","Intensity")
}
typeC=factor(rep("Biotic (alt.)",lC*length(Av3)))
data3DC=do.call(rbind,data_mC)
data3DC=data.frame(data3DC,typeC)
colnames(data3DC)=c("Scan","Mass_to_charge_ratio","Intensity","Type")
NameC=factor(rep(split_name_new[indexC],each=length(Av3)))

NameS=c(NameA, NameB, NameC)
data3D=rbind(data3DA, data3DB, data3DC)
#dataframe with scan#, m/z, intensity coordinates, sample type, and name for the 48 ranked
## importance variables
data3D=cbind(data3D,NameS)
ind_all=which(data3D[,3]>0)

#########################################################################################
#Split the scan# and m/z values for the ranked importance variables into two separate
# components for different number of features.
#Determine the proportion of samples that contain these features with color coding.
xx3=list()
yy3=list()
Av2=list()
group_col2=list()
INT_b = list()

DIMM=c(20, 60,150)
	for(j in 1:3){
		vari2 = which(variables2[,2]%in% feature_rank5[1:DIMM[j],1])
		#Ranked importance variables average across the MC samples
		average_imp2=variables2[vari2,] 
		# Highest importance variables in the original sample
		base2=as.character(unlist(dff[1:DIMM[j],1]))
		#Most important importance variables in the original sample that are also among
		#the lowest ranked importance variables across the MC samples
		Av2[[j]]=intersect(base2,average_imp2[,1]) 
		
		#Splitting a vector string
		xx3[[j]] = as.numeric(unlist(data_split(Av2[[j]])[[1]]))
		yy3[[j]] = as.numeric(unlist(data_split(Av2[[j]])[[2]]))
		#Calculate the proportion of samples that have the importance variables sorted
		#by abiotic, biotic contemporary, and altered biotic 
		#Intensity values for the importance variables in the vector Av2 for each species
		data_stat2 = calc_proportion(Av2[[j]])
            group_col = numeric(length(Av2[[j]]))
		for (k in 1:length(Av2[[j]])){
			if (data_stat2[k,1] > 0.5){
			group_col[k]= 1
			}else if (data_stat2[k,2] > 0.5 && data_stat2[k,3] <= 0.5) {
			group_col[k]= 2
			}else if (data_stat2[k,3] > 0.5 && data_stat2[k,2] <= 0.5) { 
			group_col[k]= 3
			}else if (data_stat2[k,2] > 0.5 && data_stat2[k,3] > 0.5) {
			group_col[k]= 4
			} else {
			group_col[k]= 5
			}
		}

		group_col2[[j]]  = c(group_col ,rep(6,(K-1)-length(Av2[[j]])))
		INT_b[[j]] = which(as.character(unlist(dff[,1])) %in% Av2[[j]])
}

#######################################################################################
#Color coding for Figure 3A
data_stat3=calc_proportion(datadf)

group_col3 = numeric(length(datadf))
for (i in 1:length(datadf)){
	if (data_stat3[i,1] > data_stat3[i,2] && data_stat3[i,1] > data_stat3[i,3]) {
	group_col3[i]= 1
	}else if (data_stat3[i,2] > data_stat3[i,3] && data_stat3[i,2] > data_stat3[i,1]) {
	group_col3[i]= 2
	} else { 
	group_col3[i]= 3
	}
}


#####Graphs###
#######################################################################################
#Scatterplots of scan#:m/z:intensity for importance variables in ranked order
#Figure 3
#Plotting legend outside plot is based on code found in
#"Plot a legend outside of the plotting area in base graphics?":
# "https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics , 
#stack overflow modified

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 7, 0, 0),ps=8, family="sans", new=TRUE)
  on.exit(par(opar))
  plot.new()
  legend(...)
}

par(mfrow=c(2,2),xpd=T,family="sans", ps=8, mai = c(0.5, 0.5, 0.3, 0.3), 
oma = c(1.3, 1.5, 0, 0))

plot(xx, yy, col=c("#B15928","#33A02C","#1F78B4")[group_col3],xlab="",ylab="",
     pch=c(19, 19,19)[group_col3],main="",cex=0.57)
text(100,200,label=substitute(paste(bold(('a')))), col="black")

j=1
plot(c(xx3[[j]],xx[-INT_b[[j]]])[1:300],c(yy3[[j]],yy[-INT_b[[j]]])[1:300],
     col=c("#B15928","#33A02C","#1F78B4","#6A3D9A","#CAB2D6","grey")[group_col2[[j]]][1:300],
     xlab="",ylab="",pch=c(19,19,19,19,19,1)[group_col2[[j]]][1:300],main="")
text(100,200,label=substitute(paste(bold(('b')))), col="black")
legend(c(0,200),c("Biotic (cont.)","Biotic (alt.)","Biotic (cont. and alt.)",
                  expression(""<="50% of the samples"),"IV"),
       col=c("#33A02C","#1F78B4","#6A3D9A","#CAB2D6","grey"), pch=c(19,19,19,19,1),bty = "n")

j=2
plot(c(xx3[[j]],xx[-INT_b[[j]]])[1:300],c(yy3[[j]],yy[-INT_b[[j]]])[1:300],
     col=c("#B15928","#33A02C","#1F78B4","#6A3D9A","#CAB2D6","grey")[group_col2[[j]]][1:300],
     xlab="",ylab="",pch=c(19, 19,19,19,19,1)[group_col2[[j]]][1:300],main="")
text(100,200,label=substitute(paste(bold(('c')))), col="black")
legend(c(0,200),c("Biotic (cont.)","Biotic (alt.)","Biotic (cont. and alt.)",
                  expression(""<="50% of the samples"),"IV"),
       col=c("#33A02C","#1F78B4","#6A3D9A","#CAB2D6","grey"), pch=c(19,19,19,19,1),bty = "n")

j=3
plot(c(xx3[[j]],xx[-INT_b[[j]]])[1:300],c(yy3[[j]],yy[-INT_b[[j]]])[1:300],
     col=c("#B15928","#33A02C","#1F78B4","#6A3D9A","#CAB2D6","grey")[group_col2[[j]]][1:300],
     xlab="",ylab="",pch=c(19, 19,19,19,19,1)[group_col2[[j]]][1:300],main="")
text(100,200,label=substitute(paste(bold(('d')))), col="black")
legend(c(0,200),c("Abiotic","Biotic (cont.)","Biotic (alt.)","Biotic (cont. and alt.)",
                  expression(""<="50% of the samples"),"IV"),
       col=c("#B15928","#33A02C","#1F78B4","#6A3D9A","#CAB2D6","grey"),
       pch=c(19,19,19,19,19,1),bty = "n")

mtext("Scan #", side = 1,ps=8, outer = TRUE, line = 0)
mtext(expression(italic("m/z")), side = 2,ps=8, outer = TRUE, line = -1)

add_legend("topleft",legend=c("Abiotic","Biotic (cont.)","Biotic (alt.)"),
           col=c("#B15928","#33A02C","#1F78B4"), pch=c(19,19,19),horiz=TRUE,bty = "n")



#######################################################################################
#Figure 4
##3D-PCA
cols=c("#B15928","#33A02C","#1F78B4")
Cols=function(vec){
cols=c("#B15928","#33A02C","#1F78B4")
return(cols[as.numeric(as.factor(vec))])
}

windowsFonts(A = windowsFont("sans")) 
plot3d(pr.out3$x[,1:3], col = Cols(y2), type = 's', radius = .2  )
text3d(pr.out3$x[,1:3],
       texts=split_name_new,
       ps= 0.6, pos=3, family="A")
 
bbox3d(color = c("grey", "black"), emission = "grey", 
         specular = "grey", shininess = 5, alpha = 0.8)

########################################################################################
#Histograms over distribution of proportion of samples
#Figure 4S

par(mfrow=c(2,2),xpd=T,family="sans", ps=8, mai = c(0.5, 0.5, 0.3, 0.3), 
oma = c(1.3, 1.5, 0, 0))

hist(prop_A,col="#B15928",xlab="",ylab="", main="",xlim=c(0,1), ylim=c(0,3300))
text(-0.15,3300,label=substitute(paste(bold(('a')))), col="black",xpd=NA)
title(xlab="Proportion of abiotic samples", line=2)

hist(prop_B,col="#6A3D9A",xlab="",ylab="", main="",xlim=c(0,1), ylim=c(0,3300))
text(-0.15,3300,label=substitute(paste(bold(('b')))), col="black",xpd=NA)
title(xlab="Proportion of biotic samples", line=2)

hist(prop_B_Cont,col="#33A02C",xlab="",ylab="", main="",xlim=c(0,1), ylim=c(0,3300))
text(-0.15,3300,label=substitute(paste(bold(('c')))), col="black",xpd=NA)
title(xlab="Proportion of contemporary biotic samples", line=2)

hist(prop_B_Alt,col="#1F78B4",xlab="",ylab="", main="", xlim=c(0,1), ylim=c(0,3300))
text(-0.15,3300,label=substitute(paste(bold(('d')))), col="black",xpd=NA)
title(xlab="Proportion of altered biotic samples", line=2)

mtext("Number of features", side = 2, outer = TRUE, line = -0.5)


##########################################################################################
#Barplot of the lowest ranked importance variables in ranked order 
#Figure 5 (a-d)

par(mfcol = c(4, 1), mar = numeric(4),oma = c(6, 4, .5, .5), mai = c(0.05, 0.2, 0.3, 0.2),
    mgp = c(2, .6, 0),family="sans")

barplot(D_A_array,beside=TRUE,ylim=c(0,1),col="#B15928" ,las=2,axes=FALSE,
        main="Abiotic samples")
text(-2.2,1.12,label=substitute(paste(bold(('a')))),col="black",xpd=NA)
axis(2L)
box()

barplot(D_B_array,beside=TRUE,ylim=c(0,1),col="#6A3D9A",las=2,axes=FALSE,
        main="Biotic samples")
text(-2.2,1.12,label=substitute(paste(bold(('b')))), col="black",xpd=NA)
axis(2L)
box()

barplot(D_B_Cont_array,beside=TRUE,ylim=c(0,1),col="#33A02C" ,las=2,axes=FALSE,
        main="Contemporary biotic samples")
text(-2.2,1.12,label=substitute(paste(bold(('c')))), col="black",xpd=NA)
axis(2L)
box()

barplot(D_B_Alt_array,beside=TRUE,names.arg=Av[1:DIM],ylim=c(0,1),col="#1F78B4",
        las=2,axes=FALSE,main="Altered biotic samples")
text(-2.2,1.12,label=substitute(paste(bold(('d')))), col="black",xpd=NA)
axis(2L)
box()

mtext("Proportion of samples", side = 2, outer = TRUE, line = 2.2,ps=8)
mtext("Importance variables", side = 1, outer = TRUE, line = 4.8,ps=8)

#########################################################################################
#Boxplot of the normalized intensity values for samples containing the first
#48 ranked importance features
#Figure 6

add_legend <- function(...) {
  opar = par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
  mar=c(0, 7, 0, 0),ps=8, family="sans", new=TRUE)
  on.exit(par(opar))
  plot.new()
  legend(...)
}

par(mfrow = c(7, 7), mar = numeric(4),oma = c(5, 4, .5, .5), mai = c(0.05, 0.2, 0.3, 0.2),
    mgp = c(2, .6, 0),family="sans", ps=8)

for(i in 1:42){
  boxplot(Dbar[,i]~Type,data=Dbar,ylab="",main="",xlab="",
          col=c("#B15928","#6A3D9A","#1F78B4", "#33A02C"),axes=FALSE)
  title(main=paste( Av[i]))
  axis(2L)
  box()
}
for(i in 43:48){
  boxplot(Dbar[,i]~Type,data=Dbar,ylab="",main="",xlab="",
          col=c("#B15928","#6A3D9A","#1F78B4", "#33A02C"),las=2)
  title(main=paste( Av[i]))
  box()
}

mtext("Normalized intensity", side = 2, outer = TRUE, line =1)

add_legend("topleft",legend=c("A: abiotic","B: biotic", "B(alt.): altered biotic",
                              "B(cont.): contemporary biotic"),
           col=c("#B15928","#6A3D9A","#1F78B4","#33A02C"),horiz=TRUE,bty = "n",
           fill=c("#B15928","#6A3D9A","#1F78B4","#33A02C"),text.width=c(0.2,0.2,0.2,0.2))

################################################################################
#Plot scan#:m/z:intensity for the 48 ranked importance variables with sample names for
#each point
#Figure 7
fig=plot_ly(data=data3D[ind_all,], x=~Scan,y=~Mass_to_charge_ratio,
            z=~Intensity, type="scatter3d", mode="markers",marker=list(size = 3), color=~Type,colors=c("#B15928","#33A02C","#1F78B4"),text = ~paste("", NameS))
fig=fig %>% layout(legend = list(x = 0.4, y = 0.95,orientation = 'h',
                                 itemsizing='constant',bgcolor ="ghostwhite"))
fig





