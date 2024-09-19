## "Creating graphs for gas chromatography-mass spectrometry data using
## logistic regression with elastic net penalty and the Benjamini-Hochberg procedure "

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
## In this file we create graphs regarding the significant variables from
## logistic regression with elastic net penalty and the Benjamini-Hochberg procedure.


library(dplyr)         # for dataframe computation
library(MALDIquant)    # for chemometrics processing
library(caret)         # for machine learning
library(mlr3)          # for machine learning
library("mlr3verse")   # for machine learning
library("mlr3learners")# for machine learning
library("mlr3tuning")  # for machine learning
library("data.table")  # for rbindlist
library(glmnet)        # for logistic regression with elastic net penalty 
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
## by Dahna, A., datascience+, 03/09/2019
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

#####################################################################################
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
        log_int=ifelse(M_new[[i]][,j+1]>0, log10(M_new[[i]][,j+1]),0)
	       S[[j]] = createMassSpectrum(mass=seq(1,MM,1), log_int,
	       metaData=list(name="Chrom"))  
	 }
   chrom = smoothIntensity(S, method="MovingAverage",halfWindowSize=5)#smooth the data
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

#Scan# and mass spectrum for each sample
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


  
#########################################################################################
#Tune alpha in the elastic net model on the 150 training data

set.seed(99)
resampling_cv10 = rsmp("cv", folds = 10)  #For 10-fold cross validation
resampling_cv10$instantiate(task)

#Create the learner
learner = lrn("classif.cv_glmnet", id = "glmnet", predict_type="prob")

graph = po("removeconstants")  %>>% po("scale")  %>>% learner
plot(graph)

graph_learner = as_learner(graph)
#graph_learner$param_set$ids()

#glmnet
graph_learner$param_set$values$glmnet.alpha = to_tune(p_dbl(0,1))

graph_learner$id = "graph_learner"

RS = tnr("random_search")

future::plan("multisession")

#using stratified random sampling inside each fold
instance = tune(
  tuner = RS,
  task =  task,
  learner = graph_learner,
  resampling = resampling_cv10,
  measure = msr("classif.ce"),
  term_evals = 20
)

instance$result_y
instance$result

######################################################################################
#Tune the overall penalty parameter in the elastic net model
set.seed(99)
dim_en = dim(training_transformed)[2]-1
cv_lambda=cv.glmnet(as.matrix(training_transformed[,1:dim_en]),y,
                    alpha=instance$result$glmnet.alpha,family="binomial")
plot(cv_lambda)
best_lambda=cv_lambda$lambda.min

elastic_mod=glmnet(as.matrix(training_transformed[,1:dim_en]),y,
                   alpha=instance$result$glmnet.alpha,family="binomial")
elastic_coef = predict(elastic_mod, type = "coefficients",
                       alpha=instance$result$glmnet.alpha, 
                       s = best_lambda,family="binomial")

ind_p = which(c(elastic_coef[,1] > 0))
ind_n = which(c(elastic_coef[,1] < 0))

var_p=colnames(training_transformed[,ind_p - 1])
var_n = colnames(training_transformed[,ind_n - 1])
Av=c(var_p,var_n) #non-zero coefficients


K=(dim(training_transformed)[2]-1)

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

#######Calculate the proportion of samples that have the significant variables sorted
#######by abiotic, biotic contemporary, and altered biotic 
calc_proportion=function(Av) {
#Intensity values for the variables in the vector Av for each species
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

####################################################################################
#Split the vector Xscan#.m/z into the scan# component and the m/z value component

#K is the number of features 
#Omit Perc_nonzero feature
datadf=colnames(training_transformed[,1:(dim(training_transformed)[2]-2)])

#Splitting a vector string based on code found in
#"Splitting Strings in R programming – strsplit() method":
# https://www.geeksforgeeks.org/splitting-strings-in-r-programming-strsplit-method/
datadf_st = strsplit(datadf, split = "[.]+")
datadf_st2=list()
for (i in 1:length(datadf_st)){
datadf_st2[[i]] = strsplit(datadf_st[[i]],split = '""')
}

#Get the first element of a list based on code found in
#"R list get first item of each element":
#https://stackoverflow.com/questions/44176908/r-list-get-first-item-of-each-element , 
#stack overflow
datadf_st21= unlist(sapply(datadf_st2, function(x) x[1]))

#gsub based on code found in "Remove Character From String in R":
#https://sparkbyexamples.com/r-programming/remove-character-from-string-in-r/ 
#by Nelamali,N., March 27, 2024
xx = as.numeric(gsub('[X]','',datadf_st21)) #Scan numbers
yy= as.numeric(unlist(sapply(datadf_st2, function(x) x[2]))) #m/z values

######################################################################################
#Split the scan# and m/z values for the significant variables into two separate
#components for different number of features.

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

##################################################################################################
#Benjamini-Hochberg (BH) method to determine significant different variables between
#biotic and abiotic species.
#We used code (modified) for BH found in 
#James et al., "An introduction to statistical learning" 2ed, Springer, 2023
m = ncol(training_transformed)-1
x1 = training_transformed  %>% filter(y=="A")
x2 = training_transformed %>% filter(y=="B")

n1=nrow(x1)
n2=nrow(x2)

B=1000
set.seed(99)

#Using the permutation test to create the null distribution for the t-statistic
 TT = rep(NA, m)
 TT.star = matrix(NA, ncol = m, nrow = B)
 for(j in 1:m){
  	TT[j] = t.test(x1[,j],x2[,j],var.equal = FALSE)$statistic
 	for (b in 1:B){
 		sam = sample(c(x1[,j], x2[,j]))
		 TT.star[b,j] = t.test(sam[1:n1],sam[(n1 + 1):(n1 + n2)], var.equal = FALSE)$statistic
 	}
 }

c = sort(abs(TT))
FDR <- RR <- VV <- rep(NA , m )
 for (j in 1:m){
    R = sum(abs(TT) >= c[j])
    V = sum(abs(TT.star) >= c[j])/B
    RR[j] = R
    VV[j] = V
    FDR[j] = V/R
 }

max(RR[FDR <= 0.001])
index=abs(TT) >= min (c[FDR < 0.001]) 

head(training_transformed[,index]) #Significant different features

Av2=colnames(training_transformed[,index])
#######Calculate the proportion of samples that have the significant variables sorted

data_stat2=calc_proportion(Av2)
#Split the scan# and m/z values for the significant variables into two
#separate components for different number of features.

data_split2=function(Av2){
datadf3=Av2

datadf_st3 = strsplit(datadf3, split = "[.]+")
datadf_st23=list()
for (i in 1:length(datadf_st3)){
datadf_st23[[i]] = strsplit(datadf_st3[[i]],split = '""')
}

datadf_st213= unlist(sapply(datadf_st23, function(x) x[1]))

xx3 = as.numeric(gsub('[X]','',datadf_st213))
yy3= as.numeric(unlist(sapply(datadf_st23, function(x) x[2])))
return(list(xx3=xx3,yy3=yy3))
}
xx3=as.numeric(unlist(data_split2(Av2)[1]))
yy3=as.numeric(unlist(data_split2(Av2)[2]))

##########################################################################################
#Determine the proportion of samples that contain the significant features with color coding.
          
group.col = numeric(length(Av))
		for (k in 1:length(Av)){
			if (data_stat[k,1]>0.5  && data_stat[k,3]<=0.5){
			group.col[k]= 1
			}else if (data_stat[k,2]>0.5 && data_stat[k,3]<=0.5) {
			group.col[k]= 2
			}else if (data_stat[k,3]>0.5 && data_stat[k,2]<=0.5 && data_stat[k,1]<=0.5){ 
			group.col[k]= 3
			}else if (data_stat[k,2]>0.5 && data_stat[k,3]>0.5) {
			group.col[k]= 4
                  }else if (data_stat[k,1]>0.5 && data_stat[k,3]>0.5) {
			group.col[k]= 5
                  }else {
			group.col[k]= 6
			}	
}


group.col2 = numeric(length(Av2))
		for (k in 1:length(Av2)){
			if (data_stat2[k,1]>0.5  && data_stat2[k,3]<=0.5){
			group.col2[k]= 1
			}else if (data_stat2[k,2]>0.5 && data_stat2[k,3]<=0.5) {
			group.col2[k]= 2
			}else if (data_stat2[k,3]>0.5 && data_stat2[k,2]<=0.5 && data_stat2[k,1]<=0.5){ 
			group.col2[k]= 3
			}else if (data_stat2[k,2]>0.5 && data_stat2[k,3]>0.5) {
			group.col2[k]= 4
                  }else if (data_stat2[k,1]>0.5 && data_stat2[k,3]>0.5) {
			group.col2[k]= 5
                  }else {
		  group.col2[k]= 6
			}	
}

 

####################
data_stat3=calc_proportion(datadf)

group.col3 = numeric(length(datadf))
for (i in 1:length(datadf)){
	if (data_stat3[i,1]>data_stat3[i,2]&& data_stat3[i,1]>data_stat3[i,3]) {
	group.col3[i]= 1
	}else if (data_stat3[i,2]> data_stat3[i,3]&& data_stat3[i,2]>data_stat3[i,1]) {
	group.col3[i]= 2
	} else {
	group.col3[i]= 3
	}
}


#Scatterplots of scan#:m/z:intensity for the significant different variables 

#Plotting legend outside plot is based on code found in
#"Plot a legend outside of the plotting area in base graphics?: "https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics ,
#stack overflow
#modified
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 7, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot.new()
  legend(...)
}


par(mfrow=c(2,2),xpd=T,family="sans", ps=9, mai = c(0.5, 0.5, 0.3, 0.3), 
oma = c(1.3, 1.5, 0, 0))

plot(xx, yy, col=c("#B15928","#33A02C","#1F78B4")[group.col3],xlab="",ylab="",
     pch=c(19, 19,19)[group.col3],main="",cex=0.57)
text(100,200,label=substitute(paste(bold(('a')))), col="black")

plot(xx2,yy2,col=c("#B15928","#33A02C","#1F78B4","#6A3D9A","#A6CEE3","#CAB2D6")[group.col],xlab="",ylab="",
     pch=c(19, 19,19,19,19,19)[group.col],main="", cex=0.85)
text(100,200,label=substitute(paste(bold(('b')))), col="black")
legend(c(0,200),c("Abiotic","Biotic (cont.)","Biotic (alt.)","Biotic (cont. and alt.)",
                  "Abiotic and biotic (alt.)",expression(""<="50% of the samples")),col=c("#B15928","#33A02C","#1F78B4","#6A3D9A","#A6CEE3","#CAB2D6"),
       cex=0.85, pch=c(19,19,19,19,19,19),bty = "n")

plot(xx3,yy3,col=c("#B15928","#33A02C","#1F78B4","#6A3D9A","#A6CEE3","#CAB2D6")[group.col2],xlab="",ylab="",
     pch=c(19,19,19,19,19,19)[group.col2],main="",cex=0.85)
text(202,198,label=substitute(paste(bold(('c')))), col="black")
legend(c(0,197),c("Biotic (cont.)","Biotic (alt.)","Biotic (cont. and alt.)",
                  expression(""<="50% of the samples")),col=c("#33A02C","#1F78B4","#6A3D9A","#CAB2D6"),cex=0.85,
       pch=c(19,19,19,19),bty = "n")
mtext("Scan #", side = 1,ps=8, outer = TRUE, line = 0)
mtext(expression(italic("m/z")), side = 2,ps=8, outer = TRUE, line = -1)

add_legend("topleft",legend=c("Abiotic","Biotic (cont.)","Biotic (alt.)"),
           col=c("#B15928","#33A02C","#1F78B4"),cex=0.8, pch=c(19,19,19),horiz=TRUE,bty = "n")


 




