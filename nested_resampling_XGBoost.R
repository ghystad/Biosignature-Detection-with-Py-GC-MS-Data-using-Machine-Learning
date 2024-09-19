## "Nested resampling for gas chromatography-mass spectrometry data analyzed
## by extreme gradient boosting classification method"

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
## EXtreme Gradient Boosting (XGBoost) is used as the machine learning method.
## We apply nested resampling to obtain the generalized prediction performance.


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

##########################################################################################
#Nested resampling, XGBoost

set.seed(99)
inner_loop = rsmp("cv", folds = 10) # inner loop for nested CV
outer_loop = rsmp("cv", folds = 5)  # outer loop for nested CV
outer_loop$instantiate(task)

#Create the learner
learner = lrn("classif.xgboost", id = "xboost", predict_type="prob")

graph = po("removeconstants")  %>>% learner
plot(graph)

graph_learner = as_learner(graph)
#graph_learner$param_set$ids()

#XGboost
graph_learner$param_set$values$xboost.nrounds = 
  to_tune(p_int(80,1000))
graph_learner$param_set$values$xboost.min_child_weight = 
  to_tune(p_int(1,6))
graph_learner$param_set$values$xboost.gamma =
  to_tune(p_dbl(1e-5, 1e0, logscale = TRUE))
graph_learner$param_set$values$xboost.subsample = 
  to_tune(p_dbl(0.8,1))
graph_learner$param_set$values$xboost.colsample_bytree  =
  to_tune(p_dbl(0.5,1))
graph_learner$param_set$values$xboost.alpha = 
  to_tune(p_dbl(0.5,1))

#graph_learner$param_set$values$xboost.lambda =
#to_tune(p_dbl(0,1))

graph_learner$id = "graph_learner"

future::plan("multisession")

at = AutoTuner$new(
  learner = graph_learner,
  resampling = inner_loop,
  measure = msr("classif.ce"),
  terminator = trm("evals", n_evals = 20),
  tuner = tnr("random_search"),
  store_models = TRUE)


rr = resample(task = task, learner = at, resampling = outer_loop, 
              store_models = TRUE)
#extract_inner_tuning_results(rr)

#Unbiased prediction performance (testing error)
rr$aggregate()

rr$aggregate(msr("classif.auc"))





