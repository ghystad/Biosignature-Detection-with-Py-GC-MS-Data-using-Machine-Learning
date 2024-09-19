#Graphs for the preprocessing steps for chromatograms

library(dplyr)         # for dataframe computation
library(MALDIquant)    # for chemometrics processing


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

#Scan number:m/z:intensity values for 134 samples from previous file
M=readRDS("C:/Desktop/Templeton_grant_research/Second_paper_material/Paper/Excel_publication/M.RDS")

#Scan number:m/z:intensity values for 16 samples from previous file
M2=readRDS("C:/Desktop/Templeton_grant_research/Second_paper_material/Paper/Excel_publication/M2.RDS")


M_new=c(M,M2)
N=length(M_new)
species_n2=rep(2,16)
species_names_new=c(species_names,species_names2)

y=c(species_n,species_n2)
y=factor(y,labels=c("A","B"))


MZ=151
i=52
#Grass
S=list()
     	for(j in 1:MZ){
	     S[[j]] = createMassSpectrum(mass=seq(1,MM,1), intensity=unlist(M_new[[i]][,j+1]),
	  metaData=list(name="Chrom")) 
	 }

par(mfrow=c(2,2),family="sans", ps=8, mai = c(0.5, 0.5, 0.3, 0.3), 
oma = c(1.3, 1.5, 0, 0))

#"Chromatogram of Grass (B)"
#m/z=77 and m/z=117
plot(S[[28]],col="#1F78B4", main="Chromatogram of grass ",xlab="",ylab="",
     xlim=c(890,1255),ylim=c(0,70000),xaxt='n',lwd=1.5)
axis(side=1,at=seq(900,1250,50))
text(900,70000,label=substitute(paste(bold(('a')))), col="black")
title(ylab="Intensity", line=2)
lines(S[[68]],col="#6A3D9A",lwd=1.5)
legend(1168,73500,legend=c("m/z=77","m/z = 117"),
       col=c("#1F78B4","#6A3D9A"),fill=c("#1F78B4","#6A3D9A"), bty = "n")

chrom = transformIntensity(S,method="sqrt") #stabilize the variance

chrom = smoothIntensity(chrom, method="MovingAverage",halfWindowSize=5)#smooth the data
plot(chrom[[28]],col="#1F78B4",xlim=c(890,1250),ylim=c(0,200),
     main="SQRT transform and smoothed",xlab="",ylab="",xaxt='n',lwd=1.5)
axis(side=1,at=seq(900,1250,50))
text(900,200,label=substitute(paste(bold(('b')))), col="black")
title(ylab="Intensity", line=2)
legend("topright",legend=c("m/z=77"),col=c("#1F78B4"),fill=c("#1F78B4"),bty = "n")

chrom = removeBaseline(chrom, method="SNIP") #remove the baseline

  #Put the processed chromatograms back into a dataframe
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


noise=estimateNoise(S2[[28]],method=c("MAD"))
noise_SNR=noise[,2]*4

plot(S2[[28]],col="#1F78B4",xlim=c(890,1250),ylim=c(0,0.35),
     main="Baseline removed, peak detection",
     xlab="",ylab="",xaxt='n',lwd=1.5)
points(peaks[[28]],col="#B15928")
lines(noise_SNR,col="#A6CEE3", lwd=1.5)
axis(side=1,at=seq(900,1250,50))
text(900,0.35,label=substitute(paste(bold(('c')))), col="black")
title(ylab="Normalized intensity", line=2)
legend("topright",legend=c("m/z=77"),col=c("#1F78B4"),fill=c("#1F78B4"),bty = "n")



#######################################################################################
i=107
#Peat (B)
S=list()
  for(j in 1:MZ){
	   S[[j]] = createMassSpectrum(mass=seq(1,MM,1), 
	                               intensity=unlist(M_new[[i]][,j+1]),
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

S3=list()
	 for(t in 1:MZ){
	    S3[[t]] = createMassSpectrum(mass=seq(1,MM,1), 
	                                 intensity= processed_mass_dataframe[,t],
	    metaData=list(name="Chrom_normalized"))  
	  }

    peaks3 = detectPeaks(S3, method="MAD", halfWindowSize=20, SNR=4)
#######################################################################################
i=112
#Pumpkin (B)
S=list()
   for(j in 1:MZ){
	    S[[j]] = createMassSpectrum(mass=seq(1,MM,1),
	                                intensity=unlist(M_new[[i]][,j+1]),
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
S4=list()
   for(t in 1:MZ){
	     S4[[t]] = createMassSpectrum(mass=seq(1,MM,1),
	                                  intensity= processed_mass_dataframe [,t],
	     metaData=list(name="Chrom_normalized"))  
	 }
peaks4 = detectPeaks(S4, method="MAD", halfWindowSize=20, SNR=4)


#######################################################################################

plot(S2[[28]],col="#1F78B4",xlim=c(980,1060),ylim=c(0,0.35),
     main="Chromatograms, m/z=77",xlab="",ylab="",xaxt='n',lwd=1.5)
lines(S3[[28]],col="#33A02C",lwd=1.5)
lines(S4[[28]],col="#6A3D9A",lwd=1.5)
axis(side=1,at=seq(980,1060,20))
text(982,0.35,label=substitute(paste(bold(('d')))), col="black")
title(ylab="Normalized Intensity", line=2)

points(peaks[[28]],col="#B15928")
points(peaks3[[28]],col="#B15928")
points(peaks4[[28]],col="#B15928")

legend("topright",legend=c("Grass","Peat","Pumpkin"),
       col=c("#1F78B4","#33A02C","#6A3D9A"),
       fill=c("#1F78B4","#33A02C","#6A3D9A"),bty = "n")
mtext("Scan #", side = 1,ps=8, outer = TRUE, line = 0)





