# ------------ General information  ---------------------
# The package "rawDiag" might cause issues during installation (see below), this script has been tested and runs under R version 4.1.1 
# For details see https://github.com/fgcz/rawDiag
#
# Naming of raw files should be as follows for the script to run without problems
# The file name should contain Standard for HeLa standard samples
# For the HeLa standard samples, on the second position in the file name there should be the peptide amount in ng separated by "_" 
# i.e. "210130_10_Standard_TIC_Norm_FMS-45-65_R1.raw" for 10 ng

# ------------ User Interaction needed ---------------------
# Please adjust below how you prepared your Evotips for TIC normalization
amount_prepared <- 150  ## How much volume [ul] have you actually prepared when loading tips (Sample amount + Buffer A as diluent)
sample_taken <- 3      ## How much digest did you take per sample [ul] 
amount_loaded <- 100    ## How much volume of the prepared dilution did you load onto the tips

# Specify the path to the directory containing the raw files, the results will be saved in the same folder
# default is ".", same directory as this script/R project
Path = "."

# More interaction needed in line 82-89


# Check if  Package 'rawDiag' needs installation
# load required packages
if(!any(grepl("rawDiag", installed.packages() )) ) {
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("rawDiag")

  }

library("rawDiag")
library("parallel")
library("dplyr")
library("data.table")
library("ggplot2")

# Select all the Files that are the HELA standard in the Path, if your Raw File naming is not consistent this needs to be adjusted
filenames <- list.files(path=Path, full.names = T)
standardFiles <- filenames[grepl("Standard",filenames, ignore.case = T )&grepl(".raw",filenames)]


# load the raw files into a dataframe
RAW <- mclapply(standardFiles, read.raw, mc.cores=1)
RAW <- as.data.table(plyr::rbind.fill(RAW))


# assign the correct amount of HELA peptides as factor labels based on order of files in standardFiles

RAW$peptideAmount = sapply( RAW$filename, function(x){
  as.numeric(strsplit(x, "_")[[1]][2])
})

# calculate the total TIC per Raw File, logged and unlogged
HeLa_Dilution_2CV <- RAW[ , .(TotalTIC=sum(TIC),TotalTIC_logged=log10(sum(TIC))), by = .(peptideAmount)]



#########
#      Plot to initially show how the regression curve looks like to decide which Data points to use below 
#########

ggplot(HeLa_Dilution_2CV,
       aes(x=TotalTIC_logged, y=log10(peptideAmount)))+
  geom_point(size=4)+
  geom_line(linewidth=1, lty="dashed")+
  scale_x_continuous("\nlog10(Total TIC)")+
  scale_y_continuous("log10(HeLa peptide amount [ng])\n")+
  ggtitle("Calibration curve - first iteration")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))




###      decide which data points to take
#    put to FALSE if you want to take one out
setorder(HeLa_Dilution_2CV, peptideAmount)

Datapoints <- c(FALSE, #10 ng
                TRUE, #31.25 ng
                TRUE, #62.5 ng
                TRUE, #125 ng
                TRUE, #250 ng
                TRUE, #500 ng
                TRUE, #750 ng
                TRUE  #1000 ng
)



# actual regression line on which calculations are based on
HeLa_Dilution_2CV$symbols = ifelse(Datapoints, "1", "2")

regression = ggplot(HeLa_Dilution_2CV, #[Datapoints],
       aes(x=TotalTIC_logged, y=log10(peptideAmount), shape=symbols))+
  geom_point(size=4)+
  geom_line(linewidth=1,lty="dashed")+
  scale_shape_manual(values=c(19,1))+
  scale_x_continuous("\nlog10(Total TIC)")+
  scale_y_continuous("log10(HeLa peptide amount [ng])\n")+
  ggtitle("Calibration curve - second iteration")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        legend.position = "none")


lm <- lm(log10(HeLa_Dilution_2CV$peptideAmount)[Datapoints]~HeLa_Dilution_2CV$TotalTIC_logged[Datapoints])

regression = regression + 
  annotate("text", 
         x = min(HeLa_Dilution_2CV$TotalTIC_logged)+0.5, 
         y = max( log10(HeLa_Dilution_2CV$peptideAmount))-0.5, 
         label = paste0("intercept = ",round(lm$coefficients[1],2),"\n",
                       "slope = ",round(lm$coefficients[2],2),"\n",
                       "R^2 = ",round(summary(lm)$r.squared,4))
         )+
  geom_abline(slope=lm$coefficients[2], intercept=lm$coefficients[1], col="#0065bd", linewidth=1.2)

regression



##########
## Processing of actual samples 

# Select all the Files that are the Samples in the Path
# This always needs to be adjusted to a pattern that matches/identifies your samples

TICNormSamples <- filenames[!grepl("Standard",filenames)&grepl(".raw",filenames)]

#load Raw Files
RAW <- mclapply(TICNormSamples, read.raw, mc.cores=1)
RAW <-as.data.table(plyr::rbind.fill(RAW))


# calculate Total TIC per Sample
Results <-RAW[ ,.(TotalTIC = sum(TIC),loggedTotalTIC = log10(sum(TIC))), by=.(filename)]


# concentration calculation
Results$`Concentration [ng/ul]` <- round((amount_prepared/amount_loaded)*(round(10^((Results$loggedTotalTIC*lm$coefficients[2])+lm$coefficients[1]),2))/sample_taken,2)
Results$Injected_Amount = 10^((Results$loggedTotalTIC*lm$coefficients[2])+lm$coefficients[1])

# pipetting scheme calculation
Results$`Sample [ul]` <- round(((600/100)/Results$Conc)*220,1)
Results$`Buffer_A [ul]` <- round(220-44-Results$Sample,1)
Results$`Procal [ul] (at 10 fmol/ul)` <- 44
Results$`Total [ul]` <- 220


# write the Results and pipetting scheme as a txt  
fwrite(Results,file = file.path(Path,"TIC_Norm_Results.txt"), sep="\t")




### ------------------------------------------------------------------------------------------------------------------------------------------------
### make a PDF with the calibration curve, regression line and min and max values and all TIC chromatograms of all samples to check quality manually
pdf(file=file.path(Path, "Regression and chromatograms.pdf"),width=20, height=10)


# plot the final calibration curve with values of the linear regression 
# and vertical lines showing the maximum and minimum values of actual samples 
# to see if they are still in the linear range or need to be redone with higher/lower dilution

regression +
geom_vline(xintercept = min(Results$loggedTotalTIC), col="#D95117", linewidth=1.2)+
geom_vline(xintercept = max(Results$loggedTotalTIC), col="#D95117", linewidth=1.2)


# plot all TIC chromatograms
for (file in unique(RAW$filename)){
  
  sub <- RAW[filename == file,.(StartTime,TIC)]
  print(ggplot(sub,aes(x=StartTime, y=TIC))+
          geom_line()+
          ggtitle(file)+
          scale_y_continuous("TIC \n")+  
          scale_x_continuous("\nRetention time [min]", breaks = seq(0,12,1))+
          theme_classic()
  )
}
dev.off()

