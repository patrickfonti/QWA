#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
# GENERAL DESCRIPTION:
# 
# This scripts summarizes the annual stats data into overall annual summary
# files. In addition, specific files for each combination of parameter-aggregation 
# approach (mean, median, 25th quantile, 75th quantile)-resolution are created.
#    
# INCLUDED PARAMETERS: 
# 
# MRW - ring width (integrating tangential fluctuations) [µm]
# LA - lumen area [µm2]
# TCA - total cell area (lumen area + cell wall area) [µm2]
# DRAD - radial lumen diameter [µm]
# DTAN - tangential lumen diameter [µm]
# CWA - cell wall area (unit: cell)  [µm2]
# CWAACC - accumulated cell wall area (unit: sum of each intra-ring band) [µm2]
# CWTALL - average cell wall thickness (mean of tangential and radial) [µm]
# CWTALL.ADJ - adjusted average cell wall thickness (earlywood: tangential, latewood: mean of tangential and radial) [µm]
# CWTTAN - tangential cell wall thickness [µm]
# CWTRAD - radial cell wall thickness [µm]
# RTSR - Mork's index (radial thickness-to-span ratio: 4*CWTTAN/DRAD) [-]
# CTSR - adjusted Mork's index (circular thickness-to-span ratio: 4*CWTALL/D, where D is the tracheid diameter assuming a circular lumen) [-]
# DCWT - anatomical density based on cWT (assuming a circular lumen area with homogenous cell wall thickness around it; for latewood-like cells, taking CWTALL, for earlywood-like cells, taking CWTTAN to avoid pit artefacts) [-]
# DCWA - anatomical density based on CWA (CWA/(CWA+LA)) [-]
# DCWT2 - special anatomical density based on CWT (CWTRAD/DRAD) [-]
# TB2 - cell wall reinforcement index ((2*CWTTAN/DTAN)^2) [-]
# KH - theoretical hydraulic conductivity as approximated by Poiseuille's law and adjusted to elliptical tubes [m^4*s^-1*MPa^-1]
# DH - mean hydraulic diameter ((S(D^4)/n)^0.25, where S is the sum, D is the hydraulic cell diameter, n the number of cells) [µm]
# N.BAND - number of considered cells in intra-ring band; some of the cells in the last band in each year are in common with those in the second last band [#]
# N.TRUE - number of considered cells in intra-ring band; no overlapping / shared cells between last and second last band in each year [#]
#
# 
# v1.4, 23 August 2019
# 
# (c) Georg von Arx
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#


### 1. Clean and load library ####
rm(list=ls()) # clean desk 
library(reshape2)


### 2. Define directory containing all the summary data files (only on top level, no sub-directories!) ####
# setwd("D:/TrainingSchools_QWA/NADEF_QWA_2019_Cody/Analysis/R Analysis/Series")
setwd("G:/____DensRecalc20190822/Series")
homedir <- getwd()


### 3. Get file names and for each file the wood piece identifier and the used intra-annual resolution ####
(stats.files <- grep("AnnualStats", list.files(path=homedir), value=TRUE))   #find the annual stats files
(stats.files <- grep("Summary", stats.files, value=TRUE, invert=TRUE))   #find the annual stats files
woodid.files <- data.frame(MU=rep(NA,length(stats.files)))
reso.files <- data.frame(MU=rep(NA,length(stats.files)))
aggr.files <- data.frame(MU=rep(NA,length(stats.files)))

for (i in c(1:length(stats.files)))
{
  ipos <- lapply(strsplit(stats.files[i], ''), function(x) which(x=='_'))   #get position of all "_"
  # ipos <- unlist(ipos)[(length(unlist(ipos))-1)]   #get position of last "_"    
  ipos <- unlist(ipos)[(length(unlist(ipos))-2)]   #get position of second last "_"    
  woodid.files[i,1] <- substr(stats.files[i],1 , ipos-1)   #extract string after last "_"
  
  ipos <- lapply(strsplit(stats.files[i], ''), function(x) which(x=='_'))   #get position of all "_"
  ipos <- unlist(ipos)[(length(unlist(ipos)))]   #get position of last "_"    
  reso.files[i,1] <- substr(stats.files[i],ipos+1, nchar(stats.files[i]))   #extract string after last "_"
  reso.files[i,1] <- sub("mu.txt", "", reso.files[i,1])   #extract full image code  
  
  ipos <- lapply(strsplit(stats.files[i], ''), function(x) which(x=='_'))   #get position of all "_"
  ipos <- unlist(ipos)[(length(unlist(ipos)))]   #get position of second last "_"    
  aggr.files[i,1] <- substr(stats.files[i],1 , ipos-1)   #extract string after last "_"
  ipos <- lapply(strsplit(aggr.files[i,1], ''), function(x) which(x=='_'))   #get position of all "_"
  ipos <- unlist(ipos)[(length(unlist(ipos)))]   #get position of last "_"    
  aggr.files[i,1] <- substr(aggr.files[i,1], ipos+1, nchar(aggr.files[i,1]))   #extract string after last "_"
}


t <- Sys.time()
stats <- NULL

### 4. Process data from each combination of intra-ring position - band aggregation method ####
### 4.1 Loop through each intra-ring resolution (µm) ####
for (r in c(1:length(unique(reso.files$MU))))   
{
  print(paste("(", r, "/", length(unique(reso.files$MU)), ") - ", "Processing resolution: ", unique(reso.files$MU)[r], " µm", sep=""))  
  
  ### 4.1.1 Loop through each band-aggregation method (mean, median, q25, q75) ####
  for (s in c(1:length(unique(aggr.files$MU))))   
  {
    print(paste("  (", LETTERS[s], ") ", "Processing aggregation method: ", unique(aggr.files$MU)[s], sep=""))
    
    ### 4.1.1.1 Loop through each annual statistics file ####
    for (i in c(1:length(stats.files)))   #for each file
    {  
      if (aggr.files$MU[i]==unique(aggr.files$MU)[s])
      {
        if (reso.files$MU[i]==unique(reso.files$MU)[r])
        { 
          ### 4.1.1.1.1 Create overall summary file #### 
          print(paste("     (", make.unique(rep(letters, length.out=100), sep='')[i], ") ", "Processing wood piece: ", woodid.files[i,1], sep=""))
          if (is.null(stats))
          {
            stats <- read.table(stats.files[i], header=T)   #annual stats data 
            stats <- stats[!is.na(stats$YEAR),]
          } else
          {
            df <- NULL
            df <- read.table(stats.files[i], header=T)   #annual stats data
            df <- df[!is.na(df$YEAR),]        
            stats <- rbind(stats, df)
          }
        }   
      }
    }
  
    # library(dplyr)
    # stats2 <- distinct(stats, WOODID, YEAR)   #eliminate last-band information from lower part (with original bands)
    # d1 <- stats[,c(1:4)]
    # d2 <- stats2[,c(1:4)]
    # d3 <- full_join(d1, d2)
    # duplicated(stats) | duplicated(stats[nrow(stats):1, ])[nrow(stats):1]
    # d4 <- which(duplicated(stats) | duplicated(stats, fromLast = TRUE)==TRUE)
    # 
    # library(data.table)
    # DT <- data.table(stats, key = c("WOODID", "YEAR"))
    # DT[unique(DT[duplicated(DT)]),which=T]
    # DT[,count := .N,by = list(WOODID,YEAR)][count>1, which=T]
    
    write.table(stats, file=paste("Summary_AnnualStats_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
  
    ### 4.1.1.1.2 Create wide formate files ####  
    ### 4.1.1.1.2.1 MRW ####
    ### MRW
    df <- stats[,c("WOODID", "YEAR" ,"MRW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MRW")
    write.table(df, file=paste("MRW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### EWW
    df <- stats[,c("WOODID", "YEAR" ,"EWW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="EWW")
    write.table(df, file=paste("EWW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### LWW
    df <- stats[,c("WOODID", "YEAR" ,"LWW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="LWW")
    write.table(df, file=paste("LWW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### 4.1.1.1.2.2 N.BAND ####
    ### N.BAND
    df <- stats[,c("WOODID", "YEAR" ,"N.BAND")]
    df <- dcast(df, YEAR ~ WOODID, value.var="N.BAND")
    write.table(df, file=paste("N.BAND_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### N.BAND.EW
    df <- stats[,c("WOODID", "YEAR" ,"N.BAND.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="N.BAND.EW")
    write.table(df, file=paste("N.BAND.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### N.BAND.LW
    df <- stats[,c("WOODID", "YEAR" ,"N.BAND.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="N.BAND.LW")
    write.table(df, file=paste("N.BAND.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXN.BAND
    df <- stats[,c("WOODID", "YEAR" ,"MXN.BAND")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXN.BAND")
    write.table(df, file=paste("MXN.BAND_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
   
    ### MNN.BAND
    df <- stats[,c("WOODID", "YEAR" ,"MXN.BAND")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXN.BAND")
    write.table(df, file=paste("MXN.BAND_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
  
    ### 4.1.1.1.2.3 N.TRUE ####
    ### N.TRUE
    df <- stats[,c("WOODID", "YEAR" ,"N.TRUE")]
    df <- dcast(df, YEAR ~ WOODID, value.var="N.TRUE")
    write.table(df, file=paste("N.TRUE_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### N.TRUE.EW
    df <- stats[,c("WOODID", "YEAR" ,"N.TRUE.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="N.TRUE.EW")
    write.table(df, file=paste("N.TRUE.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### N.TRUE.LW
    df <- stats[,c("WOODID", "YEAR" ,"N.TRUE.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="N.TRUE.LW")
    write.table(df, file=paste("N.TRUE.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXN.TRUE
    df <- stats[,c("WOODID", "YEAR" ,"MXN.TRUE")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXN.TRUE")
    write.table(df, file=paste("MXN.TRUE_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNN.TRUE
    df <- stats[,c("WOODID", "YEAR" ,"MXN.TRUE")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXN.TRUE")
    write.table(df, file=paste("MXN.TRUE_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### 4.1.1.1.2.4 LA ####
    ### MLA
    df <- stats[,c("WOODID", "YEAR" ,"MLA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MLA")
    write.table(df, file=paste("MLA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### MLA.EW
    df <- stats[,c("WOODID", "YEAR" ,"MLA.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MLA.EW")
    write.table(df, file=paste("MLA.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### MLA.LW
    df <- stats[,c("WOODID", "YEAR" ,"MLA.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MLA.LW")
    write.table(df, file=paste("MLA.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXLA
    df <- stats[,c("WOODID", "YEAR" ,"MXLA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXLA")
    write.table(df, file=paste("MXLA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNLA
    df <- stats[,c("WOODID", "YEAR" ,"MNLA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNLA")
    write.table(df, file=paste("MNLA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### 4.1.1.1.2.5 DRAD ####
    ### DRAD
    df <- stats[,c("WOODID", "YEAR" ,"DRAD")]
    df <- dcast(df, YEAR ~ WOODID, value.var="DRAD")
    write.table(df, file=paste("DRAD_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### DRAD.EW
    df <- stats[,c("WOODID", "YEAR" ,"DRAD.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="DRAD.EW")
    write.table(df, file=paste("DRAD.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### DRAD.LW
    df <- stats[,c("WOODID", "YEAR" ,"DRAD.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="DRAD.LW")
    write.table(df, file=paste("DRAD.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXDRAD
    df <- stats[,c("WOODID", "YEAR" ,"MXDRAD")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXDRAD")
    write.table(df, file=paste("MXDRAD_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNDRAD
    df <- stats[,c("WOODID", "YEAR" ,"MNDRAD")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNDRAD")
    write.table(df, file=paste("MNDRAD_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### 4.1.1.1.2.6 DTAN ####
    ### DTAN
    df <- stats[,c("WOODID", "YEAR" ,"DTAN")]
    df <- dcast(df, YEAR ~ WOODID, value.var="DTAN")
    write.table(df, file=paste("DTAN_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### DTAN.EW
    df <- stats[,c("WOODID", "YEAR" ,"DTAN.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="DTAN.EW")
    write.table(df, file=paste("DTAN.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### DTAN.LW
    df <- stats[,c("WOODID", "YEAR" ,"DTAN.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="DTAN.LW")
    write.table(df, file=paste("DTAN.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXDTAN
    df <- stats[,c("WOODID", "YEAR" ,"MXDTAN")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXDTAN")
    write.table(df, file=paste("MXDTAN_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNDTAN
    df <- stats[,c("WOODID", "YEAR" ,"MNDTAN")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNDTAN")
    write.table(df, file=paste("MNDTAN_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### 4.1.1.1.2.7 MTCA ####
    ### MTCA
    df <- stats[,c("WOODID", "YEAR" ,"MTCA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MTCA")
    write.table(df, file=paste("MTCA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### MTCA.EW
    df <- stats[,c("WOODID", "YEAR" ,"MTCA.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MTCA.EW")
    write.table(df, file=paste("MTCA.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### MTCA.LW
    df <- stats[,c("WOODID", "YEAR" ,"MTCA.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MTCA.LW")
    write.table(df, file=paste("MTCA.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXTCA
    df <- stats[,c("WOODID", "YEAR" ,"MXTCA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXTCA")
    write.table(df, file=paste("MXTCA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNTCA
    df <- stats[,c("WOODID", "YEAR" ,"MNTCA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNTCA")
    write.table(df, file=paste("MNTCA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### 4.1.1.1.2.8 MCWA ####
    ### MCWA
    df <- stats[,c("WOODID", "YEAR" ,"MCWA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWA")
    write.table(df, file=paste("MCWA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MCWA.EW
    df <- stats[,c("WOODID", "YEAR" ,"MCWA.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWA.EW")
    write.table(df, file=paste("MCWA.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### MCWA.LW
    df <- stats[,c("WOODID", "YEAR" ,"MCWA.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWA.LW")
    write.table(df, file=paste("MCWA.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXCWA
    df <- stats[,c("WOODID", "YEAR" ,"MXCWA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXCWA")
    write.table(df, file=paste("MXCWA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNCWA
    df <- stats[,c("WOODID", "YEAR" ,"MNCWA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNCWA")
    write.table(df, file=paste("MNCWA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### 4.1.1.1.2.9 CWAACC ####
    ### CWAACC
    df <- stats[,c("WOODID", "YEAR" ,"CWAACC")]
    df <- dcast(df, YEAR ~ WOODID, value.var="CWAACC")
    write.table(df, file=paste("CWAACC_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### CWAACC.EW
    df <- stats[,c("WOODID", "YEAR" ,"CWAACC.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="CWAACC.EW")
    write.table(df, file=paste("CWAACC.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### CWAACC.LW
    df <- stats[,c("WOODID", "YEAR" ,"CWAACC.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="CWAACC.LW")
    write.table(df, file=paste("CWAACC.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXCWAACC
    df <- stats[,c("WOODID", "YEAR" ,"MXCWAACC")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXCWAACC")
    write.table(df, file=paste("MXCWAACC_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNCWAACC
    df <- stats[,c("WOODID", "YEAR" ,"MNCWAACC")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNCWAACC")
    write.table(df, file=paste("MNCWAACC_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### 4.1.1.1.2.10 CWTALL ####
    ### MCWTALL
    df <- stats[,c("WOODID", "YEAR" ,"MCWTALL")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTALL")
    write.table(df, file=paste("MCWTALL_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MCWTALL.EW
    df <- stats[,c("WOODID", "YEAR" ,"MCWTALL.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTALL.EW")
    write.table(df, file=paste("MCWTALL.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MCWTALL.LW
    df <- stats[,c("WOODID", "YEAR" ,"MCWTALL.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTALL.LW")
    write.table(df, file=paste("MCWTALL.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MXCWTALL
    df <- stats[,c("WOODID", "YEAR" ,"MXCWTALL")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXCWTALL")
    write.table(df, file=paste("MXCWTALL_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNCWTALL
    df <- stats[,c("WOODID", "YEAR" ,"MNCWTALL")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNCWTALL")
    write.table(df, file=paste("MNCWTALL_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### 4.1.1.1.2.11 CWTALL.ADJ ####
    ### MCWTALL.ADJ
    df <- stats[,c("WOODID", "YEAR" ,"MCWTALL.ADJ")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTALL.ADJ")
    write.table(df, file=paste("MCWTALL.ADJ_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MCWTALL.ADJ.EW
    df <- stats[,c("WOODID", "YEAR" ,"MCWTALL.ADJ.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTALL.ADJ.EW")
    write.table(df, file=paste("MCWTALL.ADJ.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MCWTALL.ADJ.LW
    df <- stats[,c("WOODID", "YEAR" ,"MCWTALL.ADJ.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTALL.ADJ.LW")
    write.table(df, file=paste("MCWTALL.ADJ.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MXCWTALL.ADJ
    df <- stats[,c("WOODID", "YEAR" ,"MXCWTALL.ADJ")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXCWTALL.ADJ")
    write.table(df, file=paste("MXCWTALL.ADJ_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNCWTALL.ADJ
    df <- stats[,c("WOODID", "YEAR" ,"MNCWTALL.ADJ")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNCWTALL.ADJ")
    write.table(df, file=paste("MNCWTALL.ADJ_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### 4.1.1.1.2.12 CWTTAN ####
    ### MCWTTAN
    df <- stats[,c("WOODID", "YEAR" ,"MCWTTAN")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTTAN")
    write.table(df, file=paste("MCWTTAN_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MCWTTAN.EW
    df <- stats[,c("WOODID", "YEAR" ,"MCWTTAN.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTTAN.EW")
    write.table(df, file=paste("MCWTTAN.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MCWTTAN.LW
    df <- stats[,c("WOODID", "YEAR" ,"MCWTTAN.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTTAN.LW")
    write.table(df, file=paste("MCWTTAN.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MXCWTTAN
    df <- stats[,c("WOODID", "YEAR" ,"MXCWTTAN")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXCWTTAN")
    write.table(df, file=paste("MXCWTTAN_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNCWTTAN
    df <- stats[,c("WOODID", "YEAR" ,"MNCWTTAN")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNCWTTAN")
    write.table(df, file=paste("MNCWTTAN_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### 4.1.1.1.2.13 CWTRAD ####
    ### MCWTRAD
    df <- stats[,c("WOODID", "YEAR" ,"MCWTRAD")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTRAD")
    write.table(df, file=paste("MCWTRAD_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MCWTRAD.EW
    df <- stats[,c("WOODID", "YEAR" ,"MCWTRAD.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTRAD.EW")
    write.table(df, file=paste("MCWTRAD.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MCWTRAD.LW
    df <- stats[,c("WOODID", "YEAR" ,"MCWTRAD.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCWTRAD.LW")
    write.table(df, file=paste("MCWTRAD.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXCWTRAD
    df <- stats[,c("WOODID", "YEAR" ,"MXCWTRAD")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXCWTRAD")
    write.table(df, file=paste("MXCWTRAD_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNCWTRAD
    df <- stats[,c("WOODID", "YEAR" ,"MNCWTRAD")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNCWTRAD")
    write.table(df, file=paste("MNCWTRAD_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### 4.1.1.1.2.14 RTSR ####
    ### MRTSR
    df <- stats[,c("WOODID", "YEAR" ,"MRTSR")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MRTSR")
    write.table(df, file=paste("MRTSR_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MRTSR.EW
    df <- stats[,c("WOODID", "YEAR" ,"MRTSR.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MRTSR.EW")
    write.table(df, file=paste("MRTSR.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MRTSR.LW
    df <- stats[,c("WOODID", "YEAR" ,"MRTSR.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MRTSR.LW")
    write.table(df, file=paste("MRTSR.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXRTSR
    df <- stats[,c("WOODID", "YEAR" ,"MXRTSR")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXRTSR")
    write.table(df, file=paste("MXRTSR_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNRTSR
    df <- stats[,c("WOODID", "YEAR" ,"MNRTSR")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNRTSR")
    write.table(df, file=paste("MNRTSR_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### 4.1.1.1.2.15 CTSR ####
    ### MCTSR
    df <- stats[,c("WOODID", "YEAR" ,"MCTSR")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCTSR")
    write.table(df, file=paste("MCTSR_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MCTSR.EW
    df <- stats[,c("WOODID", "YEAR" ,"MCTSR.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCTSR.EW")
    write.table(df, file=paste("MCTSR.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MCTSR.LW
    df <- stats[,c("WOODID", "YEAR" ,"MCTSR.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MCTSR.LW")
    write.table(df, file=paste("MCTSR.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXCTSR
    df <- stats[,c("WOODID", "YEAR" ,"MXCTSR")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXCTSR")
    write.table(df, file=paste("MXCTSR_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNCTSR
    df <- stats[,c("WOODID", "YEAR" ,"MNCTSR")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNCTSR")
    write.table(df, file=paste("MNCTSR_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### 4.1.1.1.2.16 DCWT ####
    ### MDCWT
    df <- stats[,c("WOODID", "YEAR" ,"MDCWT")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDCWT")
    write.table(df, file=paste("MDCWT_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MDCWT.EW
    df <- stats[,c("WOODID", "YEAR" ,"MDCWT.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDCWT.EW")
    write.table(df, file=paste("MDCWT.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### MDCWT.LW
    df <- stats[,c("WOODID", "YEAR" ,"MDCWT.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDCWT.LW")
    write.table(df, file=paste("MDCWT.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXDCWT
    df <- stats[,c("WOODID", "YEAR" ,"MXDCWT")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXDCWT")
    write.table(df, file=paste("MXDCWT_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNDCWT
    df <- stats[,c("WOODID", "YEAR" ,"MNDCWT")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNDCWT")
    write.table(df, file=paste("MNDCWT_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### 4.1.1.1.2.17 DCWT2 ####
    ### MDCWT2
    df <- stats[,c("WOODID", "YEAR" ,"MDCWT2")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDCWT2")
    write.table(df, file=paste("MDCWT2_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)
    
    ### MDCWT2.EW
    df <- stats[,c("WOODID", "YEAR" ,"MDCWT2.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDCWT2.EW")
    write.table(df, file=paste("MDCWT2.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)
    
    ### MDCWT2.LW
    df <- stats[,c("WOODID", "YEAR" ,"MDCWT2.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDCWT2.LW")
    write.table(df, file=paste("MDCWT2.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)
    
    ### MXDCWT2
    df <- stats[,c("WOODID", "YEAR" ,"MXDCWT2")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXDCWT2")
    write.table(df, file=paste("MXDCWT2_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)
    
    ### MNDCWT2
    df <- stats[,c("WOODID", "YEAR" ,"MNDCWT2")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNDCWT2")
    write.table(df, file=paste("MNDCWT2_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)
    
    ### 4.1.1.1.2.18 DCWA #### 
    ### MDCWA
    df <- stats[,c("WOODID", "YEAR" ,"MDCWA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDCWA")
    write.table(df, file=paste("MDCWA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### MDCWA.EW
    df <- stats[,c("WOODID", "YEAR" ,"MDCWA.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDCWA.EW")
    write.table(df, file=paste("MDCWA.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### MDCWA.LW
    df <- stats[,c("WOODID", "YEAR" ,"MDCWA.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDCWA.LW")
    write.table(df, file=paste("MDCWA.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXDCWA
    df <- stats[,c("WOODID", "YEAR" ,"MXDCWA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXDCWA")
    write.table(df, file=paste("MXDCWA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNDCWA
    df <- stats[,c("WOODID", "YEAR" ,"MNDCWA")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNDCWA")
    write.table(df, file=paste("MNDCWA_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### 4.1.1.1.2.19 TB2 ####    
    ### MTB2
    df <- stats[,c("WOODID", "YEAR" ,"MTB2")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MTB2")
    write.table(df, file=paste("MTB2_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MTB2.EW
    df <- stats[,c("WOODID", "YEAR" ,"MTB2.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MTB2.EW")
    write.table(df, file=paste("MTB2.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MTB2.LW
    df <- stats[,c("WOODID", "YEAR" ,"MTB2.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MTB2.LW")
    write.table(df, file=paste("MTB2.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXTB2
    df <- stats[,c("WOODID", "YEAR" ,"MXTB2")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXTB2")
    write.table(df, file=paste("MXTB2_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNTB2
    df <- stats[,c("WOODID", "YEAR" ,"MNTB2")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNTB2")
    write.table(df, file=paste("MNTB2_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### 4.1.1.1.2.20 KH ####
    ### MKH
    df <- stats[,c("WOODID", "YEAR" ,"MKH")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MKH")
    write.table(df, file=paste("MKH_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   
    
    ### MKH.EW
    df <- stats[,c("WOODID", "YEAR" ,"MKH.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MKH.EW")
    write.table(df, file=paste("MKH.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MKH.LW
    df <- stats[,c("WOODID", "YEAR" ,"MKH.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MKH.LW")
    write.table(df, file=paste("MKH.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXKH
    df <- stats[,c("WOODID", "YEAR" ,"MXKH")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXKH")
    write.table(df, file=paste("MXKH_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNKH
    df <- stats[,c("WOODID", "YEAR" ,"MNKH")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNKH")
    write.table(df, file=paste("MNKH_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### 4.1.1.1.2.21 DH ####
    ### MDH
    df <- stats[,c("WOODID", "YEAR" ,"MDH")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDH")
    write.table(df, file=paste("MDH_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### MDH.EW
    df <- stats[,c("WOODID", "YEAR" ,"MDH.EW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDH.EW")
    write.table(df, file=paste("MDH.EW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)  
    
    ### MDH.LW
    df <- stats[,c("WOODID", "YEAR" ,"MDH.LW")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MDH.LW")
    write.table(df, file=paste("MDH.LW_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MXDH
    df <- stats[,c("WOODID", "YEAR" ,"MXDH")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MXDH")
    write.table(df, file=paste("MXDH_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE) 
    
    ### MNDH
    df <- stats[,c("WOODID", "YEAR" ,"MNDH")]
    df <- dcast(df, YEAR ~ WOODID, value.var="MNDH")
    write.table(df, file=paste("MNDH_", unique(aggr.files$MU)[s], "_", unique(reso.files$MU)[r], "mu.txt", sep=""), row.names = FALSE)   

    stats <- NULL   #reset dataframe    
  }
}

Sys.time() - t


### 5. clean up ####
rm(list=ls(all=TRUE))

