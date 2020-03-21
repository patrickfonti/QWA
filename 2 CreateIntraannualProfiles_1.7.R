#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
# GENERAL DESCRIPTION:
# 
# This scripts creates intra-annual profiles of anatomical parameters based on 
# summary files (one per radial wood piece) of ROXAS text file output for rings  
# and cells. No data is shared across ring borders, but restricted to the 
# respective target ring.Profile data is plotted and saved to file. In addition, 
# files for means (overall, earlywood, latewood), minimum and maximum values per 
# ring are created.
#
# 
# SOME DETAILS:
# 
# - Measured distances from the ring border for each cell (distance of center 
#   of mass) are re-scaled to account for tangential fluctuations in ring width
#   by the factor (RWi/MRW), where RWi is the ring width at the cell's tangential 
#   position and MRW is the average ring width. 
# - The last tangential intra-ring band of each ring is forced to have the 
#   same width as specified, which usually means that some cells of the second
#   last band will be re-used for aggregations in the last band.
# - For each tangential intra-ring band 4 different types of aggregating are
#   available: mean, median, 25th quantile, 75th quantile. Data and plots are 
#   created for each combination of data aggregation type and resolution.
# 
# 
# CUSTOMIZABLE OPTIONS:
# 
# YR.START, YR.END - target period
# STATS - aggregation method(s) for intra-ring tangential bands (mean, median, 25th and 75th quantile)
# RESO - resolution(s) / band width(s) (µm) of the intra-ring profiles
# SMOOTHER - smoothing spline for the intra-ring profiles
# PLOT.IMG.ID - alternating white-grey background in intra-ring profile plots to indicate when data from a next image is used 
# LINE.COL - color gradient to visualize RESO in intra-ring profile plots 
# MORK.EWLW - threshold of Mork's index to indicate transition from earlywood to latewood
# Several variables to define quality control (QC) features added to the intra-annual lumen area (LA) profiles
# 
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
# TO DO:
#   
# Suppress warnings (suppressWarnings()) when using min/max function and getting message "No non-missing values found in at least one group. Returning '-Inf' for such groups to be consistent with base"
#   
#     
# v1.7, 30 October 2019
# 
# (c) Georg von Arx
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#


### 1. Clean and load libraries ####
rm(list=ls()) # clean desk 
library(plyr)
library(data.table)
library(zoo)
library(dplyr)


### 2. Let user customize some options ####
### 2.1 General options ####
YR.START <- 1978   #minimum potential start year of profiles
YR.END <- 2017   #maximum potentitial ending year of profiles

# STATS <- c("mean","median","q25","q75")   #different methods for aggregation within each intra-ring band  
STATS <- c("median")   #different methods for aggregation within each intra-ring band  

# RESO <- c(2,5,10,20,30,40,50,60,80,100,120,160,200,300,500)   #spatial resolution of profiles (micrometers), i.e. width of intra-ring bands
# RESO <- c(20)   #spatial resolution of profiles (micrometers), i.e. width of intra-ring bands
RESO <- c(20)   #spatial resolution of profiles (micrometers), i.e. width of intra-ring bands

SMOOTHER <- 1   #uneven numbers only! Defines the width of the rolling mean filter with respect to RESO; 1 means no smoothing

PLOT.IMG.ID <- TRUE   #when this option is activated, the background of the intra-annual profiles shows periods for which data is taken from a specific image

LINE.COL <- colorRampPalette(c("red","green","blue"))

MORK.EWLW <- 1


### 2.2 Quality control (QC) options applied to intra-ring lumen area (LA) profiles ####
NCELL.MIN	<- 5   #minimum number of cells per band; if below, the band will be flagged

LA.REF	<- 0.9   #lumen area quantile per image (segments) for horizontal reference line
LA.THR <- 0.9   #maximum allowed deviation in lumen area quantile between neighboring image (segments); if deviating more, the reference quantile line will be highlighted

DUP.RW	<- 0.05   #maximum deviation in ring width of neighboring rings to consider them in duplicate identification   
DUP.KS	<- 0.98   #minimum p-value of Kolmogorov-Smirnoff test between lumen area profiles of neighboring rings to highlight them as potential duplicates (cross-dating mistake!) 

RM.KERNEL	<- 9   #kernel size of running mean smoothing; the resulting smoothed LA profiles serves as a reference to identify (large) outliers 
OUT.THR	<- 1.5   #threshold factor by which LA at a specific point should be larger then the smoothed LA reference profile to flag it as a potential outlier
OUT.QU	<- 0.75   #minimum quantile of smoothed LA reference profile to allow flagging potential outliers


### 3. Define directory containing all the summary data files (only on top level, no sub-directories!) ####
# setwd("D:/TrainingSchools_QWA/NADEF_QWA_2019_Cody/Analysis/R Analysis/Series")
setwd("C:/Users/ferriz/Desktop/roxas1/output")
homedir <- getwd()


### 4. Process data from each radial wood piece ####
(cell.files <- grep("Output_Cells.txt", list.files(path=homedir), value=TRUE))   #find the summary cell files
(ring.files <- grep("Output_Rings.txt", list.files(path=homedir), value=TRUE))   #find the summary ring files

### Get file with information about years/rings to be excluded (simple text file with 2 columns: TREE, YEAR)
exclude <- read.table("_YearsToExclude.txt", header=T)

### Define some functions
RobustMax <- function(x) {if (length(x)>0) max(x) else -Inf}   #to avoid applying max to numeric of length 0 argument
RobustMin <- function(x) {if (length(x)>0) min(x) else -Inf}   #to avoid applying max to numeric of length 0 argument

t <- Sys.time()   #initiate keeping track of processing speed

### Start main loop through all files
for (i in c(2:length(cell.files)))
{
  ### 4.1 Load data ####    
  woodid <- sub("_Output_Cells.txt", "", cell.files[i])     
  print(paste("(", i, "/", length(cell.files), ") - ", "Processing wood piece: ", woodid, sep=""))

  cells <- NULL
  rings <- NULL
  cells <- read.table(cell.files[i], header=T)   #cell data
  rings <- read.table(ring.files[i], header=T)   #ring width data: used for standardization
  cells <- cells[cells$YEAR>min(rings$YEAR,na.rm=T),]   #drop innermost/first year that is not complete
  rings <- rings[rings$YEAR>min(rings$YEAR,na.rm=T),]   #drop innermost/first year that is not complete
  
  rings <- rings[rings$YEAR%in%intersect(unique(rings$YEAR), unique(cells$YEAR)),]   #exclude data from years that are not present in both cells and rings dataframes
  cells <- cells[cells$YEAR%in%intersect(unique(rings$YEAR), unique(cells$YEAR)),]   #exclude data from years that are not present in both cells and rings dataframes
  
  ### Create dataframe with WOODID - IMAGE - YEAR structure to re-install time series integrity if some rings are missing/dropped later on
  ringTemplate <- rings
  ringTemplate$WOODID <- woodid
  colnames(ringTemplate)[1] <- "IMAGE"
  ringTemplate <- ringTemplate[,c("YEAR", "WOODID","IMAGE")]
  ringTemplate$WOODID <- as.factor(ringTemplate$WOODID)
  
  dt <- NULL
  dt <- cells[!duplicated(cells[,c("ID","YEAR")]),c("ID","YEAR")]   #create dataframe with IMAGE-YEAR from cells file
  rownames(dt) <- NULL 
  
  ringTemplate <- merge(ringTemplate,dt,"YEAR",all=TRUE)
  ringTemplate$IMAGE <- ifelse(!is.na(ringTemplate$ID),as.character(ringTemplate$ID),as.character(ringTemplate$IMAGE))  #take IMAGE/ID info from cells in case there is such info
  ringTemplate$IMAGE <- as.factor(ringTemplate$IMAGE)
  ringTemplate$ID <- NULL
  
  
  ### 4.2 Prepare data frame for analysis ####
  cells <- left_join(cells, rings[,c("YEAR", "MRW")], by="YEAR")   #merge the two data files
  cells[] <- lapply(cells, function(x){if (is.numeric(x)) replace(x, x < 0, NA) else x})   #replace error codes (negative values) by NA   
  cells$RW.CELL <- cells$RADDISTR / cells$RRADDISTR * 100   #get ring width at tangential position of each cell
  cells$RADDISTR.ST <- cells$RADDISTR / cells$RW.CELL * cells$MRW   #standardized absolute distance from ring border (accounting for wedging, curvatures and other tangential irregularities of ring width)
  cells <- cells[order(cells$YEAR, cells$RADDISTR.ST),]   #order cells by YEAR and standardized absolute distance from ring border
  if("BEND" %in% colnames(cells)){colnames(cells)[colnames(cells)=="BEND"] <- "TB2"}   #for backwards compatibility: replace BEND by TB2
  cells$WOODID <- as.factor(woodid)   #make sure WOODID is a factor
  

  ### 4.3 Add CWT-based density: assume a circular lumen area with homogenous cell wall thickness around it; for latewood-like cells, take overall average CWT, for earlywood-like cells, only consider CWTTAN, which avoids pit artefacts ####
  cells$RADIUS <- sqrt(cells$LA / pi)
  cells$WA <- ifelse(cells$RTSR < 1, 
                     (cells$RADIUS + cells$CWTTAN)^2 * pi - cells$LA,
                     (cells$RADIUS + cells$CWTALL)^2 * pi - cells$LA)
  cells$DCWT <- cells$WA / (cells$LA + cells$WA)
  

  ### 4.4 Add mean CWT: mean of radial and tangential CWT if Mork index latewood-like, in earlywood-like cells take CWTTAN ####
  cells$CWTALL.ADJ <- ifelse(cells$RTSR < 1, cells$CWTTAN, cells$CWTALL)
  

  ### 4.5 Remove rows where the cell center falls exactly on ring border (mistake/inaccurate!) and correct small interpolation mistakes where relative intra-ring positions >100% to avoid calculations problems later on ####
  cells <- cells[!(cells$RADDISTR==0),]    
  cells$RRADDISTR <- ifelse(cells$RRADDISTR>100, 100, cells$RRADDISTR)  
  
  
  ### 4.6 Loop for each aggregation method ("mean","median","q25","q75") ####
  for (s in c(1:length(STATS)))     
  {  
    print(paste("  (", make.unique(rep(LETTERS, length.out=100), sep='')[s], ") ", STATS[s], ": ", sep=""))  
    
    
    ### 4.6.1 Process data for each intra-ring resolution step ####
    for (r in c(1:length(RESO)))
    {
      print(paste("     (", make.unique(rep(letters, length.out=100), sep='')[r], ") ", "Resolution: ", RESO[r], " microns", sep=""))
      
      
      ### 4.6.1.1 Add column with resolution (micrometers) for each ring; insert row if one band is missing ####
      cells$RADDISTR.BAND <- RESO[r]* round((cells$RADDISTR.ST-RESO[r]/2.01)/RESO[r] , 0)   #add column with resolution (micrometers) for each ring; insert row if one band is missing
  
      
      ### 4.6.1.2 Calculate aggregates for each intra-ring resolution band (standardized absolute distance from ring border) #### 
      if (tolower(STATS[s])=="mean")   #aggregation method = "mean"
        {  
        dt <- NULL
        dt <- data.table(cells)
        setkey(dt, "YEAR", "WOODID", "ID", "RADDISTR.BAND")
        dt <- dt[, list(RRADDISTR=mean(RRADDISTR, na.rm=TRUE)
                        , MRW=RobustMax(MRW)
                        , LA=as.numeric(mean(LA, na.rm=TRUE))
                        , TCA=mean(TCA, na.rm=TRUE)
                        , DRAD=mean(DRAD, na.rm=TRUE)
                        , DTAN=mean(DTAN, na.rm=TRUE)
                        , CWA=mean(CWA, na.rm=TRUE)
                        , CWAACC=sum(CWA, na.rm=TRUE)
                        , CWTALL=mean(CWTALL, na.rm=TRUE)
                        , CWTALL.ADJ=mean(CWTALL.ADJ, na.rm=TRUE)
                        , CWTTAN=mean(CWTTAN, na.rm=TRUE)
                        , CWTRAD=mean(CWTRAD, na.rm=TRUE) 
                        , RTSR=mean(RTSR, na.rm=TRUE)
                        , CTSR=mean(CTSR, na.rm=TRUE) 
                        , DCWT=mean(DCWT, na.rm=TRUE)
                        , DCWA=mean(RWD, na.rm=TRUE)
                        , DCWT2=mean(RWD2, na.rm=TRUE)
                        , TB2=mean(TB2, na.rm=TRUE)
                        , KH=mean(KH, na.rm=TRUE)
                        , N.BAND=sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)
                        , DH=(sum(DH^4)/(sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)))^0.25), 
                 by=key(dt)]
        df0 <- as.data.frame(dt)
      }
      
      if (tolower(STATS[s])=="median")   #aggregation method = "median"
        {  
        dt <- NULL
        dt <- data.table(cells)
        setkey(dt, "YEAR", "WOODID", "ID", "RADDISTR.BAND")
        dt <- dt[, list(RRADDISTR=median(RRADDISTR, na.rm=TRUE)
                        , MRW=RobustMax(MRW)
                        , LA=as.numeric(median(LA, na.rm=TRUE))
                        , TCA=median(TCA, na.rm=TRUE)
                        , DRAD=median(DRAD, na.rm=TRUE)
                        , DTAN=mean(DTAN, na.rm=TRUE)
                        , CWA=median(CWA, na.rm=TRUE)
                        , CWAACC=sum(CWA, na.rm=TRUE)
                        , CWTALL=median(CWTALL, na.rm=TRUE)
                        , CWTALL.ADJ=median(CWTALL.ADJ, na.rm=TRUE)
                        , CWTTAN=median(CWTTAN, na.rm=TRUE)
                        , CWTRAD=median(CWTRAD, na.rm=TRUE) 
                        , RTSR=median(RTSR, na.rm=TRUE)
                        , CTSR=median(CTSR, na.rm=TRUE) 
                        , DCWT=median(DCWT, na.rm=TRUE)
                        , DCWA=median(RWD, na.rm=TRUE)
                        , DCWT2=median(RWD2, na.rm=TRUE)
                        , TB2=median(TB2, na.rm=TRUE)
                        , KH=median(KH, na.rm=TRUE)
                        , N.BAND=sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)
                        , DH=(sum(DH^4)/(sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)))^0.25), 
                 by=key(dt)]
        df0 <- as.data.frame(dt)
      }
      
      if (tolower(STATS[s])=="q25")   #aggregation method = "25th quantile"
        {  
        dt <- NULL
        dt <- data.table(cells)
        setkey(dt, "YEAR", "WOODID", "ID", "RADDISTR.BAND")
        dt <- dt[, list(RRADDISTR=quantile(RRADDISTR, probs=0.25, na.rm=TRUE)
                        , MRW=RobustMax(MRW)
                        , LA=as.numeric(quantile(LA, probs=0.25, na.rm=TRUE))
                        , TCA=quantile(TCA, probs=0.25, na.rm=TRUE)
                        , DRAD=quantile(DRAD, probs=0.25, na.rm=TRUE)
                        , DTAN=quantile(DTAN, probs=0.25, na.rm=TRUE)
                        , CWA=quantile(CWA, probs=0.25, na.rm=TRUE)
                        , CWAACC=sum(CWA, na.rm=TRUE)
                        , CWTALL=quantile(CWTALL, probs=0.25, na.rm=TRUE)
                        , CWTALL.ADJ=quantile(CWTALL.ADJ, probs=0.25, na.rm=TRUE)
                        , CWTTAN=quantile(CWTTAN, probs=0.25, na.rm=TRUE)
                        , CWTRAD=quantile(CWTRAD, probs=0.25, na.rm=TRUE) 
                        , RTSR=quantile(RTSR, probs=0.25, na.rm=TRUE)
                        , CTSR=quantile(CTSR, probs=0.25, na.rm=TRUE) 
                        , DCWT=quantile(DCWT, probs=0.25, na.rm=TRUE)
                        , DCWA=quantile(RWD, probs=0.25, na.rm=TRUE)
                        , DCWT2=quantile(RWD2, probs=0.25, na.rm=TRUE)
                        , TB2=quantile(TB2, probs=0.25, na.rm=TRUE)
                        , KH=quantile(KH, probs=0.25, na.rm=TRUE)
                        , N.BAND=sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)
                        , DH=(sum(DH^4)/(sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)))^0.25), 
                 by=key(dt)]
        df0 <- as.data.frame(dt)
      }
      
      if (tolower(STATS[s])=="q75")   #aggregation method = "75th quantile"
        {  
        dt <- NULL
        dt <- data.table(cells)
        setkey(dt, "YEAR", "WOODID", "ID", "RADDISTR.BAND")
        dt <- dt[, list(RRADDISTR=quantile(RRADDISTR, probs=0.75, na.rm=TRUE)
                        , MRW=RobustMax(MRW)
                        , LA=as.numeric(quantile(LA, probs=0.75, na.rm=TRUE))
                        , TCA=quantile(TCA, probs=0.75, na.rm=TRUE)
                        , DRAD=quantile(DRAD, probs=0.75, na.rm=TRUE)
                        , DTAN=quantile(DTAN, probs=0.75, na.rm=TRUE)
                        , CWA=quantile(CWA, probs=0.75, na.rm=TRUE)
                        , CWAACC=sum(CWA, na.rm=TRUE)
                        , CWTALL=quantile(CWTALL, probs=0.75, na.rm=TRUE)
                        , CWTALL.ADJ=quantile(CWTALL.ADJ, probs=0.75, na.rm=TRUE)
                        , CWTTAN=quantile(CWTTAN, probs=0.75, na.rm=TRUE)
                        , CWTRAD=quantile(CWTRAD, probs=0.75, na.rm=TRUE) 
                        , RTSR=quantile(RTSR, probs=0.75, na.rm=TRUE)
                        , CTSR=quantile(CTSR, probs=0.75, na.rm=TRUE) 
                        , DCWT=quantile(DCWT, probs=0.75, na.rm=TRUE)
                        , DCWA=quantile(RWD, probs=0.75, na.rm=TRUE)
                        , DCWT2=quantile(RWD2, probs=0.75, na.rm=TRUE)
                        , TB2=quantile(TB2, probs=0.75, na.rm=TRUE)
                        , KH=quantile(KH, probs=0.75, na.rm=TRUE)
                        , N.BAND=sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)
                        , DH=(sum(DH^4)/(sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)))^0.25), 
                 by=key(dt)]
        df0 <- as.data.frame(dt)
      }
      
      
      ### 4.6.1.3 Correct for potentially reduced last band compared to the specified resolution ####
      ### Find out for each ring the max(RADDISTR.ST)
      dt <- NULL
      dt <- data.table(cells)
      setkey(dt, "YEAR", "WOODID", "ID")
      dt <- dt[, list(RADDISTR.ST.LAST=max(RADDISTR.ST, na.rm=TRUE), RADDISTR.BAND.MAX=max(RADDISTR.BAND, na.rm=TRUE)), by=key(dt)]
      df01 <- as.data.frame(dt)
      
      ### Define for each ring the RADDISTR.ST that is max(RADDISTR.ST)-RESO; the minimum possible value is 0 (in very narrow rings and/or wide bands)
      df01$RADDISTR.ST.LAST <- df01$RADDISTR.ST-RESO[r]
      df01$RADDISTR.ST.LAST <- ifelse(df01$RADDISTR.ST.LAST<0, 0, df01$RADDISTR.ST.LAST)
      
      ### Add columns to cells dataframe with information about max. abs. distance and last band 
      cells <- full_join(cells, df01[,c("YEAR","WOODID","ID","RADDISTR.ST.LAST","RADDISTR.BAND.MAX")], by=c("ID", "YEAR", "WOODID"))
      
      ### Calculate a summary of all relevant parameters for this updated last band
      cells1 <- cells[cells$RADDISTR.ST>=cells$RADDISTR.ST.LAST,]   #subset dataframe with only data used for summarizing last updated band
      cells1 <- cells1[!is.na(cells1$YEAR),]   #for some reason, the previous line creates NA lines...
      cells$RADDISTR.ST.LAST <- NULL   #reset
      cells$RADDISTR.BAND.MAX <- NULL   #reset
      
      if (tolower(STATS[s])=="mean")   #aggregation method = "mean"
      {  
        dt <- NULL
        dt <- data.table(cells1)
        setkey(dt, "YEAR", "WOODID", "ID", "RADDISTR.BAND.MAX")
        dt <- dt[, list(RRADDISTR=mean(RRADDISTR, na.rm=TRUE)
                        , MRW=RobustMax(MRW)
                        , LA=as.numeric(mean(LA, na.rm=TRUE))
                        , TCA=mean(TCA, na.rm=TRUE)
                        , DRAD=mean(DRAD, na.rm=TRUE)
                        , DTAN=mean(DTAN, na.rm=TRUE)
                        , CWA=mean(CWA, na.rm=TRUE)
                        , CWAACC=sum(CWA, na.rm=TRUE)
                        , CWTALL=mean(CWTALL, na.rm=TRUE)
                        , CWTALL.ADJ=mean(CWTALL.ADJ, na.rm=TRUE)
                        , CWTTAN=mean(CWTTAN, na.rm=TRUE)
                        , CWTRAD=mean(CWTRAD, na.rm=TRUE) 
                        , RTSR=mean(RTSR, na.rm=TRUE)
                        , CTSR=mean(CTSR, na.rm=TRUE) 
                        , DCWT=mean(DCWT, na.rm=TRUE)
                        , DCWA=mean(RWD, na.rm=TRUE)
                        , DCWT2=mean(RWD2, na.rm=TRUE)
                        , TB2=mean(TB2, na.rm=TRUE)
                        , KH=mean(KH, na.rm=TRUE)
                        , N.BAND=sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)
                        , DH=(sum(DH^4)/(sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)))^0.25), 
                 by=key(dt)]
        df02 <- as.data.frame(dt)
      }
      
      if (tolower(STATS[s])=="median")   #aggregation method = "median"
      {  
        dt <- NULL
        dt <- data.table(cells1)
        setkey(dt, "YEAR", "WOODID", "ID", "RADDISTR.BAND.MAX")
        dt <- dt[, list(RRADDISTR=median(RRADDISTR, na.rm=TRUE)
                        , MRW=RobustMax(MRW)
                        , LA=as.numeric(median(LA, na.rm=TRUE))
                        , TCA=median(TCA, na.rm=TRUE)
                        , DRAD=median(DRAD, na.rm=TRUE)
                        , DTAN=mean(DTAN, na.rm=TRUE)
                        , CWA=median(CWA, na.rm=TRUE)
                        , CWAACC=sum(CWA, na.rm=TRUE)
                        , CWTALL=median(CWTALL, na.rm=TRUE)
                        , CWTALL.ADJ=median(CWTALL.ADJ, na.rm=TRUE)
                        , CWTTAN=median(CWTTAN, na.rm=TRUE)
                        , CWTRAD=median(CWTRAD, na.rm=TRUE) 
                        , RTSR=median(RTSR, na.rm=TRUE)
                        , CTSR=median(CTSR, na.rm=TRUE) 
                        , DCWT=median(DCWT, na.rm=TRUE)
                        , DCWA=median(RWD, na.rm=TRUE)
                        , DCWT2=median(RWD2, na.rm=TRUE)
                        , TB2=median(TB2, na.rm=TRUE)
                        , KH=median(KH, na.rm=TRUE)
                        , N.BAND=sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)
                        , DH=(sum(DH^4)/(sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)))^0.25), 
                 by=key(dt)]
        df02 <- as.data.frame(dt)
      }
      
      if (tolower(STATS[s])=="q25")   #aggregation method = "25th quantile"
      {  
        dt <- NULL
        dt <- data.table(cells1)
        setkey(dt, "YEAR", "WOODID", "ID", "RADDISTR.BAND.MAX")
        dt <- dt[, list(RRADDISTR=quantile(RRADDISTR, probs=0.25, na.rm=TRUE)
                        , MRW=RobustMax(MRW)
                        , LA=as.numeric(quantile(LA, probs=0.25, na.rm=TRUE))
                        , TCA=quantile(TCA, probs=0.25, na.rm=TRUE)
                        , DRAD=quantile(DRAD, probs=0.25, na.rm=TRUE)
                        , DTAN=quantile(DTAN, probs=0.25, na.rm=TRUE)
                        , CWA=quantile(CWA, probs=0.25, na.rm=TRUE)
                        , CWAACC=sum(CWA, na.rm=TRUE)
                        , CWTALL=quantile(CWTALL, probs=0.25, na.rm=TRUE)
                        , CWTALL.ADJ=quantile(CWTALL.ADJ, probs=0.25, na.rm=TRUE)
                        , CWTTAN=quantile(CWTTAN, probs=0.25, na.rm=TRUE)
                        , CWTRAD=quantile(CWTRAD, probs=0.25, na.rm=TRUE) 
                        , RTSR=quantile(RTSR, probs=0.25, na.rm=TRUE)
                        , CTSR=quantile(CTSR, probs=0.25, na.rm=TRUE) 
                        , DCWT=quantile(DCWT, probs=0.25, na.rm=TRUE)
                        , DCWA=quantile(RWD, probs=0.25, na.rm=TRUE)
                        , DCWT2=quantile(RWD2, probs=0.25, na.rm=TRUE)
                        , TB2=quantile(TB2, probs=0.25, na.rm=TRUE)
                        , KH=quantile(KH, probs=0.25, na.rm=TRUE)
                        , N.BAND=sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)
                        , DH=(sum(DH^4)/(sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)))^0.25), 
                 by=key(dt)]
        df02 <- as.data.frame(dt)
      }
      
      if (tolower(STATS[s])=="q75")   #aggregation method = "75th quantile"
      {  
        dt <- NULL
        dt <- data.table(cells1)
        setkey(dt, "YEAR", "WOODID", "ID", "RADDISTR.BAND.MAX")
        dt <- dt[, list(RRADDISTR=quantile(RRADDISTR, probs=0.75, na.rm=TRUE)
                        , MRW=RobustMax(MRW)
                        , LA=as.numeric(quantile(LA, probs=0.75, na.rm=TRUE))
                        , TCA=quantile(TCA, probs=0.75, na.rm=TRUE)
                        , DRAD=quantile(DRAD, probs=0.75, na.rm=TRUE)
                        , DTAN=quantile(DTAN, probs=0.75, na.rm=TRUE)
                        , CWA=quantile(CWA, probs=0.75, na.rm=TRUE)
                        , CWAACC=sum(CWA, na.rm=TRUE)
                        , CWTALL=quantile(CWTALL, probs=0.75, na.rm=TRUE)
                        , CWTALL.ADJ=quantile(CWTALL.ADJ, probs=0.75, na.rm=TRUE)
                        , CWTTAN=quantile(CWTTAN, probs=0.75, na.rm=TRUE)
                        , CWTRAD=quantile(CWTRAD, probs=0.75, na.rm=TRUE) 
                        , RTSR=quantile(RTSR, probs=0.75, na.rm=TRUE)
                        , CTSR=quantile(CTSR, probs=0.75, na.rm=TRUE) 
                        , DCWT=quantile(DCWT, probs=0.75, na.rm=TRUE)
                        , DCWA=quantile(RWD, probs=0.75, na.rm=TRUE)
                        , DCWT2=quantile(RWD2, probs=0.75, na.rm=TRUE)
                        , TB2=quantile(TB2, probs=0.75, na.rm=TRUE)
                        , KH=quantile(KH, probs=0.75, na.rm=TRUE)
                        , N.BAND=sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)
                        , DH=(sum(DH^4)/(sum(MRW, na.rm=TRUE)/mean(MRW, na.rm=TRUE)))^0.25), 
                 by=key(dt)]
        df02 <- as.data.frame(dt)
      }
      
      colnames(df02)[4] <- "RADDISTR.BAND"   #rename column
      cells1 <- NULL   #reset
      
      ### Replace the current last band by this new one
      df0 <- rbind(df02, df0)   #merge the the two dataframes with the one holding the updated last-band information on top
      
      dt <- NULL
      dt <- data.table(df0)
      dt[, N.TRUE := RobustMin(N.BAND), by = c("YEAR", "WOODID", "ID", "RADDISTR.BAND")]   #add column with minimum (original) number for each band -> the same as n, except for the last band, for which the count before updating is given
      df0 <- as.data.frame(dt)
      df0 <- df0[!duplicated(df0[,c("YEAR","WOODID","ID","RADDISTR.BAND")]),]   #eliminate last-band information from lower part (with original bands)
      df0 <- arrange(df0, YEAR,RADDISTR.BAND)   #bring data frame into original order: year-band
  
      
      #### THE FOLLOWING BLOCK COMMENTED AS IT IS THE OLD WAY OF CALCULATING EWW & LWWW 
      # ### 4.6.1.4 Calculate EWW, LWW and position of EW-LW transition ###
      # ### Create smoothed Mork's profile when using high resolution to locate EW-LW transition more confidentially
      # smooth <- ifelse(RESO[r]<=2, 13, ifelse(RESO[r]<=5, 9, ifelse(RESO[r]<=10, 5, ifelse(RESO[r]<=20, 3, 1))))   #only smooth with resolution 2, 5, 10 and 20 microns
      # 
      # ### CONSIDER TO PERFORM THIS FOR EACH RING INDIVIDUALLY!
      # if (smooth>1)
      # {
      #   margin <- c(mean(df0$RTSR.EW[1:((smooth-1)/2)], na.rm=TRUE), NA ,mean(df0$RTSR.EW[c((ln<-length(df0$RTSR.EW)):(ln-((smooth-1)/2)))], na.rm=TRUE))
      #   df0$RTSR.S <- rollapply(df0$RTSR.EW, smooth, fill=margin, mean, na.rm=TRUE)      
      # } else
      # {  
      #   df0$RTSR.S <- df0$RTSR.EW 
      # }
      
      # ### If the intra-ring profile does not contain any band with Mork > 1 after smoothing, replace the values by the original values
      # dt <- NULL
      # dt <- data.table(df0)
      # # dt <- data.table(ddply(dt, "YEAR", mutate, MAXMORK=RobustMax(RTSR.S)))   #find for each year the max. Mork
      # dt <- data.table(ddply(dt, "YEAR", mutate, MAXMORK=max(RTSR.S,na.rm=T)))   #find for each year the max. Mork
      # 
      # dt$RTSR.S <- ifelse(dt$MAXMORK<1,dt$RTSR.EW,dt$RTSR.S)   #Replace RTSR.S by RTSR if max MAXMORK < 1
      # dt$MAXMORK <- NULL
      # df0 <- as.data.frame(dt)
      
      # ### Find band with maximum Mork index excluding the very first part of the earlywood
      # dt <- NULL
      # # dt <- data.table(df0[df0$RRADDISTR>=50,])
      # dt <- data.table(df0[(df0$RRADDISTR/100*df0$MRW>=50 | df0$RRADDISTR>33),])   #only consider bands that are at least 50 micrometer away from ring border; note that RRADDISTR is based on individual cells and only approximate
      # setkey(dt, "YEAR")
      # #     dt <- dt[, list(BANDMXMORK=RADDISTR.BAND[which(RTSR.S==RobustMax(RTSR.S))]), by=key(dt)]   #for some reason this sometimes skipped some years... 
      # suppressWarnings(dt <- dt[, list(BANDMXMORK=RADDISTR.BAND[which(RTSR.S==max(RTSR.S,na.rm=T))]), by=key(dt)])
      # ewlw <- as.data.frame(dt)
      # 
      # # ewlw <- distinct(ewlw, YEAR)   #to make sure only the first of equal max. values within a year is taken!
      # ewlw <- ewlw[!duplicated(ewlw[,c("YEAR")]),]   #replace distinct like this for compatibility reasons
      # 
      # df0 <- left_join(df0, ewlw, by="YEAR")  
      # ewlw <- NULL
      
      # ### Loop from position of maximum smoothed Mork to beginning of ring and exit when Mork threshold is crossed
      # df0.sub <- NULL
      # df0.sub <- df0[!is.na(df0$RADDISTR.BAND) & df0$RADDISTR.BAND<=df0$BANDMXMORK,]   #only get data before maximum Mork
      # 
      # dt <- data.table(df0.sub) 
      # # dt <- data.table(ddply(dt, "YEAR", mutate, MAXMORK=RobustMax(RTSR.S)))   #add for each year the max. Mork
      # dt <- data.table(ddply(dt, "YEAR", mutate, MAXMORK=max(RTSR.S,na.rm=T)))   #add for each year the max. Mork
      # dt <- dt[dt$MAXMORK>=1,]   #extract years with max. Mork >=1; this is required to have a LW-EW transition    
      # dt <- setorder(dt, YEAR, -RADDISTR.BAND)
      # dt <- dt[RTSR.S<1, .SD[1], by="YEAR"]   #extract band with LW-EW transition  
      # dt$MAXMORK <- NULL    
      
      # ### Calculate some statistics per year: mean Mork, position of max. and min. 
      # dt1 <- data.table(df0.sub, key="YEAR", "RADDISTR.BAND")
      # dt1 <- dt1[, list(MRTSR=mean(RTSR, na.rm=TRUE)
      #                   , MINBAND=RobustMin(RADDISTR.BAND)
      #                   , MAXRPOS=RobustMax(RRADDISTR)
      #                   , MINRPOS=RobustMin(RRADDISTR))
      #             , by=key(dt1)]
      # dt1 <- dt1[!is.na(dt1$YEAR),]
      
      # ### Add information about band of max. Mork index to dt1; relevant for large bands  
      # dt2 <- data.table(df0.sub, key="YEAR")
      # dt2 <- dt2[, list(BANDMXMORK=RADDISTR.BAND[which(RTSR.S==RobustMax(RTSR.S))], MRW=RobustMax(MRW)), by=key(dt2)] 
      # dt2 <- dt2[!is.na(dt2$YEAR),]
      # dt1 <- merge(dt1, dt2, "YEAR", all=TRUE)
      # dt2 <- NULL
      # 
      # dt$BANDMXMORK <- NULL   #remove column to avoid doubling during merge in next block 
      # dt$MRW <- NULL   #remove column to avoid doubling during merge in next block    
      
      # ### Add statistics to main data table and replace any NA (if no EW-LW transition found!) by sensible values
      # ewlw <- as.data.frame(merge(dt, dt1, "YEAR", all=TRUE))
      # # ewlw$V2 <- NULL    #remove unused column created during data.table manipulation
      # # ewlw$RADDISTR.BAND <- ifelse(!is.na(ewlw$RADDISTR.BAND),ewlw$RADDISTR.BAND,ifelse(ewlw$MRTSR<1,ewlw$MINBAND,ewlw$BANDMXMORK))   #if no transition found, assign all to EW or LW depending on mean Mork
      # # ewlw$RRADDISTR <- ifelse(!is.na(ewlw$RRADDISTR),ewlw$RRADDISTR,ifelse(ewlw$MRTSR<1,ewlw$MINRPOS,ewlw$MAXRPOS))   #if no transition found, assign all to EW or LW depending on mean Mork
      # ewlw$RADDISTR.BAND <- ifelse(!is.na(ewlw$RADDISTR.BAND),ewlw$RADDISTR.BAND,ifelse(ewlw$MRTSR<1,ewlw$BANDMXMORK,ewlw$MINBAND))   #if no transition found, assign all to EW or LW depending on mean Mork
      # ewlw$RRADDISTR <- ifelse(!is.na(ewlw$RRADDISTR),ewlw$RRADDISTR,ifelse(ewlw$MRTSR<1,ewlw$MAXRPOS,ewlw$MINRPOS))   #if no transition found, assign all to EW or LW depending on mean Mork
      
      # ### Calculate EWW and LWW and the position of the EW-LW transition
      # ewlw <- ewlw[,c("YEAR","RADDISTR.BAND","RRADDISTR","MRW","BANDMXMORK")]
      # ewlw$EWW <- round(ewlw$MRW * (ewlw$RRADDISTR/100))
      # ewlw$LWW <- ewlw$MRW - ewlw$EWW
      # ewlw <- ewlw[,c("YEAR","EWW","LWW","RADDISTR.BAND","RRADDISTR")]
      # colnames(ewlw)[4] <- "EWLW.BAND"
      # colnames(ewlw)[5] <- "EWLW.RPOS"
      
      # ### Merge newly created data with overall dataframe
      # df0 <- left_join(df0, ewlw, by="YEAR")
      # # df0$EWLW.ID <- as.factor(ifelse(df0$RADDISTR.BAND<=df0$EWLW.BAND,"ew","lw"))
      # # df0$EWLW.ID <- as.factor(ifelse(df0$MRW >= RESO[r],   # 
      # #                                 ifelse(df0$RADDISTR.BAND<=df0$EWLW.BAND,"ew","lw"),
      # #                                 ifelse(df0$EWW>=df0$LWW,"ew","lw")))    
      # # df0$EWLW.ID <- as.factor(ifelse(df0$MRW >= RESO[r],   # 
      # #                                 ifelse(df0$EWLW.BAND==0, ifelse(df0$RTSR.S<1,"ew","lw"),
      # #                                        ifelse(df0$RADDISTR.BAND<=df0$EWLW.BAND, "ew","lw")),
      # #                                 ifelse(df0$RTSR.S<1,"ew","lw"))) 
      # df0$EWLW.ID <- as.factor(ifelse(df0$MRW >= RESO[r],   # 
      #                                 ifelse(df0$EWLW.BAND==0, ifelse(df0$RTSR.S<1,"ew","lw"),
      #                                        ifelse(df0$RADDISTR.BAND<=df0$EWLW.BAND, "ew",
      #                                               ifelse(df0$RTSR.S<1,"ew","lw"))),
      #                                 ifelse(df0$RTSR.S<1,"ew","lw"))) 
      # 
      # df0 <- df0[!is.na(df0$YEAR),]   #for some reasons there can be years with NA
      # # d <- df0[,c(1:6,22:28)]
      # # fix(d)
      #### THE PREVIOUS BLOCK COMMENTED AS IT IS THE OLD WAY OF CALCULATING EWW & LWWW 
      
      
      ### 4.6.1.4 Calculate EWW, LWW and position of EW-LW transition ####
      ### Only calculate EWW and LWW on first run per sample as they are independent from aggregation method (STATS) and intra-ring resolution (RESO)
      if (s==1 & r==1)  
      {
        
        ### 4.6.1.4.1 Create data frame with median Mork's index (RTSR) in 10-µm steps   
        dt <- NULL
        dt <- data.table(cells[,c("YEAR", "MRW", "RRADDISTR", "RADDISTR.ST", "RTSR")])
        dt$YEAR <- as.factor(dt$YEAR)
        dt$RADDISTR.BAND <- 10 * round((dt$RADDISTR.ST- 10 / 2.01) / 10 , 0)   #add column where each cell is assigned to 10-µm band; insert row if one band is missing
        
        setkey(dt, "YEAR", "RADDISTR.BAND")
        dt <- dt[, list(RRADDISTR=median(RRADDISTR, na.rm=TRUE)
                        , MRW=median(MRW, na.rm=TRUE)
                        , RTSR=median(RTSR, na.rm=TRUE))
                    , by=key(dt)]
        df01 <- as.data.frame(dt)
        
        ### 4.6.1.4.2 Smooth Mork's profile for each ring separately; fill marginal values with averages (n = (smooth-1)/2)    
        smooth <- 5   #define smoothing kernel size
  
        if (smooth>1)
        {      
          for (y in unique(df01$YEAR))
          {
            df02 <- df01[df01$YEAR==y,]
            margin <- c(mean(df02$RTSR[1:((smooth-1)/2)], na.rm=TRUE), NA ,mean(df02$RTSR[c((ln<-length(df02$RTSR)):(ln-((smooth-1)/2)+1))], na.rm=TRUE))   #define the values at lower margin (first value), if NA (second value), and at upper margin (third value)
            df02$RTSR.S <- rollapply(df02$RTSR, smooth, fill=margin, median, na.rm=TRUE)   #apply rolling mean      
            df01$RTSR.S[df01$YEAR==y] <- df02$RTSR.S
          }  
        } else
        {  
          df01$RTSR.S <- df01$RTSR 
        }
        
        ### 4.6.1.4.3 If the intra-ring profile does not contain any band with Mork > MORK.EWLW after smoothing, replace the values by the original values
        dt <- NULL
        dt <- data.table(df01)
        dt <- data.table(ddply(dt, "YEAR", mutate, MAXMORK=max(RTSR.S,na.rm=T)))   #find for each year the max. Mork
        
        dt$RTSR.S <- ifelse(dt$MAXMORK<MORK.EWLW, dt$RTSR, dt$RTSR.S)   #Replace RTSR.S by RTSR if max MAXMORK < 1
        dt$MAXMORK <- NULL
        df01 <- as.data.frame(dt)
        
        ### 4.6.1.4.4 Find band with maximum Mork index excluding the very first part of the earlywood
        dt <- NULL
        dt <- data.table(df01[(df01$RRADDISTR/100*df01$MRW>=50 | df01$RRADDISTR>33),])   #only consider bands that are at least 50 micrometer away from ring border OR later than relative intra-ring position 33%
        # setkey(dt, "YEAR")
        # suppressWarnings(dt <- dt[, list(BANDMXMORK=RADDISTR.BAND[which(RTSR.S==max(RTSR.S,na.rm=T))]), by=key(dt)])
        dt <- dt[, list(BANDMXMORK=RADDISTR.BAND[which(RTSR.S==max(RTSR.S,na.rm=T))]), by="YEAR"]
        ewlw <- as.data.frame(dt)
        
        ewlw <- ewlw[!duplicated(ewlw$YEAR),]   #to make sure only the first of equal max. values within a year is taken!
        
        df01 <- left_join(df01, ewlw, by="YEAR")  
        ewlw <- NULL      
  
        ### 4.6.1.4.5 Loop from position of maximum smoothed Mork to beginning of ring and exit when Mork threshold is crossed
        df01.sub <- NULL
        df01.sub <- df01[!is.na(df01$RADDISTR.BAND) & df01$RADDISTR.BAND<=df01$BANDMXMORK,]   #only get data before maximum Mork
        
        dt <- data.table(df01.sub) 
        dt <- data.table(ddply(dt, "YEAR", mutate, MAXMORK=max(RTSR.S,na.rm=T), MINMORK=min(RTSR.S,na.rm=T)))   #add for each year the max. & min. Mork
        dt <- dt[dt$MAXMORK>=MORK.EWLW & dt$MINMORK<MORK.EWLW,]   #extract years with max. Mork >= MORK.EWLW & min. Mork < MORK.EWLW; this is required to have a LW-EW transition    
        dt$MAXMORK <- NULL   #remove unused column          
        dt$MINMORK <- NULL   #remove unused column          
        dt <- setorder(dt, YEAR, -RADDISTR.BAND)
        dt <- dt[RTSR.S<MORK.EWLW, .SD[1], by=YEAR]   #extract band with LW-EW transition  
        
        ### 4.6.1.4.6 Calculate some statistics per year: mean Mork, position of max. and min. 
        # dt1 <- data.table(df01.sub, key="YEAR", "RADDISTR.BAND")   #=> BETTER USE DF01, NOT SUBSET?
        dt1 <- data.table(df01, key="YEAR")   #=> BETTER USE DF01, NOT SUBSET?

        dt1 <- dt1[, list(MRTSR=mean(RTSR, na.rm=TRUE))
                          # , MINBAND=RobustMin(RADDISTR.BAND)
                          # , MAXRPOS=RobustMax(RRADDISTR)
                          # , MINRPOS=RobustMin(RRADDISTR))
                   , by=key(dt1)]
        dt1 <- dt1[!is.na(dt1$YEAR),]      
        
        ### 4.6.1.4.7 Add information about band of max. Mork index to dt1; relevant if large intra-ring resolution was opted for 
        # dt2 <- data.table(df01.sub)
        dt2 <- data.table(df01.sub[(df01.sub$RRADDISTR/100*df01.sub$MRW>=50 | df01.sub$RRADDISTR>33),])   #only consider bands that are at least 50 micrometer away from ring border OR later than relative intra-ring position 33%
        dt2 <- dt2[, list(BANDMXMORK=RADDISTR.BAND[which(RTSR.S==RobustMax(RTSR.S))], MRW=RobustMax(MRW)), by="YEAR"] 
        dt2 <- dt2[!is.na(dt2$YEAR),]
        dt1 <- merge(dt1, dt2, "YEAR", all=TRUE)
        dt2 <- NULL
        
        dt$BANDMXMORK <- NULL   #remove column to avoid doubling during merge in next block 
        dt$MRW <- NULL   #remove column to avoid doubling during merge in next block          
        
        ### 4.6.1.4.8 Add statistics to main ewlw data table and replace any NA (if no EW-LW transition found!) by sensible values
        ewlw <- as.data.frame(merge(dt, dt1, "YEAR", all=TRUE))
        # ewlw$RADDISTR.BAND <- ifelse(!is.na(ewlw$RADDISTR.BAND),ewlw$RADDISTR.BAND,ifelse(ewlw$MRTSR<MORK.EWLW,ewlw$BANDMXMORK,ewlw$MINBAND))   #if no transition found, assign all to EW or LW depending on mean Mork
        ewlw$RADDISTR.BAND <- ifelse(!is.na(ewlw$RADDISTR.BAND),ewlw$RADDISTR.BAND,ifelse(ewlw$MRTSR<MORK.EWLW,ewlw$BANDMXMORK,0))   #if no transition found, assign all to EW or LW depending on mean Mork
        # ewlw$RRADDISTR <- ifelse(!is.na(ewlw$RRADDISTR),ewlw$RRADDISTR,ifelse(ewlw$MRTSR<MORK.EWLW,ewlw$MAXRPOS,ewlw$MINRPOS))   #if no transition found, assign all to EW or LW depending on mean Mork
        ewlw$RRADDISTR <- ifelse(!is.na(ewlw$RRADDISTR),ewlw$RRADDISTR,ifelse(ewlw$MRTSR<MORK.EWLW,100,0))   #if no transition found, assign all to EW or LW depending on mean Mork
        
        ### 4.6.1.4.9 Calculate EWW and LWW and the position of the EW-LW transition
        ewlw <- ewlw[,c("YEAR","RADDISTR.BAND","RRADDISTR","MRW","BANDMXMORK")]
        ewlw$EWW <- round(ewlw$MRW * (ewlw$RRADDISTR/100))
        ewlw$LWW <- ewlw$MRW - ewlw$EWW
        ewlw <- ewlw[,c("YEAR","EWW","LWW","RADDISTR.BAND","RRADDISTR")]
        colnames(ewlw)[4] <- "EWLW.BAND"
        colnames(ewlw)[5] <- "EWLW.RPOS"
      }   #end of if (s==1 & r==1)
      
      ### 4.6.1.4.10 Merge newly created data with overall dataframe
      df0$YEAR <- as.factor(df0$YEAR)
      df0 <- left_join(df0, ewlw, by="YEAR")
      
      #### 4.6.1.4.11 Express last EW band within each ring in target resolution
      # df0$EWLW.BAND <- round(df0$EWLW.BAND / RESO[r]) * RESO[r]
      df0$EWLW.BAND <- round((df0$EWLW.BAND-(RESO[r]/2-1)) / RESO[r]) * RESO[r]  
      
      ### 4.6.1.4.12 Add column indicating for each band, whether it is earlywood or latewood 
      df0$EWLW.ID <- as.factor(ifelse(df0$MRW >= RESO[r],   # 
                                      ifelse(df0$EWLW.BAND==0, ifelse(df0$RTSR<MORK.EWLW,"ew","lw"),
                                             ifelse(df0$RADDISTR.BAND<=df0$EWLW.BAND, "ew",
                                             ifelse(df0$RTSR<MORK.EWLW,"ew","lw"))),
                                             # ifelse(df0$RADDISTR.BAND<=df0$EWLW.BAND, "ew", "lw")),
                                      ifelse(df0$RTSR<MORK.EWLW,"ew","lw"))) 
      
      df0 <- df0[!is.na(df0$YEAR),]   #for some reasons there can be years with NA
  
      ### 4.6.1.4.13 If a ring has no latewood band based on the above (likely due to low intra-ring resolution), 
      ###            check whether last band(s) have a Mork's index > MORK.EWLW
      if (length(unique(df0$YEAR)[!unique(df0$YEAR) %in% unique(df0$YEAR[df0$EWLW.ID=="lw"])]) > 0)   #check for any ring without minimum 1 latewood band
      {  
        for (y in unique(df0$YEAR)[!unique(df0$YEAR) %in% unique(df0$YEAR[df0$EWLW.ID=="lw"])])   #loop through all found rings
        {
          # if (!is.na(unique(df0$EWLW.BAND[df0$YEAR==y][1])) & !is.na(df0$RADDISTR.BAND[df0$YEAR==y]))   #exclude years with Mork's index = NA (e.g., when ring was not finished/completed at the end of the image -> stitching mistake!)
          if (length(is.na(df0$EWLW.BAND[df0$YEAR==y]))==0 & length(is.na(df0$RADDISTR.BAND[df0$YEAR==y]))==0)   #exclude years with Mork's index = NA (e.g., when ring was not finished/completed at the end of the image -> stitching mistake!)
          {  
            for (j in seq(unique(df0$EWLW.BAND[df0$YEAR==y]),max(df0$RADDISTR.BAND[df0$YEAR==y], na.rm=TRUE),by=RESO[r]))   #loop through all bands starting from latest earlywood band towards ring border
            {
              if (!is.na(df0$RTSR[df0$YEAR==y & df0$RADDISTR.BAND==j]))   #exclude bands with Mork's index = NA
              {  
                if (df0$RTSR[df0$YEAR==y & df0$RADDISTR.BAND==j] >= MORK.EWLW)   #is the mean Mork's index for the target band >=  MORK.EWLW?
                {
                  df0$EWLW.ID[df0$YEAR==y & df0$RADDISTR.BAND==j] <- "lw"   #... if yes, assign it to latewood
                }  
              }  
            }
          }  
        }  
      }  
      
      # ### Optional: evaluate how many years have no latewood and plot it
      # d <- df0[df0$YEAR %in% unique(df0$YEAR)[!unique(df0$YEAR) %in% unique(df0$YEAR[df0$EWLW.ID=="lw"])],]
      # d <- d[,c("YEAR","WOODID","ID","RADDISTR.BAND","MRW","RTSR","EWW","LWW","EWLW.BAND","EWLW.ID")]
      # fix(d)
      
      ### CONSIDER TO ADD HERE CODE TO FORCE THE LAST BAND OF EACH RING TO BE LATEWOOD
      
      
      ### 4.6.1.5 Add missing bands (bands containing no cells) ####
      ### Create data.table with number of bands for each year
      dt <- NULL
      dt <- data.table(df0)
      setkey(dt, "YEAR")
      dt <- dt[, list(MRW=RobustMax(MRW)), by=key(dt)]
      
      dt$MOD <- (dt$MRW-(RESO[r]/2)) %% RESO[r]   #get number of bands to be used considering the ring width of target ring
      dt$NUMBAND <- ifelse(dt$MRW-(RESO[r]/2)>0, (dt$MRW-(RESO[r]/2)-dt$MOD), 0)
      dt <- dt[!is.na(dt$YEAR),]   #in some cases there is a leading year with NA
      
      ### 4.6.1.6 Create data frame with completed series of intra-annual bands for each calendar year ####
      allyears <- lapply(1:length(unique(dt$YEAR)),function(y){seq(0, dt$NUMBAND[y], RESO[r])})
      bandlist <- data.table(BANDS=unlist(allyears))
      bandlist$YEAR <- NA
      bandlist$YEAR[which(bandlist$BANDS==0)] <- as.character(dt$YEAR)   #write calendar year at first intra-annual bands
      bandlist$YEAR <- na.locf(bandlist$YEAR)   #last observation carried forward
      bandlist$YEAR <- as.factor(bandlist$YEAR)
      
      df1 <- as.data.frame(merge(bandlist,dt,"YEAR"))
      colnames(df1)[2] <- "RADDISTR.BAND"
      df1[3:5] <- list(NULL)
      
      # ### Remove first bands that are not included in the first ring (partial ring) #THIS IS PROBABLY NO MORE NEEDED!
      # # min.year <- min(df0$YEAR, na.rm=TRUE)
      # first.row <- df0$RADDISTR.BAND[1] / RESO[r] + 1   
      # df1 <- df1[c(first.row:nrow(df1)),]
  
      ### 4.6.1.7 Merge with master data frame to fill gaps of bands #### 
      df2 <- full_join(df0,df1, by=c("YEAR", "RADDISTR.BAND"))
      df2 <- arrange(df2, YEAR, RADDISTR.BAND)   #bring data frame into original order: year-band
  
      ### Some dataframe housekeeping
      setDT(df2)[, MRW:= MRW[!is.na(MRW)][1L], by=YEAR]   #if NA-line, copy MRW from first MRW-value of respective YEAR
      setDT(df2)[, WOODID:=WOODID[!is.na(WOODID)][1L], by=YEAR]   #if NA-line, copy WOODID from first WOODID-value of respective YEAR
      setDT(df2)[, ID:=ID[!is.na(ID)][1L], by=YEAR]   #if NA-line, copy ID from first ID-value of respective YEAR
      df2$N.BAND <- ifelse(is.na(df2$N.BAND), 0, df2$N.BAND)   #replace NA by 0 for number of cells in target band 
      df2$N.TRUE <- ifelse(is.na(df2$N.TRUE), 0, df2$N.TRUE)   #replace NA by 0 for number of cells in target band 
      df2 <- df2[!is.na(df2$RADDISTR.BAND),]  #remove rows with NA at RADDISTR.BAND (possible if no overlapping rings and ring at lower image edge not finished)
      
      ### Create first version of intra-annual profiles dataframe 
      iap <- as.data.frame(df2)
      
      iap$EWLW.ID <- as.factor(iap$EWLW.ID)
      iap$EWLW.ID[iap$EWLW.ID==""] <- NA
      iap$YEAR  <- as.numeric(as.character(iap$YEAR))
      colnames(iap)[3] <- "IMAGE"
  
  
      ### 4.6.1.8 Add column with the relative position of each band in each year (for plotting keeping constant ring widths) ####
      iap$YR.RRADDISTR.BAND <- iap$RADDISTR.BAND / iap$MRW * 100   #express each band as a percentage of entire ring width
      iap$YR.RRADDISTR.BAND <- ifelse(iap$YR.RRADDISTR.BAND >= 100, 99.99, iap$YR.RRADDISTR.BAND)   #replace any 100 positions by 99.99 (should not be required any more)
      iap$YR.RRADDISTR.BAND <- round(100*iap$YR.RRADDISTR.BAND, 0)   
      iap$YR.RRADDISTR.BAND <- sprintf("%04d", iap$YR.RRADDISTR.BAND)   #convert to text
      iap$YR.RRADDISTR.BAND <- as.numeric(paste(as.character(iap$YEAR), as.character(iap$YR.RRADDISTR.BAND), sep="."))   #merge year and rel. position information
      
      
      ### 4.6.1.9 Add column with continuous absolute distance based on standardized bands (for plotting considering different ring widths) ####
      ### If using the "special" last band, make sure it is done correctly without creating artifacts!
      df5 <- NULL
      offset <- RESO[r]
      
      for (y in unique(iap$YEAR))   
      {
        df4 <- iap$RADDISTR.BAND[iap$YEAR==y]   #extract data for target year 
        df4[length(df4)] <- iap$MRW[iap$YEAR==y][1] - RESO[r]
        
        # combine data frames from all years
        if (y==iap$YEAR[1])
        {
          df5 <- df4
        } else
        {
          df5 <- c(df5, max(df5)+offset+df4)  
        }
      }
      
      df5 <- df5 + RESO[r]/2   #center in band
      iap$RADDIST.CONT <- df5   
      
      
      ### 4.6.1.20 Add label for x-axis (only at beginning of each calendar year) ####
      iap$X.LABEL <- ifelse(iap$RADDISTR.BAND==0, iap$YEAR, NA)
      
      
      ### 4.6.1.21 Insert line with NA where there are missing rings ####
      repeat
      {  
        missingrow <- FALSE
        for (h in 2:(nrow(iap)-1))
        {
          if (iap$YEAR[h] > iap$YEAR[h-1] + 1)
          {
            newrow <- as.data.frame(t(rep(as.numeric(NA),ncol(iap))))
            colnames(newrow) <- colnames(iap)
            newrow$YEAR <- as.numeric(iap$YEAR[h-1] + 1)
            newrow$WOODID <-  as.factor(as.character(iap$WOODID[h-1]))
            newrow$IMAGE <-  as.factor(as.character(iap$IMAGE[h-1]))        
            
            iap <- rbind(iap[1:(h-1),], newrow, iap[-(1:(h-1)),])
            newrow <- NULL
            missingrow <- TRUE
            break
          } 
        }
        if (missingrow == FALSE)
        {
          break
        }  
      }
  
      row.names(iap) <- NULL   
      iap$RADDISTR.BAND <- as.numeric(iap$RADDISTR.BAND)
  
  
      ### 4.6.1.22 Remove years with uncertain dating based on file _YearsToExclude.txt ####
      if (woodid %in% exclude$WOODID)
      {   
        iap[iap$YEAR %in% exclude$YEAR[exclude$WOODID==woodid]==TRUE, c(4:(ncol(iap)-3))] <- NA        
      }
  
      write.table(iap, file=paste(woodid, "_IntraannualProfiles_", STATS[s], "_", RESO[r], "mu.txt", sep=""), row.names = FALSE)
  
  
      ### 4.6.1.23 Create new dataframe excluding rows that have NA for the continuous radial distance ####
      # iap2 <- iap[!is.na(iap$RADDIST.CONT),]
      iap2 <- iap
      
      
      ### 4.6.1.24 Smooth the data by a rolling mean ####
      if (SMOOTHER > 1)
      {  
        margin <- c(mean(iap2$LA[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$LA[c((ln<-length(iap2$LA)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$LA <- rollapply(iap2$LA, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$DRAD[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$DRAD[c((ln<-length(iap2$DRAD)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$DRAD <- rollapply(iap2$DRAD, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$DTAN[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$DTAN[c((ln<-length(iap2$DTAN)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$DTAN <- rollapply(iap2$DTAN, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$TCA[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$TCA[c((ln<-length(iap2$TCA)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$TCA <- rollapply(iap2$TCA, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$CWA[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$CWA[c((ln<-length(iap2$CWA)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$CWA <- rollapply(iap2$CWA, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$CWAACC[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$CWAACC[c((ln<-length(iap2$CWAACC)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$CWAACC <- rollapply(iap2$CWAACC, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$CWTALL[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$CWTALL[c((ln<-length(iap2$CWTALL)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$CWTALL <- rollapply(iap2$CWTALL, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$CWTALL.ADJ[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$CWTALL.ADJ[c((ln<-length(iap2$CWTALL.ADJ)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$CWTALL.ADJ <- rollapply(iap2$CWTALL.ADJ, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$CWTTAN[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$CWTTAN[c((ln<-length(iap2$CWTTAN)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$CWTTAN <- rollapply(iap2$CWTTAN, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$CWTRAD[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$CWTRAD[c((ln<-length(iap2$CWTRAD)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$CWTRAD <- rollapply(iap2$CWTRAD, SMOOTHER, fill=margin, mean, na.rm=TRUE)      
        iap2$RTSR <- rollapply(iap2$RTSR, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$CTSR[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$CTSR[c((ln<-length(iap2$CTSR)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$CTSR <- rollapply(iap2$CTSR, SMOOTHER, fill=margin, mean, na.rm=TRUE)            
        margin <- c(mean(iap2$DCWT2[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$DCWT2[c((ln<-length(iap2$DCWT2)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$DCWT2 <- rollapply(iap2$DCWT2, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$DCWT[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$DCWT[c((ln<-length(iap2$DCWT)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$DCWT <- rollapply(iap2$DCWT, SMOOTHER, fill=margin, mean, na.rm=TRUE)
        margin <- c(mean(iap2$DCWA[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$DCWA[c((ln<-length(iap2$DCWA)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$DCWA <- rollapply(iap2$DCWA, SMOOTHER, fill=margin, mean, na.rm=TRUE)     
        margin <- c(mean(iap2$TB2[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$TB2[c((ln<-length(iap2$TB2)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$TB2 <- rollapply(iap2$TB2, SMOOTHER, fill=margin, mean, na.rm=TRUE)      
        margin <- c(mean(iap2$KH[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$KH[c((ln<-length(iap2$KH)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$KH <- rollapply(iap2$KH, SMOOTHER, fill=margin, mean, na.rm=TRUE)    
        margin <- c(mean(iap2$DH[1:((SMOOTHER-1)/2)], na.rm=TRUE), NA ,mean(iap2$DH[c((ln<-length(iap2$DH)):(ln-((SMOOTHER-1)/2)))], na.rm=TRUE))
        iap2$DH <- rollapply(iap2$DH, SMOOTHER, fill=margin, mean, na.rm=TRUE)    
      }
      
  
      ### 4.6.1.25 Define some general options for plotting intra-annual profiles ####
      ### 4.6.1.25.1 Image dimensions ####
      im.width <- round(0.3*sum(rings$MRW))   #scale image width by overall length in absolute units: 0.3 pixels per micrometer width
      im.height <- 1200
      
      ### 4.6.1.25.2 If Quality check mode is activated, perform several quality control operations ####
      if (PLOT.IMG.ID==TRUE)   
      {  
        ### 4.6.1.25.2.1 Define position of grey-white background rectangles indicating changes in data source from one image to the next, and get image caption of each data block ####
        iap3 <- iap2[,c("YEAR","IMAGE","RADDIST.CONT","X.LABEL")]   #subset intra-annual dataframe
        iap3 <- iap3[!is.na(iap3$X.LABEL),]   #extract rows at the beginning of each ring
        iap3$BLOCK <- 1   #initialize grouping column
        for (n in c(2:nrow(iap3)))   #give same index to rows getting data from the same image; if data source flips between two images where they overlap, treat each segment as independent group
        {
          if (iap3$IMAGE[n]!=iap3$IMAGE[n-1]) 
          {
            iap3$BLOCK[n] <- iap3$BLOCK[n-1] + 1
          }else
          {
            iap3$BLOCK[n] <- iap3$BLOCK[n-1]
          }  
        }  
        iap3$BLOCK <- as.factor(iap3$BLOCK)   #convert grouping index to factor
        iap3 <- iap3[!duplicated(iap3[c("BLOCK")]),]   #extract the first row of each group
        if (nrow(iap3)>1)   #if there is more than one group, define corner coordinates of grouping background rectangles
        {  
          topleft <- iap3$RADDIST.CONT[seq(2,nrow(iap3),2)]   #coordinate of every other group at top left corner
          bottomright <- iap3$RADDIST.CONT[seq(3,nrow(iap3)+1,2)]   #coordinate of every other group at bottom right corner
          if (is.na(bottomright[length(bottomright)]))   #if number of groups is even, the last coordinate at the bottom right corner is NA, but needs to get the outmost coordinate in the series
          {
            bottomright[length(bottomright)] <- max(iap2$RADDIST.CONT,na.rm=TRUE)
          }  
        }else   #if there is only one group, set coordinates to 0
        {
          topleft <- 0 
          bottomright <- 0
        } 
        
        ### 4.6.1.25.2.2 Check for too few cells in the bands ####
        iap4 <- NULL
        if (length(which(iap2$N.BAND < NCELL.MIN & !is.na(iap2$LA))) > 0)   #only extract bands if there are any with too few cells that are not NA
        {  
          iap4 <- iap2[which(iap2$N.BAND < NCELL.MIN & !is.na(iap2$LA)),c("RADDIST.CONT","N.BAND","LA")]
        }
        
        ### 4.6.1.25.2.3 Calculate nth (e.g., 90th) percentile of lumen area per image or contiguous image segment, respectively ####
        iap5 <- NULL
        dt <- NULL
        dt <- data.table(iap2)
        dt$BLOCK <- as.numeric(match(dt$IMAGE, levels(dt$IMAGE)))   #assign numeric index of unique source image to grouping column
        setDT(dt)[ , BLOCK2:=cumsum(c(1L, (BLOCK!=shift(BLOCK,1L))[-1L]))]   #create second grouping column that holds the image; if data source flips between two images where they overlap, treat each segment as independent group
        dt$BLOCK2 <- as.factor(dt$BLOCK2)
        setkey(dt, "BLOCK2")
        dt <- dt[, list(IMAGE=unique(IMAGE)   #aggregate per image segment and calculate lumen area quantile
                        , LA90=quantile(LA, LA.REF, na.rm=TRUE)
                        , X.START=min(RADDIST.CONT, na.rm=TRUE)
                        , X.END=max(RADDIST.CONT, na.rm=TRUE)+RESO[r]),
                 by=key(dt)]
        iap5 <- as.data.frame(setorder(dt, -BLOCK2))        
        
        ### 4.6.1.25.2.4 Find first value after a (block) of NAs to identify gaps in lumen area profile ####
        iap7 <- NULL
        if (sum(is.na(iap2$LA)) > 0)   #check whether there is any NA
        {
          iap7 <- iap2[which(is.na(iap2$LA)), c("RADDIST.CONT","LA")]   #extract rows with NA for lumen area
          iap7$ID <- as.numeric(row.names(iap7))
          if (nrow(iap7) > 1)   #check whether there is more than 1 NA
          {  
            for (n in 1:c(nrow(iap7)-1))
            {
              if (iap7$ID[n]+1 == iap7$ID[n+1])   #assign -999 if NA-row is following another directly adjacent NA-row
              {
                iap7$ID[n] <- -999
              }
            }  
          } 
          iap7 <- iap7[which(iap7$ID>0), "ID"] + 1   #extract NA-rows that are not flagged -999
          iap7 <- iap2$RADDIST.CONT[iap7]   #extract radial distance of first bands after a gap
          iap7 <- iap7[!is.na(iap7)]
        }
        
        ### 4.6.1.25.2.5 Check for potential duplicated rings due to wrong visual cross-dating ####
        iap8 <- NULL
        k <- 0
        for (n in c(2:length(unique(iap2$YEAR))))   #loop through all years
        {
          if (sum(!is.na(iap2$LA[iap2$YEAR==unique(iap2$YEAR)[n]]),na.rm=TRUE) > 0 & sum(!is.na(iap2$LA[iap2$YEAR==unique(iap2$YEAR)[n-1]]),na.rm=TRUE) > 0)   #exclude years with only NA
          {
            if (abs(1-(unique(iap2$MRW[iap2$YEAR==unique(iap2$YEAR)[n]]) / unique(iap2$MRW[iap2$YEAR==unique(iap2$YEAR)[n-1]]))) < DUP.RW)   #only evaluate if ring width of neighboring rings differ by <n% (e.g., 5%)
            {  
              if (ks.test(iap2$LA[iap2$YEAR==unique(iap2$YEAR)[n]], iap2$LA[iap2$YEAR==unique(iap2$YEAR)[n-1]])$p.value > DUP.KS)   #perform Kolmogorov-Smirnoff test and compare p-value against threshold
              {
                k <- k+1
                iap8[k] <- unique(iap2$YEAR)[n] - 1   #add calendar year to list
                k <- k+1
                iap8[k] <- unique(iap2$YEAR)[n]   #add calendar year to list
              }
                
            }        
          }
        }  
        # iap8 <- unique(iap2$YEAR)[c(which(iap8 > DUP.KS)-1, which(iap8 > DUP.KS))]   #assume potential duplicate ring if similar ring width and similar shape (p-vaue of ks > n (e.g., 0.95))
        
        ### 4.6.1.25.2.6 Calculate smoothed profile serving as a reference to identify (large) outliers ####
        iap9 <- iap2[,c("RADDIST.CONT","LA")]
        margin <- c(mean(iap9$LA[1:((RM.KERNEL-1)/2)], na.rm=TRUE), NA ,mean(iap9$LA[c((ln<-length(iap9$LA)):(ln-((RM.KERNEL-1)/2)))], na.rm=TRUE))
        iap9$LA <- rollapply(iap9$LA, RM.KERNEL, fill=margin, mean, na.rm=TRUE)        
      }
      
      
      ### 4.6.1.26 Plot and save intra-annual profile of LA (lumen area); add additional QC features if opted for ####
      png(file=paste(woodid, "_Intraprofile_LA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$LA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Lumen Area (LA) (µm2)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(25,1500))
           cex.axis=3, cex.lab=3, ylim=c(25,max(iap2$LA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      
      ### add QC features
      if (PLOT.IMG.ID==TRUE)
      {  
        ### Grey-white background and image captions
        lim <- par("usr")
        if (length(bottomright)>1)
        {  
          rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray95", col = "gray95")
        }
        text(x=iap3$RADDIST.CONT,y=max(iap2$LA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T),iap3$IMAGE,adj=c(0,0),cex=3)
        
        ### Highlight image-segment level lumen area quantiles if deviating strongly from neighboring segments
        if (nrow(iap5) >=2)   #only apply if minimum 2 segments
        {  
          # if (!is.na(iap5$LA90[1]) & !is.na(iap5$LA90[2]) & (iap5$LA90[2]/iap5$LA90[1] < (1-2*(1-LA.THR)) | iap5$LA90[1]/iap5$LA90[2] < (1-2*(1-LA.THR))))  #at beginning: if 1st segments is much smaller than 2nd segment; verify also that lumen area quantiles are not NA 
          if (!is.na(iap5$LA90[1]) & !is.na(iap5$LA90[2]) & iap5$LA90[1]/iap5$LA90[2] < (1-2*(1-LA.THR)))  #at beginning: if 1st segments is much smaller than 2nd segment; verify also that lumen area quantiles are not NA 
          {
            lines(x=c(iap5$X.START[1], iap5$X.END[1]), y=c(iap5$LA90[1], iap5$LA90[1]), lty=1, lwd=60, col="yellow")
          }  
          # if (!is.na(iap5$LA90[nrow(iap5)]) & !is.na(iap5$LA90[nrow(iap5)-1]) & (iap5$LA90[nrow(iap5)-1]/iap5$LA90[nrow(iap5)] < (1-2*(1-LA.THR)) | iap5$LA90[nrow(iap5)]/iap5$LA90[nrow(iap5)-1] < (1-2*(1-LA.THR))))  #at end: if last segments is much smaller than 2nd last segment; verify also that lumen area quantiles are not NA 
          if (!is.na(iap5$LA90[nrow(iap5)]) & !is.na(iap5$LA90[nrow(iap5)-1]) & iap5$LA90[nrow(iap5)]/iap5$LA90[nrow(iap5)-1] < (1-2*(1-LA.THR)))  #at end: if last segments is much smaller than 2nd last segment; verify also that lumen area quantiles are not NA
          {
            lines(x=c(iap5$X.START[nrow(iap5)], iap5$X.END[nrow(iap5)]), y=c(iap5$LA90[nrow(iap5)], iap5$LA90[nrow(iap5)]), lty=1, lwd=60, col="yellow")
          }
        }
        if (nrow(iap5) >= 3)   #for the non-marginal segments
        {  
          for (n in c(2:(nrow(iap5)-1)))
          { 
            if(!is.na(iap5$LA90[n]) & !is.na(iap5$LA90[n-1]) & !is.na(iap5$LA90[n+1]))   #verify that the considered lumen area quantiles are not NA
            {  
              # if (((2 * iap5$LA90[n]) / (iap5$LA90[n-1] + iap5$LA90[n+1])) < LA.THR | ((iap5$LA90[n-1] + iap5$LA90[n+1]) / (2 * iap5$LA90[n])) < LA.THR)   #highlight segment if considerable smaller or larger lumen area quantile than for the average of its neighboring segments 
              if ((2 * iap5$LA90[n] / (iap5$LA90[n-1] + iap5$LA90[n+1])) < LA.THR)   #highlight segment if considerable smaller or larger lumen area quantile than for the average of its neighboring segments 
              {  
                lines(x=c(iap5$X.START[n], iap5$X.END[n]), y=c(iap5$LA90[n], iap5$LA90[n]), lty=1, lwd=60, col="yellow")
              }
            }
          }
        }
        
        ### Add horizontal lines for all image-segments with their average lumen area percentile
        for (n in c(1:nrow(iap5)))
        {  
          lines(x=c(iap5$X.START[n], iap5$X.END[n]), y=c(iap5$LA90[n], iap5$LA90[n]), lty=2, lwd=6)
        }
        
        ### Add yellow circle for bands based on very few individual cells, positioned at y = corresponding lumen area 
        if (length(iap4) > 0)
        {  
          points(iap4$RADDIST.CONT, iap4$LA, pch=21, col="black", bg="yellow", lwd=5, cex=40)
        }
        
        ### Add red circle at beginning of previously interrupted (gap) profile at y = 0
        if (length(iap7) > 0)
        {  
          points(iap7, rep(lim[3]+10,length(iap7)), pch=21, col="black", bg="red", lwd=5, cex=40)
        }
        
        ### Highlight neighboring rings that are potential duplicated due to a visual cross-dating mistake
        if (length(iap8) > 0)
        {  
          for (n in c(1:length(iap8)))
          {  
            lines(y=iap2$LA[iap2$YEAR==iap8[n]], x=iap2$RADDIST.CONT[iap2$YEAR==iap8[n]], lwd=60, col="yellow")
          }
        }
        
        ### Add green triangles to mark potential (large) outliers
        lines(y=iap9$LA, x=iap9$RADDIST.CONT, lwd=3)   #add smoothed reference line
        # iap10 <- which(iap2$LA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]/iap9$LA > OUT.THR & iap9$LA > quantile(iap9$LA, OUT.QU, na.rm=TRUE))   #extract (large) outliers
        # iap10 <- which(iap2$LA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]-iap9$LA > max(iap9$LA,na.rm=TRUE)/3 & iap9$LA > quantile(iap9$LA, OUT.QU, na.rm=TRUE))   #extract (large) outliers
        iap10 <- which((iap2$LA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]/iap9$LA > OUT.THR | iap2$LA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]-iap9$LA > quantile(iap9$LA,0.99,na.rm=TRUE)/2) & iap9$LA > quantile(iap9$LA, OUT.QU, na.rm=TRUE))   #extract (large) outliers
        points(x=iap2$RADDIST.CONT[iap10], y=iap2$LA[iap10], pch=25, col="black", bg="green", lwd=5, cex=30)   #plot outliers
      }
      
      lines(y=iap2$LA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Lumen Area (LA) (µm2)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()

      
      ### 4.6.1.27 Plot and save intra-annual profile of DRAD (radial lumen diameter) ####
      png(file=paste(woodid, "_Intraprofile_DRAD_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$DRAD[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Radial lumen diameter (DRAD) (µm)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(25,1500))
           cex.axis=3, cex.lab=3, ylim=c(0,max(iap2$DRAD[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$DRAD[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Radial lumen diameter (DRAD) (µm)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()

      
      ### 4.6.1.28 Plot and save intra-annual profile of DTAN (tangential lumen diameter) ####
      png(file=paste(woodid, "_Intraprofile_DTAN_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$DTAN[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Tangential lumen diameter (DTAN) (µm)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(25,1500))
           cex.axis=3, cex.lab=3, ylim=c(0,max(iap2$DTAN[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$DTAN[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Tangential lumen diameter (DTAN) (µm)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()

      
      ### 4.6.1.29 Plot and save intra-annual profile of TCA (total cell area) ####
      png(file=paste(woodid, "_Intraprofile_TCA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)        
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$TCA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Cell Area (TCA) (µm2)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0,4000))
           cex.axis=3, cex.lab=3, ylim=c(25,max(iap2$TCA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$TCA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Cell Area (TCA) (µm2)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()    
      
      
      ### 4.6.1.30 Plot and save intra-annual profile of CWA (cell wall area) ####
      png(file=paste(woodid, "_Intraprofile_CWA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$CWA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Cell Wall Area (CWA) (µm2)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0,2000))
           cex.axis=3, cex.lab=3, ylim=c(25,max(iap2$CWA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$CWA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Cell Wall Area (CWA) (µm2)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()      
      
      
      ### 4.6.1.31 Plot and save intra-annual profile of CWAACC (accumulated cell wall area per band) ####
      png(file=paste(woodid, "_Intraprofile_CWAACC_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$CWAACC[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Accumulated Cell Wall Area (CWAACC) (µm2)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0,2000))
           cex.axis=3, cex.lab=3, ylim=c(25,max(iap2$CWAACC[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$CWAACC[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Accumulated Cell Wall Area (CWAACC) (µm2)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()  
      
      
      ### 4.6.1.32 Plot and save intra-annual profile of CWTALL (mean cell wall thickness) ####
      png(file=paste(woodid, "_Intraprofile_CWTALL_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$CWTALL[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Overall Cell Wall thickness (CWTALL) (µm)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(1.5,10.0))
           cex.axis=3, cex.lab=3, ylim=c(1.5,max(iap2$CWTALL[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$CWTALL[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Overall Cell Wall thickness (CWTALL) (µm)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()        
      
      
      ### 4.6.1.33 Plot and save intra-annual profile of CWTALL.ADJ (adjusted mean cell wall thickness) ####
      png(file=paste(woodid, "_Intraprofile_CWTALL.ADJ_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$CWTALL.ADJ[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Overall adjusted Cell Wall thickness (CWTALL.ADJ) (µm)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(1.5,10.0))
           cex.axis=3, cex.lab=3, ylim=c(1.5,max(iap2$CWTALL.ADJ[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$CWTALL.ADJ.ADJ[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Overall adjusted Cell Wall thickness (CWTALL.ADJ) (µm)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()         
 
      
      ### 4.6.1.34 Plot and save intra-annual profile of CWTRAD (cell wall thickness of radial wall) ####
      png(file=paste(woodid, "_Intraprofile_CWTRAD_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$CWTRAD[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Radial Cell Wall thickness (CWTRAD) (µm)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(1.5,10.0))
           cex.axis=3, cex.lab=3, ylim=c(1.5,max(iap2$CWTRAD[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$CWTRAD[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Radial Cell Wall thickness (CWTRAD) (µm)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()    
      
      
      ### 4.6.1.35 Plot and save intra-annual profile of CWTTAN (cell wall thickness of tangential wall) ####
      png(file=paste(woodid, "_Intraprofile_CWTTAN_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$CWTTAN[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Tangential Cell Wall thickness (CWTTAN) (µm)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(1.5,9.0))
           cex.axis=3, cex.lab=3, ylim=c(1.5,max(iap2$CWTRAD[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$CWTTAN[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Tangential Cell Wall thickness (CWTTAN) (µm)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()

      
      ### 4.6.1.36 Plot and save intra-annual profile of RTSR (Mork's index) ####
      png(file=paste(woodid, "_Intraprofile_Mork_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$RTSR[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Mork's index", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0.2,7))
           cex.axis=3, cex.lab=3, ylim=c(0.2,max(iap2$RTSR[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$RTSR[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Mork's index", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()
      
      
      ### 4.6.1.37 Plot and save intra-annual profile of CTSR (adjusted Mork's index) ####
      png(file=paste(woodid, "_Intraprofile_Mork.ADJ_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$CTSR[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Adjusted Mork's index (4*CWTALL/D)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0.2,7))
           cex.axis=3, cex.lab=3, ylim=c(0.2,max(iap2$CTSR[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$CTSR[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Adjusted Mork's index (4*CWTALL/D)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()
      
      
      ### 4.6.1.38 Plot and save intra-annual profile of DCWT (CWT-based relative anatomical density) ####
      png(file=paste(woodid, "_Intraprofile_DCWT_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$DCWT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Rel. Anatomical Density (CWT based)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0.2,0.9))
           cex.axis=3, cex.lab=3, ylim=c(0.1,max(iap2$DCWT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$DCWT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Rel. Anatomical Density (CWT based)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()
      
      
      ### 4.6.1.39 Plot and save intra-annual profile of DCWT2 (special relative anatomical density: CWTRAD/DRAD) ####
      png(file=paste(woodid, "_Intraprofile_DCWT2_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$DCWT2[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],
           xaxt="n", yaxt="n", xlab="", ylab="Rel. Anatomical Density (CWTRAD/DRAD)", type="l", lwd=0, col=LINE.COL(length(RESO))[r],
           # cex.axis=3, cex.lab=3, ylim=c(0.2,0.9))
           cex.axis=3, cex.lab=3, ylim=c(0,max(iap2$DCWT2[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$DCWT2[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],
            xaxt="n", yaxt="n", xlab="", ylab="Rel. Anatomical Density (CWTRAD/DRAD)", type="l", lwd=4, col=LINE.COL(length(RESO))[r],
            cex.axis=3, cex.lab=3)
      par(xaxt="s")
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()
      
      
      ### 4.6.1.40 Plot and save intra-annual profile of DCWA (CWA-based relative anatomical density=best estimate!) ####
      png(file=paste(woodid, "_Intraprofile_DCWA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$DCWA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Rel. Anatomical Density (CWA based)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0.1,0.95))
           cex.axis=3, cex.lab=3, ylim=c(0.1,max(iap2$DCWA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$DCWA[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Rel. Anatomical Density (CWA based)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()    
      
      
      ### 4.6.1.41 Plot and save intra-annual profile of TB2 (cell wall reinforcement index) ####
      png(file=paste(woodid, "_Intraprofile_TB2_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$TB2[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Cell wall reinforcement index", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0.1,0.95))
           cex.axis=3, cex.lab=3, ylim=c(0,max(iap2$TB2[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$TB2[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Cell wall reinforcement index", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()        
      
      
      ### 4.6.1.42 Plot and save intra-annual profile of KH (theoretical hydraulic conductivity) ####
      png(file=paste(woodid, "_Intraprofile_KH_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$KH[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Theor. hydraul. conductivity (m4*s-1*MPa-1)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0.1,0.95))
           cex.axis=3, cex.lab=3, ylim=c(0,max(iap2$KH[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$KH[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Theor. hydraul. conductivity (m4*s-1*MPa-1)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()        
      
      
      ### 4.6.1.43 Plot and save intra-annual profile of DH (mean hydraulic diameter) ####
      png(file=paste(woodid, "_Intraprofile_DH_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$DH[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Mean hydraulic diameter (µm)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0.1,0.95))
           cex.axis=3, cex.lab=3, ylim=c(0,max(iap2$DH[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$DH[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Mean hydraulic diameter (µm)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()        
      
      
      ### 4.6.1.44 Plot and save intra-annual profile of N.BAND (Number of cells per band, allowing duplicates) ####
      png(file=paste(woodid, "_Intraprofile_N.BAND_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$N.BAND[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Number of cells per band (N.BAND)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0.1,0.95))
           cex.axis=3, cex.lab=3, ylim=c(0,max(iap2$N.BAND[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$N.BAND[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Number of cells per band (N.BAND)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()        
      
      
      ### 4.6.1.45 Plot and save intra-annual profile of N.TRUE (Number of cells per band, no duplicates for last band) ####
      png(file=paste(woodid, "_Intraprofile_N.TRUE_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=im.width,height=im.height)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=iap2$N.TRUE[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
           xaxt="n", yaxt="n", xlab="", ylab="Number of cells per band (N.TRUE)", type="l", lwd=0, col=LINE.COL(length(RESO))[r], 
           # cex.axis=3, cex.lab=3, ylim=c(0.1,0.95))
           cex.axis=3, cex.lab=3, ylim=c(0,max(iap2$N.TRUE[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END],na.rm=T)))
      if (PLOT.IMG.ID==TRUE)
      {  
        lim <- par("usr")
        rect(bottomright, lim[3]-1, topleft, lim[4]+1, border = "gray90", col = "gray90")
      }
      lines(y=iap2$N.TRUE[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            x=iap2$RADDIST.CONT[iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], 
            xaxt="n", yaxt="n", xlab="", ylab="Number of cells per band (N.TRUE)", type="l", lwd=4, col=LINE.COL(length(RESO))[r], 
            cex.axis=3, cex.lab=3)
      par(xaxt="s") 
      axis(1, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, cex.axis=2, labels=iap2$X.LABEL[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END])
      axis(3, at=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END]     
           , tck=0.015, lwd.ticks=2, labels=FALSE)
      axis(2, tck=0.015, lwd.ticks=2, labels=TRUE, cex.axis=2)
      axis(4, tck=0.015, lwd.ticks=2, labels=FALSE)
      abline(v=iap2$RADDIST.CONT[!is.na(iap2$X.LABEL) & iap2$YEAR>=YR.START & iap2$YEAR<=YR.END], lwd=2, lty=3)
      dev.off()        
      
      
      ### 4.6.1.45 Prepare calculations of annual statistics ####
      ### Make sure NaNs are replaced by NAs
      iap2 <- data.frame(lapply(iap2, function(x){replace(x, is.infinite(x), NA)}))
      iap2 <- data.frame(lapply(iap2, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
  
  
      ### Add column that gives information whether ring includes only earlwood or only latewood bands
      iap2 <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), mutate, ONLYEW=RobustMax(RTSR))
      iap2 <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), mutate, ONLYLW=RobustMin(RTSR))
      iap2$EWLW.ID2 <- ifelse(iap2$ONLYEW<1,"only.ew", ifelse(iap2$ONLYLW>1,"only.lw","NA"))
      iap2$ONLYEW <- NULL
      iap2$ONLYLW <- NULL
  
      
      ### 4.6.1.46 Calculate annual ring width (MRW), earlywood width (EWW), latewood width(LWW) ####
      mrw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MRW=mean(MRW, na.rm=T)) 
      eww <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, EWW=mean(EWW, na.rm=T)) 
      lww <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, LWW=mean(LWW, na.rm=T)) 


      ### 4.6.1.47 Calculate annual cell lumen area (LA) statistics ####
      mla <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MLA=mean(LA, na.rm=TRUE))    
      mla <- data.frame(lapply(mla, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mla.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MLA.EW=mean(LA[EWLW.ID=="ew"], na.rm=TRUE))    
      mla.ew <- data.frame(lapply(mla.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mla.ew$MLA.EW <- ifelse(is.na(mla.ew$MLA.EW),mla$MLA,mla.ew$MLA.EW)
      
      mla.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MLA.LW=mean(LA[EWLW.ID=="lw"], na.rm=TRUE))    
      mla.lw <- data.frame(lapply(mla.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mla.lw$MLA.LW <- ifelse(is.na(mla.lw$MLA.LW),mla$MLA,mla.lw$MLA.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$LA),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXLA=RobustMax(LA), 
                      RPOS.MXLA=RRADDISTR[which.max(LA)], 
                      APOS.MXLA=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(LA)], 
                      BAND.MXLA=round(RADDISTR.BAND[which.max(LA)])),              
               by=key(dt)]   
      mxla <- as.data.frame(dt)   
      mxla <- merge(mxla, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="lw"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.ew",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$LA),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNLA=RobustMin(LA), 
                      RPOS.MNLA=RRADDISTR[which.min(LA)], 
                      APOS.MNLA=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(LA)], 
                      BAND.MNLA=round(RADDISTR.BAND[which.min(LA)])),  
               by=key(dt)]    
      mnla <- as.data.frame(dt)   
      mnla <- merge(mnla, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA

      
      ### 4.6.1.48 Calculate annual total cell area (TCA) statistics ####
      mtca <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MTCA=mean(TCA, na.rm=TRUE))    
      mtca <- data.frame(lapply(mtca, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mtca.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MTCA.EW=mean(TCA[EWLW.ID=="ew"], na.rm=TRUE))    
      mtca.ew <- data.frame(lapply(mtca.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mtca.ew$MTCA.EW <- ifelse(is.na(mtca.ew$MTCA.EW),mtca$MTCA, mtca.ew$MTCA.EW)
      
      mtca.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MTCA.LW=mean(TCA[EWLW.ID=="lw"], na.rm=TRUE))    
      mtca.lw <- data.frame(lapply(mtca.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mtca.lw$MTCA.LW <- ifelse(is.na(mtca.lw$MTCA.LW),mtca$MTCA, mtca.lw$MTCA.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$TCA),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXTCA=RobustMax(TCA), 
                      RPOS.MXTCA=RRADDISTR[which.max(TCA)], 
                      APOS.MXTCA=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(TCA)], 
                      BAND.MXTCA=round(RADDISTR.BAND[which.max(TCA)])),  
               by=key(dt)]     
      mxtca <- as.data.frame(dt)   
      mxtca <- merge(mxtca, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="lw"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.ew",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$TCA),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNTCA=RobustMin(TCA), 
                      RPOS.MNTCA=RRADDISTR[which.min(TCA)], 
                      APOS.MNTCA=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(TCA)], 
                      BAND.MNTCA=round(RADDISTR.BAND[which.min(TCA)])),  
               by=key(dt)]    
      mntca <- as.data.frame(dt)   
      mntca <- merge(mntca, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      
      ### 4.6.1.49 Calculate annual radial cell lumen diameter (DRAD) statistics ####
      mdrad <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, DRAD=mean(DRAD, na.rm=TRUE))    
      mdrad <- data.frame(lapply(mdrad, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mdrad.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, DRAD.EW=mean(DRAD[EWLW.ID=="ew"], na.rm=TRUE))    
      mdrad.ew <- data.frame(lapply(mdrad.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mdrad.ew$DRAD.EW <- ifelse(is.na(mdrad.ew$DRAD.EW),mdrad$DRAD,mdrad.ew$DRAD.EW)
      
      mdrad.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, DRAD.LW=mean(DRAD[EWLW.ID=="lw"], na.rm=TRUE))    
      mdrad.lw <- data.frame(lapply(mdrad.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mdrad.lw$DRAD.LW <- ifelse(is.na(mdrad.lw$DRAD.LW),mdrad$DRAD,mdrad.lw$DRAD.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$DRAD),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXDRAD=RobustMax(DRAD), 
                      RPOS.MXDRAD=RRADDISTR[which.max(DRAD)], 
                      APOS.MXDRAD=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(DRAD)], 
                      BAND.MXDRAD=round(RADDISTR.BAND[which.max(DRAD)])),              
               by=key(dt)]   
      mxdrad <- as.data.frame(dt)   
      mxdrad <- merge(mxdrad, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="lw"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.ew",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$DRAD),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNDRAD=RobustMin(DRAD), 
                      RPOS.MNDRAD=RRADDISTR[which.min(DRAD)], 
                      APOS.MNDRAD=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(DRAD)], 
                      BAND.MNDRAD=round(RADDISTR.BAND[which.min(DRAD)])),  
               by=key(dt)]    
      mndrad <- as.data.frame(dt)   
      mndrad <- merge(mndrad, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA  
      
      
      ### 4.6.1.50 Calculate annual tangential cell lumen diameter (DTAN) statistics ####
      mdtan <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, DTAN=mean(DTAN, na.rm=TRUE))    
      mdtan <- data.frame(lapply(mdtan, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mdtan.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, DTAN.EW=mean(DTAN[EWLW.ID=="ew"], na.rm=TRUE))    
      mdtan.ew <- data.frame(lapply(mdtan.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mdtan.ew$DTAN.EW <- ifelse(is.na(mdtan.ew$DTAN.EW),mdtan$DTAN,mdtan.ew$DTAN.EW)
      
      mdtan.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, DTAN.LW=mean(DTAN[EWLW.ID=="lw"], na.rm=TRUE))    
      mdtan.lw <- data.frame(lapply(mdtan.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mdtan.lw$DTAN.LW <- ifelse(is.na(mdtan.lw$DTAN.LW),mdtan$DTAN,mdtan.lw$DTAN.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$DTAN),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXDTAN=RobustMax(DTAN), 
                      RPOS.MXDTAN=RRADDISTR[which.max(DTAN)], 
                      APOS.MXDTAN=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(DTAN)], 
                      BAND.MXDTAN=round(RADDISTR.BAND[which.max(DTAN)])),              
               by=key(dt)]   
      mxdtan <- as.data.frame(dt)   
      mxdtan <- merge(mxdtan, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="lw"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.ew",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$DTAN),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNDTAN=RobustMin(DTAN), 
                      RPOS.MNDTAN=RRADDISTR[which.min(DTAN)], 
                      APOS.MNDTAN=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(DTAN)], 
                      BAND.MNDTAN=round(RADDISTR.BAND[which.min(DTAN)])),  
               by=key(dt)]    
      mndtan <- as.data.frame(dt)   
      mndtan <- merge(mndtan, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA  

      
      ### 4.6.1.51 Calculate annual cell wall area (CWA) statistics ####
      mcwa <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWA=mean(CWA, na.rm=TRUE))    
      mcwa <- data.frame(lapply(mcwa, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mcwa.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWA.EW=mean(CWA[EWLW.ID=="ew"], na.rm=TRUE))    
      mcwa.ew <- data.frame(lapply(mcwa.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mcwa.ew$MCWA.EW <- ifelse(is.na(mcwa.ew$MCWA.EW),mcwa$MCWA, mcwa.ew$MCWA.EW)
      
      mcwa.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWA.LW=mean(CWA[EWLW.ID=="lw"], na.rm=TRUE))    
      mcwa.lw <- data.frame(lapply(mcwa.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mcwa.lw$MCWA.LW <- ifelse(is.na(mcwa.lw$MCWA.LW),mcwa$MCWA, mcwa.lw$MCWA.LW)
      
      dt <- NULL
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWA),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXCWA=RobustMax(CWA), 
                      RPOS.MXCWA=RRADDISTR[which.max(CWA)], 
                      APOS.MXCWA=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(CWA)], 
                      BAND.MXCWA=round(RADDISTR.BAND[which.max(CWA)])),  
               by=key(dt)]  
      mxcwa <- as.data.frame(dt)   
      mxcwa <- merge(mxcwa, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWA),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNCWA=RobustMin(CWA), 
                      RPOS.MNCWA=RRADDISTR[which.min(CWA)], 
                      APOS.MNCWA=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(CWA)], 
                      BAND.MNCWA=round(RADDISTR.BAND[which.min(CWA)])),  
               by=key(dt)]    
      mncwa <- as.data.frame(dt)   
      mncwa <- merge(mncwa, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA

      
      ### 4.6.1.52 Calculate annual accumulated cell wall area (CWAACC) statistics ####
      cwaacc <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, CWAACC=sum(CWAACC, na.rm=TRUE))    
      cwaacc <- data.frame(lapply(cwaacc, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      cwaacc.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, CWAACC.EW=sum(CWAACC[EWLW.ID=="ew"], na.rm=TRUE))    
      cwaacc.ew <- data.frame(lapply(cwaacc.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      cwaacc.ew$CWAACC.EW <- ifelse(is.na(cwaacc.ew$CWAACC.EW),cwaacc$CWAACC,cwaacc.ew$CWAACC.EW)
      
      cwaacc.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, CWAACC.LW=sum(CWAACC[EWLW.ID=="lw"], na.rm=TRUE))    
      cwaacc.lw <- data.frame(lapply(cwaacc.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      cwaacc.lw$CWAACC.LW <- ifelse(cwaacc.lw$CWAACC.LW==0,NA,cwaacc.lw$CWAACC.LW)
      # cwaacc.lw$CWAACC.LW <- ifelse(is.na(cwaacc.lw$CWAACC.LW),cwaacc$CWAACC,cwaacc.lw$CWAACC.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWAACC),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXCWAACC=RobustMax(CWAACC), 
                      RPOS.MXCWAACC=RRADDISTR[which.max(CWAACC)], 
                      APOS.MXCWAACC=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(CWAACC)], 
                      BAND.MXCWAACC=round(RADDISTR.BAND[which.max(CWAACC)])),              
               by=key(dt)]   
      mxcwaacc <- as.data.frame(dt)   
      mxcwaacc <- merge(mxcwaacc, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="lw"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.ew",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWAACC),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNCWAACC=RobustMin(CWAACC), 
                      RPOS.MNCWAACC=RRADDISTR[which.min(CWAACC)], 
                      APOS.MNCWAACC=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(CWAACC)], 
                      BAND.MNCWAACC=round(RADDISTR.BAND[which.min(CWAACC)])),  
               by=key(dt)]    
      mncwaacc <- as.data.frame(dt)   
      mncwaacc <- merge(mncwaacc, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA    
      
      
      ### 4.6.1.53 Calculate annual mean cell wall thickness (CWTALL) statistics ####
      mcwtall <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTALL=mean(CWTALL, na.rm=TRUE))    
      mcwtall <- data.frame(lapply(mcwtall, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mcwtall.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTALL.EW=mean(CWTALL[EWLW.ID=="ew"], na.rm=TRUE))    
      mcwtall.ew <- data.frame(lapply(mcwtall.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mcwtall.ew$MCWTALL.EW <- ifelse(is.na(mcwtall.ew$MCWTALL.EW),mcwtall$MCWTALL, mcwtall.ew$MCWTALL.EW)
      
      mcwtall.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTALL.LW=mean(CWTALL[EWLW.ID=="lw"], na.rm=TRUE))    
      mcwtall.lw <- data.frame(lapply(mcwtall.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mcwtall.lw$MCWTALL.LW <- ifelse(is.na(mcwtall.lw$MCWTALL.LW),mcwtall$MCWTALL, mcwtall.lw$MCWTALL.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWTALL),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXCWTALL=RobustMax(CWTALL), 
                      RPOS.MXCWTALL=RRADDISTR[which.max(CWTALL)], 
                      APOS.MXCWTALL=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(CWTALL)], 
                      BAND.MXCWTALL=round(RADDISTR.BAND[which.max(CWTALL)])),  
               by=key(dt)]  
      mxcwtall <- as.data.frame(dt)   
      mxcwtall <- merge(mxcwtall, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWTALL),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNCWTALL=RobustMin(CWTALL), 
                      RPOS.MNCWTALL=RRADDISTR[which.min(CWTALL)], 
                      APOS.MNCWTALL=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(CWTALL)], 
                      BAND.MNCWTALL=round(RADDISTR.BAND[which.min(CWTALL)])), 
               by=key(dt)]    
      mncwtall <- as.data.frame(dt)   
      mncwtall <- merge(mncwtall, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      
      ### 4.6.1.54 Calcuate annual mean cell wall thickness adjusted (CWTALL.ADJ) statistics ####
      mcwtalladj <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTALL.ADJ=mean(CWTALL.ADJ, na.rm=TRUE))    
      mcwtalladj <- data.frame(lapply(mcwtalladj, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mcwtalladj.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTALL.ADJ.EW=mean(CWTALL.ADJ[EWLW.ID=="ew"], na.rm=TRUE))    
      mcwtalladj.ew <- data.frame(lapply(mcwtalladj.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mcwtalladj.ew$MCWTALL.ADJ.EW <- ifelse(is.na(mcwtalladj.ew$MCWTALL.ADJ.EW),mcwtalladj$MCWTALL.ADJ, mcwtalladj.ew$MCWTALL.ADJ.EW)
      
      mcwtalladj.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTALL.ADJ.LW=mean(CWTALL.ADJ[EWLW.ID=="lw"], na.rm=TRUE))    
      mcwtalladj.lw <- data.frame(lapply(mcwtalladj.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mcwtalladj.lw$MCWTALL.ADJ.LW <- ifelse(is.na(mcwtalladj.lw$MCWTALL.ADJ.LW),mcwtalladj$MCWTALL.ADJ, mcwtalladj.lw$MCWTALL.ADJ.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWTALL.ADJ),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXCWTALL.ADJ=RobustMax(CWTALL.ADJ), 
                      RPOS.MXCWTALL.ADJ=RRADDISTR[which.max(CWTALL.ADJ)], 
                      APOS.MXCWTALL.ADJ=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(CWTALL.ADJ)], 
                      BAND.MXCWTALL.ADJ=round(RADDISTR.BAND[which.max(CWTALL.ADJ)])),
               by=key(dt)]    
      mxcwtalladj <- as.data.frame(dt)   
      mxcwtalladj <- merge(mxcwtalladj, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWTALL.ADJ),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNCWTALL.ADJ=RobustMin(CWTALL.ADJ), 
                      RPOS.MNCWTALL.ADJ=RRADDISTR[which.min(CWTALL.ADJ)], 
                      APOS.MNCWTALL.ADJ=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(CWTALL.ADJ)], 
                      BAND.MNCWTALL.ADJ=round(RADDISTR.BAND[which.min(CWTALL.ADJ)])), 
               by=key(dt)]    
      mncwtalladj <- as.data.frame(dt)   
      mncwtalladj <- merge(mncwtalladj, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      
      ### 4.6.1.55 Calculate annual tangential cell wall thickness (CWTTAN) statistics ####
      mcwttan <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTTAN=mean(CWTTAN, na.rm=TRUE))    
      mcwttan <- data.frame(lapply(mcwttan, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mcwttan.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTTAN.EW=mean(CWTTAN[EWLW.ID=="ew"], na.rm=TRUE))    
      mcwttan.ew <- data.frame(lapply(mcwttan.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mcwttan.ew$MCWTTAN.EW <- ifelse(is.na(mcwttan.ew$MCWTTAN.EW),mcwttan$MCWTTAN, mcwttan.ew$MCWTTAN.EW)
      
      mcwttan.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTTAN.LW=mean(CWTTAN[EWLW.ID=="lw"], na.rm=TRUE))    
      mcwttan.lw <- data.frame(lapply(mcwttan.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mcwttan.lw$MCWTTAN.LW <- ifelse(is.na(mcwttan.lw$MCWTTAN.LW),mcwttan$MCWTTAN, mcwttan.lw$MCWTTAN.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWTTAN),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXCWTTAN=RobustMax(CWTTAN), 
                      RPOS.MXCWTTAN=RRADDISTR[which.max(CWTTAN)], 
                      APOS.MXCWTTAN=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(CWTTAN)], 
                      BAND.MXCWTTAN=round(RADDISTR.BAND[which.max(CWTTAN)])), 
               by=key(dt)]   
      mxcwttan <- as.data.frame(dt)   
      mxcwttan <- merge(mxcwttan, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWTTAN),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNCWTTAN=RobustMin(CWTTAN), 
                      RPOS.MNCWTTAN=RRADDISTR[which.min(CWTTAN)], 
                      APOS.MNCWTTAN=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(CWTTAN)], 
                      BAND.MNCWTTAN=round(RADDISTR.BAND[which.min(CWTTAN)])), 
               by=key(dt)]    
      mncwttan <- as.data.frame(dt)   
      mncwttan <- merge(mncwttan, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      
      ### 4.6.1.56 Calculate annual radial cell wall thickness (CWTRAD) statistics ####
      mcwtrad <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTRAD=mean(CWTRAD, na.rm=TRUE))    
      mcwtrad <- data.frame(lapply(mcwtrad, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mcwtrad.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTRAD.EW=mean(CWTRAD[EWLW.ID=="ew"], na.rm=TRUE))    
      mcwtrad.ew <- data.frame(lapply(mcwtrad.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mcwtrad.ew$MCWTRAD.EW <- ifelse(is.na(mcwtrad.ew$MCWTRAD.EW),mcwtrad$MCWTRAD, mcwtrad.ew$MCWTRAD.EW)
      
      mcwtrad.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCWTRAD.LW=mean(CWTRAD[EWLW.ID=="lw"], na.rm=TRUE))    
      mcwtrad.lw <- data.frame(lapply(mcwtrad.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mcwtrad.lw$MCWTRAD.LW <- ifelse(is.na(mcwtrad.lw$MCWTRAD.LW),mcwtrad$MCWTRAD, mcwtrad.lw$MCWTRAD.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWTRAD),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXCWTRAD=RobustMax(CWTRAD), 
                      RPOS.MXCWTRAD=RRADDISTR[which.max(CWTRAD)], 
                      APOS.MXCWTRAD=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(CWTRAD)], 
                      BAND.MXCWTRAD=round(RADDISTR.BAND[which.max(CWTRAD)])), 
               by=key(dt)]   
      mxcwtrad <- as.data.frame(dt)   
      mxcwtrad <- merge(mxcwtrad, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CWTRAD),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNCWTRAD=RobustMin(CWTRAD), 
                      RPOS.MNCWTRAD=RRADDISTR[which.min(CWTRAD)], 
                      APOS.MNCWTRAD=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(CWTRAD)], 
                      BAND.MNCWTRAD=round(RADDISTR.BAND[which.min(CWTRAD)])), 
               by=key(dt)]   
      mncwtrad <- as.data.frame(dt)   
      mncwtrad <- merge(mncwtrad, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      
      ### 4.6.1.57 Calculate annual Mork's index (RTSR) statistics ####
      mmork <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MRTSR=mean(RTSR, na.rm=TRUE))    
      mmork <- data.frame(lapply(mmork, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA
      
      mmork.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MRTSR.EW=mean(RTSR[EWLW.ID=="ew"], na.rm=TRUE))    
      mmork.ew <- data.frame(lapply(mmork.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mmork.ew$MRTSR.EW <- ifelse(is.na(mmork.ew$MRTSR.EW),mmork$MRTSR, mmork.ew$MRTSR.EW)
      
      mmork.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MRTSR.LW=mean(RTSR[EWLW.ID=="lw"], na.rm=TRUE))    
      mmork.lw <- data.frame(lapply(mmork.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mmork.lw$MRTSR.LW <- ifelse(is.na(mmork.lw$MRTSR.LW),mmork$MRTSR, mmork.lw$MRTSR.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)      
      dt <- dt[!is.na(dt$RTSR),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXRTSR=RobustMax(RTSR), 
                      RPOS.MXRTSR=RRADDISTR[which.max(RTSR)], 
                      APOS.MXRTSR=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(RTSR)], 
                      BAND.MXRTSR=round(RADDISTR.BAND[which.max(RTSR)])), 
               by=key(dt)]    
      mxmork <- as.data.frame(dt)   
      mxmork <- merge(mxmork, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$RTSR),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNRTSR=RobustMin(RTSR), 
                      RPOS.MNRTSR=RRADDISTR[which.min(RTSR)], 
                      APOS.MNRTSR=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(RTSR)], 
                      BAND.MNRTSR=round(RADDISTR.BAND[which.min(RTSR)])), 
               by=key(dt)]    
      mnmork <- as.data.frame(dt)   
      mnmork <- merge(mnmork, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      
      ### 4.6.1.58 Calculate annual circular thickness-to-span ratio (CTSR) statistics ####
      mctsr <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCTSR=mean(CTSR, na.rm=TRUE))    
      mctsr <- data.frame(lapply(mctsr, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mctsr.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCTSR.EW=mean(CTSR[EWLW.ID=="ew"], na.rm=TRUE))    
      mctsr.ew <- data.frame(lapply(mctsr.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mctsr.ew$MCTSR.EW <- ifelse(is.na(mctsr.ew$MCTSR.EW),mctsr$MCTSR, mctsr.ew$MCTSR.EW)
      
      mctsr.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MCTSR.LW=mean(CTSR[EWLW.ID=="lw"], na.rm=TRUE))    
      mctsr.lw <- data.frame(lapply(mctsr.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mctsr.lw$MCTSR.LW <- ifelse(is.na(mctsr.lw$MCTSR.LW),mctsr$MCTSR, mctsr.lw$MCTSR.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CTSR),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXCTSR=RobustMax(CTSR), 
                      RPOS.MXCTSR=RRADDISTR[which.max(CTSR)], 
                      APOS.MXCTSR=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(CTSR)], 
                      BAND.MXCTSR=round(RADDISTR.BAND[which.max(CTSR)])), 
               by=key(dt)]   
      mxctsr <- as.data.frame(dt)   
      mxctsr <- merge(mxctsr, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$CTSR),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNCTSR=RobustMin(CTSR), 
                      RPOS.MNCTSR=RRADDISTR[which.min(CTSR)], 
                      APOS.MNCTSR=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(CTSR)], 
                      BAND.MNCTSR=round(RADDISTR.BAND[which.min(CTSR)])), 
               by=key(dt)]    
      mnctsr <- as.data.frame(dt)   
      mnctsr <- merge(mnctsr, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      

      ### 4.6.1.59 Caculate annual relative anatomical density based on CWT (DCWT) statistics ####
      mdcwt <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDCWT=mean(DCWT, na.rm=TRUE))    
      mdcwt <- data.frame(lapply(mdcwt, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mdcwt.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDCWT.EW=mean(DCWT[EWLW.ID=="ew"], na.rm=TRUE))    
      mdcwt.ew <- data.frame(lapply(mdcwt.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mdcwt.ew$MDCWT.EW <- ifelse(is.na(mdcwt.ew$MDCWT.EW),mdcwt$MDCWT, mdcwt.ew$MDCWT.EW)
  
      mdcwt.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDCWT.LW=mean(DCWT[EWLW.ID=="lw"], na.rm=TRUE))    
      mdcwt.lw <- data.frame(lapply(mdcwt.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mdcwt.lw$MDCWT.LW <- ifelse(is.na(mdcwt.lw$MDCWT.LW),mdcwt$MDCWT, mdcwt.lw$MDCWT.LW)
  
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$DCWT),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXDCWT=RobustMax(DCWT), 
                      RPOS.MXDCWT=RRADDISTR[which.max(DCWT)], 
                      APOS.MXDCWT=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(DCWT)], 
                      BAND.MXDCWT=round(RADDISTR.BAND[which.max(DCWT)])), 
               by=key(dt)]    
      mxdcwt <- as.data.frame(dt)   
      mxdcwt <- merge(mxdcwt, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
     
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$DCWT),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNDCWT=RobustMin(DCWT), 
                      RPOS.MNDCWT=RRADDISTR[which.min(DCWT)], 
                      APOS.MNDCWT=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(DCWT)], 
                      BAND.MNDCWT=round(RADDISTR.BAND[which.min(DCWT)])), 
               by=key(dt)]   
      mndcwt <- as.data.frame(dt)   
      mndcwt <- merge(mndcwt, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
  
      
      ### 4.6.1.60 Calculate annual special relative anatomical density (DCWT2=CWTRAD/DRAD) statistics ####
      mdcwt2 <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDCWT2=mean(DCWT2, na.rm=TRUE))
      mdcwt2 <- data.frame(lapply(mdcwt2, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA

      mdcwt2.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDCWT2.EW=mean(DCWT2[EWLW.ID=="ew"], na.rm=TRUE))
      mdcwt2.ew <- data.frame(lapply(mdcwt2.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA
      mdcwt2.ew$MDCWT2.EW <- ifelse(is.na(mdcwt2.ew$MDCWT2.EW),mdcwt2$MDCWT2, mdcwt2.ew$MDCWT2.EW)

      mdcwt2.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDCWT2.LW=mean(DCWT2[EWLW.ID=="lw"], na.rm=TRUE))
      mdcwt2.lw <- data.frame(lapply(mdcwt2.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA
      # mdcwt2.lw$MDCWT2.LW <- ifelse(is.na(mdcwt2.lw$MDCWT2.LW),mdcwt2$MDCWT2, mdcwt2.lw$MDCWT2.LW)

      dt <- NULL
      dt <- data.table(iap2[iap2$EWLW.ID=="lw"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.ew",])
      dt <- dt[!is.na(dt$DCWT2),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXDCWT2=RobustMax(DCWT2),
                      RPOS.MXDCWT2=RRADDISTR[which.max(DCWT2)],
                      APOS.MXDCWT2=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(DCWT2)],
                      BAND.MXDCWT2=round(RADDISTR.BAND[which.max(DCWT2)])),
               by=key(dt)]
      mxdcwt2 <- as.data.frame(dt)
      mxdcwt2 <- merge(mxdcwt2, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA

      dt <- NULL
      dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- dt[!is.na(dt$DCWT2),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNDCWT2=RobustMin(DCWT2),
                      RPOS.MNDCWT2=RRADDISTR[which.min(DCWT2)],
                      APOS.MNDCWT2=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(DCWT2)],
                      BAND.MNDCWT2=round(RADDISTR.BAND[which.min(DCWT2)])),
               by=key(dt)]
      mndcwt2 <- as.data.frame(dt)
      mndcwt2 <- merge(mndcwt2, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
  
      ### 4.6.1.61 Calculate annual relative anatomical density based on CWA (DCWA) statistics ####
      mdcwa <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDCWA=mean(DCWA, na.rm=TRUE))    
      mdcwa <- data.frame(lapply(mdcwa, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
  
      mdcwa.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDCWA.EW=mean(DCWA[EWLW.ID=="ew"], na.rm=TRUE))    
      mdcwa.ew <- data.frame(lapply(mdcwa.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mdcwa.ew$MDCWA.EW <- ifelse(is.na(mdcwa.ew$MDCWA.EW),mdcwa$MDCWA,mdcwa.ew$MDCWA.EW)
  
      mdcwa.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDCWA.LW=mean(DCWA[EWLW.ID=="lw"], na.rm=TRUE))    
      mdcwa.lw <- data.frame(lapply(mdcwa.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mdcwa.lw$MDCWA.LW <- ifelse(is.na(mdcwa.lw$MDCWA.LW),mdcwa$MDCWA,mdcwa.lw$MDCWA.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$DCWA),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXDCWA=RobustMax(DCWA), 
                      RPOS.MXDCWA=RRADDISTR[which.max(DCWA)], 
                      APOS.MXDCWA=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(DCWA)], 
                      BAND.MXDCWA=round(RADDISTR.BAND[which.max(DCWA)])), 
               by=key(dt)]   
      mxdcwa <- as.data.frame(dt)   
      mxdcwa <- merge(mxdcwa, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$DCWA),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNDCWA=RobustMin(DCWA), 
                      RPOS.MNDCWA=RRADDISTR[which.min(DCWA)], 
                      APOS.MNDCWA=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(DCWA)], 
                      BAND.MNDCWA=round(RADDISTR.BAND[which.min(DCWA)])), 
               by=key(dt)]   
      mndcwa <- as.data.frame(dt)   
      mndcwa <- merge(mndcwa, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
  
  
      ### 4.6.1.62 Calculate annual cell reinforcement index, (t/b)2 (TB2) statistics ####
      mtb2 <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MTB2=mean(TB2, na.rm=TRUE))    
      mtb2 <- data.frame(lapply(mtb2, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mtb2.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MTB2.EW=mean(TB2[EWLW.ID=="ew"], na.rm=TRUE))    
      mtb2.ew <- data.frame(lapply(mtb2.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mtb2.ew$MTB2.EW <- ifelse(is.na(mtb2.ew$MTB2.EW),mtb2$MTB2,mtb2.ew$MTB2.EW)
      
      mtb2.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MTB2.LW=mean(TB2[EWLW.ID=="lw"], na.rm=TRUE))    
      mtb2.lw <- data.frame(lapply(mtb2.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mtb2.lw$MTB2.LW <- ifelse(is.na(mtb2.lw$MTB2.LW),mtb2$MTB2,mtb2.lw$MTB2.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$TB2),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXTB2=RobustMax(TB2), 
                      RPOS.MXTB2=RRADDISTR[which.max(TB2)], 
                      APOS.MXTB2=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(TB2)], 
                      BAND.MXTB2=round(RADDISTR.BAND[which.max(TB2)])), 
               by=key(dt)]  
      mxtb2 <- as.data.frame(dt)   
      mxtb2 <- merge(mxtb2, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$TB2),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNTB2=RobustMin(TB2), 
                      RPOS.MNTB2=RRADDISTR[which.min(TB2)], 
                      APOS.MNTB2=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(TB2)], 
                      BAND.MNTB2=round(RADDISTR.BAND[which.min(TB2)])), 
               by=key(dt)]    
      mntb2 <- as.data.frame(dt)   
      mntb2 <- merge(mntb2, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
  
  
      ### 4.6.1.63 Calculate annual theoretical hydraulic conductivity (KH) statistics ####
      mkh <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MKH=mean(KH, na.rm=TRUE))    
      mkh <- data.frame(lapply(mkh, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mkh.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MKH.EW=mean(KH[EWLW.ID=="ew"], na.rm=TRUE))    
      mkh.ew <- data.frame(lapply(mkh.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mkh.ew$MKH.EW <- ifelse(is.na(mkh.ew$MKH.EW),mkh$MKH,mkh.ew$MKH.EW)
      
      mkh.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MKH.LW=mean(KH[EWLW.ID=="lw"], na.rm=TRUE))    
      mkh.lw <- data.frame(lapply(mkh.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mkh.lw$MKH.LW <- ifelse(is.na(mkh.lw$MKH.LW),mkh$MKH,mkh.lw$MKH.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$KH),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXKH=RobustMax(KH), 
                      RPOS.MXKH=RRADDISTR[which.max(KH)], 
                      APOS.MXKH=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(KH)], 
                      BAND.MXKH=round(RADDISTR.BAND[which.max(KH)])), 
               by=key(dt)]  
      mxkh <- as.data.frame(dt)   
      mxkh <- merge(mxkh, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="lw"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.ew",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$KH),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNKH=RobustMin(KH), 
                      RPOS.MNKH=RRADDISTR[which.min(KH)], 
                      APOS.MNKH=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(KH)], 
                      BAND.MNKH=round(RADDISTR.BAND[which.min(KH)])),
               by=key(dt)]    
      mnkh <- as.data.frame(dt)   
      mnkh <- merge(mnkh, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
  
      
      ### 4.6.1.64 Calcuate annual mean hydraulic diameter (DH) statistics ####
      mdh <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDH=mean(DH, na.rm=TRUE))    
      mdh <- data.frame(lapply(mdh, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      mdh.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDH.EW=mean(DH[EWLW.ID=="ew"], na.rm=TRUE))    
      mdh.ew <- data.frame(lapply(mdh.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      mdh.ew$MDH.EW <- ifelse(is.na(mdh.ew$MDH.EW),mdh$MDH,mdh.ew$MDH.EW)
      
      mdh.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, MDH.LW=mean(DH[EWLW.ID=="lw"], na.rm=TRUE))    
      mdh.lw <- data.frame(lapply(mdh.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      # mdh.lw$MDH.LW <- ifelse(is.na(mdh.lw$MDH.LW),mdh$MDH,mdh.lw$MDH.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$DH),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXDH=RobustMax(DH), 
                      RPOS.MXDH=RRADDISTR[which.max(DH)], 
                      APOS.MXDH=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(DH)], 
                      BAND.MXDH=round(RADDISTR.BAND[which.max(DH)])), 
               by=key(dt)]   
      mxdh <- as.data.frame(dt)   
      mxdh <- merge(mxdh, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="lw"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.ew",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$DH),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNDH=RobustMin(DH), 
                      RPOS.MNDH=RRADDISTR[which.min(DH)], 
                      APOS.MNDH=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(DH)], 
                      BAND.MNDH=round(RADDISTR.BAND[which.min(DH)])), 
               by=key(dt)]   
      mndh <- as.data.frame(dt)   
      mndh <- merge(mndh, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA

      
      ### 4.6.1.65 Calculate annual number of cells (N.BAND) statistics ####
      # ncells.band <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, NCELLS.BAND=sum(N.BAND, na.rm=TRUE))
      # ncells.true <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, NCELLS.TRUE=sum(N.TRUE, na.rm=TRUE))
      ncells.band <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, N.BAND=sum(N.BAND, na.rm=TRUE))    
      ncells.band <- data.frame(lapply(ncells.band, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      ncells.band.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, N.BAND.EW=sum(N.BAND[EWLW.ID=="ew"], na.rm=TRUE))    
      ncells.band.ew <- data.frame(lapply(ncells.band.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      ncells.band.ew$N.BAND.EW <- ifelse(is.na(ncells.band.ew$N.BAND.EW),ncells.band$N.BAND,ncells.band.ew$N.BAND.EW)
      
      ncells.band.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, N.BAND.LW=sum(N.BAND[EWLW.ID=="lw"], na.rm=TRUE))    
      ncells.band.lw <- data.frame(lapply(ncells.band.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      ncells.band.lw$N.BAND.LW <- ifelse(ncells.band.lw$N.BAND.LW==0,NA,ncells.band.lw$N.BAND.LW)
      # ncells.band.lw$N.BAND.LW <- ifelse(is.na(ncells.band.lw$N.BAND.LW),ncells.band$N.BAND,ncells.band.lw$N.BAND.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$N.BAND),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXN.BAND=RobustMax(N.BAND), 
                      RPOS.MXN.BAND=RRADDISTR[which.max(N.BAND)], 
                      APOS.MXN.BAND=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(N.BAND)], 
                      BAND.MXN.BAND=round(RADDISTR.BAND[which.max(N.BAND)])),              
               by=key(dt)]   
      mxncells.band <- as.data.frame(dt)   
      mxncells.band <- merge(mxncells.band, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="lw"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.ew",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$N.BAND),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNN.BAND=RobustMin(N.BAND), 
                      RPOS.MNN.BAND=RRADDISTR[which.min(N.BAND)], 
                      APOS.MNN.BAND=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(N.BAND)], 
                      BAND.MNN.BAND=round(RADDISTR.BAND[which.min(N.BAND)]),
                      MRW=mean(MRW,na.rm=TRUE)),                      
               by=key(dt)]  
      
      dt$APOS.MNN.BAND <- ifelse(!is.na(dt$APOS.MNN.BAND), dt$APOS.MNN.BAND, (dt$MNN.BAND+RESO[r]/2) )   #if there minimum number of cells is 0, there is no RRADDISTR; assume center position of band instead
      dt$RPOS.MNN.BAND <- ifelse(!is.na(dt$RPOS.MNN.BAND), dt$RPOS.MNN.BAND, (dt$APOS.MNN.BAND*100/dt$MRW))   #if there minimum number of cells is 0, there is no RRADDISTR; assume center position of band instead
      dt$MRW <- NULL   #remove temporary column that was only used to troubleshoot NAs in RPOS.MNN.BAND
      
      mnncells.band <- as.data.frame(dt)   
      mnncells.band <- merge(mnncells.band, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA         
      
      
      ### 4.6.1.66 Calculate annual number of cells (N.TRUE) statistics ####
      ncells.true <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, N.TRUE=sum(N.TRUE, na.rm=TRUE))    
      ncells.true <- data.frame(lapply(ncells.true, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      ncells.true.ew <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, N.TRUE.EW=sum(N.TRUE[EWLW.ID=="ew"], na.rm=TRUE))    
      ncells.true.ew <- data.frame(lapply(ncells.true.ew, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      ncells.true.ew$N.TRUE.EW <- ifelse(is.na(ncells.true.ew$N.TRUE.EW),ncells.true$N.TRUE,ncells.true.ew$N.TRUE.EW)
      
      ncells.true.lw <- ddply(iap2, c("YEAR", "WOODID", "IMAGE"), summarise, N.TRUE.LW=sum(N.TRUE[EWLW.ID=="lw"], na.rm=TRUE))    
      ncells.true.lw <- data.frame(lapply(ncells.true.lw, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      ncells.true.lw$N.TRUE.LW <- ifelse(ncells.true.lw$N.TRUE.LW==0,NA,ncells.true.lw$N.TRUE.LW)
      # ncells.true.lw$N.TRUE.LW <- ifelse(is.na(ncells.true.lw$N.TRUE.LW),ncells.true$N.TRUE,ncells.true.lw$N.TRUE.LW)
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="ew"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.lw",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$N.TRUE),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MXN.TRUE=RobustMax(N.TRUE), 
                      RPOS.MXN.TRUE=RRADDISTR[which.max(N.TRUE)], 
                      APOS.MXN.TRUE=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.max(N.TRUE)], 
                      TRUE.MXN.TRUE=round(RADDISTR.BAND[which.max(N.TRUE)])),              
               by=key(dt)]   
      mxncells.true <- as.data.frame(dt)   
      mxncells.true <- merge(mxncells.true, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA
      
      dt <- NULL
      # dt <- data.table(iap2[iap2$EWLW.ID=="lw"|is.na(iap2$EWLW.ID)|iap2$EWLW.ID2=="only.ew",])
      dt <- data.table(iap2)
      dt <- dt[!is.na(dt$N.TRUE),]
      setkey(dt, "YEAR", "WOODID", "IMAGE")
      dt <- dt[, list(MNN.TRUE=RobustMin(N.TRUE), 
                      RPOS.MNN.TRUE=RRADDISTR[which.min(N.TRUE)], 
                      APOS.MNN.TRUE=mean(MRW,na.rm=TRUE)/100*RRADDISTR[which.min(N.TRUE)], 
                      TRUE.MNN.TRUE=round(RADDISTR.BAND[which.min(N.TRUE)]), 
                      MRW=mean(MRW,na.rm=TRUE)),
               by=key(dt)] 
      
      dt$APOS.MNN.TRUE <- ifelse(!is.na(dt$APOS.MNN.TRUE), dt$APOS.MNN.TRUE, (dt$MNN.TRUE+RESO[r]/2))   #if there minimum number of cells is 0, there is no RRADDISTR; assume center position of band instead
      dt$RPOS.MNN.TRUE <- ifelse(!is.na(dt$RPOS.MNN.TRUE), dt$RPOS.MNN.TRUE, (dt$APOS.MNN.TRUE*100/dt$MRW))   #if there minimum number of cells is 0, there is no RRADDISTR; assume center position of band instead
      dt$MRW <- NULL   #remove temporary column that was only used to troubleshoot NAs in RPOS.MNN.TRUE
      
      mnncells.true <- as.data.frame(dt)   
      mnncells.true <- merge(mnncells.true, ringTemplate, all=TRUE)   #to make sure any missing ring gets NA         
      
      
      ### 4.6.1.67 Summarize all annual statistics and save them to file ####
      df.ann <- Reduce(function(x, y) merge(x, y, all=TRUE), 
                       list(mrw, eww, lww, 
                            ncells.band, ncells.band.ew, ncells.band.lw, mxncells.band, mnncells.band, 
                            ncells.true, ncells.true.ew, ncells.true.lw, mxncells.true, mnncells.true, 
                            mla, mla.ew, mla.lw, mxla, mnla, 
                            mdrad, mdrad.ew, mdrad.lw, mxdrad, mndrad, 
                            mdtan, mdtan.ew, mdtan.lw, mxdtan, mndtan, 
                            mtca, mtca.ew, mtca.lw, mxtca, mntca,
                            mcwa, mcwa.ew, mcwa.lw, mxcwa, mncwa, 
                            cwaacc, cwaacc.ew, cwaacc.lw, mxcwaacc, mncwaacc,
                            mcwtall, mcwtall.ew, mcwtall.lw, mxcwtall, mncwtall,
                            mcwtalladj, mcwtalladj.ew, mcwtalladj.lw, mxcwtalladj, mncwtalladj, 
                            mcwttan, mcwttan.ew, mcwttan.lw, mxcwttan, mncwttan, 
                            mcwtrad, mcwtrad.ew, mcwtrad.lw, mxcwtrad, mncwtrad, 
                            mmork, mmork.ew, mmork.lw, mxmork, mnmork, 
                            mctsr, mctsr.ew, mctsr.lw, mxctsr, mnctsr, 
                            mdcwt, mdcwt.ew, mdcwt.lw, mxdcwt, mndcwt, 
                            mdcwt2, mdcwt2.ew, mdcwt2.lw, mxdcwt2, mndcwt2, 
                            mdcwa, mdcwa.ew, mdcwa.lw, mxdcwa, mndcwa, 
                            mtb2, mtb2.ew, mtb2.lw, mxtb2, mntb2, 
                            mkh, mkh.ew, mkh.lw, mxkh, mnkh, 
                            mdh, mdh.ew, mdh.lw, mxdh, mndh)) 
      
      df.ann <- data.frame(lapply(df.ann, function(x){replace(x, is.infinite(x), NA)}))   #replace error codes (negative values) by NA 
      df.ann <- data.frame(lapply(df.ann, function(x){replace(x, is.nan(x), NA)}))   #replace error codes (negative values) by NA 
      
      write.table(df.ann, file=paste(woodid, "_AnnualStats_", STATS[s], "_", RESO[r], "mu.txt", sep=""), row.names = FALSE)
  
  #     ### Just for fun...
  #     cor.test(mxdcwt$MXDCWT, mxdcwa$MXDCWA)
  #     cor.test(mndcwt$MNDCWT, mndcwa$MNDCWA)
  #     cor.test(mxdcwt$MXDCWT, mxcwtrad$MXCWTRAD)
  #     cor.test(mxcwtrad$MXCWTRAD, mxdcwa$MXDCWA)
  #     cor.test(mxcwtrad$MXCWTRAD, mxcwttan$MXCWTTAN)
      
      
      ### 4.6.1.68 Plot and save annual ring width (MRW) ####    
      png(file=paste(woodid, "_MRW_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mrw$MRW[mrw$YEAR>=YR.START & mrw$YEAR<=YR.END], 
           x=mrw$YEAR[mrw$YEAR>=YR.START & mrw$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Ring width (MRW) (µm)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #ylim=c(150,1500))
      par(xaxt="s") 
      axis(1, at=mrw$YEAR[mrw$YEAR>=YR.START & mrw$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mrw$YEAR[mrw$YEAR>=YR.START & mrw$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mrw$YEAR[mrw$YEAR>=YR.START & mrw$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.69 Plot and save annual earlywood width (EWW) ####    
      png(file=paste(woodid, "_EWW_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=eww$EWW[eww$YEAR>=YR.START & eww$YEAR<=YR.END], 
           x=eww$YEAR[eww$YEAR>=YR.START & eww$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Earlywood width (EWW) (µm)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #ylim=c(150,1500))
      par(xaxt="s") 
      axis(1, at=eww$YEAR[eww$YEAR>=YR.START & eww$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=eww$YEAR[eww$YEAR>=YR.START & eww$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=eww$YEAR[eww$YEAR>=YR.START & eww$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.70 Plot and save annual latewood width (LWW) ####    
      png(file=paste(woodid, "_LWW_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=lww$LWW[lww$YEAR>=YR.START & lww$YEAR<=YR.END], 
           x=lww$YEAR[lww$YEAR>=YR.START & lww$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Latewood width (LWW) (µm)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #ylim=c(150,1500))
      par(xaxt="s") 
      axis(1, at=lww$YEAR[lww$YEAR>=YR.START & lww$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=lww$YEAR[lww$YEAR>=YR.START & lww$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=lww$YEAR[lww$YEAR>=YR.START & lww$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
  
      ### 4.6.1.71 Plot and save annual maximum lumen area (MXLA) ####    
      png(file=paste(woodid, "_MXLA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxla$MXLA[mxla$YEAR>=YR.START & mxla$YEAR<=YR.END], 
           x=mxla$YEAR[mxla$YEAR>=YR.START & mxla$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. Lumen area (MXLA) (µm2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(150,1500))
      par(xaxt="s") 
      axis(1, at=mxla$YEAR[mxla$YEAR>=YR.START & mxla$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxla$YEAR[mxla$YEAR>=YR.START & mxla$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxla$YEAR[mxla$YEAR>=YR.START & mxla$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.72 Plot and save annual minimum lumen area (MNLA) ####    
      png(file=paste(woodid, "_MNLA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mnla$MNLA[mnla$YEAR>=YR.START & mnla$YEAR<=YR.END], 
           x=mnla$YEAR[mnla$YEAR>=YR.START & mnla$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. Lumen area (MNLA) (µm2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(5,700))
      par(xaxt="s") 
      axis(1, at=mnla$YEAR[mnla$YEAR>=YR.START & mnla$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mnla$YEAR[mnla$YEAR>=YR.START & mnla$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mnla$YEAR[mnla$YEAR>=YR.START & mnla$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.73 Plot and save annual maximum cell area (MXTCA) ####    
      png(file=paste(woodid, "_MXTCA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxtca$MXTCA[mxtca$YEAR>=YR.START & mxtca$YEAR<=YR.END], 
           x=mxtca$YEAR[mxtca$YEAR>=YR.START & mxtca$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. Cell area (MXTCA) (µm2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(200,5000))
      par(xaxt="s") 
      axis(1, at=mxtca$YEAR[mxtca$YEAR>=YR.START & mxtca$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxtca$YEAR[mxtca$YEAR>=YR.START & mxtca$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxtca$YEAR[mxtca$YEAR>=YR.START & mxtca$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.74 Plot and save annual minimum cell area (MNTCA) ####    
      png(file=paste(woodid, "_MNTCA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mntca$MNTCA[mntca$YEAR>=YR.START & mntca$YEAR<=YR.END], 
           x=mntca$YEAR[mntca$YEAR>=YR.START & mntca$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. Cell area (MNTCA) (µm2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(50,1000))
      par(xaxt="s") 
      axis(1, at=mntca$YEAR[mntca$YEAR>=YR.START & mntca$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mntca$YEAR[mntca$YEAR>=YR.START & mntca$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mntca$YEAR[mntca$YEAR>=YR.START & mntca$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()    
  
      
      ### 4.6.1.75 Plot and save annual maximum accumulated cell wall area (MXCWAACC) ####    
      png(file=paste(woodid, "_MXCWAACC_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxcwa$MXCWA[mxcwa$YEAR>=YR.START & mxcwa$YEAR<=YR.END], 
           x=mxcwa$YEAR[mxcwa$YEAR>=YR.START & mxcwa$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. Acc. Cell wall area (MXCWAACC) (µm2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           # cex.axis=1.3, cex.lab=2, ylim=c(100,4000))
           cex.axis=1.3, cex.lab=2)
      par(xaxt="s") 
      axis(1, at=mxcwa$YEAR[mxcwa$YEAR>=YR.START & mxcwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxcwa$YEAR[mxcwa$YEAR>=YR.START & mxcwa$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxcwa$YEAR[mxcwa$YEAR>=YR.START & mxcwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.76 Plot and save annual minimum accumulated cell wall area (MNCWAACC) ####    
      png(file=paste(woodid, "_MNCWAACC_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mncwa$MNCWA[mncwa$YEAR>=YR.START & mncwa$YEAR<=YR.END], 
           x=mncwa$YEAR[mncwa$YEAR>=YR.START & mncwa$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. Acc. Cell wall area (MNCWAACC) (µm2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           # cex.axis=1.3, cex.lab=2, ylim=c(20,500))
           cex.axis=1.3, cex.lab=2)
      par(xaxt="s") 
      axis(1, at=mncwa$YEAR[mncwa$YEAR>=YR.START & mncwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mncwa$YEAR[mncwa$YEAR>=YR.START & mncwa$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mncwa$YEAR[mncwa$YEAR>=YR.START & mncwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()      
      
      
      ### 4.6.1.77 Plot and save annual maximum cell wall area (MXCWA) ####    
      png(file=paste(woodid, "_MXCWA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxcwa$MXCWA[mxcwa$YEAR>=YR.START & mxcwa$YEAR<=YR.END], 
           x=mxcwa$YEAR[mxcwa$YEAR>=YR.START & mxcwa$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. Cell wall area (MXCWA) (µm2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(100,4000))
      par(xaxt="s") 
      axis(1, at=mxcwa$YEAR[mxcwa$YEAR>=YR.START & mxcwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxcwa$YEAR[mxcwa$YEAR>=YR.START & mxcwa$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxcwa$YEAR[mxcwa$YEAR>=YR.START & mxcwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.78 Plot and save annual minimum cell wall area (MNCWA) ####    
      png(file=paste(woodid, "_MNCWA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mncwa$MNCWA[mncwa$YEAR>=YR.START & mncwa$YEAR<=YR.END], 
           x=mncwa$YEAR[mncwa$YEAR>=YR.START & mncwa$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. Cell wall area (MNCWA) (µm2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(20,500))
      par(xaxt="s") 
      axis(1, at=mncwa$YEAR[mncwa$YEAR>=YR.START & mncwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mncwa$YEAR[mncwa$YEAR>=YR.START & mncwa$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mncwa$YEAR[mncwa$YEAR>=YR.START & mncwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()    
      
      
      ### 4.6.1.79 Plot and save annual maximum tangential cell wall thickness (MXCWTTAN) ####    
      png(file=paste(woodid, "_MXCWTTAN_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxcwttan$MXCWTTAN[mxcwttan$YEAR>=YR.START & mxcwttan$YEAR<=YR.END], 
           x=mxcwttan$YEAR[mxcwttan$YEAR>=YR.START & mxcwttan$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. tangential cell wall thickness (MXCWTTAN) (µm)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(1.5,10))
      par(xaxt="s") 
      axis(1, at=mxcwttan$YEAR[mxcwttan$YEAR>=YR.START & mxcwttan$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxcwttan$YEAR[mxcwttan$YEAR>=YR.START & mxcwttan$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxcwttan$YEAR[mxcwttan$YEAR>=YR.START & mxcwttan$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.80 Plot and save annual minimum tangential cell wall thickness (MNCWTTAN) ####    
      png(file=paste(woodid, "_MNCWTTAN_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mncwttan$MNCWTTAN[mncwttan$YEAR>=YR.START & mncwttan$YEAR<=YR.END], 
           x=mncwttan$YEAR[mncwttan$YEAR>=YR.START & mncwttan$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. tangential cell wall thickness (MNCWTTAN) (µm)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(1,6))
      par(xaxt="s") 
      axis(1, at=mncwttan$YEAR[mncwttan$YEAR>=YR.START & mncwttan$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mncwttan$YEAR[mncwttan$YEAR>=YR.START & mncwttan$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mncwttan$YEAR[mncwttan$YEAR>=YR.START & mncwttan$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()      
      
      
      ### 4.6.1.81 Plot and save annual maximum radial cell wall thickness (MXCWTRAD) ####    
      png(file=paste(woodid, "_MXCWTRAD_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxcwtrad$MXCWTRAD[mxcwtrad$YEAR>=YR.START & mxcwtrad$YEAR<=YR.END], 
           x=mxcwtrad$YEAR[mxcwtrad$YEAR>=YR.START & mxcwtrad$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. radial cell wall thickness (MXCWTRAD) (µm)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(1.5,10))
      par(xaxt="s") 
      axis(1, at=mxcwtrad$YEAR[mxcwtrad$YEAR>=YR.START & mxcwtrad$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxcwtrad$YEAR[mxcwtrad$YEAR>=YR.START & mxcwtrad$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxcwtrad$YEAR[mxcwtrad$YEAR>=YR.START & mxcwtrad$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.82 Plot and save annual minimum radial cell wall thickness (MNCWTRAD) ####    
      png(file=paste(woodid, "_MNCWTRAD_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mncwtrad$MNCWTRAD[mncwtrad$YEAR>=YR.START & mncwtrad$YEAR<=YR.END], 
           x=mncwtrad$YEAR[mncwtrad$YEAR>=YR.START & mncwtrad$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. radial cell wall thickness (MNCWTRAD) (µm)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(1,6))
      par(xaxt="s") 
      axis(1, at=mncwtrad$YEAR[mncwtrad$YEAR>=YR.START & mncwtrad$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mncwtrad$YEAR[mncwtrad$YEAR>=YR.START & mncwtrad$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mncwtrad$YEAR[mncwtrad$YEAR>=YR.START & mncwtrad$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()        
      
      
      ### 4.6.1.83 Plot and save annual maximum Mork's index (MXRTSR) ####    
      png(file=paste(woodid, "_MXRTSR_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxmork$MXRTSR[mxmork$YEAR>=YR.START & mxmork$YEAR<=YR.END], 
           x=mxmork$YEAR[mxmork$YEAR>=YR.START & mxmork$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. Mork's index (MXRTSR)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(1.5,10))
      par(xaxt="s") 
      axis(1, at=mxmork$YEAR[mxmork$YEAR>=YR.START & mxmork$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxmork$YEAR[mxmork$YEAR>=YR.START & mxmork$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxmork$YEAR[mxmork$YEAR>=YR.START & mxmork$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.84 Plot and save annual minimum Mork's index (MNRTSR) ####    
      png(file=paste(woodid, "_MNRTSR_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mnmork$MNRTSR[mnmork$YEAR>=YR.START & mnmork$YEAR<=YR.END], 
           x=mnmork$YEAR[mnmork$YEAR>=YR.START & mnmork$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. Mork's index (MNRTSR)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(1,6))
      par(xaxt="s") 
      axis(1, at=mnmork$YEAR[mnmork$YEAR>=YR.START & mnmork$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mnmork$YEAR[mnmork$YEAR>=YR.START & mnmork$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mnmork$YEAR[mnmork$YEAR>=YR.START & mnmork$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()        
      
      
      ### 4.6.1.85 Plot and save annual maximum relative anatomical density based on CWT (DCWT) ####    
      png(file=paste(woodid, "_MXDCWT_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxdcwt$MXDCWT[mxdcwt$YEAR>=YR.START & mxdcwt$YEAR<=YR.END], 
           x=mxdcwt$YEAR[mxdcwt$YEAR>=YR.START & mxdcwt$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. Anatomical density (MXDCWT)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.3,1))
      par(xaxt="s") 
      axis(1, at=mxdcwt$YEAR[mxdcwt$YEAR>=YR.START & mxdcwt$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxdcwt$YEAR[mxdcwt$YEAR>=YR.START & mxdcwt$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxdcwt$YEAR[mxdcwt$YEAR>=YR.START & mxdcwt$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.86 Plot and save annual minimum relative anatomical density based on CWT (DCWT) ####    
      png(file=paste(woodid, "_MNDCWT_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mndcwt$MNDCWT[mndcwt$YEAR>=YR.START & mndcwt$YEAR<=YR.END], 
           x=mndcwt$YEAR[mndcwt$YEAR>=YR.START & mndcwt$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. Anatomical density (MNDCWT)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.1,0.7))
      par(xaxt="s") 
      axis(1, at=mndcwt$YEAR[mndcwt$YEAR>=YR.START & mndcwt$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mndcwt$YEAR[mndcwt$YEAR>=YR.START & mndcwt$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mndcwt$YEAR[mndcwt$YEAR>=YR.START & mndcwt$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()  
      
      
      ### 4.6.1.87 Plot and save annual maximum special anatomical density (DCWT2): CWTRAD/DRAD ####
      png(file=paste(woodid, "_MXDCWT2_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxdcwt2$MXDCWT2[mxdcwt2$YEAR>=YR.START & mxdcwt2$YEAR<=YR.END],
           x=mxdcwt2$YEAR[mxdcwt2$YEAR>=YR.START & mxdcwt2$YEAR<=YR.END],
           xaxt="n",  yaxt="n", xlab="",ylab="Max. Anatomical density (MXDCWT2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r],
           # cex.axis=1.3, cex.lab=2, ylim=c(0.3,1))
           cex.axis=1.3, cex.lab=2)
      par(xaxt="s")
      axis(1, at=mxdcwt2$YEAR[mxdcwt2$YEAR>=YR.START & mxdcwt2$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3,
           labels=mxdcwt2$YEAR[mxdcwt2$YEAR>=YR.START & mxdcwt2$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxdcwt2$YEAR[mxdcwt2$YEAR>=YR.START & mxdcwt2$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()


      ### 4.6.1.88 Plot and save annual minimum special anatomical density (DCWT2): CWTRAD/DRAD ####
      png(file=paste(woodid, "_MNDCWT2_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mndcwt2$MNDCWT2[mndcwt2$YEAR>=YR.START & mndcwt2$YEAR<=YR.END],
           x=mndcwt2$YEAR[mndcwt2$YEAR>=YR.START & mndcwt2$YEAR<=YR.END],
           xaxt="n",  yaxt="n", xlab="",ylab="Min. Anatomical density (MNDCWT2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r],
           # cex.axis=1.3, cex.lab=2, ylim=c(0.1,0.7))
           cex.axis=1.3, cex.lab=2)
      par(xaxt="s")
      axis(1, at=mndcwt2$YEAR[mndcwt2$YEAR>=YR.START & mndcwt2$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3,
           labels=mndcwt2$YEAR[mndcwt2$YEAR>=YR.START & mndcwt2$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mndcwt2$YEAR[mndcwt2$YEAR>=YR.START & mndcwt2$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.89 Plot and save annual maximum relative anatomical density based on CWA (DCWA) ####    
      png(file=paste(woodid, "_MXDCWA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxdcwa$MXDCWA[mxdcwa$YEAR>=YR.START & mxdcwa$YEAR<=YR.END], 
           x=mxdcwa$YEAR[mxdcwa$YEAR>=YR.START & mxdcwa$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. Anatomical density (MXDCWA)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.3,1))
      par(xaxt="s") 
      axis(1, at=mxdcwa$YEAR[mxdcwa$YEAR>=YR.START & mxdcwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxdcwa$YEAR[mxdcwa$YEAR>=YR.START & mxdcwa$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxdcwa$YEAR[mxdcwa$YEAR>=YR.START & mxdcwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.90 Plot and save annual minimum relative anatomical density based on CWA (DCWA) ####    
      png(file=paste(woodid, "_MNDCWA_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mndcwa$MNDCWA[mndcwa$YEAR>=YR.START & mndcwa$YEAR<=YR.END], 
           x=mndcwa$YEAR[mndcwa$YEAR>=YR.START & mndcwa$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. Anatomical density (MNDCWA)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.1,0.7))
      par(xaxt="s") 
      axis(1, at=mndcwa$YEAR[mndcwa$YEAR>=YR.START & mndcwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mndcwa$YEAR[mndcwa$YEAR>=YR.START & mndcwa$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mndcwa$YEAR[mndcwa$YEAR>=YR.START & mndcwa$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()  
      
      
      ### 4.6.1.91 Plot and save annual maximum cell wall reinforcement index (TB2) ####    
      png(file=paste(woodid, "_MXTB2_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxtb2$MXTB2[mxtb2$YEAR>=YR.START & mxtb2$YEAR<=YR.END], 
           x=mxtb2$YEAR[mxtb2$YEAR>=YR.START & mxtb2$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. cell wall reinforcement index (MXTB2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.3,1))
      par(xaxt="s") 
      axis(1, at=mxtb2$YEAR[mxtb2$YEAR>=YR.START & mxtb2$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxtb2$YEAR[mxtb2$YEAR>=YR.START & mxtb2$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxtb2$YEAR[mxtb2$YEAR>=YR.START & mxtb2$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.92 Plot and save annual minimum cell wall reinforcement index (TB2) ####    
      png(file=paste(woodid, "_MNTB2_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mntb2$MNTB2[mntb2$YEAR>=YR.START & mntb2$YEAR<=YR.END], 
           x=mntb2$YEAR[mntb2$YEAR>=YR.START & mntb2$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. cell wall reinforcement index (MNTB2)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.1,0.7))
      par(xaxt="s") 
      axis(1, at=mntb2$YEAR[mntb2$YEAR>=YR.START & mntb2$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mntb2$YEAR[mntb2$YEAR>=YR.START & mntb2$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mntb2$YEAR[mntb2$YEAR>=YR.START & mntb2$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()   
      
      
      ### 4.6.1.93 Plot and save annual maximum mean hydraulic diameter (DH) ####    
      png(file=paste(woodid, "_MXDH_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxdh$MXDH[mxdh$YEAR>=YR.START & mxdh$YEAR<=YR.END], 
           x=mxdh$YEAR[mxdh$YEAR>=YR.START & mxdh$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. mean hydraulic diameter (MXDH) (µm)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.3,1))
      par(xaxt="s") 
      axis(1, at=mxdh$YEAR[mxdh$YEAR>=YR.START & mxdh$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxdh$YEAR[mxdh$YEAR>=YR.START & mxdh$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxdh$YEAR[mxdh$YEAR>=YR.START & mxdh$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.94 Plot and save annual minimum mean hydraulic diameter (DH) ####    
      png(file=paste(woodid, "_MNDH_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mndh$MNDH[mndh$YEAR>=YR.START & mndh$YEAR<=YR.END], 
           x=mndh$YEAR[mndh$YEAR>=YR.START & mndh$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. mean hydraulic diameter (MNDH) (µm)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.1,0.7))
      par(xaxt="s") 
      axis(1, at=mndh$YEAR[mndh$YEAR>=YR.START & mndh$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mndh$YEAR[mndh$YEAR>=YR.START & mndh$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mndh$YEAR[mndh$YEAR>=YR.START & mndh$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off() 
      
      
      ### 4.6.1.95 Plot and save annual maximum number of cells (N.BAND) ####    
      png(file=paste(woodid, "_MXN.BAND_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxncells.band$MXN.BAND[mxncells.band$YEAR>=YR.START & mxncells.band$YEAR<=YR.END], 
           x=mxncells.band$YEAR[mxncells.band$YEAR>=YR.START & mxncells.band$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. number of cells (MXN.BAND)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.3,1))
      par(xaxt="s") 
      axis(1, at=mxncells.band$YEAR[mxncells.band$YEAR>=YR.START & mxncells.band$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxncells.band$YEAR[mxncells.band$YEAR>=YR.START & mxncells.band$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxncells.band$YEAR[mxncells.band$YEAR>=YR.START & mxncells.band$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.96 Plot and save annual minimum number of cells (N.BAND) ####    
      png(file=paste(woodid, "_MNN.BAND_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mnncells.band$MNN.BAND[mnncells.band$YEAR>=YR.START & mnncells.band$YEAR<=YR.END], 
           x=mnncells.band$YEAR[mnncells.band$YEAR>=YR.START & mnncells.band$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. number of cells (MNN.BAND)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.1,0.7))
      par(xaxt="s") 
      axis(1, at=mnncells.band$YEAR[mnncells.band$YEAR>=YR.START & mnncells.band$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mnncells.band$YEAR[mnncells.band$YEAR>=YR.START & mnncells.band$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mnncells.band$YEAR[mnncells.band$YEAR>=YR.START & mnncells.band$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()    
      
      
      ### 4.6.1.97 Plot and save annual maximum number of cells (N.TRUE) ####    
      png(file=paste(woodid, "_MXN.TRUE_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mxncells.true$MXN.TRUE[mxncells.true$YEAR>=YR.START & mxncells.true$YEAR<=YR.END], 
           x=mxncells.true$YEAR[mxncells.true$YEAR>=YR.START & mxncells.true$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Max. number of cells (MXN.TRUE)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.3,1))
      par(xaxt="s") 
      axis(1, at=mxncells.true$YEAR[mxncells.true$YEAR>=YR.START & mxncells.true$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mxncells.true$YEAR[mxncells.true$YEAR>=YR.START & mxncells.true$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mxncells.true$YEAR[mxncells.true$YEAR>=YR.START & mxncells.true$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()
      
      
      ### 4.6.1.98 Plot and save annual minimum number of cells (N.TRUE) ####    
      png(file=paste(woodid, "_MNN.TRUE_", STATS[s], "_", RESO[r], "mu_", SMOOTHER, "rm", ".png", sep=""),width=3000,height=600)
      par(mar=c(3,5.5,2.1,2.1))   #c(bottom, left, top, right)
      plot(y=mnncells.true$MNN.TRUE[mnncells.true$YEAR>=YR.START & mnncells.true$YEAR<=YR.END], 
           x=mnncells.true$YEAR[mnncells.true$YEAR>=YR.START & mnncells.true$YEAR<=YR.END], 
           xaxt="n",  yaxt="n", xlab="",ylab="Min. number of cells (MNN.TRUE)", type="l", lwd=3, col=LINE.COL(length(RESO))[r], 
           cex.axis=1.3, cex.lab=2) #, ylim=c(0.1,0.7))
      par(xaxt="s") 
      axis(1, at=mnncells.true$YEAR[mnncells.true$YEAR>=YR.START & mnncells.true$YEAR<=YR.END], tck=0.02, lwd.ticks=2, cex.axis=1.3, 
           labels=mnncells.true$YEAR[mnncells.true$YEAR>=YR.START & mnncells.true$YEAR<=YR.END])
      axis(2, tck=0.02, lwd.ticks=2, cex.axis=1.3, labels=TRUE)
      axis(3, at=mnncells.true$YEAR[mnncells.true$YEAR>=YR.START & mnncells.true$YEAR<=YR.END], tck=0.02, lwd.ticks=2, labels=FALSE)
      axis(4, tck=0.02, lwd.ticks=2, labels=FALSE)
      dev.off()        
    }
  }
} 

Sys.time() - t


### 5. clean up ####
rm(list=ls(all=TRUE))
