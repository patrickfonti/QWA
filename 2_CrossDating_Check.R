#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
# This script check the crossdating of roxas tree-ring width time-series at tree level
# The comparison is: 
#   1) between roxas and traditional ring width measurement, tree by tree
#   2) within the Roxas tree-ring time-series 
# Input files are outcoming from script 1_SummarizeRadialWedgeData (Only Ring files), 
#
# IMPORTANT: 
# Labels of Roxas trees is according to this MODEL SITE(3)_TREE(4)_SLIDE(2)_IMAGE(2): example YAM_107_1_1
# Labels of Reference RWL file is according to this MODEL SITE(1)_YEAR(5): example 
# Y_103    -20    25    65    38    75    78    58    73    23    30    20
# Y_103    -10    10     3    13    55    58    48    43    35    20    25
# Y_103      0    33    28    25    48    53    45    13    50    30    50
# Y_103     10    43    15    18    23    48    50    18    93    18    38
# Y_103     20    43    45    23    33     5    33    15     8    20     3
# Y_103     30    38    28    38    43    63    98    23    75    53    80
# Y_103     40    70    43    65    60    63    18     5    33    35     5
# Y_103     50     8    23     5    33    35    35    28    35    53    10
# Y_103     60    63    33    18    43    50    53    35    60    35    88
# Y_103     70   999
#
# v1.0, 21 March 2020
# 
# (c) Patrick Fonti
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#

t <- Sys.time()

### 0. Functions ###
Extract_names <- function(x) {
  Names <-matrix("NA",length(x),6) 
  colnames(Names) <- c("plot","tree","slide","image","treecode","imagecode")
  for (i in 1:length(x)) {
    ipos <-  unlist(lapply(strsplit(x[i], ''), function(a) which(a == '_')))  # get position of all "_"
    plot <- substr(x[i],1,ipos[1]-1)
    tree <- substr(x[i],ipos[1]+1,ipos[2]-1)
    slide <- substr(x[i],ipos[2]+1,ipos[3]-1)
    image <- substr(x[i],ipos[3]+1,nchar(x[i])) 
    treecode <- substr(x[i],1,ipos[2]-1)   
    imagecode <- x[i]
    Names[i,] <-c(plot,tree,slide,image,treecode,imagecode)
  } 
  return(as.data.frame(Names))
}


### 1. install required packages ####
list.of.packages <- c("dplR","tidyr", "plyr", "zoo")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
library(dplR)
library(tidyr)
library(plyr)
library(zoo)


### 2. Up-load reference rwl file ####
# tucson file with long = TRUE, 
# i.e.; only the first 7 characters can be used for series IDs. 
# This allow to have negative years (before CE).
# Tucson file precision should be with 999 (for RW unit like 112) or -9999 (for RW unit like 1120) in TRW

RWL <- "Reference.rwl" # Provide here the path and name of your reference file!!!!
REFERENCE.RW <- read.tucson(RWL, long=TRUE) 
# plot(REFERENCE.RW, plot.type="spag", zfac=1)
# write.tucson(REFERENCE.RW, "new.Reference.rwl", header = NULL, append = FALSE, long.names = TRUE,mapping.fname = "new.IDs")
# write.csv(REFERENCE.RW[,], "test.txt", row.names = TRUE)
# A<-read.rwl("test.txt", format = "csv")
# write.tucson(A, "newALL.txt", header = NULL, append = FALSE, long.names = TRUE)


### 3. Define Inoput and output directories ####
topfolder <- getwd()
WD_Input <- paste0(topfolder,"/Roxas_tree_series")                     # directory with the full series of ring data for the entire wood pieces
WD_Output <- paste0(topfolder,"/Roxas_tree_series/Crossdate_check")    # directory where to insert the output of the cross-dating analyses
dir.create(WD_Output, showWarnings = TRUE)


### 4. Get file names of ring output files and write them into a dataframe in long and wide format ####
rf <- list.files(path=WD_Input, pattern="Output_Rings.txt", full.names=TRUE, recursive=TRUE, ignore.case=TRUE, include.dirs=TRUE, no.. = FALSE)

## Long
ROXAS.RW.long <- NULL
for (i in 1:length(rf)) {
df.ring <- read.table(rf[i], header=TRUE, sep=" ", na.strings=c("NA",""))   #ring data
df.ring <- df.ring[,c("ID","YEAR","MRW")]
ROXAS.RW.long <- rbind(ROXAS.RW.long,df.ring)
}
## Label corrections for YAMAL file
ROXAS.RW.long$ID <- gsub("Yamal_new","YAM",ROXAS.RW.long$ID)
ROXAS.RW.long$ID <- gsub("TP","P",ROXAS.RW.long$ID)

ROXAS.RW.long <- cbind(Extract_names(as.character(ROXAS.RW.long$ID))[,"treecode"],ROXAS.RW.long)
colnames(ROXAS.RW.long)[1]<-"Treecode"


End_slide <- ddply(ROXAS.RW.long, .(Treecode, ID), summarize, 
                   end_slide = max(YEAR),
                   x_position = mean(YEAR))

## Wide
ROXAS.RW.wide <- spread(ROXAS.RW.long[,c("Treecode","YEAR","MRW")], key = Treecode, value = MRW)


### 5. Check consistency of TRW files ####
colnames(REFERENCE.RW) <- gsub("Y","YAM",colnames(REFERENCE.RW))
REFERENCE.RW$YEAR<-as.numeric(as.character(rownames(REFERENCE.RW)))
Match<-colnames(ROXAS.RW.wide[,-1]) %in% colnames(REFERENCE.RW)
for (i in 1:length(Match)) {
  if(Match[i]==FALSE) (print(paste0(colnames(ROXAS.RW.wide[,-1][i]), " does not have a match in the reference file!")))
}

### 6. Export ROXAS as tucson format ###
ROXAS.rwl <- merge(ROXAS.RW.wide, REFERENCE.RW, by="YEAR")[,1:ncol(ROXAS.RW.wide)] # to ensure that YEAR is a continous timeseries
colnames(ROXAS.rwl) <- gsub(".x","",colnames(ROXAS.rwl))
rownames(ROXAS.rwl) <- as.numeric((ROXAS.rwl[,"YEAR"]))
as.rwl(ROXAS.rwl)
ROXAS.rwl <- ROXAS.rwl[,-1]/1000
plot(as.rwl(ROXAS.rwl), plot.type="spag", zfac=1)
write.tucson(ROXAS.rwl, paste0(WD_Output,"/Roxas.rwl"), header = NULL, append = FALSE, long.names = TRUE,mapping.fname = "ROXAS.RWL_new.IDs")


### 7. Compare raw series by raw series and export in a PDF file ####
pdf(paste0(WD_Output,"/1.ROXAS_REFERENCE_CHECK.pdf"))
plot(as.rwl(ROXAS.rwl), plot.type="spag", zfac=1, col="red")
require(zoo)
Detrend<-function (x) {
  detrend(as.data.frame(cbind(1:length(x),na.approx(x,na.rm=FALSE))), method="Spline",nyrs=30)
}
ROXAS.detChrono <- chron(Detrend(as.rwl(ROXAS.rwl)), prefix="YAM")
ROXAS.detChrono$YEAR <- as.numeric(rownames(ROXAS.rwl))
plot(ROXAS.detChrono$YEAR,ROXAS.detChrono$samp.depth, type="l", lwd=2, col="red", ylab = "Sample depth", xlab = "Calendar Year", bty="n")
abline(h=10, col="grey", lty=2)

# plot(as.rwl(ROXAS.rwl), plot.type="seg")
Crossdate.df<-NULL
for (i in 1:length(colnames(ROXAS.RW.wide[,-1]))) {
  Series.code <- colnames(ROXAS.RW.wide[,-1])[i]
  print(Series.code)
  MERGED.serie <- merge(ROXAS.RW.wide[,c("YEAR",Series.code)], REFERENCE.RW[,c("YEAR",Series.code)], by="YEAR" )
  xmin<-MERGED.serie$YEAR[min(which(!is.na(MERGED.serie[,2])),na.rm=TRUE)]
  xmax<-MERGED.serie$YEAR[max(which(!is.na(MERGED.serie[,2])),na.rm=TRUE)]
  ymax<-max(MERGED.serie[,2]/1000,na.rm = TRUE)
  COR <- round(cor(MERGED.serie[,2:3], use = "pairwise.complete.obs" )[1,2],2)
  plot(MERGED.serie$YEAR,MERGED.serie[,2]/1000, type="l", lwd=2, bty="n", xlim=c(xmin, xmax),
       main=Series.code, col="red", ylab = "Ring width [mm]", xlab = "Calendar Year", sub = paste("r = ",COR))
  lines(MERGED.serie$YEAR,MERGED.serie[,3])
  legend("right", c("Roxas", "Reference"), col=c("red", "black"), lwd=2, bty="n")
  Crossdate.df<-rbind(Crossdate.df,c(Series.code,COR, xmin, xmax, round(mean(MERGED.serie[,2]/1000, na.rm = TRUE),2)))
  # Add info on images 
  Image_label <- End_slide[which(End_slide$Treecode==Series.code),]
  abline(v=Image_label$end_slide, col="grey")
  label <- paste(Extract_names(as.character(Image_label$ID))[,c("slide")],Extract_names(as.character(Image_label$ID))[,c("image")], sep="_")
  text(Image_label$x_position,ymax,label, cex=0.5, col="grey")
}
dev.off()

colnames(Crossdate.df) <- c("Series", "r", "Start", "End", "Average RW")
Crossdate.df<-as.data.frame(Crossdate.df)
write.table(Crossdate.df,paste0(WD_Output,"/1.ROXAS_REFERENCE_table.txt"),row.names = FALSE, sep="\t")


### 8. Compare Roxas series ####
pdf(paste0(WD_Output,"/2.ROXAS_CROSSDATE_CHECK.pdf"))
rownames(ROXAS.rwl)<-as.numeric(rownames(ROXAS.rwl))
for (i in 1:length(colnames(ROXAS.rwl))) {
  tryCatch({ 
    ccf.series.rwl(ROXAS.rwl, colnames(ROXAS.rwl)[i], seg.length = 20, bin.floor = 20, main = colnames(ROXAS.rwl)[i])
  }, error = function(e) {print(paste0("Could not make a plot for ",colnames(ROXAS.rwl)[i]))
    }
  )
  
#series.rwl.plot(ROXAS.rwl, colnames(ROXAS.rwl)[i], seg.length = 20, bin.floor = 20)
}
dev.off()

Sys.time() - t

