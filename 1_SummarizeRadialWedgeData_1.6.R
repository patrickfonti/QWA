#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
# This script creates summary cell and ring data files for each radial wood piece. 
# For the overlapping years, data is taken from the input file with more cells. 
# Input files from the same wood piece do not need to be in a chronological order, 
# but output is in ascending year and relative intra-annual position (cells) order.
#
#
# v1.6, 23 August 2019
# 
# (c) Georg von Arx
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#


### 1. Clean up ####
rm(list=ls()) # clean desk 


### 2. Define top directory containing all the data in whatever hierarchical structure ####
WD_Input <- "Roxas_files"   #directory with the cell and ring data from individual images
WD_output <- "Roxas_tree_series"   #directory with the full series of cell and ring data for the entire wood pieces

setwd(WD_Input)
topfolder <- getwd()


### 3. Get file names of cell and ring output files and write them into a metadata table ####
### cells
cf <- list.files(path=topfolder, pattern="Output_Cells.txt", full.names=TRUE, recursive=TRUE, ignore.case=TRUE, include.dirs=TRUE, no.. = FALSE)
len <- length(cf)
metadata <- data.frame(tree=rep(NA,len), slide=rep(NA,len), image=rep(NA,len), imagecode=rep(NA,len), cell.file=rep(NA,len), ring.file=rep(NA,len))
metadata$cell.file <- as.matrix(cf)

### rings
rf <- list.files(path=topfolder, pattern="Output_Rings.txt", full.names=TRUE, recursive=TRUE, ignore.case=TRUE, include.dirs=TRUE, no.. = FALSE)
metadata$ring.file <- as.matrix(rf)


### 4. Check consistency of cell and ring files within each row ####
for (i in c(1:len))
{
  ### Extract full image code of cell file (without path and any other text)
  ipos <- lapply(strsplit(cf[i], ''), function(x) which(x=='/'))   #get position of all "/"
  ipos <- unlist(ipos)[(length(unlist(ipos)))]   #get position of last "/"
  cf.code <- substr(cf[i],ipos+1, nchar(cf[i]))   #extract string after last "/"
  cf.code <- sub("_Output_Cells.txt", "", cf.code)   #extract full image code   
  
  ### Extract full image code of ring file (without path and any other text)
  ipos <- lapply(strsplit(rf[i], ''), function(x) which(x=='/'))   #get position of all "/"
  ipos <- unlist(ipos)[(length(unlist(ipos)))]   #get position of last "/"
  rf.code <- substr(rf[i],ipos+1, nchar(rf[i]))   #extract string after last "/"
  rf.code <- sub("_Output_Rings.txt", "", rf.code)   #extract full image code   
  
  if (cf.code!=rf.code)
  {
    stop(paste("cell and ring files are not consistent at line ", i, "!", sep=""))
  }  
}  


### 5. Extract plot, tree, slide and image information from the file name ####
for (i in c(1:len))
{
  ### Extract full image code (without path and any other text)
  ipos <- lapply(strsplit(cf[i], ''), function(x) which(x=='/'))   #get position of all "/"
  ipos <- unlist(ipos)[(length(unlist(ipos)))]   #get position of last "/"
  txt <- substr(cf[i],ipos+1, nchar(cf[i]))   #extract string after last "/"
  txt <- sub("_Output_Cells.txt", "", txt)   #extract full image code  
  
  ### Extract information from image code and write it to metadata dataframe
  ipos <- lapply(strsplit(txt, ''), function(x) which(x=='_'))   #get position of all "_"
  ipos2 <- unlist(ipos)    #get position of last "_"
  metadata$plot[i] <- substr(txt,1,ipos2[1]-1)
  metadata$tree[i] <- substr(txt,ipos2[1]+1,ipos2[2]-1)
  metadata$slide[i] <- substr(txt,ipos2[2]+1,ipos2[3]-1)
  metadata$image[i] <- substr(txt,ipos2[3]+1,nchar(txt)) 
  metadata$treecode[i] <- substr(txt,1,ipos2[2]-1)   
  metadata$imagecode[i] <- txt
}   

metadata$treecode <- as.factor(metadata$treecode)


### 6. Remove any summary files from previous runs ####
for (i in c(nrow(metadata):1))
{
  if (is.na(metadata$treecode[i]))
  {
    metadata <- metadata[-i,]
  }  
}
  

### 7. Summarize the data ####
setwd(WD_output)

t <- Sys.time()

for (i in c(1:length(levels(metadata$treecode))))
{
  # print(levels(metadata$treecode)[i])   #loop control
  print(paste("(", i, "/", length(levels(metadata$treecode)), ") - ", "Processing wood piece: ", metadata$treecode[i], sep=""))
  
  ### 7.1 Get subset data for target wedge
  mdf <- metadata[metadata$treecode==levels(metadata$treecode)[i],]
  
  ### 7.2 Get number of output files for target wedge
  numfiles <- nrow(mdf)
  
  ### 7.3 Loop through each output file and extract information
  for (f in (1:numfiles))
  {  
    # print(mdf$cell.file[f])
    print(paste("   (", make.unique(rep(LETTERS, length.out=100), sep='')[f], ") ", mdf$cell.file[f], ": ", sep=""))  
    if (f==1)   #on first run, initialize cells and rings dataframes
    {
      cells <- NULL
      cells <- read.table(mdf$cell.file[f], header=TRUE, sep="\t", na.strings=c("NA",""))   #cell data
      if("BEND" %in% colnames(cells)){colnames(cells)[colnames(cells)=="BEND"] <- "TB2"}   #for backwards compatibility
      if("CRI" %in% colnames(cells)){colnames(cells)[colnames(cells)=="CRI"] <- "TB2"}   #for backwards compatibility
      if("CA" %in% colnames(cells)){colnames(cells)[colnames(cells)=="CA"] <- "LA"}   #for backwards compatibility
      cells <- cells[, c("ID", "YEAR", "CID", "RADDIST", "RADDISTR", "RRADDISTR", "LA",  
                         "ASP", "MAJAX", "KH", "CWTPI", "CWTBA", "CWTLE", "CWTRI", "CWTTAN", 
                         "CWTRAD", "CWTALL", "RTSR", "CTSR", "DRAD", "DTAN", "TB2",   
                         "CWA", "RWD")]
      
      rings <- NULL
      rings <- read.table(mdf$ring.file[f], header=TRUE, sep="\t", na.strings=c("NA",""))   #ring data
      rings <- rings[,c("ID","YEAR","MRW")]
      if (rings$YEAR[1] %in% levels(as.factor(cells$YEAR))==FALSE)   #in exceptional cases the innermost ring contains no cells and should therefore be discarded
      {
        rings <- rings[c(2:nrow(rings)),]
      }  
    } else
    {
      df <- NULL
      df <- read.table(mdf$cell.file[f], header=TRUE, sep="\t", na.strings=c("NA",""))   #cell data
      if("BEND" %in% colnames(df)){colnames(df)[colnames(df)=="BEND"] <- "TB2"}   #for backwards compatibility
      if("CRI" %in% colnames(df)){colnames(df)[colnames(df)=="CRI"] <- "TB2"}   #for backwards compatibility
      if("CA" %in% colnames(df)){colnames(df)[colnames(df)=="CA"] <- "LA"}   #for backwards compatibility
      df <- df[, c("ID", "YEAR", "CID", "RADDIST", "RADDISTR", "RRADDISTR", "LA",  
                   "ASP", "MAJAX", "KH", "CWTPI", "CWTBA", "CWTLE", "CWTRI", "CWTTAN", 
                   "CWTRAD", "CWTALL", "RTSR", "CTSR", "DRAD", "DTAN", "TB2",   
                   "CWA", "RWD")]
      df2 <- NULL
      df2 <- read.table(mdf$ring.file[f], header=TRUE, sep="\t", na.strings=c("NA",""))   #ring data
      df2 <- df2[,c("ID","YEAR","MRW")]
      if (df2$YEAR[1] %in% levels(as.factor(df$YEAR))==FALSE)   #in exceptional cases the innermost ring contains no cells and should therefore be discarded
      {
          df2 <- df2[c(2:nrow(df2)),]
      }  

      ### Find overlapping years and select data from ring with more cells    
      overlap.yr1 <- NULL
      overlap.yr2 <- NULL
      
      if (length(which(levels(as.factor(df$YEAR)) %in% levels(as.factor(rings$YEAR))==TRUE))>0)   #check whether there is any overlapping year: cell within ring data
      {
        overlap.yr1 <- levels(as.factor(df$YEAR))[which(levels(as.factor(df$YEAR)) %in% levels(as.factor(cells$YEAR))==TRUE)]   #cell data
        overlap.yr2 <- levels(as.factor(df2$YEAR))[which(levels(as.factor(df2$YEAR)) %in% levels(as.factor(rings$YEAR))==TRUE)]   #ring data
        overlap.yr <- unique(sort(c(overlap.yr1, overlap.yr2)))   #combine information of overlapping years from cell and ring data
        df.double <- NULL
        df.double2 <- NULL
        for (y in 1:length(overlap.yr))   #create dataframe of overlapping years with data from file with more cells 
        {
          # print(paste("overlapping years = ",overlap.yr[y], sep=""))
          print(paste("       overlapping years = ",overlap.yr[y], sep=""))
          A <- length(which(df$YEAR==as.numeric(overlap.yr[y])))   #how many rows (=cells) in target year of new file
          B <- length(which(cells$YEAR==as.numeric(overlap.yr[y])))   #how many rows (=cells) in target year of previously added file?
          ifelse(A < B, df.double <- rbind(df.double, cells[cells$YEAR==overlap.yr[y],]), df.double <- rbind(df.double, df[df$YEAR==overlap.yr[y],])) #for overlapping rings, create new cell data frame with data taken from ring with more cells (either previously added [B] or newly targeted file [A])          
          ifelse(A < B, df.double2 <- rbind(df.double2, rings[rings$YEAR==overlap.yr[y],]), df.double2<-rbind(df.double2, df2[df2$YEAR==overlap.yr[y],])) #for overlapping rings, create new ring data frame with data taken from ring with more cells (either previously added [B] or newly targeted file [A])          
        } 
        cells <- cells[-which(cells$YEAR %in% as.numeric(overlap.yr)), ]    #remove data of overlapping years from summary file
        cells <- rbind(cells, df.double)   #add updated data of overlapping years
        if (sum(which(df$YEAR %in% as.numeric(overlap.yr)))>0)   #in exceptional cases it is possible that there are no cells in the overlapping rings!
        {  
          df <- df[-which(df$YEAR %in% as.numeric(overlap.yr)), ]   #remove data of overlapping years from new file 
        }
        cells <- rbind(cells, df)   #add new data from non-overlapping years to summary file
        if (nrow(df.double2)>0)   #in some exceptional cases there is no ring data in the overlapping year
        {
          rings <- rings[-which(rings$YEAR %in% as.numeric(overlap.yr)), ]   #remove data of overlapping years from summary file
          rings <- rbind(rings, df.double2)   #add updated data of overlapping years
          if (length(which(df2$YEAR %in% as.numeric(overlap.yr)))>0)    #remove data of overlapping years from new file, but only if there are overlapping years
          {
            df2 <- df2[-which(df2$YEAR %in% as.numeric(overlap.yr)), ]
          }  
        }
        rings <- rbind(rings, df2)   #add new data from non-overlapping years to summary file
      } else
      {  
        cells <- rbind(cells, df)   #if no overlapping years, simply append new data
        rings <- rbind(rings, df2)   #if no overlapping years, simply append new data  
        # rings <- arrange(rings, YEAR,-MRW)   #order by ring and decreasing MRW (only relevant if duplicate ring, one without cells, so not detected above)
        rings <- rings[order(rings$YEAR, -rings$MRW),]   #order by ring and decreasing MRW (only relevant if duplicate ring, one without cells, so not detected above)
        rings <- rings[!duplicated(rings[,c("YEAR")]),]   #delete any duplicate ring, if any (the second smaller based on the previous ordering)
      }        
    }   
  }    
  
  ### 7.4 Correct for occasionally wrong CWA measurements (cells at edge, with large lateral wall integration)
  cells$MAX.CWT <- ifelse(!is.na(cells$CWTPI), cells$CWTPI, 0)
  cells$MAX.CWT <- ifelse(!is.na(cells$CWTBA), ifelse(cells$CWTBA > cells$MAX.CWT, cells$CWTBA, cells$MAX.CWT), cells$MAX.CWT)
  cells$MAX.CWT <- ifelse(!is.na(cells$CWTLE), ifelse(cells$CWTLE > cells$MAX.CWT, cells$CWTLE, cells$MAX.CWT), cells$MAX.CWT)
  cells$MAX.CWT <- ifelse(!is.na(cells$CWTRI), ifelse(cells$CWTRI > cells$MAX.CWT, cells$CWTRI, cells$MAX.CWT), cells$MAX.CWT)
  
  cells$LR <- sqrt(cells$LA / pi)
  cells$MAX.CWA <- 2 * (pi * (cells$LR + cells$MAX.CWT) ^ 2 - cells$LA)   #2x theoretical value to allow some margin
  
  cells$CWA <- ifelse(abs(cells$CWA) > cells$MAX.CWA, NA, cells$CWA)
  cells$RWD <- ifelse(abs(cells$CWA) > cells$MAX.CWA, NA, cells$RWD)

  cells$MAX.CWT <- NULL
  cells$LR <- NULL
  cells$MAX.CWA <- NULL
      
  ### 7.5 Add cell area 
  cells$TCA <- ifelse(cells$CWA>0, cells$LA + cells$CWA, ifelse(is.na(cells$CWA), NA, -(cells$LA-cells$CWA)))
  
  ### 7.6 Add hydraulic diameter of each cell
  cells$a <- 2*sqrt(cells$ASP*cells$LA/pi)
  cells$b <- cells$a/cells$ASP
  cells$DH <- sqrt((2*cells$a^2*cells$b^2)/(cells$a^2+cells$b^2))
  cells$a <- NULL
  cells$b <- NULL
  
  ### 7.7 Add special anatomical cell density
  cells$RWD2 <- cells$CWTRAD/cells$DRAD
  
  ### 7.8 Order data by ascending year
  cells <- cells[order(cells$YEAR, cells$RRADDISTR),]
  rings <- rings[order(rings$YEAR),]
  
  ### 7.9 Write data to files
  output <- paste(mdf$tree[1], "_Output_Cells.txt", sep="")
  write.table(cells, file=output, row.names = FALSE)
  
  output2 <- paste(mdf$tree[1], "_Output_Rings.txt", sep="")
  write.table(rings, file=output2, row.names = FALSE)  
}  
  
Sys.time() - t


### 8. Clean up ####
rm(list=ls()) # clean desk 

