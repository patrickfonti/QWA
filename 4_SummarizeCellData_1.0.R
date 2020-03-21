###########################################################################
# This scripts summarizes the cells data, extracts only the columns 
# needed for further analysis (ID, IMAGE.ID, YEAR, RRADDISTR, LA, CWTRAD),
# and writes it to a new text file.
#   
# 25 July 2019, Georg von Arx
###########################################################################


### 1. Clean and load library ####
rm(list=ls()) # clean desk 


### 2. Define directory containing all the cell files (only on top level, no sub-directories!) ####
setwd("D:/TrainingSchools_QWA/NADEF_QWA_2019_Cody/Analysis/R Analysis/Series")
homedir <- getwd()


### 3. Combine cell files into a single file ####
(cells.files <- grep("Output_Cells", list.files(path=homedir), value=TRUE))   #find the annual stats files

cells <- NULL

for (i in c(1:length(cells.files)))
{
  if (i==1)
  {
    cells <- read.table(cells.files[i], header=T)   #cell data      
  } else
  {  
    dt <- read.table(cells.files[i], header=T)   #cell data
    cells <- rbind(cells,dt)
  }
}


### 4. Check file, add core ID, and replace intra-ring position values <0 and >100 ####
head(cells)

cells <- cells[,c("ID","YEAR","RRADDISTR","LA","CWTRAD","KH")]   #crop data
head(cells)

colnames(cells)[1] <- "IMAGE.ID"   #change column label
cells$ID <- substr(cells$IMAGE.ID,1,12)   #create new column with core ID

cells <- cells[,c("ID","IMAGE.ID","YEAR","RRADDISTR","LA","CWTRAD","KH")]   #re-arrange columns in data file
cells1 <- cells
cells <- cells1
cells$RRADDISTR <- ifelse(cells$RRADDISTR < 0, cells$RRADDISTR <- 0, cells$RRADDISTR)   #correct values <0
cells$RRADDISTR <- ifelse(cells$RRADDISTR < 100, cells$RRADDISTR, cells$RRADDISTR <- 100)   #correct vaues >100

head(cells)
str(cells)


### 5. Write data to file ####
write.table(cells, "CellData.txt", sep="\t", na="NA", dec=".")


### 6. Clean up ####
rm(list=ls()) # clean desk 
