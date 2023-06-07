library(ggplot2)
library(viridis)
library(dplyr)
library(tidyverse)

mergeContacts <- function(exp, chrom, res, norm){
  
  # Set directory variables here
  dir.root <- "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/Lamin_HiC/contact_data"
  dir.1 <-paste0(dir.root,"/",exp,"/Rep1/contacts_TopDom/")
  dir.2 <-paste0(dir.root,"/",exp,"/Rep2/contacts_TopDom/")
  dir.3 <-paste0(dir.root,"/",exp,"/Rep3/contacts_TopDom/")
  dir.4 <-paste0(dir.root,"/",exp,"/Rep4/contacts_TopDom/")
  dir.5 <-paste0(dir.root,"/",exp,"/Rep5/contacts_TopDom/")
  dir.6 <-paste0(dir.root,"/",exp,"/Rep6/contacts_TopDom/")
  print(dir.1)
  
  ## Read in table
  m1 <- read.table(file.path(paste0(dir.1, chrom, "-observed_", res,"Kb_",norm,"norm.txt")), sep= "\t", header= F)
  m2 <- read.table(file.path(paste0(dir.2, chrom, "-observed_", res,"Kb_",norm,"norm.txt")), sep= "\t", header= F)
  m3 <- read.table(file.path(paste0(dir.3, chrom, "-observed_", res,"Kb_",norm,"norm.txt")), sep= "\t", header= F)
  m4 <- read.table(file.path(paste0(dir.4, chrom, "-observed_", res,"Kb_",norm,"norm.txt")), sep= "\t", header= F)
  m5 <- read.table(file.path(paste0(dir.5, chrom, "-observed_", res,"Kb_",norm,"norm.txt")), sep= "\t", header= F)
  m6 <- read.table(file.path(paste0(dir.6, chrom, "-observed_", res,"Kb_",norm,"norm.txt")), sep= "\t", header= F)
  
  ## Save positions for later
  ch3 <- m1$V1
  pos1 <- m1$V2
  pos2 <- m1$V3
  
  m1 <- m1[,-c(1:3)]
  m2 <- m2[,-c(1:3)]
  m3 <- m3[,-c(1:3)]
  m4 <- m4[,-c(1:3)]
  m5 <- m5[,-c(1:3)]
  m6 <- m6[,-c(1:3)]
  
  ## Sum matrices
  mat <- m1+m2+m3+m4+m5+m6
  
  ## Add positions back 
  mat <- cbind(ch3, pos1, pos2, mat)
  head(mat, 20)
  
  write.table(mat, file = paste0(dir.root,"/",exp,"/",chrom, "-observed_", res,"Kb_",norm,"norm.txt"), sep = "\t", row.names = F, col.names = F)
  cat(paste('finished processing ',exp,": ",chrom,'\n\n'))
}


## loop for import
chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", 
            "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
            "chr22", "chrX", "chrY")

for (chrom in 1:length(chroms)){
  
  chrom <- chroms[chrom]
  print(chrom)
  
  mergeContacts(exp="24hrAuxin", chrom, res=25, norm="VC")
}


for (chrom in 1:length(chroms)){
  
  chrom <- chroms[chrom]
  print(chrom)
  
  mergeContacts(exp="untreated", chrom, res=25, norm="VC")
}


for (chrom in 1:length(chroms)){
  
  chrom <- chroms[chrom]
  print(chrom)
  
  mergeContacts(exp="withdraw", chrom, res=25, norm="VC")
}   
