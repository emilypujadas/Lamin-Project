library(TopDom)

##------------------------------ values for variables ------------------------------##

## Change these to your condition labels
conds <- c("24hrAuxin", "untreated", "withdraw")
## Change these to the chromosome annotations you're using
chroms = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X')

windowSize = 5 ## TopDom window
res = 25 ## Change to your resolution
norm = "VC" ## Change to your normalization

##-------------------------------values for variables------------------------------##

## Function to loop through and call TADs
callTADs <- function(exp, conds, res, windowSize, norm, chroms){

for (cond in 1:length(conds)){

  cond <- conds[cond] ## Condition labels

  ## Input and output paths on Quest
  input.path = paste0('/projects/b1042/BackmanLab/HiC2/opt/juicer/work/',exp,'/contact_data/',cond,'/merged/')
  output.path = paste0('/projects/b1042/BackmanLab/HiC2/opt/juicer/work/',exp,'/contact_data/',cond,'/merged/domains/')


  for (chrom in chroms) {
    print(chrom)

    input.fname = paste('chr',chrom,'-observed_',res,'Kb_',norm,'norm.txt',sep="")
    output.fname = paste('chr',chrom,'-TADs-window',as.character(windowSize),sep="")

    topDom_file = paste(input.path,input.fname,sep="")
    topDom_outFile = paste(output.path,output.fname,sep="")

    topDom_out <- TopDom(topDom_file, window.size=windowSize, outFile=topDom_outFile)
    }
  }
}

callTADs("Lamin_HiC", conds=c("24hrAuxin", "untreated", "withdraw"), res = 25, windowSize = 5, norm = "VC", chroms = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'))
