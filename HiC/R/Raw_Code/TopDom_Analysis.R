rm(list = ls())
dev.off()

# import packages
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyverse)
library(forcats)
library(stringr)

#########################################################################################
# TopDom hist - for merged contacts
#########################################################################################

# Set directory variables here
dir.domains <- "/Volumes/My Passport/projects/Emily/Lamins_HiC/HiC_Analysis_Dir/Lamin_HiC/contact_domains/topdom_domains"

chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", 
            "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
            "chr22", "chrX", "chrY")

## Condition labels
#exps <- c("24hrAuxin","untreated", "withdraw")
exps <- c("24hrAuxin","untreated")

files <- list.files(paste0(dir.domains,"/",exps[1],"/merged/domains/"), pattern=".bed", all.files=T)
files

## loop for import

all.my.files <- c()
for (exp in 1:length(exps)){
  
  exp <- exps[exp]
  filepath <- paste0(dir.domains,"/",exp,"/merged/domains/")

    
    for (file in 1:length(files)){
      infile <- files[file]
      
      cat(paste('appending TopDom .BED files',exp,":",infile,'\n\n'))
      
      infile <- read.table(file.path(filepath, infile), sep= "\t", header= FALSE)
      
      infile <-data.frame(infile)
      infile <- infile[,c(1:4)]
      
      infile <- infile %>% 
        mutate(V1=str_trim(infile$V1, "right"))%>% # trim white space from chrom names
        mutate(name =sapply(strsplit(exp, "_"), `[`, 1)) %>%
        mutate(size= abs(infile$V2-infile$V3)) 
      
      all.my.files <- rbind(all.my.files,infile)
      
    }
    
  }
  
## Wranlge data
tads.length <- all.my.files
tads.length <- tads.length[tads.length$V4 == "domain",]
colnames(tads.length) <- c("chr","start","end","type", "name","size")

library(plyranges)
library(GenomicRanges)

## Convert to Granges object
get.names  <- tads.length  %>% 
  transform( seqnames= tads.length$chr, start = tads.length$start, end = tads.length$end)  %>% 
  as_granges()

tads.length$regions = paste0(seqnames(get.names),':',start(get.names),'-',end(get.names))

# get size summary
tads.length %>% group_by(name) %>% 
  summarize(m = mean(size), # calculates the mean
            s = sd(size),   # calculates the standard deviation
            n = n()) 

#########################################################################################
# TopDom hist - looping through the files must be done on Quest if using Reps , since the contact files
# are very large. use R analytics node to do the first part below.
#########################################################################################

# Set directory variables here
dir.domains <- "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/Lamin_HiC/contact_data"

chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", 
            "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
            "chr22", "chrX", "chrY")

## Condition labels
exps <- c("24hrAuxin","untreated", "withdraw")

reps <- c("Rep1","Rep2", "Rep3", "Rep4", "Rep5", "Rep6")

files <- list.files(paste0(dir.domains,"/",exps[1],"/",reps[1],"/contacts_TopDom/domains/"), pattern=".bed", all.files=T)
files

## loop for import

all.my.files <- c()
for (exp in 1:length(exps)){
  
  exp <- exps[exp]
  filepath <- paste0(dir.domains,"/",exp,"/")
  
  for (rep in 1:length(reps)){
    
    rep <- reps[rep]
    filepath <- paste0(dir.domains,"/",exp,"/",rep, "/contacts_TopDom/domains/")
    
    for (file in 1:length(files)){
      infile <- files[file]
      
      cat(paste('appending TopDom .BED files',exp,":",infile,'\n\n'))
      
      infile <- read.table(file.path(filepath, infile), sep= "\t", header= FALSE)
      
      infile <-data.frame(infile)
      infile <- infile[,c(1:4)]
      
      infile <- infile %>% 
        mutate(name =sapply(strsplit(exp, "_"), `[`, 1)) %>% 
        mutate(rep =sapply(strsplit(rep, "_"), `[`, 1)) %>% 
        mutate(size= abs(infile$V2-infile$V3)) 
      
      all.my.files <- rbind(all.my.files,infile)
      
    }
    
  }
  
}

##--------------------Read in table from Quest here--------------------##

# Set directory variables here
dir.domains <- "/Volumes/My Passport/projects/Emily/Lamins_HiC/HiC_Analysis_Dir/Lamin_HiC/contact_domains/topdom_domains"

## Read in table
all.my.files <- read.table(file.path(dir.domains, "lamins_TopDom_domains.csv"), sep= ",", header= T)
all.my.files <- all.my.files[,-1]

## Wranlge data
tads.length <- all.my.files
tads.length <- tads.length[tads.length$V4 == "domain",]
colnames(tads.length) <- c("chr","start","end","type", "name","rep", "size")

library(plyranges)
library(GenomicRanges)

##--------------------If getting rid of redundant domains, do--------------------##

## Convert to Granges object
tads.length  <- tads.length  %>% 
  transform( seqnames= tads.length$chr, start = tads.length$start, end = tads.length$end)  %>% 
  as_granges()

tads.length$regions = paste0(seqnames(tads.length),':',start(tads.length),'-',end(tads.length))
tads.length

## Group treatments
auxin <- as.data.frame(tads.length) %>%
  filter(name == "24hrAuxin")%>%
  distinct(regions, .keep_all = TRUE)

untreated <- as.data.frame(tads.length) %>%
  filter(name == "untreated")%>%
  distinct(regions, .keep_all = TRUE)

withdraw <- as.data.frame(tads.length) %>%
  filter(name == "withdraw")%>%
  distinct(regions, .keep_all = TRUE)

tads.length <- rbind(auxin, untreated, withdraw)

# get size summary
tads.length %>% group_by(name) %>% 
  summarize(m = mean(size), # calculates the mean
            s = sd(size),   # calculates the standard deviation
            n = n()) 

##--------------------if merging contact files to call domains on, do--------------------##

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

#########################################################################################
# Compaction analysis by loci
#########################################################################################

plottheme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                 panel.grid.major = element_line(),
                 #panel.grid.minor = element_line(colour = "tomato", size=.25, linetype = "dashed"),
                 strip.background=element_blank(), axis.text.x=element_text(family="Arial", colour="black"), axis.title.x=element_text(face="bold", size=20,family="Arial", colour="black", vjust=-2),
                 axis.text.y=element_text(family="Arial", colour="black"), 
                 axis.title.y=element_text(face="bold", size=15, family="Arial", colour="black", vjust=2),
                 text = element_text(size=10, family="Arial"), 
                 plot.title=element_text(size=20, face="bold", family="Arial", color="tomato", hjust=0.0, vjust=10, lineheight=1.5),
                 plot.subtitle=element_text(size=15, face="bold", family="Arial", color="black", hjust=0.0, lineheight=1.5),
                 legend.title = element_text(size=12, color= "tomato",face="bold"), 
                 legend.text = element_text(size=10),
                 legend.key=element_rect(fill=NA),
                 axis.ticks=element_line(colour="black"))

dir <- "/Volumes/My Passport/projects/Emily/Lamins_HiC/HiC_Analysis_Dir/Lamin_HiC/contact_domains/topdom_domains"

list.files(dir)

chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", 
            "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
            "chr22", "chrX")


DTL_Ratio <- function(infile, dir, res=25000, binsize = 2000000 ){

## bring in data
file <- read.table(file.path(dir, infile), sep= "\t", header= FALSE)

        ## Wrangle data
        index <- as.numeric(file$V2)
        file <- file[,-c(1:3)]
        rownames(file) <- index
        colnames(file) <- index
        file$index <- index
        
        ## Separate data by distal and local; local <= 3Mb, distal >= 3Mb
        distal.int <- file[index >=3000000 ,]
        local.int <- file[index  <=3000000 ,]
        distal.int <- distal.int[,-c(ncol(distal.int))]
        local.int <- local.int[,-c(ncol(local.int))]
        
        ## Gather contact frequencies
        distal <-  colSums(distal.int,na.rm=T)
        local <-  colSums(local.int, na.rm=T)
        ratio <- data.frame(distal, local,row.names = index)
        
        ## Generate log2 ratio
        ratio <- ratio  %>% 
        mutate(ratio = abs(log2((distal/local))))
        
        ## Remove non-values
        ratio[is.na(ratio$ratio),] <-0
        ratio[!is.finite(ratio$ratio),] <- 0
        
        ## Bin and reframe data
        ratio$index <- as.numeric(index)
        bin = round(((nrow(ratio)*res)/binsize)) -1
        ratio <- mutate(ratio, bins = ntile(index, bin))
        
        ## Gather binned sums
        sums = (tapply(ratio$ratio, ratio$bins, sum)/(binsize/res))
        index <- round(log10(seq((binsize),max(bin)*binsize, by = binsize)),2)
        
        ## New data frame
        ratio <- data.frame(sums, index)

return(ratio)
}

ratio.treat <- DTL_Ratio("chr3-observed_25Kb_VCnorm.txt", 
                         "/Volumes/My Passport/projects/Emily/Lamins_HiC/HiC_Analysis_Dir/Lamin_HiC/contact_domains/topdom_domains/24hrAuxin/merged")
ratio.ctrl <- DTL_Ratio("chr3-observed_25Kb_VCnorm.txt", 
                        "/Volumes/My Passport/projects/Emily/Lamins_HiC/HiC_Analysis_Dir/Lamin_HiC/contact_domains/topdom_domains/untreated/merged")

# For Chrom 5: scale_x_discrete(labels = c(6.30, 7.15, 7.41, 7.58, 7.70, 7.79, 7.87, 7.93, 7.99, 8.04, 8.09, 8.14, 8.19, 8.25),
                             #breaks = c(6.30, 7.15, 7.41, 7.58, 7.70, 7.79, 7.87, 7.93, 7.99, 8.04, 8.09, 8.14, 8.19, 8.25))+

## Prepare data for plotting
ratio.treat$group <- rep("24hrAuxin", nrow(ratio.treat))
ratio.ctrl$group <- rep("untreated", nrow(ratio.ctrl))
ratio<- rbind(ratio.treat, ratio.ctrl)
ratio$factor <- seq(1,nrow(ratio), by = 1)
colnames(ratio) <- c("ratio", "index", "group", "factor" )

library(ggplot2)
library(forcats)
library(hrbrthemes)

## Plot DTL as log values
p <- ratio %>% 
  mutate(index = fct_reorder(as.character(index), factor)) %>% 
  ggplot(aes(x= index, y= ratio, group=group)) + 
  geom_line(aes(color=group)) + 
  plottheme + scale_x_discrete(labels = c(6.30, 7.15, 7.41, 7.58, 7.70, 7.79, 7.87, 7.93, 7.99, 8.04, 8.09, 8.14, 8.19, 8.25),
                               breaks = c(6.30, 7.15, 7.41, 7.58, 7.70, 7.79, 7.87, 7.93, 7.99, 8.04, 8.09, 8.14, 8.19, 8.25))+
  ylab("Log2 Ratio")+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1, size = 15),
        axis.text.y = element_text(angle = 0,size = 20))+ xlab("Log10 BP")+
  labs(title="", subtitle=" Distal-To-Local Ratio: Chromosome 3", #subtitle="Top 1000 Variable Regions", 
       caption="")
p

#########################################################################################
# Generate histograms of contacts
#########################################################################################
library(viridis)
library(dplyr)
library(tidyverse)

## Wrangle data to plot TAD contacts by Loci
tads.hist1 <- tads.length[,c(1,2,5)]
tads.hist2 <- tads.length[,c(1,3,5)]
colnames(tads.hist1) <- c("chr","locus","name")
colnames(tads.hist2) <- c("chr","locus", "name")
tads.hist <- rbind(tads.hist1,tads.hist1)
rm(tads.hist1,tads.hist2)

## Gather loci contacts by condition
tads.hist<-tads.hist[,c(2,3)]
tads.hist <- tads.hist %>%
  gather(key="name", value="locus")

## Plot all four histograms #TopDom_TAD_freq
p.hist <- tads.hist %>%
  mutate(name = fct_reorder(name, locus)) %>%
  ggplot( aes(x=locus, color=name, fill=name)) +
  geom_histogram(alpha=0.6, bins=100) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() + labs(title="TAD loci contacts distribution.", x="Loci", y="Freq")+
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8),
    axis.text.x = element_text(angle = 90)
  ) +
  xlab("") + ylab("Frequency") + facet_wrap(~name)

p.hist

#########################################################################################
# TopDom: Tad size plots 
#########################################################################################
tads.length <- as.data.frame(tads.length)

library(ggforce)
library(PupillometryR)

# Plot tads by size boxplot #TopDom_TAD_box_all
tads.length %>%
  ggplot( aes(x=name, y=size, fill=name, col = name)) +
  geom_boxplot(notch=TRUE, notchwidth = 0.3, outlier.shape = NA, aes(color = "grey")) +
  scale_fill_manual(values=c("#C70039", "#FF5733"))+
  #scale_y_discrete(limits=c(0, 7e05))+
  theme_ipsum() + 
  scale_y_continuous(labels = c("0 kb", "100 kb", "200 kb", "300 kb", "400 kb", "500 kb", "600 kb"),breaks = c(0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5),
                     limits = quantile(tads.length$size, c(0.1, 0.9)))+
  theme(legend.position="none", plot.title = element_text(size=11)) +
  ggtitle("TADs By Size") + xlab("")

# Plot tads by size boxplot #TopDom_TAD_box_all
tads.length %>%
  ggplot( aes(x=name, y=size, fill=name)) +
  geom_flat_violin(color = NA, position = position_nudge(x = .15))+ coord_flip() +
  scale_fill_manual(values=c("#C70039", "#FF5733"))+
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) +
  theme_ipsum() + geom_boxplot(width = 0.2,notch=TRUE, notchwidth = 0.3, outlier.shape = NA, color = "white")+
  scale_y_continuous(labels = c("0 kb", "100 kb", "200 kb", "300 kb", "400 kb", "500 kb", "600 kb", "700 kb", "800 kb"),breaks = c(0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 7e5, 8e5),limits = quantile(tads.length$size, c(0.001, 0.95)))+
  theme(legend.position="none", plot.title = element_text(size=11), axis.text.x = element_text(angle = 90)) +
  ggtitle("TADs By Size") + xlab("")

## compare tads for each chromosome for each contrast
ctrl.tads <- tads.length[tads.length$name=="untreated",]

#CTRL_TAD_byChrom
ctrl.tads %>%
  group_by(chr) %>% 
  ggplot( aes(x=chr, y=size, fill=chr)) +
  geom_boxplot(notch=TRUE, notchwidth = 0.8, outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() + stat_boxplot(geom = "errorbar", width = 0.15) + 
  scale_y_continuous(labels = c("0 kb", "100 kb", "200 kb", "300 kb", "400 kb", "500 kb", "600 kb"),
  breaks = c(0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5),
  limits = quantile(tads.length$size, c(0.1, 0.9)))+
  theme(legend.position="none", plot.title = element_text(size=11), axis.text.x = element_text(angle = 90)) +
  ggtitle("Control Tads by chromosome") + xlab("")

# compare tads for each chromosome for each contrast
treat.tads <- tads.length[tads.length$name=="24hrAuxin",]

#2hrAuxin_TAD_byChrom
treat.tads %>%
  group_by(chr) %>% 
  ggplot( aes(x=chr, y=size, fill=chr)) +
  geom_boxplot(notch=TRUE, notchwidth = 0.8, outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() + stat_boxplot(geom = "errorbar", width = 0.15) + 
  scale_y_continuous(labels = c("0 kb", "100 kb", "200 kb", "300 kb", "400 kb", "500 kb", "600 kb"),
  breaks = c(0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5),
  limits = quantile(tads.length$size, c(0.1, 0.9)))+
  theme(legend.position="none", plot.title = element_text(size=11), axis.text.x = element_text(angle = 90)) +
  ggtitle("24hr Auxin Tads by chromosome") + xlab("")

##-------------------bin tads-----------------------##

res = 2500
binsize= 25000

tads.length <- tads.length[order(tads.length$size),]

## Generate bins
plot.bins <-  mutate(tads.length, bin = cut(tads.length$size, breaks=c(seq(0, (max(tads.length$size)+res), by = 50000))))

## Generate labels for bins in new column
plot.bins <- mutate(plot.bins, label = paste0(bin,"_",name))
max(tads.length$size)

## Gather frequency and names for each bin prior to plotting
x <-as.data.frame(table(plot.bins$label))
x$name =sapply(strsplit(as.character(x$Var1), "_"), `[`, 2)
x$breaks = sapply(strsplit(as.character(x$Var1), "_"), `[`, 1)

## Data wrongling
intervals <-as.data.frame(cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", x$breaks) ),
                                upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", x$breaks) )))

x$size <- intervals$upper
colnames(x) <- c("cut", "count", "name", "breaks", "size")
x <- x[order(x$size),]
x <- x[-c(169),]

# Plot large tads #TopDom_TADs_bySize
p <- x %>%
  mutate(name = fct_reorder(name, size)) %>% 
  ggplot(aes(x= size, y= count, group=name)) + geom_point(aes(color=name),size=0.1) + 
  geom_line(aes(color=name))+
  plottheme + 
  #scale_y_continuous(limits=c(0,4000), breaks=c(0,500,1000,1500,2000,2500,3000,3500,4000))+
  scale_x_continuous( labels =c("0", "1.5 MB", "3 MB","4.5 MB", "6 MB", "7.5 MB", "9 MB"), breaks = c(50000, 1500000, 3000000, 4500000, 6000000, 7500000, 9000000),  expand = c(0.025, 0.025))+
  scale_color_manual(values=c("#e88471", "#39b185", "#39b185"))+ ylab("Aggregate N")+ xlab("Aggregate Size")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size = 13))+
  labs(title="Total Number of Aggregates", subtitle="Binsize = 50kbp", #subtitle="Top 1000 Variable Regions", 
       caption="")
p

##-------------------------------------------------##

## number of tads
l.tads <- tads.length[tads.length$size >=1000000 ,] # large tads
ml.tads <- tads.length[tads.length$size >=750000 & tads.length$size <=1000000 ,]# medium large tads
m.tads <- tads.length[tads.length$size >=500000 & tads.length$size <=750000 ,] #medium tads
ms.tads <- tads.length[tads.length$size >=250000 & tads.length$size <=500000 ,]# medium small tads
s.tads <- tads.length[tads.length$size >=100000 & tads.length$size <=250000 ,]  #small tads
st.tads <- tads.length[tads.length$size <=100000,]#subtads

## Relabel TADs
labelTads<- function(data, lab,n){

    ctrl <- rep(paste0(lab," untreated"), times = nrow(data[data$name == "untreated",]))
    treatment <- rep(paste0(lab," 24hrAuxin"), times = nrow(data[data$name == "24hrAuxin",]))
    withdraw <- rep(paste0(lab," withdraw"), times = nrow(data[data$name == "withdraw",]))
    label <- c(ctrl, treatment, withdraw)
    
    data$name <- label
    
    data$factor <- rep(n,times = nrow(data))
    
    return(data)

}

l.tads <- labelTads(l.tads, "L", 1)
ml.tads <- labelTads(ml.tads, "ML", 2)
m.tads <- labelTads(m.tads, "M", 3)
ms.tads <- labelTads(ms.tads, "MS", 4)
s.tads <- labelTads(s.tads, "S", 5)
st.tads <- labelTads(st.tads, "ST", 6)

## Bind them together
tad.size <- rbind(l.tads, ml.tads, m.tads, ms.tads, s.tads, st.tads)

# Plot large tads #TopDom_TADs_bySize
p <- tad.size %>%
  mutate(name = fct_reorder(name, factor)) %>% 
  ggplot( aes(x=name, y=size, fill=name)) +
  geom_boxplot(notch=TRUE, notchwidth = 0.8, outlier.shape = NA) +
  scale_fill_manual(values=c("#99d8c9", "#43a2ca", "#99d8c9", "#43a2ca", "#99d8c9", "#43a2ca",
                             "#99d8c9", "#43a2ca", "#99d8c9", "#43a2ca", "#99d8c9", "#43a2ca")) +
  #geom_jitter(color="black", size=0.4, alpha=0.1) +
  theme_ipsum() + stat_boxplot(geom = "errorbar", width = 0.15) + scale_y_continuous(limits = c(0, 2.6e6))+
  scale_x_discrete( labels =c("Large treatment", "Large CTRL", "Medium Large treatment", "Medium Large CTRL",
                              "Medium treatment", "Medium CTRL","Medium Small treatment", "Medium Small CTRL",
                              "Small treatment", "Small CTRL","SubTAD treatment", "SubTAD CTRL"))+
  theme(legend.position="none", plot.title = element_text(size=11)) +
  ggtitle("Aggregates by Size") + ylab("Base Pairs")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
        axis.text.y = element_text(angle = 0,size = 20))+ xlab("Aggregate Group")

p

## Get summary statistics
getSummary <- function(data, label){

tad <- data %>% group_by(name) %>% 
  summarize(mean.size = mean(size), # calculates the mean
            sd = sd(size),   # calculates the standard deviation
            tad.n = n()) 
tad$group <- rep(label, nrow(tad))

return(tad)

}

## Run function
tad.all <- getSummary(tads.length, "all")
tads.l <- getSummary(l.tads, "large")
tads.ml <- getSummary(ml.tads, "medium large")
tads.m <- getSummary(m.tads, "medium")
tads.ms <- getSummary(ms.tads, "medium small")
tads.s <- getSummary(s.tads, "small")
tads.st <- getSummary(st.tads, "subTAD")

#####################################################################################
# LAD BED TAD intersection for LMBN1 and LMBN2                                    
#####################################################################################

library(plyranges)
library(GenomicRanges)

bed.dir <- "/Volumes/My Passport/projects/Emily/Lamins_RNAseq/Raw_Data/BEDs"
list.files(bed.dir)

##------------------------------Luay's Lamin B1 B2------------------------------##
## Using Luay's BED file for Lamin B1 and B2
bed.lmb1 <- read_bed(file.path(bed.dir, "LMB1_4DNFICCV71TZ.bed") , col_names = NULL, genome_info = "hg38")
seqlevels(bed.lmb1)

bed.lmb2 <- read_bed(file.path(bed.dir, "LMB2_4DNFIBQH62LX.bed") , col_names = NULL, genome_info = "hg38")
bed.lmb2

## Union of lmb1 and lmb2 files
bed.lam.union <- GenomicRanges::reduce(c(bed.lmb1, bed.lmb2))
seqlevels(bed.lam.union)

##----------------------------------Parse LADs----------------------------------##
## Name all of the peaks based on their genomic position
names(bed.lam.union) = paste0(seqnames(bed.lam.union),':',start(bed.lam.union),'-',end(bed.lam.union))

## Assign names to temp variable
region.names <- names(bed.lam.union)

## Convert to Granges object
## Wranlge data
tads.length <- all.my.files
tads.length <- tads.length[tads.length$V4 == "domain",]
colnames(tads.length) <- c("chr","start","end","type", "name","size")

tads.length  <- tads.length  %>% 
  transform( seqnames= tads.length$chr, start = tads.length$start, end = tads.length$end)  %>% 
  as_granges()

#tads.length$regions = paste0(seqnames(tads.length),':',start(tads.length),'-',end(tads.length))
#tads.length

## Find intersection between LADs and all genes
## intersect_rng <- join_overlap_intersect(query, subject) https://bioconductor.org/packages/devel/bioc/vignettes/plyranges/inst/doc/an-introduction.html
lamin.overlaps <- join_overlap_intersect(bed.lam.union, tads.length)
lamin.overlaps

## Find non-overlapping regions https://support.bioconductor.org/p/74077/ sp over() method gr1[!gr1 %over% gr2,]
non.lamins <- tads.length[!tads.length %over% lamin.overlaps,]
non.lamins

## Convert to data frame
lamin.overlaps <- as.data.frame(lamin.overlaps, row.names = seq(1:11087))
non.lamins <- as.data.frame(non.lamins, row.names = seq(1:5698))

lamin.overlaps <- lamin.overlaps[lamin.overlaps$name == "24hrAuxin" |lamin.overlaps$name == "untreated" ,]
non.lamins <- non.lamins[non.lamins$name == "24hrAuxin" |non.lamins$name == "untreated" ,]

##-------------------------Plots of all -------------------------##

## Relabel TADs
labelTads<- function(data, lab,n){
  
  ctrl <- rep(paste0(lab," untreated"), times = nrow(data[data$name == "untreated",]))
  treatment <- rep(paste0(lab," 24hrAuxin"), times = nrow(data[data$name == "24hrAuxin",]))

  label <- c(ctrl, treatment)
  
  data$name <- label
  
  data$factor <- rep(n,times = nrow(data))
  
  return(data)
  
}

overlaps <- labelTads(lamin.overlaps, "LADs: ", 1)
nonoverlaps <- labelTads(non.lamins, "Non-LADs:", 2)


## Bind them together
tad.size <- rbind(overlaps, nonoverlaps)

library(hrbrthemes)
library(ggforce)
library(PupillometryR)

# Plot tads by size boxplot #TopDom_TAD_box_all
tad.size %>%
  mutate(name = fct_reorder(name, factor)) %>% 
  ggplot( aes(x=name, y=size, fill=name, col = name)) +
  geom_boxplot(notch=TRUE, notchwidth = 0.3, outlier.shape = NA, aes(color = "grey")) +
  scale_fill_manual(values=c("#C70039", "#FF5733", "#C70039", "#FF5733"))+
  #scale_y_discrete(limits=c(0, 7e05))+
  theme_ipsum() + 
  scale_y_continuous(labels = c("0 kb", "100 kb", "200 kb", "300 kb", "400 kb", "500 kb", "600 kb"),breaks = c(0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5),
                     limits = quantile(tads.length$size, c(0.1, 0.9)))+
  theme(legend.position="none", plot.title = element_text(size=11)) +
  ggtitle("TADs By Size") + xlab("")

# Plot tads by size boxplot #TopDom_TAD_box_all
tad.size %>%
  mutate(name = fct_reorder(name, factor, .desc = F)) %>% 
  ggplot( aes(x=name, y=size, fill=name)) +
  geom_flat_violin(color = NA, position = position_nudge(x = .15))+ coord_flip() +
  scale_fill_manual(values=c("#800080","#00008B", "#C70039", "#FF5733"))+
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) +
  theme_ipsum() + geom_boxplot(width = 0.2,notch=TRUE, notchwidth = 0.3, outlier.shape = NA, color = "white")+
  scale_y_continuous(labels = c("0 kb", "100 kb", "200 kb", "300 kb", "400 kb", "500 kb", "600 kb", "700 kb", "800 kb"),breaks = c(0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 7e5, 8e5),limits = quantile(tads.length$size, c(0.001, 0.95)))+
  theme(legend.position="none", plot.title = element_text(size=11), axis.text.x = element_text(angle = 90)) +
  ggtitle("TADs By Size") + xlab("")


##-------------------------Stats on TADs within LADs-------------------------##

## number of tads
l.tads.o <- lamin.overlaps[lamin.overlaps$size >=1000000 ,] # large tads
ml.tads.o <- lamin.overlaps[lamin.overlaps$size >=750000 & lamin.overlaps$size <=1000000 ,]# medium large tads
m.tads.o <- lamin.overlaps[lamin.overlaps$size >=500000 & lamin.overlaps$size <=750000 ,] #medium tads
ms.tads.o <- lamin.overlaps[lamin.overlaps$size >=250000 & lamin.overlaps$size <=500000 ,]# medium small tads
s.tads.o <- lamin.overlaps[lamin.overlaps$size >=100000 & lamin.overlaps$size <=250000 ,]  #small tads
st.tads.o <- lamin.overlaps[lamin.overlaps$size <=100000,]#subtads

## number of tads
l.tads.n <- non.lamins[non.lamins$size >=1000000 ,] # large tads
ml.tads.n <- non.lamins[non.lamins$size >=750000 & non.lamins$size <=1000000 ,]# medium large tads
m.tads.n <- non.lamins[non.lamins$size >=500000 & non.lamins$size <=750000 ,] #medium tads
ms.tads.n <- non.lamins[non.lamins$size >=250000 & non.lamins$size <=500000 ,]# medium small tads
s.tads.n <- non.lamins[non.lamins$size >=100000 & non.lamins$size <=250000 ,]  #small tads
st.tads.n <- non.lamins[non.lamins$size <=100000,]#subtads

## Relabel TADs
labelTads<- function(data, lab,n){
  
  ctrl <- rep(paste0(lab," untreated"), times = nrow(data[data$name == "untreated",]))
  treatment <- rep(paste0(lab," 24hrAuxin"), times = nrow(data[data$name == "24hrAuxin",]))
  label <- c(ctrl, treatment)
  
  data$name <- label
  
  data$factor <- rep(n,times = nrow(data))
  
  return(data)
  
}

l.tads.o <- labelTads(l.tads.o, "L overlap", 1)
ml.tads.o <- labelTads(ml.tads.o, "ML overlap", 3)
m.tads.o <- labelTads(m.tads.o, "M overlap", 5)
ms.tads.o <- labelTads(ms.tads.o, "MS overlap", 7)
s.tads.o <- labelTads(s.tads.o, "S overlap", 9)
st.tads.o <- labelTads(st.tads.o, "ST overlap", 11)

l.tads.n <- labelTads(l.tads.n, "L non-overlap", 2)
ml.tads.n <- labelTads(ml.tads.n, "ML non-overlap", 4)
m.tads.n <- labelTads(m.tads.n, "M non-overlap", 6)
ms.tads.n <- labelTads(ms.tads.n, "MS non-overlap", 8)
s.tads.n <- labelTads(s.tads.n, "S non-overlap", 10)
st.tads.n <- labelTads(st.tads.n, "ST non-overlap", 12)

## Bind them together
tad.size <- rbind(l.tads.o, l.tads.n, ml.tads.o, ml.tads.n, m.tads.o, m.tads.n, ms.tads.o, ms.tads.n, s.tads.o, s.tads.n, st.tads.o, st.tads.n)

#5F4690,#1D6996,#38A6A5,#0F8554,#73AF48,#EDAD08,#E17C05,#CC503E,#94346E,#6F4070,#994E95,#666666
colors <- c(rep(c("#94346E","#994E95"), 12))
# Plot large tads #TopDom_TADs_bySize
p <- tad.size %>%
  mutate(name = fct_reorder(name, factor)) %>% 
  ggplot( aes(x=name, y=size, fill=name)) +
  geom_boxplot(notch=TRUE, notchwidth = 0.8, outlier.shape = NA) +
  scale_fill_manual(values=colors) +
  theme_ipsum() + stat_boxplot(geom = "errorbar", width = 0.15) + scale_y_continuous(limits = c(0, 2.6e6))+
  ggtitle("Aggregates by size within and outside LADs") + ylab("Base Pairs")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
        axis.text.y = element_text(angle = 0,size = 20))+ xlab("Aggregate Group")

p

## Get summary statistics
getSummary <- function(data, label){
  
  tad <- data %>% group_by(name) %>% 
    summarize(mean.size = mean(size), # calculates the mean
              sd = sd(size),   # calculates the standard deviation
              tad.n = n()) 
  tad$group <- rep(label, nrow(tad))
  
  return(tad)
  
}

tads.all.o <- getSummary(lamin.overlaps, "all overlap")
tads.all.n <- getSummary(non.lamins, "all non-overlap")
l.tads.o <- getSummary(l.tads.o, "L overlap")
ml.tads.o <- getSummary(ml.tads.o, "ML overlap")
m.tads.o <- getSummary(m.tads.o, "M overlap")
ms.tads.o <- getSummary(ms.tads.o, "MS overlap")
s.tads.o <- getSummary(s.tads.o, "S overlap")
st.tads.o <- getSummary(st.tads.o, "ST overlap")

l.tads.n <- getSummary(l.tads.n, "L non-overlap")
ml.tads.n <- getSummary(ml.tads.n, "ML non-overlap")
m.tads.n <- getSummary(m.tads.n, "M non-overlap")
ms.tads.n <- getSummary(ms.tads.n, "MS non-overlap")
s.tads.n <- getSummary(s.tads.n, "S non-overlap")
st.tads.n <- getSummary(st.tads.n, "ST non-overlap")

## Bind them together
number.topdom  <- rbind(tads.all.o,tads.all.n, l.tads.o, l.tads.n, ml.tads.o, ml.tads.n, m.tads.o, m.tads.n, ms.tads.o, ms.tads.n, s.tads.o, s.tads.n, st.tads.o, st.tads.n)
number.topdom$name <- rep(c("24hrAuxin", "untreated"), 14)


## Get ratio
mean <- mean(tad.size$size)
number.topdom$ratio <- (number.topdom$mean.size/mean)*3
number.topdom$sd <- (number.topdom$sd/mean)*3

#number.topdom <- number.topdom[-(1:2),]

plottheme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                 panel.grid.major = element_line(),
                 #panel.grid.minor = element_line(colour = "tomato", size=.25, linetype = "dashed"),
                 strip.background=element_blank(), axis.text.x=element_text(family="Arial", colour="black"), axis.title.x=element_text(face="bold", size=20,family="Arial", colour="black", vjust=-2),
                 axis.text.y=element_text(family="Arial", colour="black"), 
                 axis.title.y=element_text(face="bold", size=15, family="Arial", colour="black", vjust=2),
                 text = element_text(size=10, family="Arial"), 
                 plot.title=element_text(size=20, face="bold", family="Arial", color="tomato", hjust=0.0, vjust=10, lineheight=1.5),
                 plot.subtitle=element_text(size=15, face="bold", family="Arial", color="black", hjust=0.0, lineheight=1.5),
                 legend.title = element_text(size=12, color= "tomato",face="bold"), 
                 legend.text = element_text(size=10),
                 legend.key=element_rect(fill=NA),
                 axis.ticks=element_line(colour="black"))

number.topdom$ratio

## Plot number of TADs
p <- number.topdom %>% 
  ggplot(aes(x= group, y= tad.n, group=name)) + geom_point(aes(color=name),size=3) + 
  plottheme + 
  scale_y_continuous(limits=c(0,6000), breaks=c(0,500,1000,1500,2000,2500,3000,3500,4000, 5000, 6000))+
  scale_x_discrete( labels =c("Non-LAD TADs", "TADs in LADs", "L Overlapping TADs", "L Non-Overlapping TADs",
  "ML Overlapping TADs", "ML Non-Overlapping TADs", "M Overlapping TADs", "M Non-Overlapping TADs",
  "MS Overlapping TADs", "MS Non-Overlapping TADs", "S Overlapping TADs", "S Non-Overlapping TADs",
  "ST Overlapping TADs", "ST Non-Overlapping TADs"),  expand = c(0.025, 0.025))+
  scale_color_manual(values=c("#1571B1", "#259C25"))+ ylab("Aggregate N")+ xlab("Aggregate Group")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size = 13))+
  labs(title="Total Number of TADs", subtitle="Total Number of Aggregates", #subtitle="Top 1000 Variable Regions", 
       caption="")
p



