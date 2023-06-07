
# upregulated color: #D22B2B
# downregulated color: #0047AB

# Set directory variables here
#dir.dge.res.rsem <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Results"
#dir.plt.res <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Results/Plots"
#dir.meta.res <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Results/Metascape"
#dir.rsem <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Counts/Adli_Lamin_RNA-seq_RawData_RSEM_counts"

# data in backup drive ' My Passport'
dir.dge.res.rsem <- "/Volumes/My Passport/projects/Emily/Lamins_RNAseq/Results"
dir.plt.res <- "/Volumes/My Passport/projects/Emily/Lamins_RNAseq/Results/Plots"
dir.meta.res <- "/Volumes/My Passport/projects/Emily/Lamins_RNAseq/Results/Metascape"
dir.rsem <- "/Volumes/My Passport/projects/Emily/Lamins_RNAseq/Counts/Adli_Lamin_RNA-seq_RawData_RSEM_counts"


#########################################################################################
# Set Directory, pass in  counts, filter counts for DGE analysis
#########################################################################################

# import packages
library("DESeq2")
library("ggplot2")
library("dplyr")
library("magrittr")
library("tximport")


list.files(dir.rsem)
setwd(dir.rsem)

sample <- read.table(file.path(dir.rsem, "metadata.csv"), sep= ",", header= TRUE)
newcolnames <-sample$Alias

## Import files
files <- file.path(dir.rsem, paste0(sample$ID))
files

## Imports RSEM TPM files into a list
cts <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

## Move col names from sample info to counts matrix
colnames(cts$counts)<-newcolnames
colnames(cts$counts)

#Replace any transcripts with length 0 with a 1
cts$length[cts$length == 0] <- 1
head(cts$counts)

# logical test to see if ^ is true now
identical(sample$Alias, colnames(cts$counts))

# Store the file names in the sample info matrix as row names
row.names(sample) <- sample$Alias
head(sample, n=10)

dput(as.character(unique(sample$Conditions)))

#factor out dependencies
sample$Conditions <- factor(sample$Conditions, levels = c( "0hr","12hr","48hr","144hr"))
levels(sample$Conditions)

sample$Replicate <- factor(sample$Replicate, levels = c(1, 2))
levels(sample$Replicate)

sample$Lane <- factor(sample$Lane, levels = c(1, 2, 3, 4))
levels(sample$Lane)

sample$Label <- factor(sample$Label, levels = dput(as.character(unique(sample$Label))))
levels(sample$Label)

#########################################################################################
# Create DESeq data object and pply transformations to data for pre-DGE QC
#########################################################################################

# Create DDS object
dds <- DESeqDataSetFromTximport(cts, colData = sample, design = ~ Conditions)

# Collapse lanes together into single replicates (described here: https://www.biostars.org/p/260838/ and here: https://rdrr.io/bioc/DESeq2/man/collapseReplicates.html)
dds <- collapseReplicates(dds, dds$Label, dds$Lane)
dds$runsCollapsed

# check dds object
summary(dds)

# Pre-filtering the dataset
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 0, ]
nrow(dds)

# Variance stabilizing transformation (VST) and rlog to make data more homoskedastic
vst <- vst(dds, blind = FALSE)
#head(assay(vsd), 3)

head(assay(vst))

# regularized log transformation
rld <- rlog(dds, blind = FALSE)
#head(assay(rld),3)

# normalized counts transformation
ntd <-normTransform(dds, f=log2, pc=1)

##Add another column which defines the group that each sample is in
dds$group <- factor(paste0(dds$Conditions))

# Find size factors to normalize differences in sequencing depth across samples
##an DESeq2 algorithm for normalization of count data, conceptually comparable to CPM or TPM
dds <- estimateSizeFactors(dds)

#########################################################################################
# Base R PCA Plots
#########################################################################################

View(sample)

sample2 <- !is.na(sample$Replicate)

sample <- sample[sample2,]

PCAtheme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                panel.grid.major = element_line(),
                panel.grid.minor = element_line(colour = "tomato", size=.25, linetype = "dashed"),
                strip.background=element_blank(), axis.text.x=element_text(family="Arial", colour="dark blue"), axis.title.x=element_text(face="bold", size=15,family="Arial", colour="dark blue", vjust=-2),
                axis.text.y=element_text(family="Arial", colour="dark blue"), 
                axis.title.y=element_text(face="bold", size=15, family="Arial", colour="dark blue", vjust=2),
                text = element_text(size=10, family="Arial"), 
                plot.title=element_text(size=20, face="bold", family="Arial Rounded MT", color="tomato", hjust=0.0, vjust=10, lineheight=1.5),
                plot.subtitle=element_text(size=15, face="bold", family="Arial Rounded MT", color="black", hjust=0.0, lineheight=1.5),
                legend.title = element_text(size=12, color= "tomato",face="bold"), 
                legend.text = element_text(size=10),
                legend.key=element_rect(fill=NA),
                axis.ticks=element_line(colour="black"),
                #plot.margin=unit(c(1,1,1,1),"line")
)

library(matrixStats)
library(ggrepel)

## log normalize count table (only takes matrix as input)
#cts.norm <- as.matrix(round(cts$counts))
cts.norm <- as.matrix(round(counts(dds)))
scaled = DESeq2::rlog(cts.norm, blind = TRUE)

## First we transform the count data to log2 scale using the 'regularized log' transformation.
#scaled = scale(cts.norm)

## Next we calculate the per-row variance in the transformed table.
vars = rowVars(scaled)

## Then we re-order the rows so that the most variable rows are on top.
scaled = scaled[order(vars, decreasing = TRUE),]

## Here we perform PCA using `prcomp`. We transpose (`t`) the matrix and only operate on the top 1000 most variable regions. 
#pca = prcomp(t(scaled[1:1000,]))

## To perform PCA on the entire data set
pca = prcomp(t(scaled))

## We calculate the percent variance explained by each principal component by extracting these values from the summary.
percentVar = round(100 * summary(pca)$importance[2,])
cumuVar = round(100 * summary(pca)$importance[3,])

## Look at top variance
head(pca$x, n=12)

## Plot top PCAs for PCA1 and 2
prTab = data.frame(pca$x, colData(dds))

prTab$Conditions

prTab %>% #PCA1_2_AllSamples_withlabels
  ggplot(aes(PC1, PC2, color = Conditions)) +
  geom_point(size = 3, alpha = 0.75) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  scale_color_manual(values = c("steelblue3", "#d44842", "#9f2a63", "#65156e")) + 
  coord_fixed() + geom_text_repel(segment.color = 'transparent', aes(label = NA), size = 3, family = "Arial Narrow", hjust = 0)+
  labs(title="Principle Component Analysis", subtitle="PC1 & PC2: All Samples", #subtitle="Top 1000 Variable Regions", 
       caption="")+PCAtheme

## Plot top PCAs for PCA2 and 3 
prTab %>% #PCA2_3_AllSamples_NOlabels
  ggplot(aes(PC2, PC3, color = Conditions)) +
  geom_point(size = 3, alpha = 0.75) + 
  xlab(paste0("PC2: ", percentVar[2], "% variance")) + ylab(paste0("PC3: ", percentVar[3], "% variance")) + 
  scale_color_manual(values = c("steelblue3", "#d44842", "#9f2a63", "#65156e")) + 
  coord_fixed() + geom_text_repel(aes(label = Label), size = 3, family = "Arial Narrow", hjust = 0)+
  labs(title="Principle Component Analysis", subtitle="PC2 & PC3: All Samples", #subtitle="Top 1000 Variable Regions", 
       caption="")+PCAtheme

## Plot cumulative variance
b = barplot(percentVar, horiz = FALSE, las = 1, ylab = 'Variance (%)', ylim = c(0,105))
lines(b, cumuVar, type = 'b', col = 'red')

## Plot Eigens
eigenvalues = pca$sdev^2
plot(c(1:32), eigenvalues, xlab = 'Principal Components', ylab = 'eigenvalue', las = 1)

## determing which PCAs to keep based on: 
which(eigenvalues > mean(eigenvalues))

#########################################################################################
# Generate PCA plot and Heatmap for QC of count data
#########################################################################################

# set a new directory, to save CSV results and plots in
dir.plt.res <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Results/Plots"
list.files(dir.plt.res)
setwd(dir.plt.res)

# optional QC method: Plot dispersion estimates
dds <- estimateDispersions(dds)
plotDispEsts(dds)

# PCA Plot: function from DESeq2 package, for plotting VSD transform
##calculates the PCA and stores it in pcaData variable
pcaData <- plotPCA(vst, intgroup = c( "Conditions"), returnData = TRUE) ##use 'ntop=' argument to use <500 most variable genes

##variance variable required for GGplot to plot custom PCA data
percentVar <- round(100 * attr(pcaData, "percentVar"))

# stylize the GGplot using the 'theme' function
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
             panel.grid.major = element_line(),
             panel.grid.minor = element_line(colour = "tomato", 
                                             size=.25, 
                                             linetype = "dashed"),
             strip.background=element_blank(),
             axis.text.x=element_text(family="Arial", 
                                      colour="dark blue"),
             axis.title.x=element_text(face="bold", 
                                       size=15,
                                       family="Arial", 
                                       colour="dark blue",
                                       vjust=-2),
             axis.text.y=element_text(family="Arial", 
                                      colour="dark blue"), 
             axis.title.y=element_text(face="bold", 
                                       size=15,
                                       family="Arial", 
                                       colour="dark blue",
                                       vjust=2),
             text = element_text(size=10, family="Arial"), 
             plot.title=element_text(size=20, 
                                     face="bold", 
                                     family="Arial Rounded MT",
                                     color="tomato",
                                     hjust=0.0,
                                     vjust=10,
                                     lineheight=1.5),
             plot.subtitle=element_text(size=15, 
                                        face="bold", 
                                        family="Arial Rounded MT",
                                        color="black",
                                        hjust=0.0,
                                        lineheight=1.5),
             legend.title = element_text(size=12, 
                                         color= "tomato",
                                         face="bold"), 
             legend.text = element_text(size=10),
             legend.key=element_rect(fill=NA),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))

# Stores the PCA plot in pca.plot variable
pca.plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Conditions)) + 
  geom_point( size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  #xlim(-15,20) + ylim(-15,15)
  theme + scale_colour_brewer(palette = "Paired") + 
  labs(title="Principle Component Analysis", subtitle="For All Samples", 
       caption="")

# Generates PCA plot with the following additional parameters
pca.plot + scale_x_continuous(breaks=seq(-10, 20, 5), 
                              labels = sprintf("%1.0f%%", seq(-10, 20, 5))) + 
  scale_y_continuous(breaks=seq(-6, 6, 1),  labels = sprintf("%1.0f%%", seq(-6, 6, 1)))

library(pheatmap)
library(viridis)
library(dendsort)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

# Hierachal clustering map plot for sample QC
##Extract the rlog matrix from the object
rld_mat <- assay(rld)

# Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R functions
head(rld_cor)

# Plot heatmap using the correlation matrix and the metadata object
# organize hierachal clusters
mat <- rld_cor
mat_cluster_cols <- hclust(dist(t(mat)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- sort_hclust(hclust(dist(mat)))
length(mat)
colors <- colorRampPalette(c("#0047AB", "white", "#D22B2B"))(64)

# plot heat map
p<-pheatmap(mat, 
            color = colors,
            fontsize = 10,
            fontsize_row = 6, 
            fontsize_col = 6,
            cellwidth = 8,
            cellheight = 8,
            cluster_cols= mat_cluster_cols,
            cluster_rows= mat_cluster_rows,
            main = "Sample Correlation",
            annotation_legend = TRUE,
            scale="row")
p

#########################################################################################
# Differential expression analysis and 'results' matrices
#########################################################################################

## Generate DGE DESeqDataSet data frame
dds <- DESeq(dds)

## Which contrasts to use?
resultsNames(dds)

dput(as.character(unique(sample$Conditions)))

## Generate a results matrix for Treatment vs. Control with an alpha of 0.05
res.12hr <- results(dds, contrast = c("Conditions", "12hr", "0hr"), alpha = 0.05)
res.48hr <- results(dds, contrast = c("Conditions", "48hr", "0hr"), alpha = 0.05)
res.144hr <- results(dds, contrast = c("Conditions", "144hr", "0hr"), alpha = 0.05)

## Wash off results matrix
res.12hrv144hr <- results(dds, contrast = c("Conditions", "12hr", "144hr"), alpha = 0.05)
res.48hrv144hr <- results(dds, contrast = c("Conditions", "48hr", "144hr"), alpha = 0.05)

#########################################################################################
# Annotating the results using Annotation Hub and reading results out
#########################################################################################

setwd(dir.dge.res.rsem)

Annotate_results <- function(results, filename, db ="OrgDb",org="Homo sapiens",
                             key="ENSEMBL",col="SYMBOL"  ){
  
  require("AnnotationHub")
  ah <- AnnotationHub() ## Connect to Annotation Hub
  qcall <- query(ah, c(db,org)) ## Make a query to your organism/database
  ano <- ah[[names(qcall)]] ## Annotation object
  
  ## Remove .## version of gene in ensembl ID
  rownames(results) <- gsub("\\..*","",rownames(results))
  
  ## Determine how many of the differentially expressed genes can be mapped to a gene name
  table(rownames(results) %in% keys(ano, key))
  
  ## Map gene symbols onto the results matrix using Annotation Hub function mapIDS()
  results$symbol <- mapIds(ano, rownames(results), column = col, keytype = key)
  
  results$padj <- ifelse(is.na(results$padj), 1, results$padj) ## NA filtering, padj
  results$pvalue <- ifelse(is.na(results$pvalue), 1, results$pvalue) ## NA filtering, pv
  results <- na.omit(results)
  #results$symbol <- ifelse(is.na(results$symbol), "Unknown", results$symbol) ## NA filtering, symbol
  
  ## Order results and save Files as CSV function
  results <- results[order(results$padj),]
  results <- results[c("symbol", "baseMean","log2FoldChange", "lfcSE", "stat", "pvalue", "padj" )]
  #write.csv(as.data.frame(results), file=paste(filename, "csv", sep="."))
  
  return(results)
  
}

## Annotate 10nm vs control
res.12hr <- Annotate_results(res.12hr, "laminKO_12hr_vs_Control_DGE")
res.48hr <- Annotate_results(res.48hr, "laminKO_48hr_vs_Control_DGE")
res.144hr <- Annotate_results(res.144hr, "laminKO_144hr_vs_Control_DGE")
res.12hrv144hr <- Annotate_results(res.12hrv144hr, "laminKO_12hr_vs_Washout_DGE")
res.48hrv144hr <- Annotate_results(res.48hrv144hr, "laminKO_48hr_vs_Washout_DGE")

#########################################################################################
# MA Plot using functions written for DESeq2
#########################################################################################

# upregulated color: #D22B2B
# downregulated color: #0047AB

## ggPlot that takes a DESeq results object and a title and generates an MA plot
ggPlot_MA = function (results, title, filename, coef=coef, pv=0.01, dot.trans=0.5, x.label="log10(mean counts)", 
                      y.label="log2(fold change) treatment/control", de.col="#0047AB", f.col="grey", dir.plt= dir.plt.res
                        ){
  require("ggplot2")
  require("apeglm")
  
  MAtheme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(),
                   panel.grid.minor = element_line(colour = "black", size=.25, linetype = "dashed"),
                   strip.background=element_blank(), axis.text.x=element_text(family="Arial", colour="black"),
                   axis.title.x=element_text(face="bold", size=15, family="Arial", colour="black", vjust=-2),
                   axis.text.y=element_text(family="Arial", colour="black"), 
                   axis.title.y=element_text(face="bold", size=15,family="Arial", colour="black", vjust=2),
                   text = element_text(size=10, family="Arial"), 
                   plot.title=element_text(size=20, face="bold", family="Arial Rounded MT", color="black", hjust=0.0, vjust=10, lineheight=1.5),
                   plot.subtitle=element_text(size=15, face="bold", family="Arial Rounded MT", color="black", hjust=0.0, lineheight=1.5),
                   legend.title = element_text(size=12, color= "black", face="bold"), 
                   legend.text = element_text(size=10),
                   legend.key=element_rect(fill=NA),
                   axis.ticks=element_line(colour="black"))
  
  res <- lfcShrink(dds, coef=coef, type="apeglm")          
  
  res <-  as.data.frame(res) %>% 
    mutate(sig = pvalue < pv) %>% 
    ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = sig)) + 
    geom_point(alpha = dot.trans) +
    coord_cartesian(ylim = c(-10,10)) +
    scale_color_manual(values = c("TRUE" = de.col, "FALSE" =  f.col)) +
    geom_hline(yintercept = c(1,-1), linetype = 'dotted', color = 'darkred') +
    MAtheme+
    labs(title = title,
         x = x.label,
         y = y.label)
  ggsave(paste(dir.plt.res,"/",filename, ".png", sep=""), dpi=200)
  
  return(res)
  
}
resultsNames(dds)
### Call function for results
ggPlot_MA(res.12hr, "Lamin KO at 12 Hours","12hr_vs_cMAplot1", "Conditions_12hr_vs_0hr")
ggPlot_MA(res.48hr, "Lamin KO at 48 Hours","48hr_vs_cMAplot1", "Conditions_48hr_vs_0hr")
ggPlot_MA(res.144hr, "Lamin KO at 144 Hours","144hr_vs_cMAplot1", "Conditions_144hr_vs_0hr")

## Need to change the contrasts to look at these two
#ggPlot_MA(res.12hrv144hr, "Lamin KO, 12hrs vs 144hrs Hours","12hr_vs_144hrMAplot1", "Conditions_12hr_vs_144hr")
#ggPlot_MA(res.48hrv144hr, "Lamin KO, 48hrs vs 144hrs Hours","48hr_vs_144hrMAplot1", "Conditions_48hr_vs_144hr")

# upregulated color: #D22B2B
# downregulated color: #0047AB

## A second function that calls on a function from ggpubr to generate an MA plot from DESeq results object
plot_MA <- function(coef, title, results, filename, dncol="#0047AB", upcol="#D22B2B", dir.plt=
                      "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Bonini/bulk_RNAseq_2021/results/plots_rsem"){
  
  require("apeglm")
  require("ggpubr")
  
  ## Remove .## version of gene in ensembl ID
  rownames(results) <- gsub("\\..*","",rownames(results))
  rownames(dds) <- gsub("\\..*","",rownames(dds))
  
  # shrink noisy log2 fold change estimates for coefficients: DESeq2 function
  res <- lfcShrink(dds, coef=coef, type="apeglm")
  res <- res[rownames(res)%in% rownames(results),]
  
  g1<-ggmaplot(res, 
               main = title,
               fdr = 0.10, fc = 1, size = 2,
               palette = c(upcol, dncol, "darkgray"),
               legend = "top", top = -5,
               genenames = results$symbol,
               font.label = c("bold", 5, "black"),
               font.legend = c("bold", 20, "black"),
               font.main = c("bold", 30,"black", "Arial"),
               ggtheme = ggplot2::theme_minimal()) + theme(text = element_text(size = 20, color= "black", family = "Arial Bold"))
  
  ggsave(paste(dir.plt.res,"/",filename, ".png", sep=""), dpi=300)
  return(g1)
  
}

## Find coefficients for following function
resultsNames(dds)

### Call function for results
plot_MA("Conditions_12hr_vs_0hr", "Lamin KO at 12 Hours", res.12hr,"12hr_vs_cMAplot2")
plot_MA("Conditions_48hr_vs_0hr", "Lamin KO at 48 Hours", res.48hr,"48hr_vs_cMAplot2")
plot_MA("Conditions_144hr_vs_0hr", "Lamin KO at 144 Hours", res.144hr,"144hr_vs_cMAplot2")

#########################################################################################
# Base R Volacano plot of DESeq results
#########################################################################################

# upregulated color: #D22B2B
# downregulated color: #0047AB

# lfc.up = positive integer (ex. 0.58, 1, 2)
# lfc.dwm = negative integer (ex. -0.58, -1, -2)
# pvadj = positive integer (ex. 0.01, 0.05. 0.1)

DE_Vol_Plot <- function(results, plt.title, filename, lfc.up=1, lfc.dwn=-1,pvadj=0.01,
                        y.limits=c(0,150), x.limits=c(-6.5,6.5),dncol="#0047AB", upcol="#D22B2B", dir.plt=
                          "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Results/Plots") {
  require(ggrepel)
  require(ggplot2)
  
  ## Generate data frame
  genes.rna <- data.frame(symbol = results$symbol,lfc = results$log2FoldChange, padj=results$padj, base=results$baseMean)
  
  resultsDN <- as.data.frame(genes.rna[which(genes.rna$lfc < lfc.dwn),])
  resultsDN <- as.data.frame(resultsDN[which(resultsDN$padj < pvadj),])
  length(rownames(resultsDN)) ## Gives count of down regulated
  
  resultsUP <- as.data.frame(genes.rna[which(genes.rna$lfc > lfc.up),])
  resultsUP <- as.data.frame(resultsUP[which(resultsUP$padj < pvadj),])
  length(rownames(resultsUP)) ## Gives count of up regulated  
  
  ## add a column of NAs
  genes.rna$diffexpressed <- "NO"
  ## If log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  genes.rna$diffexpressed[genes.rna$lfc > lfc.up & genes.rna$padj < pvadj] <- "Upregulated"
  ## If log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  genes.rna$diffexpressed[genes.rna$lfc < lfc.dwn & genes.rna$padj < pvadj] <- "Downregulated"
  
  ## Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  genes.rna$delabel <- NA
  genes.rna$delabel[genes.rna$diffexpressed != "NO"] <- genes.rna$symbol[genes.rna$diffexpressed != "NO"]
  
  ## Jitter labeling
  pos <- position_jitter(width = 0.2, seed = 1)
  
  vplot <- genes.rna %>% mutate(mean.intensity = base/max(base)) %>% 
    ggplot(aes(x=lfc, y=-log10(padj), col=diffexpressed, size = mean.intensity, label=delabel
    )) +
    geom_point(position = pos, alpha=0.7)+
    scale_color_manual(name="Upregulated", values=c(dncol, "lightgrey", upcol))+
    geom_vline(xintercept=c(lfc.dwn, lfc.up), col="red", lty = "dashed") +
    geom_hline(yintercept=-log10(pvadj), col="red", lty = "dashed")+ 
    geom_text_repel(segment.color = 'transparent', position = pos, size = 3, xlim  = c(-6,6)) +
    labs(title=plt.title, subtitle=paste(paste("Downregulated:",length(rownames(resultsDN))), paste("| Upregulated:",length(rownames(resultsUP))))) +
    scale_y_continuous(name="-Log 10 Padj", limits=y.limits, breaks=c(0,50,100,150,200,250, 300, 350))+
    scale_x_continuous(name="Log 2 Fold Change", limits=x.limits, breaks=c(-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7))+
    theme( axis.line = element_line(color = "black", size = 1, linetype = "solid"),
           axis.text.y = element_text( color="black", size=8, family = "Arial Narrow"),
           axis.text.x = element_text(angle = 0, color="black", size=8, family = "Arial Narrow"),
           plot.title = element_text(hjust = 0.5, color="black", face = "bold", size=15, family = "Arial Rounded MT Bold"),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank()) 
  
  #ggsave(paste(dir.plt,"/",filename, ".png", sep=""), dpi=300)
  return(vplot)
  
}

### Call function for results
DE_Vol_Plot(res.12hr,"Lamin KO at 12 Hours: Differentially Expressed Genes", "res.12hr.volplot" )
DE_Vol_Plot(res.48hr,"Lamin KO at 48 Hours: Differentially Expressed Genes", "res.48hr.volplot" )
DE_Vol_Plot(res.144hr,"Lamin KO at 144 Hours: Differentially Expressed Genes", "res.144hr.volplot" )

### No labels (must comment out geom_text_repel() and label = ..)
DE_Vol_Plot(res.12hr,"Lamin KO at 12 Hours: Differentially Expressed Genes", "res.12hr.volplot_nolabels" )
DE_Vol_Plot(res.48hr,"Lamin KO at 48 Hours: Differentially Expressed Genes", "res.48hr.volplot_nolabels" )
DE_Vol_Plot(res.144hr,"Lamin KO at 144 Hours: Differentially Expressed Genes", "res.144hr.volplot_nolabels" )

#########################################################################################
# Generate heatmap of top Differential results 
#########################################################################################

#https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/
## As gene.n (number of tiles) increases, tile.width must increase, to keep image proportional
Heat_Map_Counts <- function(results, title, subtitle, caption, sample.n = 1:8, padj = 0.05, lfc.abs = 1, gene.n=gene.n, 
                            n_colors=n_colors, breaks=breaks, labels = as.character(labels), tile.width=tile.width,
                            axis.text.x = element_text( color="black", size=7, family = "Arial Narrow", angle = 90, vjust = 0.5, hjust=1),
                            axis.text.y = axis.text.y) {
  
  require(ggplot2) # ggplot() for plotting
  require(dplyr) # data reformatting
  require(tidyr) # data reformatting
  require(stringr) # string manipulation
  require(viridis) # colors package
  require(DESeq2) # DESeq2
  
  
  ## Regularized log transformation, count extraction, and normalization
  rld <- DESeq2::rlog(dds, blind = FALSE)
  assay <-t(scale(t(assay(rld)[,sample.n])))
  
  ## Remove .## version of gene in ensembl ID
  rownames(assay) <- gsub("\\..*","",rownames(assay))
  
  ## Organize genes by significance and extract
  sig.genes <- rownames(results[results$padj <= padj & abs(results$log2FoldChange) > abs(lfc.abs),])
  results <- as.data.frame(results[rownames(results) %in% sig.genes,])
  results <- results %>% 
    dplyr::distinct(symbol, .keep_all = T)
  
  ## Order genes by LFC, pAdj, and extract N genes for heatmap
  #results <- results[order(results$padj, decreasing = F),]
  #results <- as.data.frame(results[order(results$log2FoldChange, decreasing = FALSE),])[gene.n,] # order the results DOWNREG
  #results <- as.data.frame(results[order(results$log2FoldChange, decreasing = TRUE),])[gene.n,] # order the results UPREG
  results <- as.data.frame(results)[gene.n,] # unordered results for clustering
  
  ## Prepare data frame for heatmapping 
  top.heat <- as.data.frame(assay[rownames(assay)%in%rownames(results),]) %>% 
    dplyr::mutate(label = rownames(results)) %>%
    dplyr::mutate(symbol = results$symbol) %>% # make rownames a new column
    tidyr::gather(key="samples",value="value", -symbol, -label) %>% # convert data to long format
    stats::setNames(c("label", "symbol","samples","value")) %>% # rename columns
    dplyr::mutate(regions=factor(symbol)) %>% # convert factors
    dplyr::mutate(samples=factor(samples)) %>%  # convert factors
    dplyr::mutate(value=as.numeric((value))) #convert value to numeric (also converts '-' to NA, gives a warning)
  
  ## Compute a distance calculation on both dimensions of the matrix for clustering
  #distance.gene <- dist(as.data.frame(assay[rownames(results),])) ## average distance method
  distance.gene <- as.dist((1 - cor(t(as.data.frame(assay[rownames(results),]))))/2) ## correlation method
  distance.sample <- dist(t(as.data.frame(assay[rownames(results),])))
  
  # Cluster based on the distance calculations
  #cluster.gene <- hclust(distance.gene, method="average") ## average distance method
  cluster.gene=hclust(distance.gene) ## correlation method
  cluster.sample <- hclust(distance.sample, method="average")
  
  # Re-factor sample and genes for ggplot2
  top.heat$label <- factor(top.heat$label, levels=cluster.gene$labels[cluster.gene$order])
  top.heat$samples <- factor(top.heat$samples, levels=cluster.sample$labels[cluster.sample$order])
  top.heat$regions <- factor(top.heat$regions, levels=top.heat$regions[cluster.gene$order])
  
  # upregulated color: #D22B2B
  # downregulated color: #0047AB
  ## Apply viridis_pal function to generate a palette
  #palette <- colorRampPalette(c("#D22B2B", "#0047AB"))(n_colors)
  palette <- colorRampPalette(c("#0047AB", "white", "#D22B2B"))(n_colors)
  
  ## Bin heat map values using cut() and factor
  top.heat <- top.heat %>% # create a new variable from count
    mutate(countfactor=cut(value,breaks=c(breaks), labels=c(labels))) %>% # change level order
    mutate(countfactor=factor(as.character(countfactor),levels=rev(levels(countfactor)))) %>% # factor bins
    arrange(factor(top.heat$label, levels=cluster.gene$labels[cluster.gene$order])) %>% # factor labels
    arrange(factor(top.heat$regions, levels=top.heat$regions[cluster.gene$order])) # factor regions
  
  
  ## Plot heat map of top 1000
  heat.map <- top.heat %>% 
    ggplot(aes(x=samples,y=regions,fill=countfactor))+
    geom_tile(size=1)+coord_fixed(ratio=tile.width)+
    scale_fill_manual(values=c(palette),na.value = "grey90")+ ## color scale 2
    labs(title=title, subtitle=subtitle, caption=caption) +
    scale_y_discrete(name=NULL)+
    scale_x_discrete(name=NULL)+
    theme(axis.text.y = axis.text.y ,
          axis.text.x = axis.text.x ,
          plot.title = element_text(hjust = 0.5, color="black", face = "bold", size=15,
                                    family = "Arial Rounded MT Bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  
  return(heat.map)
  
}

## Call function for N=1000 and N=100
Heat_Map_Counts(res.12hr, "Top 100 Differentially Expressed Genes", "", "", gene.n=1:100, tile.width = 0.4, breaks=seq(from=-1.75, to=1.75, by=0.25),
                n_colors=14, labels=seq(from=-1.75, to=1.5, by=0.25), axis.text.y = element_text( color="black", size=4.5, family = "Arial Narrow"))
Heat_Map_Counts(res.12hr, "Top 1000 Differentially Expressed Genes", "", "", gene.n=1:673, tile.width = .04, breaks=seq(from=-2, to=1.75, by=0.25),
                n_colors=15, labels=seq(from=-2, to=1.5, by=0.25),axis.text.y =element_blank())
Heat_Map_Counts(res.12hr, "All Differentially Expressed Genes", "", "", gene.n=1:2039, tile.width = .01, breaks=seq(from=-2.25, to=2.5, by=0.25),
                n_colors=18, labels=seq(from=-2.25, to=2.25, by=0.25),axis.text.y =element_blank())


#####################################################################################
# DE gene pheatmap                                         
#####################################################################################

# upregulated color: #D22B2B
# downregulated color: #0047AB

## Set colors for annotations
q_colors =  12 # for no particular reason
v_colors =  viridis(q_colors, option = "E") # E= Civis
v_colors

Gene_pHeat_Map <- function(results, title, gene.n = gene.n, padj = 0.05, lfc.abs = 0.58, sample.n = 1:8, c.width = c.width, c.height = c.height, f.size.row=4,
                           f.size.col=10, f.size=10, cols=colorRampPalette(c("#0047AB", "white", "#D22B2B"))(30), 
                           meta.annotation=c("Conditions", "Replicate"),
                           anno.colors = list(Replicate = c("1" = "#000004FF", "2" = "#5E626EFF"), 
                                             Conditions = c("0hr" = "#00204DFF", "12hr" = "#00336FFF", "48hr" = "#39486BFF", "144hr" = "#575C6DFF"))
                           ) {
  require(pheatmap)
  require(viridis)
  require(dendsort)
  require(DESeq2)
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  
  ## Remove .## version of gene in ensembl ID
  #rownames(results) <- gsub("\\..*","",rownames(results))
  #rownames(dds) <- gsub("\\..*","",rownames(dds))
  
  
  ## Regularized log transformation, count extraction, and normalization
  rld <- DESeq2::rlog(dds, blind = FALSE)
  assay <-t(scale(t(assay(rld)[,sample.n])))
  
  ## Remove .## version of gene in ensembl ID
  rownames(assay) <- gsub("\\..*","",rownames(assay))
  
  ## Organize genes by significance and extract
  sig.genes <- rownames(results[results$padj <= padj & abs(results$log2FoldChange) > abs(lfc.abs),])
  results <- as.data.frame(results[rownames(results) %in% sig.genes,])
  results <- results %>% 
    dplyr::distinct(symbol, .keep_all = T)
  
  ## Order genes by  pAdj, and extract N genes for heatmap
  results <- results[order(results$padj, decreasing = F),]
  print(results)
  #results <- as.data.frame(results[order(results$log2FoldChange, decreasing = FALSE),])[gene.n,]
  results <- as.data.frame(results)[gene.n,]
  results
  
  ## Prepare data frame for heatmapping 
  top.heat <- as.data.frame(assay[rownames(assay)%in%rownames(results),]) 
  top.heat
  ## Compute a distance calculation on both dimensions of the matrix for clustering
  #distance.gene <- dist(as.data.frame(assay[rownames(results),])) ## average distance method
  distance.gene <- as.dist((1 - cor(t(as.data.frame(assay[rownames(results),]))))/2) ## correlation method
  distance.sample <- dist(t(as.data.frame(assay[rownames(results),])))
  
  # Cluster based on the distance calculations
  #cluster.gene <- hclust(distance.gene, method="average") ## average distance method
  cluster.gene=hclust(distance.gene) ## correlation method
  cluster.sample <- hclust(distance.sample, method="average")
  
  ## Annotation information in data frame object
  dfcol <- as.data.frame(colData(dds)[,meta.annotation])

  # plot heat map
  p<-pheatmap(top.heat, 
              color = cols,
              fontsize = f.size,
              fontsize_row = f.size.row, 
              fontsize_col = f.size.col,
              cellwidth = c.width,
              cellheight = c.height,
              cluster_cols= F,
              #cluster_cols= cluster.sample,
              cluster_rows= cluster.gene,
              main = title,
              annotation_colors = anno.colors,
              labels_row = results$symbol,
              show_rownames = T,
              border_color = NA,
              annotation_legend = TRUE,
              annotation_col = dfcol
              )
 
  return(p)
  
}

Gene_pHeat_Map(res.12hr.ov, "", gene.n=1:100, c.width = 7, c.height = 4) #12hr_inLADS_heatmap_top100
Gene_pHeat_Map(res.12hr.non, "", gene.n=1:100, c.width = 7, c.height = 4) #12hr_outsideLADS_heatmap_top100

## Use these for new plots
Gene_pHeat_Map(res.48hr.ov, "LAD overlapping", gene.n=1:100, c.width = 9, c.height = 6) #48hr_inLADS_heatmap_top100
Gene_pHeat_Map(res.48hr.non, "LAD non-overlapping", gene.n=1:100, c.width = 9, c.height = 6) #48hr_outsideLADS_heatmap_top100
##

Gene_pHeat_Map(res.144hr.ov, "", gene.n=1:100, c.width = 7, c.height = 4) #144hr_inLADS_heatmap_top100
Gene_pHeat_Map(res.144hr.non, "", gene.n=1:100, c.width = 7, c.height = 4) #144hr_outsideLADS_heatmap_top100

res <- res.144hr.non
res <-res%>%
  filter(abs(log2FoldChange) > 0.58 )%>%
  filter(padj < 0.05 )

Gene_pHeat_Map(res.12hr.ov, "", gene.n=1:500, c.width = 7, c.height = 0.75) #12hr_inLADS_heatmap_top500
Gene_pHeat_Map(res.12hr.non, "", gene.n=1:500, c.width = 7, c.height = 0.75) #12hr_outsideLADS_heatmap_top500

Gene_pHeat_Map(res.48hr.ov, "", gene.n=1:500, c.width = 7, c.height = 0.75) #48hr_inLADS_heatmap_top500
Gene_pHeat_Map(res.48hr.non, "", gene.n=1:500, c.width = 7, c.height = 0.75) #48hr_outsideLADS_heatmap_top500

Gene_pHeat_Map(res.144hr.ov, "", gene.n=1:400, c.width = 7, c.height = 0.75) #144hr_inLADS_heatmap_top400
Gene_pHeat_Map(res.144hr.non, "", gene.n=1:500, c.width = 7, c.height = 0.75) #144hr_outsideLADS_heatmap_top500

rm(assay)

#####################################################################################
# DE gene pheatmap                                         
#####################################################################################

# upregulated color: #D22B2B
# downregulated color: #0047AB

## Set colors for annotations
q_colors =  12 # for no particular reason
v_colors =  viridis(q_colors, option = "E") # E= Civis
v_colors

Gene_pHeat_Map <- function(results, title, padj = 0.01, lfc.abs = 1, sample.n = 1:8, c.width = c.width, c.height = c.height, f.size.row=5,
                           f.size.col=10, f.size=10, cols=colorRampPalette(c("#0047AB", "white", "#D22B2B"))(30), 
                           meta.annotation=c("Conditions", "Replicate"),
                           anno.colors = list(Replicate = c("1" = "#000004FF", "2" = "#5E626EFF"), 
                                              Conditions = c("0hr" = "#00204DFF", "12hr" = "#00336FFF", "48hr" = "#39486BFF", "144hr" = "#575C6DFF"))
) {
  require(pheatmap)
  require(viridis)
  require(dendsort)
  require(DESeq2)
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  
  ## Remove .## version of gene in ensembl ID
  #rownames(results) <- gsub("\\..*","",rownames(results))
  #rownames(dds) <- gsub("\\..*","",rownames(dds))
  
  
  ## Regularized log transformation, count extraction, and normalization
  rld <- DESeq2::rlog(dds, blind = FALSE)
  assay <-t(scale(t(assay(rld)[,sample.n])))
  
  ## Remove .## version of gene in ensembl ID
  rownames(assay) <- gsub("\\..*","",rownames(assay))
  
  ## ----
  res <- data.frame(res.48hr)
  sig.genes <- rownames(res[res$padj <= padj & abs(res$log2FoldChange) > abs(lfc.abs),])
  res <- as.data.frame(res[rownames(results) %in% sig.genes,])
  res <- res %>% 
    dplyr::distinct(symbol, .keep_all = T)
  
  res <- res[order(res$padj, decreasing = F),]
  
  res <- as.data.frame(res)[1:100,]
  res <-rownames(res)
  print(res)
  
  ## ----
  
  ## Organize genes by significance and extract
  sig.genes <- rownames(results[results$padj <= padj & abs(results$log2FoldChange) > abs(lfc.abs),])
  results <- as.data.frame(results[rownames(results) %in% sig.genes,])
  results <- results %>% 
    dplyr::distinct(symbol, .keep_all = T)
  
  ## Order genes by  pAdj, and extract N genes for heatmap
  results <- results[order(results$padj, decreasing = F),]
  #results <- as.data.frame(results[order(results$log2FoldChange, decreasing = FALSE),])[gene.n,]
  results <- as.data.frame(results[rownames(results) %in% res,])
  print(results)
  
  ## Prepare data frame for heatmapping 
  top.heat <- as.data.frame(assay[rownames(assay)%in%rownames(results),]) 
  top.heat
  ## Compute a distance calculation on both dimensions of the matrix for clustering
  #distance.gene <- dist(as.data.frame(assay[rownames(results),])) ## average distance method
  distance.gene <- as.dist((1 - cor(t(as.data.frame(assay[rownames(results),]))))/2) ## correlation method
  distance.sample <- dist(t(as.data.frame(assay[rownames(results),])))
  
  # Cluster based on the distance calculations
  #cluster.gene <- hclust(distance.gene, method="average") ## average distance method
  cluster.gene=hclust(distance.gene) ## correlation method
  cluster.sample <- hclust(distance.sample, method="average")
  
  ## Annotation information in data frame object
  dfcol <- as.data.frame(colData(dds)[,meta.annotation])
  
  # plot heat map
  p<-pheatmap(top.heat, 
              color = cols,
              fontsize = f.size,
              fontsize_row = f.size.row, 
              fontsize_col = f.size.col,
              cellwidth = c.width,
              cellheight = c.height,
              cluster_cols= F,
              #cluster_cols= cluster.sample,
              cluster_rows= cluster.gene,
              main = title,
              annotation_colors = anno.colors,
              labels_row = results$symbol,
              show_rownames = T,
              border_color = NA,
              annotation_legend = TRUE,
              annotation_col = dfcol
  )
  
  return(p)
  
}

## Use these for new plots
Gene_pHeat_Map(res.48hr.ov, "LAD overlapping", c.width = 9, c.height = 6) #48hr_inLADS_heatmap_top100
Gene_pHeat_Map(res.48hr.non, "LAD non-overlapping", c.width = 9, c.height = 6) #48hr_outsideLADS_heatmap_top100
##

#####################################################################################
# Gini Coefficients                                        
#####################################################################################

# upregulated color: #D22B2B
# downregulated color: #0047AB

library("DescTools")
library("ineq")
library("dplyr")
library(tidyr)
library(ggridges)
library(forcats)

##------------------------------Plot Histogram of Expression------------------------------##

## Break these plots up into deciles and see what it looks like 

plottheme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                 panel.grid.major = element_line(),
                 #panel.grid.minor = element_line(colour = "tomato", size=.25, linetype = "dashed"),
                 strip.background=element_blank(), axis.text.x=element_text(family="Arial", colour="dark blue"), axis.title.x=element_text(face="bold", size=20,family="Arial", colour="dark blue", vjust=-2),
                 axis.text.y=element_text(family="Arial", colour="dark blue"), 
                 axis.title.y=element_text(face="bold", size=15, family="Arial", colour="dark blue", vjust=2),
                 text = element_text(size=10, family="Arial"), 
                 plot.title=element_text(size=20, face="bold", family="Arial", color="tomato", hjust=0.0, vjust=10, lineheight=1.5),
                 plot.subtitle=element_text(size=15, face="bold", family="Arial", color="black", hjust=0.0, lineheight=1.5),
                 legend.title = element_text(size=12, color= "tomato",face="bold"), 
                 legend.text = element_text(size=10),
                 legend.key=element_rect(fill=NA),
                 axis.ticks=element_line(colour="black"))

## Remove .## version of gene in ensembl ID
rownames(dds) <- gsub("\\..*","",rownames(dds))

## Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Pre-filtering the dataset
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

## Get normalized counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
norm.counts <- counts(dds, normalized=TRUE)
#norm.counts<- counts(dds)

## Scale data for plotting
#norm.counts <-scale(norm.counts) ## For gini calculations, don't use this, use unscaled counts
norm.counts

## extract counts. Counts must be raw positive counts from non-normalized data
cts.0hr <- norm.counts[,1:2] # 0hr counts (control)
cts.12hr <- norm.counts[,3:4] # 12hr counts 
cts.48hr <- norm.counts[,5:6] # 48hr counts 
cts.144hr <- norm.counts[,7:8] # 144hr counts (washout)

## Get means of each group across replicates
mean.0hr <- as.data.frame(apply(cts.0hr, 1, function (x) mean(x)))
mean.12hr <-as.data.frame(apply(cts.12hr, 1, function (x) mean(x)))
mean.48hr <-  as.data.frame(apply(cts.48hr, 1, function (x) mean(x)))
mean.144hr <-  as.data.frame(apply(cts.144hr, 1, function (x) mean(x)))
  
## Prepare data for plotting
mean.data <- cbind(mean.0hr, mean.12hr, mean.48hr, mean.144hr)
colnames(mean.data) <- c("0Hours", "12 Hours", "48 Hours", "144 Hours")

## Factor levels
factor1 <- rep(1, times = 1, length.out = nrow(mean.0hr), each = 1)
factor2 <-rep(2, times = 1, length.out = nrow(mean.12hr), each = 1)
factor3 <-rep(3, times = 1, length.out = nrow(mean.48hr), each = 1)
factor4 <-rep(4, times = 1, length.out = nrow(mean.144hr), each = 1)
factor <- cbind(factor1, factor2, factor3, factor4)
factor <- as.data.frame(factor) %>% gather(Factor, value, na.rm = T)

## Rename with gene names and gather data into long format
mean.data$Genes <- rownames(norm.counts)
rownames(mean.data) <- make.names(rownames(norm.counts), unique=TRUE)
mean.data <- as.data.frame(mean.data) %>% gather(Time, Expression, -Genes, na.rm = T)
mean.data$Expression <- round(mean.data$Expression, 4)
mean.data$factor <- factor$value

## Set colors for individual plots
v_colors =  viridis(5, option = "E") # E= Civis
v_colors
max(mean.data$Expression)
mean(mean.data$Expression)
min(mean.data$Expression)

## Plot
p <- mean.data %>% 
  mutate(Time = fct_reorder(Time, factor)) %>% 
  ggplot(aes(x= Time, y= Expression, fill = Time)) + 
  geom_violin(trim=F, size = 1) + plottheme + 
  scale_y_continuous(limits = c(0, 3))+ #c(3, 130)) for higher expressed genes
  #stat_summary(fun.data=data_summary, geom="pointrange", color="tomato")+
  scale_fill_manual(values=c("#414D6BFF", "#7C7B78FF", "#BCAF6FFF", "#f7d13d")) +
  labs(title="Scaled Expression", subtitle="", #subtitle="Top 1000 Variable Regions", 
       caption="")
p

##------------------------------Gini Coefficients------------------------------##

## Remember to use the unscaled counts, for all positive count values
## Function to filter count data for Gini calculation
filter_data <- function(data, filter){ #filter = lamin bed data, data = cts
  
  require(dplyr)
  
  rownames(data) <- gsub("\\..*","",rownames(data))
  data <- as.data.frame(data)
  data <-data %>%
    filter(rownames(data) %in% filter$ensembl_gene_id) 
  data <- as.data.frame(apply(data, 1, function (x) mean(x)))
  
  return(data)
}



## Filter data for each group
cts.0hr.non <- filter_data(cts.0hr, non.lamins)
cts.12hr.non <- filter_data(cts.12hr, non.lamins)
cts.48hr.non <- filter_data(cts.48hr, non.lamins)
cts.144hr.non <- filter_data(cts.144hr, non.lamins)

cts.0hr.ov <- filter_data(cts.0hr, lamin.overlaps)
cts.12hr.ov <- filter_data(cts.12hr, lamin.overlaps)
cts.48hr.ov <- filter_data(cts.48hr, lamin.overlaps)
cts.144hr.ov <- filter_data(cts.144hr, lamin.overlaps)

cts.0hr <- filter_data(cts.0hr, gene.pos)
cts.12hr <- filter_data(cts.12hr, gene.pos)
cts.48hr <- filter_data(cts.48hr, gene.pos)
cts.144hr <- filter_data(cts.144hr, gene.pos)

## Calculate Gini index for each treatment group
list = c(cts.0hr, cts.12hr, cts.48hr, cts.144hr, cts.0hr.ov,cts.12hr.ov, cts.48hr.ov, cts.144hr.ov, 
         cts.0hr.non, cts.12hr.non, cts.48hr.non, cts.144hr.non )

#list = c(cts.0hr, cts.48hr, cts.144hr, cts.0hr.ov,cts.48hr.ov, cts.144hr.ov, 
         #cts.0hr.non, cts.48hr.non, cts.144hr.non )

factor.list <- c()
gini.list <- c()
for (i in 1:length(list)) {
  this <- as.data.frame(list[i])
  this <- as.data.frame(apply(this, 2, function (x) ineq::Gini(x)))
  factor <- rep(i, times = 1, length.out = 1, each = 1)
  gini.list <-rbind(gini.list,this)
  factor.list <- c(factor.list, factor)
  
  
}

## Prepare data for plotting
rownames(gini.list) <- c("0hrs", "12hrs", "48hrs", "144hrs",
                         "0hrs in LADs", "12hrs in LADs", "48hrs in LADs", "144hrs in LADs",
                         "0hrs outside LADs", "12hrs outside LADs", "48hrs outside LADs", "144hrs outside LADs" )

#rownames(gini.list) <- c("0hrs", "48hrs", "144hrs",
                         #"0hrs in LADs", "48hrs in LADs", "144hrs in LADs",
                         #"0hrs outside LADs", "48hrs outside LADs", "144hrs outside LADs" )

gini.list <- cbind(gini.list, factor.list)
colnames(gini.list)  <- c("Gini", "factor")
gini.list$Time <- c("0hrs", "12hrs", "48hrs", "6 Day Washoff", "0hrs","12hrs", "48hrs", "6 Day Washoff", "0hrs","12hrs", "48hrs", "6 Day Washoff")
#gini.list$Time <- c("0hrs", "48hrs", "6 Day Washoff", "0hrs", "48hrs", "6 Day Washoff", "0hrs", "48hrs", "6 Day Washoff")

gini.list$factor <- c(1,2,3,4,1,2,3,4,1,2,3,4)
#gini.list$factor <- c(1,2,3,1,2,3,1,2,3)
gini.list$name <- c("All genes", "All genes", "All genes", "All genes", 
                    "Genes in LADs", "Genes in LADs", "Genes in LADs", "Genes in LADs",
                    "Genes outside LADs", "Genes outside LADs", "Genes outside LADs", "Genes outside LADs")
#gini.list$name <- c("All genes", "All genes", "All genes", 
                    #"Genes in LADs", "Genes in LADs", "Genes in LADs",
                    #"Genes outside LADs", "Genes outside LADs", "Genes outside LADs")

library(viridis)

## Set colors for individual plots
q_colors =  12 # for no particular reason
v_colors =  viridis(q_colors, option = "E") # E= Civis
v_colors


p <- gini.list %>% 
  mutate(Time = fct_reorder(Time, factor)) %>%
  ggplot(aes(x= Time, y= Gini, group=name)) + 
  geom_line(aes(color=name)) + geom_point(aes(color=name), size=4) + 
  #plottheme +
  scale_y_continuous( limits=c(0.8,1), breaks=c(0.8, 0.825,0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0))+
  scale_x_discrete(limits = c("0hrs","12hrs", "48hrs","6 Day Washoff"), expand = c(0.025, 0.025)) +
  scale_color_manual(values=c("#00204DFF", "#EE82EE", "#48D1CC"))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size = 13))+
  labs(title="Gini Coefficients", subtitle="Gini Coefficients", #subtitle="Top 1000 Variable Regions", 
       caption="")
p

# upregulated color: #D22B2B
# downregulated color: #0047AB

dev.off()

## Source code to calculate cumulutive sum (lorenz curve) taken from 'ineq' package
Lorenz <- function(x, n = rep(1, length(x)), plot = FALSE)
{
  ina <- !is.na(x)    
  n <- n[ina]
  x <- as.numeric(x)[ina]
  k <- length(x)
  o <- order(x)
  x <- x[o]
  n <- n[o]
  x <- n*x
  p <- cumsum(n)/sum(n) # if everything was distributed evenly
  L <- cumsum(x)/sum(x) # the distribution of the counts
  p <- c(0,p)
  L <- c(0,L)
  L2 <- L * mean(x)/mean(n)
  Lc <- list(p,L,L2)
  names(Lc) <- c("p", "L", "L.general")
  class(Lc) <- "Lc"
  if(plot) plot(Lc)
  Lc
}

## Get mean of each gene
treatment <-apply(cts.12hr, 2, function (x) mean(x))
control <-apply(cts.0hr, 2, function (x) mean(x))

## Apply lorenz function to gene averages for plotting
treat <-Lorenz(treatment)
treat$L
ctrl <-Lorenz(control)
ctrl$L

gene.color <- "#5302a3"
name <-"12 Hour"

## Plot in Base R
plot(treat, col = gene.color, main="",
     lty = 2, lwd = 2, xlab = "Replicate", 
     ylab = "Cumulative Expression") 
lines(ctrl, col = "#febb81", lty = 1, lwd = 2)

title(paste0("Lorenz curves, ", name, " , By Treatment"))

legend("topleft", lty = 1:2, lwd = 2, cex = 1.2, legend = 
         c("Perfect Distribution", name,
           "Control"),  
       col = c("black", gene.color, "#febb81")) 

#####################################################################################
# P10:P10                                 
#####################################################################################
nrow(mean.data)
mean.data <- mean.data[mean.data$Expression > 0.01, ]
nrow(mean.data)

p10.0hr <- mean.data[mean.data$Time == "0Hours",]
p10.0hr <- p10.0hr %>% mutate(Bin = ntile(Expression, n=10))
bottom.0hr <- p10.0hr[p10.0hr$Bin==1,]
top.0hr <- p10.0hr[p10.0hr$Bin==10,]
a <-mean(top.0hr$Expression)
b <- mean(bottom.0hr$Expression)
ratio1=abs(a/b)

p10.12hr <- mean.data[mean.data$Time == "12 Hours",]
p10.12hr <- p10.12hr %>% mutate(Bin = ntile(Expression, n=10))
bottom.12hr <- p10.12hr[p10.12hr$Bin==1,]
top.12hr <- p10.12hr[p10.12hr$Bin==10,]
a <-mean(top.12hr$Expression)
b <- mean(bottom.12hr$Expression)
ratio2=(a/b)

p10.48hr <- mean.data[mean.data$Time == "48 Hours",]
p10.48hr <- p10.48hr %>% mutate(Bin = ntile(Expression, n=10))
bottom.48hr <- p10.48hr[p10.48hr$Bin==1,]
top.48hr <- p10.48hr[p10.48hr$Bin==10,]
a <-mean(top.48hr$Expression)
b <- mean(bottom.48hr$Expression)
ratio3=abs(a/b)

p10.144hr <- mean.data[mean.data$Time == "144 Hours",]
p10.144hr <- p10.144hr %>% mutate(Bin = ntile(Expression, n=10))
bottom.144hr <- p10.144hr[p10.144hr$Bin==1,]
top.144hr <- p10.144hr[p10.144hr$Bin==10,]
a <-mean(top.144hr$Expression)
b <- mean(bottom.144hr$Expression)
ratio4=(a/b)

#####################################################################################
# LAD BED DEG intersection                                     
#####################################################################################

library(plyranges)
library(GenomicRanges)
library(tidyverse) 

#bed.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Raw_Data/BEDs"
bed.dir <- "/Volumes/My Passport/projects/Emily/Lamins_RNAseq/Raw_Data/BEDs"
list.files(bed.dir)

## Get chrom names and sizes (formatted for ensemble)
#chrom.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Raw_Data/Annotations"
chrom.dir <- "/Volumes/My Passport/projects/Emily/Lamins_RNAseq/Raw_Data/Annotations"
list.files(chrom.dir)

## Genomic coordinates for all coding genes in hg38
gene.pos <- read.table(file.path(chrom.dir, "hg38_gene_positions.csv"), sep= ",", header= T)

## Convert to Genomic Ranges data
gr.gene.pos <- as.data.frame(gene.pos)
gr.gene.pos <- gr.gene.pos %>% 
  transform( seqnames= paste0("chr",gr.gene.pos$chromosome_name), start = gr.gene.pos$start_position, end = gr.gene.pos$end_position)  %>% 
  as_granges()
gr.gene.pos

##-----------------------------------BED reps-----------------------------------##
## Using BED files from https://data.4dnucleome.org/experiments-damid/4DNEXMRJPX7X/ for Labmin B1
bed.rep1 <- read_bed(file.path(bed.dir, "4DNFI2BGIZ5F.bed") , col_names = NULL, genome_info = "hg38")
bed.rep1

bed.rep2 <- read_bed(file.path(bed.dir, "4DNFIC2T8L7T.bed") , col_names = NULL, genome_info = "hg38")
bed.rep2

bed.rep.union <- GenomicRanges::reduce(c(bed.rep1, bed.rep2))

##------------------------------Luay's Lamin B1 B2------------------------------##
## Using Luay's BED file for Lamin B1 and B2
bed.lmb1 <- read_bed(file.path(bed.dir, "LMB1_4DNFICCV71TZ.bed") , col_names = NULL, genome_info = "hg38")
bed.lmb1

bed.lmb2 <- read_bed(file.path(bed.dir, "LMB2_4DNFIBQH62LX.bed") , col_names = NULL, genome_info = "hg38")
bed.lmb2

## Union of lmb1 and lmb2 files
bed.lam.union <- GenomicRanges::reduce(c(bed.lmb1, bed.lmb2))
bed.lam.union

##----------------------------------Parse LADs----------------------------------##
## Name all of the peaks based on their genomic position
names(bed.lam.union) = paste0(seqnames(bed.lam.union),':',start(bed.lam.union),'-',end(bed.lam.union))

## Assign names to temp variable
region.names <- names(bed.lam.union)

## Find intersection between LADs and all genes
## intersect_rng <- join_overlap_intersect(query, subject) https://bioconductor.org/packages/devel/bioc/vignettes/plyranges/inst/doc/an-introduction.html
lamin.overlaps <- join_overlap_intersect(bed.lam.union, gr.gene.pos)

## Find non-overlapping regions https://support.bioconductor.org/p/74077/ sp over() method gr1[!gr1 %over% gr2,]
non.lamins <- gr.gene.pos[!gr.gene.pos %over% lamin.overlaps,]

##--------------------------------------------------##
##--If analyzing Lamin 1 and 2 overlaps seperately--##
lamin1.overlaps <- join_overlap_intersect(bed.lmb1, gr.gene.pos)
lamin2.overlaps <- join_overlap_intersect(bed.lmb2, gr.gene.pos)


##---------------------------Overlaps: Data Wrangling---------------------------##

## Convert to data frame
lamin.overlaps <- as.data.frame(lamin.overlaps, row.names = seq(1:29268))
##--------------------------------------------------##
##--If analyzing Lamin 1 and 2 overlaps seperately--##
#lamin1.overlaps <- as.data.frame(lamin1.overlaps, row.names = seq(1:22140))
#lamin2.overlaps <- as.data.frame(lamin1.overlaps, row.names = seq(1:28045))

## Remove non-unique Ensembl IDs
lamin.overlaps <-lamin.overlaps %>%  
  distinct(ensembl_gene_id,.keep_all = T ) 
rownames(lamin.overlaps) <-lamin.overlaps$ensembl_gene_id

## Convert results  to DF
res.12hr.df <- as.data.frame(res.12hr) 
res.48hr.df <- as.data.frame(res.48hr)
res.144hr.df <- as.data.frame(res.144hr) 

## Filter by ensembl ID and lfc
res.12hr.ov <-res.12hr.df %>%
  filter(rownames(res.12hr.df) %in% lamin.overlaps$ensembl_gene_id) 
res.48hr.ov <-res.48hr.df %>%
  filter(rownames(res.48hr.df) %in% lamin.overlaps$ensembl_gene_id) 
res.144hr.ov <-res.144hr.df %>%
  filter(rownames(res.144hr.df) %in% lamin.overlaps$ensembl_gene_id) 

## Volcano plots of lamin intersctions
DE_Vol_Plot(res.12hr.ov,"Lamin KO at 12 Hours: Within LADs", "res.12hr.LADs.volplot" )
DE_Vol_Plot(res.48hr.ov,"Lamin KO at 48 Hours: Within LADs", "res.48hr.LADs.volplot" )
DE_Vol_Plot(res.144hr.ov,"Lamin KO at 144 Hours: Within LADs", "res.144hr.LADsvolplot" )

## Volcano plots of lamin intersctions with p value =0.1 and lfc = 0.58
DE_Vol_Plot(res.12hr.ov,"Lamin KO at 12 Hours: Within LADs pv=0.1.", "res.12hr.LADs.pv0.1.volplot" )
DE_Vol_Plot(res.48hr.ov,"Lamin KO at 48 Hours: Within LADs pv=0.1.", "res.48hr.LADs.pv0.1.volplot" )
DE_Vol_Plot(res.144hr.ov,"Lamin KO at 144 Hours: Within LADs pv=0.1.", "res.144hr.LADs.pv0.1.volplot" )

##-------------------------Non Overlaps: Data Wrangling-------------------------##

## Convert to data frame
non.lamins <- as.data.frame(non.lamins, row.names = seq(1:33360))

## Remove non-unique Ensembl IDs
non.lamins <-non.lamins %>%  
  distinct(ensembl_gene_id,.keep_all = T ) 
rownames(non.lamins) <-non.lamins$ensembl_gene_id

## Filter by ensembl ID and lfc
res.12hr.non <-res.12hr.df %>%
  filter(rownames(res.12hr.df) %in% non.lamins$ensembl_gene_id) 
res.48hr.non <-res.48hr.df %>%
  filter(rownames(res.48hr.df) %in% non.lamins$ensembl_gene_id) 
res.144hr.non <-res.144hr.df %>%
  filter(rownames(res.144hr.df) %in% non.lamins$ensembl_gene_id)

## Volcano plots of non-lamin exluded genes
DE_Vol_Plot(res.12hr.non,"Lamin KO at 12 Hours: Outside LADs", "res.12hr.nonLADs.volplot" )
DE_Vol_Plot(res.48hr.non,"Lamin KO at 48 Hours: Outside LADs", "res.48hr.nonLADs.volplot" )
DE_Vol_Plot(res.144hr.non,"Lamin KO at 144 Hours: Outside LADs", "res.144hr.nonLADsvolplot" )

## Volcano plots of lamin intersctions with p value =0.1 and lfc = 0.58
DE_Vol_Plot(res.12hr.non,"Lamin KO at 12 Hours: Outside LADs pv=0.1.", "res.12hr.nonLADs.pv0.1.volplot" )
DE_Vol_Plot(res.48hr.non,"Lamin KO at 48 Hours: Outside LADs pv=0.1.", "res.48hr.nonLADs.pv0.1.volplot" )
DE_Vol_Plot(res.144hr.non,"Lamin KO at 144 Hours: Outside LADs pv=0.1.", "res.144hr.nonLADs.pv0.1.volplot" )

##-------------------------Box Plot: Lamin intersections-------------------------##

## Merge data
res.12hr.non <-res.12hr.non %>%
  mutate(group =rep("12hr nonoverlapping", nrow(res.12hr.non)))%>%
  mutate(factor =rep(6, nrow(res.12hr.non)))
res.12hr.ov <-res.12hr.ov %>%
  mutate(group =rep("12hr overlapping", nrow(res.12hr.ov)))%>%
  mutate(factor =rep(5, nrow(res.12hr.ov)))

res.48hr.non <-res.48hr.non %>%
  mutate(group =rep("48hr nonoverlapping", nrow(res.48hr.non)))%>%
  mutate(factor =rep(4, nrow(res.48hr.non)))
res.48hr.ov <-res.48hr.ov %>%
  mutate(group =rep("48hr overlapping", nrow(res.48hr.ov)))%>%
  mutate(factor =rep(3, nrow(res.48hr.ov)))

res.144hr.non <-res.144hr.non %>%
  mutate(group =rep("144hr nonoverlapping", nrow(res.144hr.non)))%>%
  mutate(factor =rep(2, nrow(res.144hr.non)))
res.144hr.ov <-res.144hr.ov %>%
  mutate(group =rep("144hr overlapping", nrow(res.144hr.ov)))%>%
  mutate(factor =rep(1, nrow(res.144hr.ov)))

## Plot all three conditions
plot.all <- rbind(res.12hr.non, res.12hr.ov, res.48hr.non, res.48hr.ov, res.144hr.non, res.144hr.ov)

## Plot just 12 and 144 hour conditions
plot.all <- rbind(res.12hr.non, res.12hr.ov, res.144hr.non, res.144hr.ov)

plot.all <-plot.all%>%
  filter(abs(log2FoldChange) > 1 ) %>%
  filter(padj < 0.01 )

plot.all <-plot.all %>%
  mutate(log2FoldChange=abs(log2FoldChange))

plot.all <-plot.all %>%
  mutate(padj=abs(log(padj)))

library(ggplot2)
library(forcats)

plottheme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                 panel.grid.major = element_line(),
                 panel.grid.minor = element_line(),
                 strip.background=element_blank(), axis.text.x=element_text(family="Arial", colour="dark blue"), axis.title.x=element_text(face="bold", size=20,family="Arial", colour="dark blue", vjust=-2),
                 axis.text.y=element_text(family="Arial", colour="dark blue"), 
                 axis.title.y=element_text(face="bold", size=15, family="Arial", colour="dark blue", vjust=2),
                 text = element_text(size=10, family="Arial"), 
                 plot.title=element_text(size=20, face="bold", family="Arial", color="tomato", hjust=0.0, vjust=10, lineheight=1.5),
                 plot.subtitle=element_text(size=15, face="bold", family="Aria", color="black", hjust=0.0, lineheight=1.5),
                 legend.title = element_text(size=12, color= "tomato",face="bold"), 
                 legend.text = element_text(size=10),
                 legend.key=element_rect(fill=NA),
                 axis.ticks=element_line(colour="black"))

## Box plot of absolute LFC for each group "#00204DFF", "#EE82EE", "#48D1CC"
plot.all%>% #Padj_Filter0.05_LADS_LFC_BoxPlot 12hrs_144hrs_Padj0.05_LADS_LFC_BoxPlot
mutate(group = fct_reorder(group, desc(factor))) %>%
ggplot(aes(x=as.factor(group), y=log2FoldChange, fill=group)) + 
  geom_boxplot(alpha=0.75, outlier.shape = NA) + plottheme + scale_fill_manual(values = c("#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE"))+
  xlab("") +  scale_y_continuous(name="Absolute LFC", limits=c(0.5,4), breaks=c(0.5,1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0))+
  labs(title="No P Adjusted Filter") + scale_x_discrete(labels=c("12Hrs outside LADs", "12Hrs inside LADs","48Hrs outside LADs", "48Hrs inside LADs","6 Days outside LADs", "6 Days inside LADs")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## Box plot of log P Adjusted for each group
plot.all %>% #Padj_Filter0.05_LADS_Padj_BoxPlot 12hrs_144hrs_Padj0.05_LADS_Padj_BoxPlot
  mutate(group = fct_reorder(group, desc(factor))) %>%
  ggplot(aes(x=as.factor(group), y=padj, fill=group)) + 
  geom_boxplot(alpha=0.75) + plottheme + scale_fill_manual(values = c("#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE"))+
  scale_y_continuous(name="Log Padj", limits=c(0,125), breaks=c(0,25, 50, 75, 100, 125))+
  labs(title="P Adjusted < 0.05") + scale_x_discrete(labels=c("12Hrs outside LADs", "12Hrs inside LADs","48Hrs outside LADs", "48Hrs inside LADs","6 Days outside LADs", "6 Days inside LADs")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


## Box plot of SEfor each group
plot.all %>% #Padj_Filter0.05_LADS_Padj_BoxPlot LADS_SE_BoxPlot
  mutate(group = fct_reorder(group, desc(factor))) %>%
  ggplot(aes(x=as.factor(group), y=lfcSE, fill=group)) + 
  geom_boxplot(alpha=0.75, outlier.shape = NA) + plottheme + scale_fill_manual(values = c("#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE"))+
  scale_y_continuous(name="Standard Error", limits=c(0,5), breaks=c(0,1,2,3,4,5))+
  labs(title="P Adjusted < 0.05") + scale_x_discrete(labels=c("12Hrs outside LADs", "12Hrs inside LADs","48Hrs outside LADs", "48Hrs inside LADs","6 Days outside LADs", "6 Days inside LADs")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


##------------------------Bar Plot: DEGs within/out LADs------------------------##

## Get DEGs for plotting relative fraction of DEGs in each group
res <- as.data.frame(res.144hr.non)
res <-res%>%
  filter(abs(log2FoldChange) > 1 )%>%
  filter(padj < 0.01 )

## Padj: 0.05 | LFC: 0.58
#dge <- c(2372,1813,551,4563,3471,1075,1902,1486,402)

## Padj: 0.01 | LFC: 1
dge <- c(559,144,413,1717,374,1331,412,104,296)
day <- c(rep("12 hours" , 3) , rep("48 hours" , 3) , rep("6 day washoff" , 3) )
cond <- rep(c("All DEG" ,"DEG in LADs", "DEG outside LADs") , 3)
factor <- c(9,8,7,6,5,4,3,2,1)
plot.df <- data.frame(cond, day, dge, factor)

## Plot box of DEGs
plot.df %>%
  mutate(day = fct_reorder(day, factor, .desc=T)) %>%
  ggplot(aes(y=dge, x=day, color = cond)) + plottheme+
  geom_bar(position=position_dodge(0.9), stat="identity", width=0.8, fill="white") + scale_color_manual(values = c("black", "#0047AB", "#D22B2B"))+
  geom_text(aes(label=dge), position = position_dodge(width = 0.875),hjust = 0.5, vjust=-0.5, size=3.5) + xlab("Day") + ylab("Differentially Expressed Genes")


# upregulated color: #D22B2B
# downregulated color: #0047AB

##----------------------Write Out LAD Intersection Results----------------------##

## Get DEGs for plotting relative fraction of DEGs in each group
res <- res.12hr.ov
res <-res%>%
  filter(abs(log2FoldChange) > 0.58 )%>%
  filter(padj < 0.05 )

dir = "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Results/"
filename = "12hrDEGs_WithinLADs"

write.csv(res, file=paste0(dir,filename, ".csv"), row.names=FALSE)

##----------------------Trajectory Plot----------------------##


## Organize genes by significance and extract
get_Ratio <- function(results.ov, results.non){
  degs.ov <- rownames(results.ov[results.ov$padj <= 0.01 & abs(results.ov$log2FoldChange) > abs(1),])
  degs.non <- rownames(results.non[results.non$padj <= 0.01 & abs(results.non$log2FoldChange) > abs(1),])
  ratio = as.numeric(round(length(degs.ov)/length(degs.non), 4))
  return(ratio)
}

hr.12 <- get_Ratio(res.12hr.ov, res.12hr.non)
hr.48 <- get_Ratio(res.48hr.ov, res.48hr.non)
hr.144 <- get_Ratio(res.144hr.ov, res.144hr.non)

hr.12.ov <- get_Ratio(res.12hr.ov, res.12hr)
hr.48.ov <- get_Ratio(res.48hr.ov, res.48hr)
hr.144.ov <- get_Ratio(res.144hr.ov, res.144hr)

hr.12.non <- get_Ratio(res.12hr.non, res.12hr)
hr.48.non <- get_Ratio(res.48hr.non, res.48hr)
hr.144.non <- get_Ratio(res.144hr.non, res.144hr)

list.non <- rbind(hr.12.non, hr.48.non, hr.144.non)
list.ov <- rbind(hr.12.ov, hr.48.ov, hr.144.ov)
list.all <- rbind(hr.12, hr.48, hr.144)

name <- c(rep(c("LAD DEGs/Non-LAD DEGs"),3), rep(c("LAD DEGs"),3),rep(c("Non-LAD DEGs"),3))
time <- rep(c("12hrs", "48hrs","6 Day Washoff"), 3)
factor <- rep(c(1,2,3), 3)
list <- rbind(list.all,list.ov, list.non)
list <- as.data.frame(cbind(list, name, time, factor))
colnames(list) <- c("ratio","name","time","factor")

#f3e79b,#fac484,#f8a07e,#eb7f86,#ce6693,#a059a0,#5c53a5

p <- list %>% 
  mutate(time = fct_reorder(time, factor)) %>%
  ggplot(aes(x= time, y= as.numeric(ratio), group = name)) + 
  geom_line(aes(color=name)) + geom_point(aes(color=name), size=8) + 
  plottheme +
  scale_y_continuous(limits = c(0,1),breaks=c(0.0, 0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8,0.9, 1.0))+
  scale_x_discrete(limits = c("12hrs", "48hrs","6 Day Washoff"), expand = c(0.025, 0.025)) +
 scale_color_manual(values=c("#f8a07e","#ce6693","#5c53a5"))+
 theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size = 15))+ ylab("DEG Ratio") + xlab("Timepoint") +
  ggtitle("Trajectory Plot")
p



#########################################################################################
# GO using Cluster Profiler
#########################################################################################

Enrich_GO <- function(results, ont, title, filename, n.results=10, x.axis="count", plot.col="qvalue", 
                      dir.plt="/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Results/Plots/LAD_plots/Cluster_Profiler"){
  
  require(clusterProfiler)
  require(org.Hs.eg.db)
  
  GO.res.dn <- as.data.frame(results[which(results$log2FoldChange <= -0.58 & results$padj <= 0.05),])
  GO.res.up <- as.data.frame(results[which(results$log2FoldChange >= 0.58 & results$padj <= 0.05),])

  egoup <- enrichGO( gene          = rownames(GO.res.up),
                     #universe      = names(results$symbol),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = ont,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable = TRUE)
  print(egoup)
  dot.p.up <- dotplot(egoup, x=x.axis, showCategory=n.results, color=plot.col ) + ggtitle(paste("UP",title))
  ggsave(paste(dir.plt,"/",filename,"_UP",".png", sep=""), dpi=300)
  print(dot.p.up)
  
  egodn <- enrichGO(   gene          = rownames(GO.res.dn),
                       #universe      = names(results$symbol),
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = ont,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE)
  
  dot.p.dn <- dotplot(egodn, x=x.axis, showCategory=n.results, color=plot.col ) + ggtitle(paste("DN",title))
  print(dot.p.dn)
  ggsave(paste(dir.plt,"/",filename,"_DN", ".png", sep=""), dpi=300)
  print(dot.p.dn)
  
}

Enrich_GO(res.12hr, "MF", " GO Enrichment MF: 12 Hours","12hr_vs_cDN_GO_MF_dotplot")
Enrich_GO(res.12hr, "CC", " GO Enrichment CC: 12 Hours","12hr_vs_cDN_GO_CC_dotplot")
Enrich_GO(res.12hr, "BP", " GO Enrichment BP: 12 Hours","12hr_vs_cDN_GO_BP_dotplot")

Enrich_GO(res.48hr, "MF", " GO Enrichment MF: 48 Hours","48hr_vs_cDN_GO_MF_dotplot")
Enrich_GO(res.48hr, "CC", " GO Enrichment CC: 48 Hours","48hr_vs_cDN_GO_CC_dotplot")
Enrich_GO(res.48hr, "BP", " GO Enrichment BP: 48 Hours","48hr_vs_cDN_GO_BP_dotplot")

Enrich_GO(res.12hr.ov, "MF", " GO Enrichment MF in LADs: 12 Hours","12hr_vs_cDN_GO_MF_inLADs_dotplot")
Enrich_GO(res.12hr.ov, "CC", " GO Enrichment CC in LADs: 12 Hours","12hr_vs_cDN_GO_CC_inLADs_dotplot")
Enrich_GO(res.12hr.ov, "BP", "GO Enrichment BP in LADs: 12 Hours","12hr_vs_cDN_GO_BP_inLADs_dotplot")

Enrich_GO(res.12hr.non, "MF", " GO Enrichment MF outside of LADs: 12 Hours","12hr_vs_cDN_GO_MF_outLADs_dotplot")
Enrich_GO(res.12hr.non, "CC", " GO Enrichment CC outside of LADs: 12 Hours","12hr_vs_cDN_GO_CC_outLADs_dotplot")
Enrich_GO(res.12hr.non, "BP", "GO Enrichment BP outside of LADs: 12 Hours","12hr_vs_cDN_GO_BP_outLADs_dotplot")

#####################################################################################
# LAD BED DEG intersection for LMBN1 and LMBN2                                    
#####################################################################################

library(plyranges)
library(GenomicRanges)
library(tidyverse) 

bed.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Raw_Data/BEDs"
list.files(bed.dir)

## Get chrom names and sizes (formatted for ensemble)
chrom.dir <- "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Raw_Data/Annotations"
list.files(chrom.dir)

## Genomic coordinates for all coding genes in hg38
gene.pos <- read.table(file.path(chrom.dir, "hg38_gene_positions.csv"), sep= ",", header= T)

## Convert to Genomic Ranges data
gr.gene.pos <- as.data.frame(gene.pos)
gr.gene.pos <- gr.gene.pos %>% 
  transform( seqnames= paste0("chr",gr.gene.pos$chromosome_name), start = gr.gene.pos$start_position, end = gr.gene.pos$end_position)  %>% 
  as_granges()
gr.gene.pos

##------------------------------Luay's Lamin B1 B2------------------------------##
## Using Luay's BED file for Lamin B1 and B2
bed.lmb1 <- read_bed(file.path(bed.dir, "LMB1_4DNFICCV71TZ.bed") , col_names = NULL, genome_info = "hg38")
bed.lmb1

bed.lmb2 <- read_bed(file.path(bed.dir, "LMB2_4DNFIBQH62LX.bed") , col_names = NULL, genome_info = "hg38")
bed.lmb2

## Get unique regions in lmb1 and lmb2 #filter_by_non_overlaps(query, subject)
lmb1.unique <- bed.lmb1[!bed.lmb1 %over% bed.lmb2,]
lmb2.unique <- bed.lmb2[!bed.lmb2 %over% bed.lmb1,]

## Second way to get unique regions (plyranges)
lmb1.unique <- filter_by_non_overlaps(bed.lmb1, bed.lmb2)
lmb2.unique <- filter_by_non_overlaps(bed.lmb2, bed.lmb1)

##----------------------------------Parse LADs----------------------------------##
## Name all of the peaks based on their genomic position
names(lmb1.unique) = paste0(seqnames(lmb1.unique),':',start(lmb1.unique),'-',end(lmb1.unique))
names(lmb2.unique) = paste0(seqnames(lmb2.unique),':',start(lmb2.unique),'-',end(lmb2.unique))


## Assign names to temp variable
region.names.l1 <- names(lmb1.unique)
region.names.l2 <- names(lmb2.unique)

## Find intersection between LADs and all genes
## intersect_rng <- join_overlap_intersect(query, subject) https://bioconductor.org/packages/devel/bioc/vignettes/plyranges/inst/doc/an-introduction.html
lamin1.overlaps <- join_overlap_intersect(lmb1.unique, gr.gene.pos)
lamin2.overlaps <- join_overlap_intersect(lmb2.unique, gr.gene.pos)

##---------------------------Overlaps: Data Wrangling---------------------------##

## Convert to data frame
lamin1.overlaps <- as.data.frame(lamin1.overlaps, row.names = seq(1:627))
lamin2.overlaps <- as.data.frame(lamin2.overlaps, row.names = seq(1:439))

## Remove non-unique Ensembl IDs
lamin1.overlaps <-lamin1.overlaps %>%  
  distinct(ensembl_gene_id,.keep_all = T ) 
rownames(lamin1.overlaps) <-lamin1.overlaps$ensembl_gene_id

## Remove non-unique Ensembl IDs
lamin2.overlaps <-lamin2.overlaps %>%  
  distinct(ensembl_gene_id,.keep_all = T ) 
rownames(lamin2.overlaps) <-lamin2.overlaps$ensembl_gene_id

## Convert results  to DF
res.12hr.df <- as.data.frame(res.12hr) 
res.48hr.df <- as.data.frame(res.48hr)
res.144hr.df <- as.data.frame(res.144hr) 

## Filter by ensembl ID and lfc
res.12hr.l1 <-res.12hr.df %>%
  filter(rownames(res.12hr.df) %in% lamin1.overlaps$ensembl_gene_id) 
res.48hr.l1 <-res.48hr.df %>%
  filter(rownames(res.48hr.df) %in% lamin1.overlaps$ensembl_gene_id) 
res.144hr.l1 <-res.144hr.df %>%
  filter(rownames(res.144hr.df) %in% lamin1.overlaps$ensembl_gene_id) 

## Filter by ensembl ID and lfc
res.12hr.l2 <-res.12hr.df %>%
  filter(rownames(res.12hr.df) %in% lamin2.overlaps$ensembl_gene_id) 
res.48hr.l2 <-res.48hr.df %>%
  filter(rownames(res.48hr.df) %in% lamin2.overlaps$ensembl_gene_id) 
res.144hr.l2 <-res.144hr.df %>%
  filter(rownames(res.144hr.df) %in% lamin2.overlaps$ensembl_gene_id) 


##-------------------------Box Plot: Lamin intersections-------------------------##

## Merge data
res.12hr.l1 <-res.12hr.l1 %>%
  mutate(group =rep("12hr lamin 1", nrow(res.12hr.l1)))%>%
  mutate(factor =rep(6, nrow(res.12hr.l1)))
res.12hr.l2 <-res.12hr.l2 %>%
  mutate(group =rep("12hr lamin 2", nrow(res.12hr.l2)))%>%
  mutate(factor =rep(5, nrow(res.12hr.l2)))

res.48hr.l1 <-res.48hr.l1 %>%
  mutate(group =rep("48hr kamin 1", nrow(res.48hr.l1)))%>%
  mutate(factor =rep(4, nrow(res.48hr.l1)))
res.48hr.l2 <-res.48hr.l2 %>%
  mutate(group =rep("48hr lamin 2", nrow(res.48hr.l2)))%>%
  mutate(factor =rep(3, nrow(res.48hr.l2)))

res.144hr.l1 <-res.144hr.l1 %>%
  mutate(group =rep("144hr lamin 1", nrow(res.144hr.l1)))%>%
  mutate(factor =rep(2, nrow(res.144hr.l1)))
res.144hr.l2 <-res.144hr.l2 %>%
  mutate(group =rep("144hr lamin 2", nrow(res.144hr.l2)))%>%
  mutate(factor =rep(1, nrow(res.144hr.l2)))

## Plot all three conditions
plot.all <- rbind(res.12hr.l1, res.12hr.l2, res.48hr.l1, res.48hr.l2, res.144hr.l1, res.144hr.l2)

plot.all <-plot.all%>%
  filter(abs(log2FoldChange) > 0.58 ) %>%
  filter(padj < 0.1 )

plot.all <-plot.all %>%
  mutate(log2FoldChange=abs(log2FoldChange))

plot.all <-plot.all %>%
  mutate(padj=abs(log(padj)))

library(ggplot2)
library(forcats)

plottheme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                 panel.grid.major = element_line(),
                 panel.grid.minor = element_line(),
                 strip.background=element_blank(), axis.text.x=element_text(family="Arial",colour="dark blue"), axis.title.x=element_text(face="bold", size=20,family="Arial", colour="dark blue", vjust=-2),
                 axis.text.y=element_text(family="Arial",colour="dark blue"), 
                 axis.title.y=element_text(face="bold", size=15, family="Arial", colour="dark blue", vjust=2),
                 text = element_text(size=10, family="Arial"), 
                 plot.title=element_text(size=20, face="bold", family="Arial", color="tomato", hjust=0.0, vjust=10, lineheight=1.5),
                 plot.subtitle=element_text(size=15, face="bold", family="Aria", color="black", hjust=0.0, lineheight=1.5),
                 legend.title = element_text(size=12, color= "tomato",face="bold"), 
                 legend.text = element_text(size=10),
                 legend.key=element_rect(fill=NA),
                 axis.ticks=element_line(colour="black"))

## Box plot of absolute LFC for each group "#00204DFF", "#EE82EE", "#48D1CC"
plot.all%>% #Padj_Filter0.05_LADS_LFC_BoxPlot 12hrs_144hrs_Padj0.05_LADS_LFC_BoxPlot Padj0.1lamin1and2_LADS_LogPadj_BoxPlot
  mutate(group = fct_reorder(group, desc(factor))) %>%
  ggplot(aes(x=as.factor(group), y=log2FoldChange, fill=group)) + 
  geom_boxplot(alpha=0.75) + plottheme + scale_fill_manual(values = c("#48D1CC", "#EE82EE","#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE"))+
  xlab("") +  scale_y_continuous(name="Absolute LFC", limits=c(0.5,4), breaks=c(0.5,1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0))+
  labs(title="No P Adjusted Filter") + scale_x_discrete(labels=c("12Hrs Lamin1 LADs", "12Hrs Lamin2 LADs","48Hrs Lamin1 LADs", "48Hrs Lamin2 LADs","6 Days Lamin1 LADs", "6 Days Lamin2 LADs")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
## Box plot of log P Adjusted for each group
plot.all %>% #Padj_Filter0.05_LADS_Padj_BoxPlot 12hrs_144hrs_Padj0.05_LADS_Padj_BoxPlot
  mutate(group = fct_reorder(group, desc(factor))) %>%
  ggplot(aes(x=as.factor(group), y=padj, fill=group)) + 
  geom_boxplot(alpha=0.75) + plottheme + scale_fill_manual(values = c("#48D1CC", "#EE82EE", "#48D1CC", "#EE82EE","#48D1CC", "#EE82EE"))+
  scale_y_continuous(name="Log Padj", limits=c(0,125), breaks=c(0,25, 50, 75, 100, 125))+
  labs(title="P Adjusted < 0.05") + scale_x_discrete(labels=c("12Hrs Lamin1 LADs", "12Hrs Lamin2 LADs","48Hrs Lamin1 LADs", "48Hrs Lamin2 LADs","6 Days Lamin1 LADs", "6 Days Lamin2 LADs")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


##------------------------Bar Plot: DEGs within/out LADs------------------------##

## Get DEGs for plotting relative fraction of DEGs in each group
res <- res.144hr.l2
res <-res%>%
  filter(abs(log2FoldChange) > 0.58 ) %>%
  filter(padj < 0.1 )

dge <- c(125,159,11,13,17,29,9,14)
day <- c( "All genes - Lamin1 LADs", "All genes - Lamin2 LADs", "12 hrs DEGs - Lamin1 LADs", "12 hrs DEGs - Lamin2 LADs", 
          "48 hrs DEGs - Lamin1 LADs", "48 hrs DEGs - Lamin2 LADs", "144 hrs DEGs - Lamin1 LADs", "144 hrs DEGs - Lamin2 LADs")

cond <- rep(c("Lamin1 LADs" , "Lamin2 LADs") , 4)
factor <- c(8,7,6,5,4,3,2,1)
plot.df <- data.frame(cond, day, dge, factor)

## Plot box of DEGs #DEGs_LADs_L1_L2
plot.df %>%
  mutate(day = fct_reorder(day, factor, .desc=T)) %>%
  ggplot(aes(y=dge, x=day, color = cond)) + plottheme+
  geom_bar(position=position_dodge(0.9), stat="identity", width=0.8, fill="white") + 
  scale_color_manual(values = c( "black","#0047AB", "#D22B2B"))+ theme(axis.text.x = element_text(angle = 45,  hjust=1) )+
  geom_text(aes(label=dge), position = position_dodge(width = 0.875),hjust = 0.5, vjust=-0.5, size=3.5) + xlab("Day") + ylab("Differentially Expressed Genes")

# upregulated color: #D22B2B
# downregulated color: #0047AB

##----------------------Write Out LAD Intersection Results----------------------##

## Get DEGs for plotting relative fraction of DEGs in each group
res <- res.144hr.l2
res <-res%>%
  filter(abs(log2FoldChange) > 0.58 )%>%
  filter(padj < 0.05 )

dir = "/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_RNAseq/Results/"
filename = "144hrDEGs_Within_LMB2_LADs"

write.csv(res, file=paste0(dir,filename, ".csv"), row.names=FALSE)

#####################################################################################
# Plot metascape results                                  
#####################################################################################

library(viridis)
library(forcats)

#OFvOP_RNAseq_Up_Stringent

dir_meta <-"/Volumes/My Passport/projects/Emily/Lamins_RNAseq/Results/Metascape"
list.files(dir_meta)
setwd(dir_meta)

## Read in Metascape results
meta.lamin <- read.table(file.path(dir_meta, "Metascape_Results_Summary.csv"), sep= ",", header= TRUE, fill = TRUE ) #OFvOP_RNAseq_Dn_Stringent

##  Build new data frame
meta.lamin <- data.frame(meta.lamin)

LAD = "in"
#LAD = "out"
Type = "cancer"
#Type = "structural"

meta.lamin <- meta.lamin[meta.lamin$LAD == LAD & meta.lamin$Type == Type,]

# upregulated color: #D22B2B
# downregulated color: #0047AB

subset = "12 Hours within LADs "
#subset = "12 Hours outside LADs "

## Plot categories # DisGeNET_within_Cancer
p.meta<-meta.lamin %>% 
  mutate(Name = fct_reorder(Name, pvalue)) %>% 
  ggplot(aes(x = pvalue, y = Name))+
  geom_point(aes(colour = as.numeric(pvalue),size = as.numeric(pvalue)))+
  scale_size_continuous(name="-Log 10 P Value")+
  scale_color_continuous(low = "#0047AB", high = "#D22B2B", name="P Value")+
  ggtitle(paste0("DisGeNET, ",subset," : ", Type))+ xlab("-Log 10 P Value")+ ylab("Annotation")+
  scale_x_continuous( limits = c(0,7), breaks = c(0,1,2,3,4,5,6,7))+
  theme(
    axis.line = element_line(color = "darkblue", 
                             size = 1, linetype = "solid"),
    axis.text.y = element_text( color="#008080", size=14, family = "Arial Narrow"),
    axis.text.x = element_text(angle = 0, color="black", size=12, family = "Arial Narrow"),
    plot.title = element_text(hjust = 0.5, color="black", face = "bold", size=15,
                              family = "Arial Rounded MT Bold"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank()) 
p.meta 




