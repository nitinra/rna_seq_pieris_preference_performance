##tutorials: https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html 
## https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html
## http://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html
## https://rnnh.github.io/bioinfo-notebook/docs/DE_analysis_edgeR_script.html
## https://research.stowers.org/cws/CompGenomics/Projects/edgeR.html
## https://bioconductor.org/packages/devel/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html
##https://bioconductor.org/packages/devel/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html
##https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#test-go-at-sepcific-level
##https://www.spandidos-publications.com/10.3892/etm.2018.6884
##https://ycl6.github.io/GO-Enrichment-Analysis-Demo/
##https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0##
##https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html##

setwd("/Users/nitinr/Library/CloudStorage/OneDrive-UniversityofSouthCarolina/USC/Lab/Chapter3_RNA-seq/Data")

#description-of-the-biological-experiment
##packages needed for analysis##
library(palettetown)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(GenomicFeatures)
library(mixOmics)
library(HTSFilter)
library(topGO)
library(dplyr)
library(tidyverse)
library(magrittr)
library(WGCNA)
library(reshape2)
library(CorLevelPlot)
##Script to create header=

## read files##
allsample<-read.table("featurecount_results_2.txt", sep="\t",header=T, row.names=1,check.names=F)
allsample <- allsample[,-c(2,3,4,5,25)] ## only take those columns that have info about sample gene counts and leave all gene annotation.
design<-read.table("design_larvae_remout.csv", header=TRUE, sep=",", row.names=1) #design of treatments.
dgeFull <- DGEList(allsample[,-1],group = factor(design$condition),genes=allsample[,1,drop=FALSE]) ## prepare a DGE object; remove outlier sample
options(digits=3)

## Extract pseudo-counts (ie log2(K+1))
pseudoCounts <- log2(dgeFull$counts+1)
head(pseudoCounts)

##Histogram for pseudo-counts (for one sample)
hist(pseudoCounts[,"23-C3_R.transcriptAlignedfinalsortedbycoord.bam"])

##Boxplot for pseudo-counts
boxplot(pseudoCounts, col="gray", las=3)

##MDS for pseudo-counts (using limma package)
plotMDS(pseudoCounts)

##heatmap for pseudo-counts (using mixOmics package)
sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists
pheatmap(sampleDists)

##remove genes with zero counts for all samples
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)
head(dgeFull$counts)

##estimate the normalization factors
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples

#**Important: using calcNormFactors does not change the counts: it just updates the column norm.factors in $samples. It is therefore recommanded that you use the same name (dgeFull) to save the result of this function:
head(dgeFull$counts)

##From the normalization factors and the original count table, find the normalized counts and use the log2-transformation to inspect them with boxplots and a MDS. Normalized counts can be extracted from dgeFull using the function cpm:
eff.lib.size <- dgeFull$samples$lib.size*dgeFull$samples$norm.factors
normCounts <- cpm(dgeFull)
pseudoNormCounts <- log2(normCounts + 1)
write.csv(pseudoNormCounts, "pseudonormcount_input_wgcna.csv")
boxplot(pseudoNormCounts, col="gray", las=3)

##MDS for pseudo-counts_normalized (using limma package)
plotMDS(pseudoNormCounts)

##estimate common and tagwise dispersion
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeFull

##perform an exact test for the difference in expression between the two conditions “WT” and “Mt”
dgeTest <- exactTest(dgeFull)
dgeTest

##remove low expressed genes
filtData <- HTSFilter(dgeFull)$filteredData
dgeTestFilt <- exactTest(filtData)
dgeTestFilt

##plot an histogram of unadjusted p-values
hist(dgeTest$table[,"PValue"], breaks=50)

##histogram of unadjusted p-values after filtering
hist(dgeTestFilt$table[,"PValue"], breaks=50)

##extract a summary of the differential expression statistics
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))
head(resNoFilt)

resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))
head(resFilt)

##compare the number of differentially expressed genes with and without filtering (risk: 1%)
sum(resNoFilt$table$FDR < 0.01) ## before independent filtering
sum(resFilt$table$FDR < 0.01) ## after independent filtering

## write the DGE to a csv file
write.csv(resFilt, file="resFilt_07-20_2023.csv")

##Find upregulated and downregulated genes separately
sigDownReg <- resFilt$table[resFilt$table$FDR<0.01,]
sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]

##write the results in csv files

write.csv(sigDownReg, file="sigDownReg_larvae_07-20_2023.csv")
write.csv(sigUpReg, file="sigUpReg_larvae_07-20-2023.csv")

##volcano plot
volcanoData <- cbind(resFilt$table$logFC, -log10(resFilt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)
plot(volcanoData, pch=19)

plot(resFilt$table$logFC, -10*log10(resFilt$table$FDR),xlab="M", ylab="-10*log(P-val)")
thlaspi.enriched<- resFilt$table$logFC > 0 & resFilt$table$FDR <0.01
thlaspi.depleted<- resFilt$table$logFC < 0 & resFilt$table$FDR < 0.01
points(resFilt$table$logFC[thlaspi.enriched],-10*log10(resFilt$table$FDR[thlaspi.enriched]), col='red' )
points(resFilt$table$logFC[thlaspi.depleted],-10*log10(resFilt$table$FDR[thlaspi.depleted]), col='blue' )


de<-as.data.frame(resFilt$table)
de$diffexpressed <- "NO"
de$diffexpressed[de$logFC > 1.5 & de$FDR < 0.01] <- "UP"
de$diffexpressed[de$logFC < -1.5 & de$FDR <0.01] <- "DOWN"

p <- ggplot(data=de, aes(x=logFC, y=-log10(FDR), col=diffexpressed, size=diffexpressed, alpha=diffexpressed)) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-1.5, 1.5), col="#4DA167", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.01), col="#4DA167", linetype = "dashed")
p2
p3<-p2+scale_color_manual(values=c("#4F7CAC", "#0C0C0C", "#D00000"))
p3
sizes <- c("UP" = 2, "DOWN" = 2, "NO" = 1) 
p4<-p3+scale_size_manual(values=sizes)
p4
alphas <- c("UP" = 1, "DOWN" = 1, "NO" = 0.5)
p5<-p4+scale_alpha_manual(values = alphas)
p5

##Heatmap
y <- cpm(dgeFull, log=TRUE, prior.count = 1)
selY <- y[rownames(resFilt$table)[resFilt$table$FDR<0.01 & 
                                    abs(resFilt$table$logFC)>1.5],]
head(selY)
dim(selY)
selY
## heatmap with logCPM
selY2 <- y[rownames(resFilt$table)[resFilt$table$FDR<0.01 & 
                                    abs(resFilt$table$logCPM)],]
write.csv(selY2, file="heatmap_matrix_CPM_2023-08-24.csv")
hm2<-pheatmap(selY2, scale="row")
plot(hm2$tree_row)
plot(hm2$tree_col)
heatmap2<-read.csv()
geneClust_row <- cutree(as.hclust(hm2$tree_row), k=3)
write.csv(geneClust_row, file="cluster_cpm_k3_2023-08-13.csv")

## write the heatmap into a matrix
write.csv(selY, file="heatmap_matrix_2023-08-13.csv") ## write FC>1.5 results into a heatmap matrix
hm <- pheatmap(selY, scale="row" )
plot(hm$tree_row)
plot(hm$tree_col)

##clustering of genes based on k-means
heatmap<-read.csv("heatmap_matrix_2023-08-13_2.csv")
heatmap<-as.matrix(heatmap)
scaledata <- t(scale(t(heatmap)))
wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(scaledata,
                                     centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

library(cluster)
sil <- rep(0, 20)

#repeat k-means for 1:20 and extract silhouette:
for(i in 2:20){
  k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(scaledata))
  sil[i] <- mean(ss[, 3])
}

# Plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)

## Plot a new heatmap based on the #of gene clusters
geneClust_row <- cutree(as.hclust(hm$tree_row), k=3)
write.csv(geneClust_row, file="cluster_k3_2024-01-05.csv")
geneClust_row <- as.data.frame(geneClust_row)
geneClust_col <- cutree(as.hclust(hm$tree_col), k=2)
geneClust_col<- data.frame(cluster = ifelse(test = geneClust_col == 1, yes = "Cardamine", no = "Thlaspi"))

###extract each cluster and save into a different file.
geneClust <- cutree(as.hclust(hm$tree_row), k=3)
geneClust<-sort(geneClust)
head(geneClust)
length(unique(geneClust))
cluster1<-as.data.frame(names(which(geneClust==1)))
colnames(cluster1)<-c('locus_id')
fulllist<-read.csv("resFilt_06-08_2023.csv", h=T)
cluster1list<-merge(cluster1, fulllist, by ='locus_id')
write.csv(cluster1list, "cluster1_2023-08-20.csv")

cluster2<-as.data.frame(names(which(geneClust==2)))
colnames(cluster2)<-c('locus_id')
cluster2list<-merge(cluster2, fulllist, by ='locus_id')
write.csv(cluster2list, "cluster2_2023-08-20.csv")

cluster3<-as.data.frame(names(which(geneClust==3)))
colnames(cluster3)<-c('locus_id')
cluster3list<-merge(cluster3, fulllist, by ='locus_id')
write.csv(cluster3list, "cluster3_2023-08-20.csv")


color<-list(
  class =pokepal(pokemon=196, spread=5)
)
par(mar = c(1, 1, 1, 1))
pdf("heatmap.pdf", width=8, height=45,onefile = T)
hm_2 <- pheatmap(selY, scale="row",annotation_row =geneClust_row, 
                 annotation_col =geneClust_col, cutree_rows = 5, cutree_cols = 2,
                 annotation_colors =color, cellheight=2.5,show_rownames=F)
dev.off() 

##Plot CPM for each cluster
d1<-read.csv("heatmap_matrix_2023-08-13.csv", h=T)  ## read heatmap matrix file
d2<-read.csv("final_list_of_genes_07-20-2023.csv", h=T) ## read list of gene flies
d3<- merge(d1,d2, by="locus_id") ## merge the heatmap matrix and genelist by locus id
write.csv(d3, "protid_heatmap_matrix_2023-08-20.csv") ## write the new file sorted by locus id which contains protein id
d4<-read.csv("protid_heatmap_matrix_2023-08-20.csv", h=T) ## read the file with protein id
d5<-read.csv("genome_features.csv", h=T) ## genome features file (gff file)
d6<- merge(d4,d5, by="protein_id") ## merge the two file by protein id
write.csv(d6, "geneid_heatmap_matrix_2023-08-20.csv") ## write the new file with gene id. Remove columns that are not needed.

d2<- read.csv("cluster_heatmap_protid_geneid_2023-08-20.csv")
d3<-melt(d2, id.vars=c("protein_id", "cluster")) ## convert to long format
colnames(d3)<-c("protein_id", "cluster", "individual", "logFC") ## rename columns
Cardamine<-c('C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6', 'C_7', 'C_8', 'C_9', 'C_10')
Thlaspi<-c('T_1', 'T_2', 'T_3', 'T_4', 'T_5', 'T_6', 'T_7', 'T_8', 'T_9', 'T_10')

d4<-d3 %>%
  mutate(individual=as.character(individual)) %>%
  mutate(treatment=
           case_when(
             individual %in% Cardamine ~"cardamine",
             individual %in% Thlaspi ~"thlaspi",
             TRUE ~individual
           ))
str(d4)
d4$treatment<-as.factor(d4$treatment)
d4$cluster<-as.factor(d4$cluster)
ggplot(d4, aes(x=cluster, y=logFC,fill=treatment))+
  scale_fill_manual(values=c("#8C7250", "#8CB568")) +
  theme_minimal()+
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black", lwd=0.2,alpha = 0.6,show.legend = T, width=0.75)+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, stroke=0,dodge.width = .75, alpha=.6,show.legend = F)


ggplot(d4, aes(x=individual, y=logFC,fill=treatment))+
  scale_fill_manual(values=c("#8C7250", "#8CB568")) + facet_wrap(~cluster)+
  theme_minimal()+
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black", lwd=0.2,alpha = 0.6,show.legend = T, width=0.75)+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, stroke=0,dodge.width = .75, alpha=.6,show.legend = F)

ggplot(d4, aes(x=individual, y=logFC,fill=treatment))+
  scale_fill_manual(values=c("#8C7250", "#8CB568")) + facet_wrap(~cluster)+
  theme_minimal()+
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black", lwd=0.2,alpha = 0.6,show.legend = T, width=0.75)

library(rstatix)
stat.test1<-d4 %>%
  group_by(cluster) %>%
  t_test(logFC~treatment) %>%
  adjust_pvalue(method="BH") %>%
  add_significance()
stat.test1
write.csv(stat.test1, "t_test_clusterCPMExpression_2023-08-20.csv")

##
de1<-read.csv("detoxenzyme_cpm_2023-08-24.csv", h=T)
de2<-melt(de1, id.vars=c("cluster", "PFAMs")) ## convert to long format
colnames(de2)<-c("cluster", "PFAMs", "individual", "logCPM") ## rename columns
Cardamine<-c('C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6', 'C_7', 'C_8', 'C_9', 'C_10')
Thlaspi<-c('T_1', 'T_2', 'T_3', 'T_4', 'T_5', 'T_6', 'T_7', 'T_8', 'T_9', 'T_10')
de3<-de2 %>%
  mutate(individual=as.character(individual)) %>%
  mutate(treatment=
           case_when(
             individual %in% Cardamine ~"cardamine",
             individual %in% Thlaspi ~"thlaspi",
             TRUE ~individual
           ))
str(de3)



de3$treatment<-as.factor(de3$treatment)
de3$cluster<-as.factor(de3$cluster)
de3$PFAMs<-as.factor(de3$PFAMs)
ggplot(de3, aes(x=PFAMs, y=logCPM,fill=treatment))+
  scale_fill_manual(values=c("#8C7250", "#8CB568")) +
  theme_minimal()+
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black", lwd=0.2,alpha = 0.6,show.legend = T, width=0.75)+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, stroke=0,dodge.width = .75, alpha=.6,show.legend = F)

stat.test2<-de3 %>%
  group_by(PFAMs) %>%
  t_test(logCPM~treatment) %>%
  adjust_pvalue(method="BH") %>%
  add_significance()
stat.test2
write.csv(stat.test2, "t_test_detoxExpression_2023-08-24.csv")


###

##GO analysis##

## match and merge fasta_id/locus_id with protein id
file1<-read.table ("match.txt", header=T) # protein id file from protein.fa
file2<-read.table ("ella.txt", header=T, sep='\t') ## file containing fasta,locus & protein id, extracted from cds.fa
listofgenes<-merge(file1, file2, by="protein_id")
write.csv(listofgenes, file="final_list_of_genes_07-20-2023.csv") ## final list of genes

##match and merge with deg
listofgenes<-read.csv("final_list_of_genes_07-20-2023.csv", h=T)
downreg<-read.csv("sigDownReg_larvae_08-13_2023.csv", h=T)
upreg<-read.csv("sigUpReg_larvae_08-13-2023.csv", h=T)
merge_listofgenes_downreg<-merge(listofgenes, downreg, by="locus_id")
merge_listofgenes_upreg<-merge(listofgenes, upreg, by="locus_id")
downreg<-merge_listofgenes_downreg[,-c(1,2,4)]
colnames(downreg)[1]<- "geneid"
upreg<-merge_listofgenes_upreg[,-c(1,2,4)]
colnames(upreg)[1]<- "geneid"

##Create a topgo object
geneID2GO <- readMappings(file ='topgo_input.txt', sep='\t') ## open gene annotations into a top go object.
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% upreg$geneid))
names(geneList) <- geneNames

##Create a GO object
GOdata <- new("topGOdata", ontology = "MF",allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
allGO=usedGO(GOdata)
length(allGO)
##fisher test for gene enrichment
results.fisher <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher")
##results.weight<- runTest(GOdata, algorithm = "weight01", statistic = "fisher") ## if you want to use another statistic
##get summary of the results
geneData(results.fisher)

## Create a table of results of the top 100 GO hits
goEnrichment <- GenTable(
  GOdata,
  pvalue = results.fisher,
  orderBy= "results.fisher",
  topNodes = 92,
  numChar = 99)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] # filter terms for KS p<0.05
write.csv(goEnrichment, file="goEnrichment_upreg_MF_2023-08-13.csv")


##Cluster wise GO enrichment.
# Define a function to process cluster data
process_cluster <- function(cluster_list, output_filename) {
  # Merge the listofgenes with the cluster list
  merge_list <- merge(listofgenes, cluster_list, by="locus_id")
  
  # Remove unwanted columns
  cluster <- merge_list[,-c(1, 2, 4)]
  colnames(cluster)[1] <- "geneid"
  
  # Create a topGO object
  geneID2GO <- readMappings(file ='topgo_input.txt', sep='\t')
  geneNames <- names(geneID2GO)
  geneList <- factor(as.integer(geneNames %in% cluster$geneid))
  names(geneList) <- geneNames
  
  # Create a GO object
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  # Run Fisher test for gene enrichment
  results.fisher <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher")
  
  # Get summary of the results
  a <- geneData(results.fisher)
  sig_terms_value <- a["SigTerms"]
  
  # Create a table of results of the top GO hits
  goEnrichment <- GenTable(
    GOdata,
    pvalue = results.fisher,
    orderBy = "results.fisher",
    topNodes = sig_terms_value,
    numChar = 99
  )
  
  # Convert pvalue to numeric
  goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
  
  # Filter terms for pvalue < 0.05
  goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,]
  
  # Write the results to a CSV file
  write.csv(goEnrichment, file = output_filename)
}

# Process c1
process_cluster(cluster1list, "goEnrichment_c1_BP_2024-01-05.csv")

# Process c2
process_cluster(cluster2list, "goEnrichment_c2_BP_2024-01-05.csv")

# Process c3
process_cluster(cluster3list, "goEnrichment_c3_BP_2024-01-05.csv")

###########

##WGCNS analysis######
##Use the pseudoNormCounts generated from the edgeR analysis. The data needs to be transposed, where rows are samples
## and columns are genes.
input_mat <-t(pseudoNormCounts)
allowWGCNAThreads() # allow multi-threading (optional)
powers = c(c(1:10), seq(from = 12, to = 20, by = 2)) # Choose a set of soft-thresholding powers
# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)
##Plot r^2 and mean connectivity plot to see which power to choose.
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 6 ## Value picked based on the r^2 analysis.
temp_cor <- cor ## A temp directory is needed to store results.
cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)


##Build WGCNA network

netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 3,
                          pamRespectsDendro = F,
                          ##detectCutHeight = 0.75,
                          minModuleSize = 40,
                          maxBlockSize = 16000,
                          corType = "pearson",
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.30,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor # Return cor function to original namespace
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
colorh1<- mergedColors[netwk$blockGenes[[1]]]

# Plot netwk# Plot the dendrogram and the module colors underneath

plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
length(unique(netwk$colors))
netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)

###Relate Module (cluster) Assignments to Treatment Groups
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
row.names(module_df) = module_df$gene_id
module_df[1:5,]
write.csv(module_df,
          file = "gene_modules.csv",
)

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
write.csv(MEs0, "gene_model_matrix.csv")
# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
design2<-read.table("design_larvae_remout.csv", header=TRUE, sep=",")
traits <- design2 %>% 
  mutate(treatment_bin = ifelse(grepl('T3', file), 1, 0))%>%
  select(1,3)
traits$treatment_bin<-as.factor(traits$treatment_bin)

#Add new column to eigengenes
MEs0$treatment = traits$file

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )
mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

####Gene significance correlations

colorlevels=unique(colorh1)

datKME=signedKME(input_mat, MEs0[1:33], outputColumnName="MM.")

GS1=as.numeric(cor(traits$treatment_bin,input_mat, use="p"))
GeneSignificance=abs(GS1)
ModuleSignificance=tapply(GeneSignificance, colorh1, mean, na.rm=T)
sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,colorh1)

head(datKME)

# Define a function to filter genes for a given module, print relevant information, and save filtered genes to a CSV file
filterAndPrintAndSaveCSV <- function(module, MM_column) {
  FilterGenes <- abs(GS1) > 0.6 & abs(datKME[[MM_column]]) > 0.80
  cat("Module:", module, "\n")
  cat("Number of filtered genes:", sum(FilterGenes), "\n")
  filtered_gene_names <- dimnames(data.frame(input_mat))[[2]][FilterGenes]
  cat("Filtered genes:", filtered_gene_names, "\n\n")
  
  # Create a data frame with filtered genes and write it to a CSV file
  filtered_genes_df <- data.frame(GeneName = filtered_gene_names)
  csv_filename <- paste("filtered_genes_", module, ".csv", sep = "")
  write.csv(filtered_genes_df, csv_filename, row.names = FALSE)
}

# Define a list of modules and their corresponding MM column names
modules <- c("black", "blue", "green", "turquoise", "magenta")
MM_columns <- c("MM.black", "MM.blue", "MM.green", "MM.turquoise", "MM.magenta")

# Loop through the modules, filter genes, print information, and save filtered genes to CSV files
for (i in seq_along(modules)) {
  filterAndPrintAndSaveCSV(modules[i], MM_columns[i])
}

####
# GO analysis for filtered genes.
performGOEnrichmentAnalysis <- function(dataframe_name, output_prefix) {
  listofgenes <- read.csv("final_list_of_genes_07-20-2023.csv", header = TRUE)
  dataframe <- read.csv(paste("filtered_genes_", dataframe_name, ".csv", sep = ""), header = TRUE)
  
  merged <- merge(listofgenes, dataframe, by = "locus_id")
  merged <- merged[, -2]
  colnames(merged)[2] <- "geneid"
  
  geneID2GO <- readMappings(file = 'topgo_input.txt', sep = '\t')
  geneNames <- names(geneID2GO)
  geneList <- factor(as.integer(geneNames %in% merged$geneid))
  names(geneList) <- geneNames
  
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  allGO <- usedGO(GOdata)
  
  results.fisher <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher")
  geneData(results.fisher)
  
  a <- geneData(results.fisher)
  sig_terms_value <- a["SigTerms"]
  
  goEnrichment <- GenTable(
    GOdata,
    pvalue = results.fisher,
    orderBy = "results.fisher",
    topNodes = sig_terms_value,
    numChar = 99
  )
  
  goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
  goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,]
  
  output_filename <- paste("goEnrichment_", output_prefix, "_2024-01-07.csv", sep = "")
  write.csv(goEnrichment, file = output_filename)
}

# Apply the function to each dataframe
dataframes <- c("black", "blue", "green", "magenta", "turquoise")

for (df_name in dataframes) {
  performGOEnrichmentAnalysis(df_name, df_name)
}

####
##VennDiagram for looking at overlap in GO between modules and DGE.
# Load the VennDiagram package
library(VennDiagram)

# Load the necessary data
dge_results <- read.csv("goEnrichment_Thlaspi_2024-01-06.csv", header = TRUE)
colors <- c("blue", "black", "green", "magenta", "turquoise")

# Loop through each color
for (color in colors) {
  # Read the corresponding WGCNA results
  wgcna_results <- read.table(paste0("goEnrichment_", color, "_2024-01-07.csv"), sep = ",", header = TRUE)
  
  # Get unique GO IDs
  dge_ids <- unique(dge_results$GO.ID)
  wgcna_ids <- unique(wgcna_results$GO.ID)
  dge_ids <- dge_ids[!is.na(dge_ids)]
  wgcna_ids <- wgcna_ids[!is.na(wgcna_ids)]
  
  # Create the Venn diagram
  venn.plot <- venn.diagram(
    x = list(DGE = dge_ids, WGCNA = wgcna_ids),
    category.names = c("DGE", paste("WGCNA", color)),
    filename = NULL
  )
  
  # Define the PDF filename
  pdf_filename <- paste("venn_", color, ".pdf", sep = "")
  
  # Save the Venn diagram as a PDF (A4 size)
  pdf(pdf_filename, width = 8.27, height = 11.69)
  grid.draw(venn.plot)
  dev.off()
}

##Find overlap and write the identity of each overlap in a separate file:
for (color in colors) {
  # Read the corresponding WGCNA results
  wgcna_results <- read.table(paste0("goEnrichment_", color, "_2024-01-07.csv"), sep = ",", header = TRUE)
  
  # Get unique GO IDs
  dge_ids <- unique(dge_results$GO.ID)
  wgcna_ids <- unique(wgcna_results$GO.ID)
  dge_ids <- dge_ids[!is.na(dge_ids)]
  wgcna_ids <- wgcna_ids[!is.na(wgcna_ids)]
  
  # Calculate overlap
  overlap_ids <- intersect(dge_ids, wgcna_ids)
  overlap_count <- length(overlap_ids)
  total_genes <- length(unique(c(dge_ids, wgcna_ids)))
  dge_genes <- length(dge_ids)
  wgcna_genes <- length(wgcna_ids)
  
  # Perform hypergeometric test
  p_value <- phyper(overlap_count - 1, m = dge_genes, n = wgcna_genes, k = total_genes - dge_genes, lower.tail = FALSE)
  
  # Print the p-value
  cat("Hypergeometric test p-value for", color, ":", p_value, "\n")
  
  # Write the list of overlapping GO.IDs to a separate file
  overlap_filename <- paste("overlap_", color, ".txt", sep = "")
  write.table(overlap_ids, file = overlap_filename, col.names = FALSE, row.names = FALSE, quote = FALSE)
}



##Test to see which of the modules have the highest difference between treatments.

d1<-read.csv("gene_model_matrix.csv", row.names=1)
des_mat <- model.matrix(~ traits$treatment_bin)
fit <- limma::lmFit(t(d1), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(d1)) %>%
  tibble::rownames_to_column("module")
head(stats_df)
write.csv(stats_df, "limma_module_stat.csv")

traits2 <- design2 %>% 
  mutate(treatment_bin = ifelse(grepl('T3', file), 1, 0))%>%
  select(1,3)
row.names(traits2)<-traits2$file

# Define numbers of genes and samples
nSamples<-nrow(pseudoNormCounts)
nGenes<-ncol(pseudoNormCounts)

##Analysis of module - trait correlation.
module.trait.corr<-cor(MEs0[1:33], traits2$treatment_bin, use="p")
module.trait.corr.pvals<-corPvalueStudent(module.trait.corr, nSamples)

mod.trait.corr.pval<-module.trait.corr.pvals%>%
  as.data.frame()%>%
  arrange(V1)

write.csv(mod.trait.corr.pval, "module_trait_corr_pvals.csv")


# visualize module-trait association as a heatmap
heatmap.data <- merge(MEs0[1:32], traits2, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[34],
             y = names(heatmap.data)[1:32],
             col = c("blue1", "skyblue", "white", "pink", "red"))

##subset modules of interest
module<- c("turquoise", "black", "magenta", "green", "blue")
modules_of_interest_magenta = c("magenta")
modules_of_interest_green = c("green")
submod_green = module_df %>%
  subset(colors %in% modules_of_interest_green)
submod_magenta = module_df %>%
  subset(colors %in% modules_of_interest_magenta)

dd_green<-merge(submod_green, pseudoNormCounts, by = 'row.names')%>%
  select(2,4:22)%>%
  view()

dd_magenta<-merge(submod_magenta, pseudoNormCounts, by = 'row.names')%>%
  select(2,4:22)

dd2_magenta<-as.data.frame(t(dd_magenta))
dd2_green<-as.data.frame(t(dd_green))

dd2m<-as.data.frame(t(dd))



write.csv(dd2_green, "subset_green.csv")
write.csv(dd2_magenta, "subset_magenta.csv")

###### Calculate correlation and p vlaues for all the subset genes in all modules:
# List of colors
colors <- c("blue", "green", "magenta", "black", "turquoise")

# Create a list to store the correlation results and p-values
correlation_results <- list()

for (color in colors) {
  # Read the subset data
  subset_data <- read.csv(paste0("subset_", color, ".csv"), row.names = 1)
  
  # Compute the correlation
  corr <- cor(subset_data, traits$treatment_bin, use = "p")
  
  # Compute p-values
  sig_pvals <- corPvalueStudent(corr, nSamples) %>% as.data.frame()
  
  # Create a new data frame with two columns and copy data and column names
  result_df <- data.frame(locus_id = rownames(sig_pvals), pvals = sig_pvals[, 1])
  
  # Save p-values to CSV file
  write.csv(result_df, paste0("subset_corr_pvals_", color, ".csv"), row.names = FALSE)
  
  # Store correlation results
  correlation_results[[color]] <- result_df
}
###########



# Calculate the gene significance for all genes and associated p-values
gene.sig<-cor(input_mat, traits$treatment_bin, use="p")
gene.sig.pvals<-corPvalueStudent(gene.sig, nSamples)

sig_genes_thlaspi<-gene.sig.pvals%>%
  as.data.frame()%>%
  arrange(V1)

write.csv(sig_genes_thlaspi, "significant_coexpression_thlaspi.csv")

##Calculate the module membership and the associated p-values
# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.
module.membership.measures<-cor(MEs0[1:33], input_mat, use="p")
module.membership.pvals<-corPvalueStudent(module.membership.measures, nSamples)

d10<-module.membership.pvals%>%
  as.data.frame()

write.csv(d10, "mod.mem.pvals.csv" )

###Plot differences between the treatments based on linear regression results of different modules.
MEbrown_df <- d1 %>%
  tibble::rownames_to_column("file") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(traits %>%
                      dplyr::select(file, treatment_bin),
                    by = c("file" = "file")
  )

str(MEbrown_df)
MEbrown_df$treatment_bin<-as.factor(MEbrown_df$treatment_bin)
str(MEbrown_df)

##Plot the significant modules:

# Define the ME colors corresponding to the modules that are of interest
ME_colors <- c("MEturquoise", "MEblack", "MEmagenta", "MEgreen", "MEblue")
# Create a function to generate the plots
generate_plot <- function(y_var) {
  p <- ggplot(
    MEbrown_df,
    aes(
      x = treatment_bin,
      y = .data[[y_var]],
      color = treatment_bin
    )
  ) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +ylim(-0.5, 0.65)+
    ggforce::geom_sina(maxwidth = 0.3) +
    theme_classic()
  return(p)
}

# Create and display the plots for each ME color
plots <- lapply(ME_colors, generate_plot)

for (i in 1:length(ME_colors)) {
  print(plots[[i]])
}


for (i in 1:length(ME_colors)) {
  p <- generate_plot(ME_colors[i])
  pdf_file <- paste0("plot_", ME_colors[i], ".pdf")
  ggsave(filename = pdf_file, plot = p, width = 8.27, height = 11.69, units = "in", device = "pdf")
}


###
# Define a vector of color names
colors_to_filter <- c("turquoise", "black", "magenta", "green", "blue")

# Initialize an empty list to store the dataframes
filtered_dfs <- list()

# Loop through each color and process
for (color in colors_to_filter) {
  # Filter the dataframe for the current color
  filtered_df <- dplyr::filter(module_df, colors == color)
  
  # Add the filtered dataframe to the list
  filtered_dfs[[color]] <- filtered_df
  
  # Construct the filename based on the color
  file_name <- paste("me_", color, ".csv", sep = "")
  
  # Write the filtered dataframe to a CSV file
  write.csv(filtered_df, file_name)
}


#### Top hub in each module, top 20
# Modify the function to have a default top_n of 10
chooseTopHubsInEachModul <- function (datExpr, colorh, omitColors = "grey", power = 2, type = "signed", top_n = 10, ...) {
  isIndex = FALSE
  modules = names(table(colorh))
  if (!is.na(omitColors)[1]) 
    modules = modules[!is.element(modules, omitColors)]
  if (is.null(colnames(datExpr))) {
    colnames(datExpr) = 1:dim(datExpr)[2]
    isIndex = TRUE
  }
  hubs = vector("list", length(modules))
  names(hubs) = modules
  for (m in modules) {
    adj = adjacency(datExpr[, colorh == m], power = power, type = type, ...)
    hub_indices = order(rowSums(adj), decreasing = TRUE)[1:top_n]
    hubs[[m]] = colnames(adj)[hub_indices]
  }
  if (isIndex) {
    hubs = lapply(hubs, as.numeric)
    names(hubs) = modules
  }
  return(hubs)
}

topmodule<-chooseTopHubsInEachModul(
  input_mat,
  mergedColors,
  omitColors = "grey",
  power = 2,
  type = "signed") 


top_genes <- chooseTopHubsInEachModule(input_mat, mergedColors, top_n = 20)

####