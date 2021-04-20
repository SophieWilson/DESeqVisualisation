
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(NMF)
library(ggplot2)
library(EnhancedVolcano)
library(tikzDevice)

setwd('C:/Users/Mischa/Documents/Uni Masters/Module 6 - Group proj/Blogpost')

# Read the data into R
seqdata <- read.delim("GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("SampleInfo_Corrected.txt", stringsAsFactors = TRUE)

head(seqdata)
head(sampleinfo)

seq_counts <- seqdata[,-c(1,2)]
rownames(seq_counts) <- seqdata[,1]
# Taking only the first 7 characters of column name
colnames(seq_counts) <- sampleinfo$SampleName
head(seq_counts)

DGEcount <- DGEList(seq_counts)
# Look at the data
head(DGEcount)

group <- factor(paste(sampleinfo$CellType, sampleinfo$Status, sep='_'))
DGEcount$samples$group <- group
head(DGEcount)

keep <- filterByExpr(DGEcount, group = group, min.total.count=10)

ggplot(DGEcount$samples, aes(x = colnames(DGEcount), y=lib.size)) + geom_bar(stat = 'identity', fill = '#00abff', color = 'black') + 
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('Sample ID') +
  ylab('Library size (millions)') +
  ggtitle('Library size')

# Obtain log counts
log_count <- cpm(DGEcount, log=TRUE)
# Save median value for the abline
abline <- median(log_count)
# Plotting the counts with an abline and data labels.
Unnormalised_boxplot <- ggplot(stack(as.data.frame(log_count)), aes(x = ind, y = values)) +
  geom_boxplot(fill = '#00abff') +
  theme(axis.text.x = element_text(angle = 90)) + 
  geom_hline(yintercept = abline, color = 'red') +
  xlab('Sample ID') +
  ylab('Counts per million (log2)') +
  ggtitle('LogCMP values') 
Unnormalised_boxplot

## Data Visualisation

## Multidimensional Scaling Plot  ##############
plotMDS(DGEcount)
# Allowing side-by-side plots of Cell type and Status
par(mfrow=c(1,2))
# Setting cell colours for each group
col_celltype <- c("#9932CC", "#FDD20EFF")[sampleinfo$CellType]
col_status <- c("#00AFBB", "#E7B800", "#FC4E07")[sampleinfo$Status]
# Plotting cell type MDS graph
plotMDS(DGEcount, col=col_celltype)
legend("topleft",fill=c("#9932CC","#FDD20EFF"),legend=levels(sampleinfo$CellType))
title("Cell type")
# Plotting sample Status MDS graph
plotMDS(DGEcount, col=col_status)
legend("top",fill=c("#00AFBB","#E7B800","#FC4E07"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

colors <- c("#00AFBB", "#E7B800", "#FC4E07")
colors <- colors[as.numeric(sampleinfo$Status)]

shapes = c(16, 17) 
shapes <- shapes[as.numeric(sampleinfo$CellType)]
par(mfrow=c(1,1))
plotMDS(DGEcount,col=colors,cex=2, pch =shapes)
legend("top",legend=levels(sampleinfo$Status), col = c("#00AFBB", "#E7B800", "#FC4E07") ,pch=16)
legend("bottom",legend=levels(sampleinfo$CellType),pch=c(16,17))

############ PCA ################
# # transforming it so that samples are on the row and chemicals are columns
pca_x <- t(seq_counts)
pca_x <- pca_x[ , which(apply(pca_x, 2, var) != 0)]
# # doing the pca
pca_res1 <- prcomp(pca_x, center = TRUE, scale = TRUE)
# making the dataframe
pc_df1 <- data.frame(pca_res1$x, Celltype = sampleinfo$CellType, Status = sampleinfo$Status)
# plotting it
library(RColorBrewer)
mycolours <- c("#00AFBB", "#E7B800", "#FC4E07")
ggplot(pc_df1, aes(PC1, PC2, shape = Celltype)) + 
  geom_point(aes(fill = Status), size = 6) + ggtitle('PCA of Cell type and Status') + scale_shape_manual(values=c(21,24)) + guides(fill = guide_legend(override.aes = list(shape = 21)))
############################



# Estimatinve variance per row in log_count matrix
variable_genes <- apply(log_count, 1, var)
# Take gene names of the 500 most variable genes
selected_var <- names(sort(variable_genes, decreasing=TRUE))[1:500]
# Subset logcounts matrix to take only the top 500 genes
highly_variable_genes <- log_count[selected_var,]
# Plot the data
heatmap(highly_variable_genes)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlGn")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","orange")[sampleinfo$CellType]
# Reducing the number of genes in the graph
highly_variable_genes_100 <- log_count[names(sort(variable_genes, decreasing=TRUE))[1:100],]
# plotting
heatmap(highly_variable_genes_100, col=rev(morecols(50)), main="Top 100 most variable genes across samples split by cell type",ColSideColors=col.cell,scale="row")

# Specify colors
ann_colors = list(CellType = c('purple', 'orange'), Status = c('#28334AFF', '#FBDE44FF', '#F65058FF'))
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Plot
par(mfrow=c(1,1))
aheatmap(highly_variable_genes_100,col=rev(morecols(50)),main="
         100 most variable genes across samples",annCol=sampleinfo[, 3:4],annColors=ann_colors,labCol=group, scale="row")


# First we must normalise for composition bias
DGEcount <- calcNormFactors(DGEcount)
# Make a design matrix with no intercept term
designmat <- model.matrix(~ 0 + group)
colnames(designmat) <- levels(group)

# Normalisation
transformed <- voom(DGEcount, designmat, plot = TRUE)
names(transformed)


par(mfrow=c(1,2))

# Test for normalisation with a boxplot
ggplot(stack(as.data.frame(transformed)), aes(x = ind, y = values)) +
  geom_boxplot(fill = '#00abff') +
  theme(axis.text.x = element_text(angle = 90)) + 
  geom_hline(yintercept = abline, color = 'red') +
  xlab('Sample ID') +
  ylab('Counts per million (log2)') +
  ggtitle('LogCMP values') 
# Plot the earlier boxplot to compare
Unnormalised_boxplot


# Fitting the linear model
lm <- lmFit(transformed)
names(lm)
# Applying make Contrasts
cont_df <- makeContrasts(Basal_lacVsPred=basal_pregnant - basal_lactate,levels=designmat)
# Viewing the contrast dataframe
cont_df

The contrast dataframe tells the limma package which columns of the design matrix to use in the comparison test.

# Apply the contrast dataframe to the linear model to get the comparison
fit_contrast <- contrasts.fit(lm, cont_df)
dim(fit_contrast)
# Use empirical Bayes shrinkage of variance to better plot
fit_contrast <- eBayes(fit_contrast)

fit_summary <- decideTests(fit_contrast)
summary(fit_summary)

plotMD(fit_contrast, coef = 1, status = fit_summary, values = c(-1, 1), hl.col=c("blue","red"))

volcanoplot(fit_contrast,coef=1,highlight=100,names=fit_contrast$genes$SYMBOL, main="Basal_lacVsPred")
volcanoplot(fit_contrast, main="Basal_lacVsPred Gene expression")

dge <- DGEList(round(seq_counts), group = sampleinfo$Status)
dds <- as.DESeqDataSet(dge)
sdd <- DESeq(dds)
res <- results(sdd)
resultsNames(sdd)
resLFC2 <- lfcShrink(sdd, coef = 'group_pregnant_vs_lactate', type='apeglm')
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')


res <- lfcShrink(sdd, coef = 'group_pregnant_vs_lactate', type='apeglm')
EnhancedVolcano(resLFC2,
                lab = rownames(resLFC2),
                x = 'log2FoldChange',
                y = 'pvalue')
