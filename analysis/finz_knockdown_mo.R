library("DESeq2")
library("apeglm")
library("RColorBrewer")
library("pheatmap")
library("tidyverse")
library("reshape2")
library("GGally")

###############################################################################
## Differential Expression Analysis
###############################################################################

load.tecounts <- function(acc) {
  filename = paste("/Users/jonwells/Projects/feschottelab/finz-znf/data/finz-knockdown/FiNZ-MO-kd/TEcount-out/", 
                   acc, 
                   ".cntTable", 
                   sep='')
  read.table(filename,
             sep = "\t",
             skip = 1,
             row.names = 1,
             header = FALSE,
             col.names = c("Gene", acc))
} 


# Load TEcounts output into raw counts matrix
accessions <- c("FM1", "FM2", "FM3", "SC1", "SC2", "SC3")
data <- lapply(accessions, load.tecounts)
cts <- as.matrix(Reduce(cbind, data))

# Check that size of matrix is correct and delete list of data.tables if so
for (i in 1:length(accessions)) {
  stopifnot(dim(data[[i]])[[1]] == dim(cts)[[1]])
}
rm(data)

# Load metadata
coldata <- DataFrame(row.names = accessions,
                     condition=factor(rep(c("treatment","control"),each=3)),
                     replicate=factor(rep(c(1, 2, 3, 1, 2, 3))))
coldata <- as.data.frame(coldata)

# Check that sample accessions are present and in correct order.
stopifnot(all(rownames(coldata) == colnames(cts)))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ replicate + condition)

# Discard genes that have fewer than 10 reads in 3 or more samples.
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "treatment","control"))
summary(res)
plotMA(res)

###############################################################################
## Quality control
###############################################################################

# Plot PCA of non-batch corrected data
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape = replicate)) +
  geom_point(size=3) +
  xlim(-13, 13) +
  ylim(-13, 13) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)

# Remove batch effects and plot
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$replicate)
pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlim(-13, 13) +
  ylim(-13, 13) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)



sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$replicate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Get top 25 genes
top25 <- rownames(res[order(res$pvalue),][1:20,])
mat <- assay(vsd)[top25,]  # get the vsd for the top 25 genes
mat <- mat - rowMeans(mat)  # plot the deviation from the mean var stab count
pheatmap(mat,annotation_col = as.data.frame(colData(dds)[,c("condition","replicate")]))


###############################################################################
## Analysis of different groups
###############################################################################

# Load normalized counts and calculate means for control and treatment
cts.df <- as.data.frame(counts(dds, normalized=TRUE)) %>% 
  rownames_to_column('gene')
cts.df <- cts.df %>%
  mutate(controlMean = rowMeans(dplyr::select(cts.df, starts_with("SC"))),
         treatmentMean = rowMeans(dplyr::select(cts.df, starts_with("FM"))))

# Load external expression data from White et al. eLife (2017)
white.df <- read.csv("~/Projects/feschottelab/finz-znf/data/finz_expression_tpm.txt", sep="\t") %>%
  mutate(finz = ifelse(finz == "True", TRUE, FALSE),
         expressed = ifelse(expressed == "True", TRUE, FALSE))
whiteCasted.df <- dcast(dplyr::filter(white.df, stageName %in% c("50pc-epiboly", "Shield", "75pc-epiboly")), 
                        gene + cluster + finz + expressed ~ stageName, 
                        value.var="tpm") %>%
  mutate(shieldLog2FoldChange=log2(`75pc-epiboly`)-log2(`Shield`))

# Generate FiNZ dataframe      
finz.df <- as.data.frame(res[startsWith(rownames(res), "g"),]) %>% 
  rownames_to_column('gene') %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up-regulated", "down-regulated"),
         genetype="finz") %>%
  merge(dplyr::select(whiteCasted.df, gene, shieldLog2FoldChange, Shield, cluster), by="gene")

# Generate TE dataframe
tes.df <- as.data.frame(res[!((startsWith(rownames(res), "ENSDARG")) | (startsWith(rownames(res), "g"))),]) %>% 
  rownames_to_column('gene') %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up-regulated", "down-regulated")) %>%
  tidyr::separate(gene, c("tename", "tefam", "teclass"), sep=":") %>%
  mutate(retroelement = ifelse(teclass == "LINE" | teclass == "LTR", TRUE, FALSE),
         genetype="TE",
         gene=tename) %>%
  merge(read.csv("~/Projects/feschottelab/drerio-tes/data/danrer11_te_summary.txt", 
                 sep="\t"),
        by=c("tename", "tefam", "teclass")) %>%
  merge(dplyr::select(whiteCasted.df, gene, shieldLog2FoldChange, Shield, cluster), by="gene")

# Generate genes dataframe
genes.df <- as.data.frame(res[startsWith(rownames(res), "ENSDARG"),]) %>% 
  rownames_to_column('gene') %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up-regulated", "down-regulated"),
         genetype = "non-finz gene") %>%
  merge(dplyr::select(whiteCasted.df, gene, shieldLog2FoldChange, Shield, `50pc-epiboly`, `75pc-epiboly`, cluster), by="gene")

# Combine dataframes and save output
dds.df <- dplyr::bind_rows(finz.df, genes.df, tes.df)
write.table(dds.df, "~/Projects/feschottelab/finz-znf/data/finz_mo_kd_summary.txt", 
            sep="\t",
            quote = FALSE, 
            row.names = FALSE)


