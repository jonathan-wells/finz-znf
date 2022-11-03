library("DESeq2")
library("apeglm")
library("RColorBrewer")
library("pheatmap")
library("tidyverse")
library("viridis")

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

load.whitecounts <- function(acc) {
  filename = paste("/Users/jonwells/Projects/feschottelab/finz-znf/data/expression/TEcount-out/", 
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

white.meta.df <- read.csv("~/Projects/feschottelab/finz-znf/data/expression/White2017/elife-30860-supp1-v1.tsv", 
                          sep="\t")

mo_accessions <- c("FM1", "FM2", "FM3", "SC1", "SC2", "SC3")
white_accessions <- filter(white.meta.df, 
                           stageName %in% c("50pc-epiboly", "Shield", "75pc-epiboly"), 
                           sequencing=="RNASeq")$accession_number
data <- c(lapply(mo_accessions, load.tecounts), lapply(white_accessions, load.whitecounts))
cts <- as.matrix(Reduce(cbind, data))


# Check that size of matrix is correct and delete list of data.tables if so
for (i in 1:length(c(mo_accessions, white_accessions))) {
  stopifnot(dim(data[[i]])[[1]] == dim(cts)[[1]])
}
rm(data)

# Load metadata
mo_coldata <- DataFrame(row.names = mo_accessions,
                        stage=factor(rep(c("Shield"), each=6)),
                        lab=factor(rep(c("Feschotte"), each=6)),
                        replicate=factor(c(1,2,3,1,2,3)),
                        condition=factor(rep(c("treatment", "control"), each=3)),
                        label=mo_accessions)
mo_coldata <- as.data.frame(mo_coldata)
white_coldata <- DataFrame(row.names = white_accessions,
                           stage=factor(rep(c("50pc_epiboly","Shield", "75pc_epiboly"),each=5)),
                           lab=factor(rep(c("White"), each=15)),
                           replicate=factor(c(4,5,6,7,8,4,5,6,7,8,4,5,6,7,8)),
                           condition=factor(rep(c("control"), each=15)),
                           label=white_accessions)
white_coldata <- as.data.frame(white_coldata)
coldata <- dplyr::bind_rows(mo_coldata, white_coldata)

stopifnot(all(rownames(coldata) == colnames(cts)))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ replicate + condition)

# Discard genes that have fewer than 10 reads in 3 or more samples.
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
dds <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=stage, shape = condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)


# Remove batch effects and plot
assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=vsd$replicate)
pcaData <- plotPCA(vsd, intgroup=c("stage", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=stage, shape=condition)) +
  geom_point(size=3) +
  xlim(-40, 40) + ylim(-40, 40) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)

res <- results(dds, contrast=c("condition", "treatment", "control"))
summary(res)
plotMA(res)

finz <- res[startsWith(rownames(res), "g"),]
plotMA(finz)

gene_dds <- dds[startsWith(rownames(res), "ENSDARG"),]
gene_vsd <- vst(gene_dds, blind=FALSE)
assay(gene_vsd) <- limma::removeBatchEffect(assay(gene_vsd), gene_vsd$replicate)
sampleDists <- dist(t(assay(gene_vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(gene_vsd$condition, gene_vsd$lab, gene_vsd$stage, gene_vsd$replicate, sep="-")
# rownames(sampleDistMatrix) <- gene_vsd$label
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "BuGn")) )(255)
pheatmap(sampleDistMatrix,
         clustering_method="ward.D",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

tes.df <- as.data.frame(res[!((startsWith(rownames(res), "ENSDARG")) | (startsWith(rownames(res), "g"))),]) %>% 
  rownames_to_column('gene') %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up-regulated", "down-regulated"),
         significant = ifelse(padj < 0.1, TRUE, FALSE)) %>%
  tidyr::separate(gene, c("tename", "tefam", "teclass"), sep=":") %>%
  mutate(retroelement = ifelse(teclass == "LINE" | teclass == "LTR", TRUE, FALSE),
         genetype=teclass,
         gene=tename) %>%
  merge(read.csv("~/Projects/feschottelab/drerio-tes/data/danrer11_te_summary.txt", 
                 sep="\t"),
        by=c("tename", "tefam", "teclass"))
ggplot(tes.df, aes(x=teclass, y=log2FoldChange)) + geom_boxplot()
ggplot(filter(tes.df, teclass=="LTR"), aes(x=baseMean, y=log2FoldChange, color=significant)) + 
  geom_hline(yintercept=0) + 
  geom_point(alpha=0.7) + 
  scale_x_log10()


###############################################################################
## Comparison of between-stage differences and treatment vs control diffs.
###############################################################################


white_cts <- as.matrix(Reduce(cbind, c(lapply(white_accessions, load.whitecounts))))

white_dds <- DESeqDataSetFromMatrix(countData = white_cts,
                                    colData = white_coldata,
                                    design = ~ stage)

# Discard genes that have fewer than 10 reads in 3 or more samples.
keep <- rowSums(counts(white_dds) >= 10) >= 3
white_dds <- white_dds[keep,]
white_dds <- DESeq(white_dds)

mo_cts <- as.matrix(Reduce(cbind, c(lapply(mo_accessions, load.tecounts))))

mo_dds <- DESeqDataSetFromMatrix(countData = mo_cts,
                                    colData = mo_coldata,
                                    design = ~ replicate + condition)

# Discard genes that have fewer than 10 reads in 3 or more samples.
keep <- rowSums(counts(mo_dds) >= 10) >= 3
mo_dds <- mo_dds[keep,]
mo_dds <- DESeq(mo_dds)


stage_tes.df <- as.data.frame(stage_res[!((startsWith(rownames(stage_res), "ENSDARG")) | (startsWith(rownames(stage_res), "g"))),]) %>% 
  rownames_to_column('gene') %>%
  tidyr::separate(gene, c("tename", "tefam", "teclass"), sep=":") %>%
  mutate(genetype="TE",
         gene=tename)

mo_tes.df <- as.data.frame(mo_res[!((startsWith(rownames(mo_res), "ENSDARG")) | (startsWith(rownames(mo_res), "g"))),]) %>% 
  rownames_to_column('gene') %>%
  tidyr::separate(gene, c("tename", "tefam", "teclass"), sep=":") %>%
  mutate(genetype="TE",
         gene=tename)

tes.df <- data.frame(mo_tes.df$gene,
                     mo_tes.df$tefam,
                     mo_tes.df$teclass,
                     mo_tes.df$log2FoldChange) %>%
  rename(gene=`mo_tes.df.gene`,
         tefam=`mo_tes.df.tefam`,
         teclass=`mo_tes.df.teclass`,
         `TE MO vs scrambled`=`mo_tes.df.log2FoldChange`) %>%
  merge(data.frame(stage_tes.df$gene, 
                   stage_tes.df$tefam,
                   stage_tes.df$teclass,
                   stage_tes.df$log2FoldChange) %>%
          rename(gene=`stage_tes.df.gene`,
                 tefam=`stage_tes.df.tefam`,
                 teclass=`stage_tes.df.teclass`,
                 `50% epiboly vs 75% epiboly`=`stage_tes.df.log2FoldChange`)) %>%
  melt(value.name = "log2 Fold Change")

stage_res <- results(white_dds, contrast=c("stage", "50pc_epiboly", "75pc_epiboly"))
stage_finz.df <- as.data.frame(stage_res[startsWith(rownames(stage_res), "g"),]) %>% 
  rownames_to_column('gene')

mo_res <- results(mo_dds, contrast=c("condition", "treatment", "control"))
mo_finz.df <- as.data.frame(mo_res[startsWith(rownames(mo_res), "g"),]) %>% 
  rownames_to_column('gene')

finz.df <- data.frame(mo_finz.df$gene, 
                      mo_finz.df$log2FoldChange) %>%
  rename(gene=`mo_finz.df.gene`,
         `FiNZ MO vs scrambled`=`mo_finz.df.log2FoldChange`) %>%
  merge(data.frame(stage_finz.df$gene, 
                   stage_finz.df$log2FoldChange) %>%
          rename(gene=`stage_finz.df.gene`,
                 `50% epiboly vs 75% epiboly`=`stage_finz.df.log2FoldChange`)) %>%
  melt(value.name = "log2 Fold Change")

head(finz.df)

ggplot(filter(tes.df, teclass=="LTR"), aes(x=variable, y=`log2 Fold Change`)) + 
  geom_boxplot() +
  xlab("")
ggplot(finz.df, aes(x=variable, y=`log2 Fold Change`)) + 
  geom_boxplot() +
  xlab("")