library(DESeq2)
library(apeglm)


load.tecounts <- function(acc) {
  filename = paste("/Users/jonwells/Projects/feschottelab/finz-znf/data/expression/TEcount-out/", acc, ".cntTable", sep='')
  read.table(filename,
             sep = "\t",
             skip = 1,
             row.names = 1,
             header = FALSE,
             col.names = c("Gene", acc))
} 

# Load count matrix
accessions <- c("ERS1079234", "ERS1079235", "ERS1079236", "ERS1079209", "ERS1079210", "ERS1079211")
data <- lapply(accessions, load.tecounts)
cts <- as.matrix(Reduce(cbind, data))

# Check that size of matrix is correct and delete list of data.tables if so
for (i in 1:length(accessions)) {
  stopifnot(dim(data[[i]])[[1]] == dim(cts)[[1]])
}
rm(data)

# Load metadata
coldata <- DataFrame(row.names = accessions,
                     condition=factor(rep(c("treated","control"),each=3)),
                     stage=factor(rep(c("1k-cell", "Dome"),each=3)))
coldata <- as.data.frame(coldata)

# Check that sample accessions are present and in correct order.
stopifnot(all(rownames(coldata) == colnames(mat)))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ stage)
dds

dds <- DESeq(dds)
res <- results(dds)
res

summary(res)
plotMA(res)

finz <- res[startsWith(rownames(res), "g"),]
tes <- res[!((startsWith(rownames(res), "ENSDARG")) | (startsWith(rownames(res), "g"))),]
dna <- res[endsWith(rownames(res), "DNA"),]
plotMA(dna)
plotMA(finz)
plotMA(tes)
