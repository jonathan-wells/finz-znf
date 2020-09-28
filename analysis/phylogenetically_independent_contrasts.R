library(ggplot2)
library(phytools)
library(dplyr)
library(ggrepel)

cypriniformes.tree <- read.tree(file="/Users/jonwells/Projects/feschottelab/finz-znf/data/species-phylogeny/supermatrix_partitions.nex.treefile")
cypriniformes.tree <- midpoint.root(cypriniformes.tree)
plot(cypriniformes.tree)

counts.df <- read.csv("/Users/jonwells/Projects/feschottelab/finz-znf/data/finz_te_counts.txt", sep='\t')

fit <- lm(counts.df$finz~counts.df$DNA)
fig1 <- ggplot(counts.df, aes(x=DNA, y=finz, label=species)) + 
  geom_point() +
  geom_text_repel() +
  geom_smooth(method='lm', formula=y~x)
fig1

pic.df <- data.frame(pic.finz = pic(counts.df$finz, cypriniformes.tree)) %>%
  mutate(pic.interspersed = pic(counts.df$interspersed, cypriniformes.tree),
         pic.DNA = pic(counts.df$DNA, cypriniformes.tree),
         pic.LTR = pic(counts.df$LTR, cypriniformes.tree),
         pic.LINE = pic(counts.df$LINE, cypriniformes.tree),
         pic.SINE = pic(counts.df$SINE, cypriniformes.tree)) 

fig2 <- ggplot(pic.df, aes(x=pic.DNA, y=pic.finz)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~x-1)
fig2

fit <- lm(formula = pic.df$pic.finz ~ pic.df$pic.DNA)
summary(fit)
cor.test(pic.df$pic.finz, pic.df$pic.DNA, method="s")
