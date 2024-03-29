library(ggplot2)
library(phytools)
library(ape)
library(geiger)
library(nlme)
library(dplyr)
library(ggrepel)
library(stringr)


cypriniformes.tree <- read.tree(file="/Users/jonwells/Projects/feschottelab/finz-znf/data/species-phylogeny/iqtree-data/busco_supermatrix_partitions.nex.treefile")
cypriniformes.tree <- drop.tip(cypriniformes.tree, c('Chanos_chanos', 
                                                     "Triplophysa_tibetana", 
                                                     "Triplophysa_dalaica", 
                                                     "Triplophysa_siluroides"))
plot(cypriniformes.tree)



coverage.df <- read.csv("/Users/jonwells/Projects/feschottelab/finz-znf/data/finz_te_coverage.txt", sep='\t')
coverage.df <- coverage.df %>% filter(species %in% cypriniformes.tree$tip.label)
coverage.df$species <- factor(coverage.df$species, levels=cypriniformes.tree$tip.label)
coverage.df <- coverage.df[order(coverage.df$species),]


fit <- lm(coverage.df$finz_znf~coverage.df$interspersed)
fig1 <- ggplot(coverage.df, aes(x=interspersed, y=finz_znf, label=species)) + 
  geom_point() +
  geom_text_repel() +
  geom_smooth(method='lm', formula=y~x)
fig1
cor.test(coverage.df$finz_znf, coverage.df$interspersed, method="s")
cor.test(coverage.df$finz_znf, coverage.df$interspersed, method="p")


pic_coverage.df <- data.frame(pic.finz_exons = pic(coverage.df$finz_exons, cypriniformes.tree)) %>%
  mutate(pic.interspersed = pic(coverage.df$interspersed, cypriniformes.tree),
         pic.finz_znf = pic(coverage.df$finz_znf, cypriniformes.tree),
         pic.DNA = pic(coverage.df$DNA, cypriniformes.tree),
         pic.LTR = pic(coverage.df$LTR, cypriniformes.tree),
         pic.LINE = pic(coverage.df$LINE, cypriniformes.tree),
         pic.SINE = pic(coverage.df$SINE, cypriniformes.tree))
write.table(x=pic_coverage.df, 
            file="/Users/jonwells/Projects/feschottelab/finz-znf/data/finz_coverage_pics.txt", 
            row.names=F,
            sep="\t", 
            quote=F)

fig2 <- ggplot(pic_coverage.df, aes(x=pic.interspersed, y=pic.finz_znf)) +
  geom_point() 
fig2
cor.test(pic_coverage.df$pic.finz_znf, pic_coverage.df$pic.interspersed, method="s")
cor.test(pic_coverage.df$pic.finz_znf, pic_coverage.df$pic.interspersed, method="p")

#########################################################################################
## Phylogenetically independent contrasts between FiNZ-ZnF count and LTR families in 
## Danioninae species.
#########################################################################################

danioninae.tree <- extract.clade(cypriniformes.tree, 50)
plot(danioninae.tree)
  
counts.df <- read.csv("/Users/jonwells/Projects/feschottelab/finz-znf/data/finz_te_counts.txt", sep='\t')
counts.df <- counts.df %>% filter(species %in% danioninae.tree$tip.label)
counts.df$species <- factor(counts.df$species, levels=danioninae.tree$tip.label)
counts.df <- counts.df[order(counts.df$species),]


fit <- lm(counts.df$finz_znf~counts.df$LTR)
fig1 <- ggplot(counts.df, aes(x=LTR, y=finz_znf, label=species)) + 
  geom_point() +
  geom_text_repel() +
  geom_smooth(method='lm', formula=y~x)
fig1
cor.test(counts.df$finz_znf, counts.df$interspersed, method="s")
cor.test(counts.df$finz_znf, counts.df$interspersed, method="p")


pic_counts.df <- data.frame(pic.finz_exons = pic(counts.df$finz_exons, danioninae.tree)) %>%
  mutate(pic.interspersed = pic(counts.df$interspersed, danioninae.tree),
         pic.finz_znf = pic(counts.df$finz_znf, danioninae.tree),
         pic.DNA = pic(counts.df$DNA, danioninae.tree),
         pic.LTR = pic(counts.df$LTR, danioninae.tree),
         pic.LINE = pic(counts.df$LINE, danioninae.tree),
         pic.SINE = pic(counts.df$SINE, danioninae.tree))
View(pic_counts.df)
write.table(x=pic_counts.df, 
            file="/Users/jonwells/Projects/feschottelab/finz-znf/data/finz_counts_pics.txt", 
            row.names=F,
            sep="\t", 
            quote=F)

fig2 <- ggplot(pic_counts.df, aes(x=pic.interspersed, y=pic.finz_znf)) +
  geom_point() 
fig2
cor.test(pic_counts.df$pic.finz_znf, pic_counts.df$pic.interspersed, method="s")
cor.test(pic_counts.df$pic.finz_znf, pic_counts.df$pic.interspersed, method="p")

?corBrownian
