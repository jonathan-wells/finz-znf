library(ggplot2)
library(phytools)
library(ape)
library(geiger)
library(nlme)
library(dplyr)
library(ggrepel)
library(stringr)


cypriniformes.tree <- read.tree(file="/Users/jonwells/Projects/feschottelab/finz-znf/data/species-phylogeny/iqtree-data/busco_supermatrix_partitions.nex.timetree.nwk")
cypriniformes.tree <- midpoint.root(cypriniformes.tree)

keepspecies <- read.csv('/Users/jonwells/Projects/feschottelab/finz-znf/data/hiqual_species.txt', header=FALSE)
currspecies = cypriniformes.tree$tip.label
for (tip in currspecies) {
  if (!(tip %in% keepspecies$V1)) {
    cypriniformes.tree <- drop.tip(cypriniformes.tree, tip)   
    }
}

plot(cypriniformes.tree)

counts.df <- read.csv("/Users/jonwells/Projects/feschottelab/finz-znf/data/finz_te_counts.txt", sep='\t')
counts.df <- counts.df %>% filter(species %in% cypriniformes.tree$tip.label)
counts.df$species <- factor(counts.df$species, levels=cypriniformes.tree$tip.label)
counts.df <- counts.df[order(counts.df$species),]


fit <- lm(counts.df$finz_znf~counts.df$interspersed)
fig1 <- ggplot(counts.df, aes(x=interspersed, y=finz_znf, label=species)) + 
  geom_point() +
  geom_text_repel() +
  geom_smooth(method='lm', formula=y~x)
fig1
cor.test(counts.df$finz_znf, counts.df$interspersed, method="s")


pic.df <- data.frame(pic.finz_exons = pic(counts.df$finz_exons, cypriniformes.tree)) %>%
  mutate(pic.interspersed = pic(counts.df$interspersed, cypriniformes.tree),
         pic.finz_znf = pic(counts.df$finz_znf, cypriniformes.tree),
         pic.DNA = pic(counts.df$DNA, cypriniformes.tree),
         pic.LTR = pic(counts.df$LTR, cypriniformes.tree),
         pic.LINE = pic(counts.df$LINE, cypriniformes.tree),
         pic.SINE = pic(counts.df$SINE, cypriniformes.tree))
write.table(x=pic.df, 
            file="/Users/jonwells/Projects/feschottelab/finz-znf/data/finz_pics.txt", 
            row.names=F,
            sep="\t", 
            quote=F)

fig2 <- ggplot(pic.df, aes(x=pic.interspersed, y=pic.finz_znf)) +
  geom_point() 
fig2
cor.test(pic.df$pic.finz_znf, pic.df$pic.interspersed, method="s")

fit <- lm(formula = pic.df$pic.finz_exons~pic.df$pic.interspersed)
summary(fit)
cor.test(pic.df$pic.finz_exons, pic.df$pic.interspersed, method="s")


  