
# for ancestry analysis and plotting
library(LEA)
library(vcfR)
library(adegenet)
library(tidyverse)
library(ggplot2)
library(vegan)

# for isolation by distance
library(tidyverse)
library(dartR)
library(MASS)
library(usedist)
library(SpatialEpi)

# for mapping
library(rnaturalearth) #for downloading shapefiles
library(maps)
library(sf)
library(scatterpie) # for ancestry coefficients pie charts
library(tess3r) # for interpolating ancestry predictions
library(ggpubr) # for ggarrange of maps
library(unikn) # for plotting labels
library(rworldmap) # for final mapping using tess3r
library(ggnewscale) # the K cluster colors won't print properly on final tess3r map without resetting scale
library(rworldxtra) # for higher quality USA map
library(devtools)

setwd("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/")

# read in and transform data
# vcf <- read.vcfR("./VCF/noTperf_65miss_50miss_FULL.recode.vcf") # 65% missing by SNP
vcf <- read.vcfR("./VCF/noTperf_80miss_50miss_FULL.recode.vcf") # 80% missing by SNP. 3568 SNP

pop_linreg <- read.csv("./VCF/noTperf_idpoplin.txt", sep="\t")
pop_linreg <- subset(pop_linreg, id %in% (vcf@gt %>% as.data.frame %>% colnames())) %>%
  tidyr::separate(lin, into=c("lin", "reg"), sep="_")

# # can leave commented if running off existing project
genlight <- vcfR2genlight(vcf) #convert to genlight

inds <- genlight$ind.names
pops <- stringr::str_split(inds, "-")
pops <- sapply(pops, `[`, 1)

pop(genlight) <- pops

# lfmm formatted data without using LEA's buggy converter
dartR::gl2geno(genlight, outfile="./VCF/gl_geno_80", outpath="~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/")

# ancestry for all populations
pc <- snmf("./VCF/gl_geno_80.lfmm", K=1:50, entropy=T, repetitions=10) # DO NOT RUN AGAIN. LOAD FROM PROJECT

#loading saved files from previous LEA run
# pc <- load.snmfProject("./VCF/gl_geno_80.snmfProject")

plot(pc, col = "blue", pch = 19, cex = 1.2, pty="s") # if this looks bad when K > pop, limit maxK = pop


# checking number of clusters
# adegenet::find.clusters(genind, max.n.clust = 50)
dapc <- adegenet::dapc.genind(x = genind, pop = genlight2@pop, n.pca = 50, n.da = 100)
scatter(dapc)
adegenet::optim.a.score(dapc)
dapc <- adegenet::dapc.genind(x = genind, pop = genlight2@pop, n.pca = 6, n.da = 100)
scatter(dapc)

# cross-entropy slope quantiles
reps <- 10
maxK <- unique(pop_linreg$pop) %>% length() %>% as.numeric() # number of populations as maxK cap
maxK = 50 # K-clusters explored in total

# determining the optimal K clusters
# code block from: https://chazhyseni.github.io/NALgen/post/determining_bestk/
ce <- list()
for(k in 1:maxK) ce[[k]] <- cross.entropy(pc, K=k)
ce.K <- c()
for(k in 1:maxK) ce.K[k] <- min(ce[[k]])
diff <- ce.K[-1] - ce.K[-maxK]
slope <- exp(-diff) - 1
best.K <- min(which(slope <= quantile(slope)[4])) %>% as.numeric()

# best.K by minimum cross entropy value
best.K = min(which(ce.K == min(ce.K)))
