
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

setwd("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/")

# read in and transform data
# vcf <- read.vcfR("./VCF/noTperf_65miss_50miss_FULL.recode.vcf") # 65% missing by SNP
vcf <- read.vcfR("./VCF/noTperf_80miss_50miss_FULL.recode.vcf") # 80% missing by SNP. 3568 SNP

pop_linreg <- read.csv("./VCF/noTperf_idpoplin.txt", sep="\t")
pop_linreg <- subset(pop_linreg, id %in% (vcf@gt %>% as.data.frame %>% colnames())) %>%
  tidyr::separate(lin, into=c("lin", "reg"), sep="_")

# gene flow estimates from Dsuite and Dsuite_filter_basicplots.R output
gene.bbaa <- read.csv("gene_flow/sig_bbaa_reg_CLDmethod.csv", stringsAsFactors = T) # no SNP missing filter
gene.bbaa <- gene.bbaa %>%
  group_by(P2,P3, lin_col) %>%
  reframe(no.sig=n(), mean.Dstat=mean(Dstatistic),
          mean.Zscore=mean(Z.score), mean.f4=mean(f4.ratio)) # , mean=mean(p.value.adj)

# write.csv(gene.bbaa, "./gene_flow/80perc/mean_Dstat.csv")

# # can leave commented if running off existing project
genlight <- vcfR2genlight(vcf) #convert to genlight

inds <- genlight$ind.names
pops <- stringr::str_split(inds, "-")
pops <- sapply(pops, `[`, 1)

pop(genlight) <- pops

#loading saved files from previous LEA run
pc <- load.snmfProject("./VCF/gl_geno_80.snmfProject")



### Isolation-by-distance as cause of the gradual decay in cross-entropy
# do additional clusters indicate IBD? common source of gradual decay as seen here
# genlight2 <- gl.read.vcf("./VCF/noTperf_65miss_50miss_FULL.recode.vcf") # 65% SNP missing
genlight2 <- gl.read.vcf("./VCF/noTperf_80miss_50miss_FULL.recode.vcf") # 80% SNP missing
pop(genlight2) <- pops
genlight2 <- gl.drop.pop(genlight2, c("NC109B", "NC109C", "NC109D", "NC90", "VA96", "VA71"))
genind <- gl2gi(genlight2)

# checking number of clusters
# adegenet::find.clusters(genind, max.n.clust = 50)
dapc <- adegenet::dapc.genind(x = genind, pop = genlight2@pop, n.pca = 50, n.da = 100)
scatter(dapc)
adegenet::optim.a.score(dapc)
dapc <- adegenet::dapc.genind(x = genind, pop = genlight2@pop, n.pca = 6, n.da = 100)
scatter(dapc)

# check for private alleles
pa <- gl.report.pa(genlight2, method="one2rest") # check that private alleles exist
pa$priv1 %>% mean()
# check for IBD within Western and Appalachian separately

### WESTERN
west_popdrop <- pop_linreg %>% filter(lin != "Western")
w <- pop_linreg %>% filter(lin == "Western")
genlight_west <- gl.drop.pop(genlight2, unique(west_popdrop$pop))

# genepop <- gl2genepop(genlight2)
genind <- gl2gi(genlight_west)
genepop <- genind2genpop(genind)

gen.w <- dist.genpop(genepop, method=2) #get genetic distance using chord distance (Cavalli-Sforza & Edwards 1967)
lin.gen.w <- gen.w/(1-gen.w)

locations <- read.csv("sNMF_structure/group_coords.csv", stringsAsFactors = T)
coords <- locations %>%
  filter(group %in% pop_linreg$pop) %>%
  dplyr::select(-lineage) %>%
  filter(group %in% w$pop & group != "NC109B" & group != "NC109C" & group != "NC109D" & group != "VA96" & group != "VA71" & group != "NC90")
rownames(coords) <- coords$group
coords <- coords %>% dplyr::select(long, lat)
coord.km <- SpatialEpi::latlong2grid(coords)
rownames(coord.km) <- rownames(coords)

dist.geo.w <- dist(coord.km)
log.dist.geo.w <- log10(dist.geo.w)

vegan::mantel(lin.gen.w, log.dist.geo.w) # 80%: R=0.3707 p=0.008

  dens <- kde2d(as.vector(log.dist.geo.w), as.vector(lin.gen.w), n=300)
  myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
  plot(log.dist.geo.w, lin.gen.w, pch=20,cex=.5)
  image(dens, col=transp(myPal(300),.7), add=TRUE)
  abline(lm(as.vector(lin.gen.w)~as.vector(log.dist.geo.w)))

### APPALACHIAN
app_popdrop <- pop_linreg %>% filter(lin != "Appalachian")
a <- pop_linreg %>% filter(lin == "Appalachian")
genlight_app <- gl.drop.pop(genlight2, unique(app_popdrop$pop))

# genepop <- gl2genepop(genlight2)
genind <- gl2gi(genlight_app)
genepop <- genind2genpop(genind)

gen.a <- dist.genpop(genepop, method=2) #get genetic distance
lin.gen.a <- gen.a/(1-gen.a)

locations <- read.csv("sNMF_structure/group_coords.csv", stringsAsFactors = T)
coords <- locations %>%
  filter(group %in% pop_linreg$pop) %>%
  dplyr::select(-lineage) %>%
  filter(group %in% a$pop & group != "NC109B" & group != "NC109C" & group != "NC109D" & group != "VA96" & group != "VA71")
rownames(coords) <- coords$group
coords <- coords %>% dplyr::select(long, lat)
coord.km <- SpatialEpi::latlong2grid(coords)
rownames(coord.km) <- rownames(coords)

dist.geo.a <- dist(coord.km)
log.dist.geo.a <- log10(dist.geo.a)

vegan::mantel(lin.gen.a, log.dist.geo.a) # 80%: R=0.4613 p=0.001

dens <- kde2d(as.vector(log.dist.geo.a), as.vector(lin.gen.a), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(log.dist.geo.a, lin.gen.a, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(lin.gen.a)~as.vector(log.dist.geo.a)))

### EASTERN
east_popdrop <- pop_linreg %>% filter(lin != "Eastern")
e <- pop_linreg %>% filter(lin == "Eastern")
genlight_east <- gl.drop.pop(genlight2, unique(east_popdrop$pop))

# genepop <- gl2genepop(genlight2)
genind <- gl2gi(genlight_east)
genepop <- genind2genpop(genind)

gen.e <- dist.genpop(genepop, method=2) #get genetic distance
lin.gen.e <- gen.e/(1-gen.e)

locations <- read.csv("sNMF_structure/group_coords.csv", stringsAsFactors = T)
coords <- locations %>%
  filter(group %in% pop_linreg$pop) %>%
  dplyr::select(-lineage) %>%
  filter(group %in% e$pop & group != "NC109B" & group != "NC109C" & group != "NC109D" & group != "VA96" & group != "VA71")
rownames(coords) <- coords$group
coords <- coords %>% dplyr::select(long, lat)
coord.km <- SpatialEpi::latlong2grid(coords)
rownames(coord.km) <- rownames(coords)

dist.geo.e <- dist(coord.km)
log.dist.geo.e <- log10(dist.geo.e)

vegan::mantel(lin.gen.e, log.dist.geo.e) # 80%: R=0.9974, p=0.167

dens <- kde2d(as.vector(log.dist.geo.e), as.vector(lin.gen.e), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(log.dist.geo.e, lin.gen.e, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(lin.gen.e)~as.vector(log.dist.geo.e)))

