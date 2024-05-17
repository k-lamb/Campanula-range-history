library(vcfR)
library(SNPfiltR )
library(adegenet)
library(ggplot2)
library(dplyr)
library(dartR)

### running into an issue where I need to decide which comes first: SNP or INDV filtering
### argument for SNP: an initial filter removes extremely bad sites which inflate indv, so indv is more meaningful
### argument for INDV: removing INDV first should boost sites retained by SNP filter

setwd("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/")

vcf <- read.vcfR("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/VCF/Cp_VCF/CpStacks_p25.vcf", convertNA=T)
genlight <- vcfR2genlight(vcf) #convert to genlight

inds <- genlight$ind.names
pops <- stringr::str_split(inds, "-")
pops <- sapply(pops, `[`, 1)

pop <- data.frame(id=inds, pop=pops)

# popmap<-popmap[popmap$id %in% colnames(vcfR@gt)[-1],]

pop_linreg <- read.csv("VCF/noTperf_idpoplin.txt", sep="\t")
linmap <- pop_linreg[,c(1,3)]
names(linmap)[2] <-"pop"

# looking at missingness
# missing_by_snp(vcf) # adjust up to 0.7 = 1459 SNPs


cutoff = 0.80
pca <- assess_missing_data_pca(vcf, popmap = pop, thresholds=c(cutoff), clustering=T) # 0.85,0.80,0.75,0.70,0.65 ... 0.65 looks best
# assess_missing_data_tsne(vcf, popmap = linmap, thresholds=c(cutoff)) # 0.85,0.80,0.75,0.70,0.65 ... 

# find real missing... not sure why it isn't working in SNPfiltR
vcf_cut <- missing_by_snp(vcf, cutoff=cutoff)
snps <- vcfR::extract.gt(vcf_cut, element = "GT", IDtoRowNames  = F, as.numeric = T, convertNA = T, return.alleles = F)
missing_sample <- colSums(is.na(snps))/dim(snps)[1]
pca <- pca[[1]]
pca$missing <- missing_sample
pca$pop_sep <- pop$pop

# and replot by missing... want to avoid high individual missing issues by removing those nearest to center
# ggplot(pca, aes(x=PC1, y=PC2, color=missing_sample))+
#   geom_point(size=3)+
#   theme_bw()

# find individuals with high missing 
missing <- pca %>%
  dplyr::group_by(rownames(.)) %>%
  dplyr::reframe(mean_miss_pop = mean(missing)) %>%
  arrange(mean_miss_pop) %>%
  filter(mean_miss_pop > 0.5) %>%
  as.data.frame()

# rejigger the vcf with the reduced set of inds
missing <- missing[,1]
sort(missing)
length(missing)

# write missing to file
write.table(missing, "VCF_filtering/Chloroplast_missing_list.txt", row.names = F, quote = F, col.names = F)

# remove missing from vcfR VCF
save <- vcf_cut

vcf <- save
vcf@gt <- subset(as.data.frame(vcf@gt), 
                 select = -which(names(as.data.frame(vcf@gt)) %in% missing)) %>% 
  as.matrix()

### futzing with cutoffs
# redo pca
vcf <- min_mac(vcf, min.mac = 1) # remove monomorphs
pop_redux <- subset(pop_linreg, id %in% (vcf@gt %>% as.data.frame %>% colnames())) # new popmap
pop_2 <- pop_redux[,c(1,3)]
names(pop_2)[2] <- "pop"

# cutoff = 0.75

# rechecking how it looks now with 50% individual missing cutoff
missing_by_snp(vcf)

# making new chloroplast pop map
popmap <- colnames(vcf@gt)
popmap <- popmap[-1] %>% as.data.frame()
names(popmap) <- "ind"

# separate out the pop map material
inds <- popmap
popmap <- popmap %>% separate(ind, into=c("pop", "ind"), sep="-")
popmap <- cbind(popmap$pop, inds)

write.table(popmap, "./VCF/CHLOROPLAST_vcf_popmap.txt", row.names=F, col.names=F, sep="\t", quote=F)

### figuring out how to rewrite to vcf using dartR because vcfR outputs are producing errors in other programs
### this section produces working but very weird VCF's that cause Moments to fail. moving the filtering scheme to VCFtools but leaving this to make the reason for departure clear
# 
# # gl2vcf direct from vcfR is failing because of number of items to replace error:
# # gl2vcf(x=genlight2, plink_path="~/Downloads/plink_mac_20231018/", outfile="noTperf_65miss_50indmiss", outpath="./VCF/deprecated/") # need to download and unzip plink from: https://www.cog-genomics.org/plink/
# 
# x <- gl2gi(genlight2) # gets vcfR object into dartR object
# pop(x) <- inds$pop # add population information
# x <- gi2gl(x) # convert back to genlight but this time in dartR format for genlight objects
# gl2vcf(x=x, plink_path="~/Downloads/plink_mac_20231018/", outfile="noTperf_65miss_50indmiss", outpath="./VCF") # need to download and unzip plink from: https://www.cog-genomics.org/plink/
