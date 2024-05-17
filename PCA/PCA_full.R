library(vcfR)
library(adegenet)
library(tidyverse)
library(SNPfiltR)
library(dartR) #dartR requires SNPRelate from BiocManager
library(FactoMineR) #for PCA tools
library(adegenet)
library(dartR)

set.seed(1066)

setwd("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/")

### DAPC

# vcf <- read.vcfR("./VCF/noTperf_65miss_50miss_FULL.recode.vcf") # 65% SNP missing
vcf <- read.vcfR("VCF/noTperf_80miss_50miss_FULL.recode.vcf") # 80% SNP missing 
vcf <- SNPfiltR::min_mac(vcf, min.mac=1) # remove invariant sites. should have 0% loss
genlight <- vcfR2genlight(vcf) #convert to genlight
pop_linreg <- read.csv("VCF/noTperf_idpoplin.txt", sep="\t")
pop_linreg <- subset(pop_linreg, id %in% (vcf@gt %>% as.data.frame %>% colnames()))
pop_groups_full <- pop_linreg[,3]

inds <- genlight$ind.names
pops <- stringr::str_split(inds, "-")
pops <- sapply(pops, `[`, 1)

pop(genlight) <- as.factor(pops)

## transform to genind -- 2001 format
genind <- dartR::gl2gi(genlight)



### PCA 
pca_obj <- 
  genlight %>%
  tab() %>%
  as.data.frame() %>%
  PCA(scale.unit = F,
      graph = F)

pca_proj_met <-
  pca_obj %>%
  .$ind %>%
  .$coord %>% 
  as.data.frame() %>%
  mutate(id = rownames(.))

pca_proj_met <- merge(pca_proj_met, pop_linreg, by="id")

pca_proj_met <- pca_proj_met %>% 
  separate(col = lin, into = c("lineage", "CZ"), sep = "_")


# PCA no labels
ggplot(pca_proj_met, aes(x=Dim.1, y=Dim.2)) +
  geom_point(aes(color = lineage, shape = CZ), size = 4) +
  scale_color_manual(values=c("green4", "skyblue3", "purple4")) +
  # scale_shape_manual(values=c(15, 19, 17, 8, 15, 8, 15, 19, 17))+
  ggtitle("Lineage and Contact Zone") +
  theme_bw()

pca_avg <- pca_proj_met %>% 
  group_by(pop, CZ, lineage) %>% 
  reframe(mean_PC1 = mean(Dim.1),
          mean_PC2 = mean(Dim.2))

# PCA labels
jpeg("PCA/PCA_ALL_80perc.jpeg", res=200, width=7, height=7, unit="in")
ggplot(pca_proj_met %>% filter(pop != "NC109B" & pop != "NC109C" & pop != "NC109D" & pop != "NC90"), 
       aes(x=Dim.1, y=Dim.2, label=pop)) +
  geom_point(aes(color = lineage, shape = CZ), size = 2, alpha = 0.5) +
  ggrepel::geom_label_repel(data=pca_avg %>% filter(pop != "NC109B" & pop != "NC109C" & pop != "NC109D" & pop != "NC90"), 
                            aes(x = mean_PC1, y = mean_PC2, color = lineage), size = 3, max.overlaps = 15) +
  scale_color_manual(values=c("green4", "skyblue3", "purple4")) +
  # ggtitle("Lineage and Contact Zone") +
  xlab("PC1")+
  ylab("PC2")+
  theme_bw()+
  theme(text = element_text(size = 18), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio=1)
dev.off()
