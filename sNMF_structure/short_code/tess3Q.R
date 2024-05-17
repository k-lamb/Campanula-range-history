
### if starting from LEA/snmf, script requires these libraries:
# # for ancestry analysis and plotting
# library(LEA)
# library(vcfR)
library(tidyverse) # can also suffice with just: library(dplyr) and library(tidyr)

### if starting from pre-written data frame, requires only these libraries:
# for mapping
library(ggplot2)
library(rnaturalearth) # for USA and world downloading shapefiles
library(sf) # for state and lake shapefiles
library(scatterpie) # for ancestry coefficients pie charts
library(tess3r) # for interpolating ancestry predictions
library(ggpubr) # for ggarrange of maps
library(ggnewscale) # the K cluster colors won't print properly on final tess3r map without resetting scale
library(ggrepel) # repels population name from pie charts so you can see them easier

setwd("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/")

### if loading data generated elsewhere:
# need a data frame with 3+Kn columns: (1) id, (2) pop, (3) lineage (lin), and all K cluster proportions
# can also start downstream with K clusters per population. if so, load locations, then skip to plotting
# row sums should equal 1 for all K clusters per individual
q_mat_all <- read.csv("~/Downloads/tess3Q_df.csv", stringsAsFactors = T) # note: region (reg) is not required for this script to run

# import 3 column data frame with coordinates: (1) pop, (2) lat, (3) long
locations <- read.csv("sNMF_structure/group_coords.csv", stringsAsFactors = T) %>% dplyr::select(-lineage)
names(locations)[1] <- "pop"

# ### if pulling data from snmf in LEA:
# # loading data
# # the VCF is here just to strip out the names of individuals from the total subset which will be used for structure
# # pop_linreg is a data frame with the col id, pop, lineage, and reg
# vcf <- read.vcfR("./VCF/noTperf_80miss_50miss_FULL.recode.vcf") # 80% missing by SNP
# pop_linreg <- read.csv("./VCF/noTperf_idpoplin.txt", sep="\t")
# pop_linreg <- subset(pop_linreg, id %in% (vcf@gt %>% as.data.frame %>% colnames())) %>%
#   tidyr::separate(lin, into=c("lin", "reg"), sep="_")
# 
# #loading saved files from previous LEA run
# best.K <- 6
# pc <- load.snmfProject("./VCF/gl_geno_80.snmfProject")
# best = which.min(cross.entropy(pc, K = best.K)) # best run within K
# 
# # oroduce a Q matrix with individuals assigned to population, lineage, and region
# # q_mat is a data frame with dimension equal to the number of individuals (y) by the number of clusters (x)
# q_mat <- LEA::Q(pc, K = best.K, run = best) %>% as.data.frame() # all row sums=1.0. can check with: q_mat %>% rowSums()
# colnames(q_mat) <- paste0("K", 1:ncol(q_mat)) # rename columns in q matrix to K1:Kn
# q_mat_all <- cbind(pop_linreg, q_mat) # cbinds information of individual ID, POPULATION, and LINEAGE to each Q matrix row


### if starting from pre-written data frame:
# need to: (1) drop unwanted populations (too few individuals); (2) group q matrix values by population. Using population mean of K-cluster
# drop populations with 3 or fewer samples
q_mat_all <- q_mat_all %>%
  dplyr::select(-id) %>%
  group_by(pop, lin) %>%
  filter(n() > 3) %>% # drop populations with fewer than 3 individuals
  dplyr::summarise(across(everything(), mean)) # find mean ancestry per population from individuals in those populations



###############################################################################################################
# PREP FOR PLOTTING

# retain only the q matrix (no id, pop, etc. columns allowed in ggtess3Q)
q_redux <- q_mat_all %>%
  ungroup() %>% # won't unselect columns without first ungrouping
  dplyr::select(-c(pop, lin))

q_mat_all <- merge(locations, q_mat_all, by="pop") # add lat/long to q matrix

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
states <- sf::st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- sf::st_as_sf(map("lakes", plot = FALSE, fill = TRUE))

# set up coordinate matrix for ggtess3Q
coord.mat <- q_mat_all %>% dplyr::select(long, lat) %>% as.matrix() 

### plot specific variables
buffer <- 0.5 # buffer for lat/long constraints of map
col.6 <- c("aquamarine", "skyblue2", "darkgreen", "purple3", "brown4", "deeppink") # need as many colors as you have K clusters
colors <- tess3r::CreatePalette(color.vector=col.6, palette.length = length(col.6)) # palette.length is the number of K clusters
res <- 1e2 # resolution for the tess3Q interpolation. higher gives smoother curves

### PLOTTING
ggtess3Q(q_redux, coord.mat, map.polygon = NULL, col.palette = colors, resolution = c(res,res))+
  geom_sf() +
  geom_sf(data = lakes, fill = "white") +
  geom_sf(data = states, fill = NA) +
  coord_sf(xlim = c(min(coord.mat[,"long"])-buffer, max(coord.mat[,"long"])+buffer), 
           ylim = c(min(coord.mat[,"lat"])-buffer, max(coord.mat[,"lat"])+buffer), expand = FALSE) +
  theme_bw()+
  xlab("Longitude")+
  ylab("Latitude")+
  ggnewscale::new_scale_fill()+
  scatterpie::geom_scatterpie(data=q_mat_all, aes(x=long, y=lat, group=lin), 
                              cols=paste0("K", 1:ncol(q_redux)), pie_scale = 1, color="black") + # pie_scale adjusts size but it's weird
  scale_fill_manual(values=col.6)+
  ggrepel::geom_label_repel(data=q_mat_all, aes(x=long, y=lat, label=pop, color=lin), size=4,
                            max.overlaps = 20, # adjust up if many populations overlap
                            # arrow = arrow(length = unit(0.02, "npc")), # adds arrow instead of line connecting label to pie chart
                            box.padding = 0.25)+
  scale_color_manual(values=c("green4", "skyblue3", "purple1"))+ # lineage colors for labels in repel
  theme(text = element_text(size = 18), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"))

