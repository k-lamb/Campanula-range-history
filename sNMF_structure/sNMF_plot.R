
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

best.K = 6 # enter value from sNMF.R
best = which.min(cross.entropy(pc, K = best.K)) # best run within K

q_mat <- LEA::Q(pc, K = best.K, run = best) %>% as.data.frame()
colnames(q_mat) <- paste0("K", 1:ncol(q_mat))
q_mat_all <- cbind(pop_linreg, q_mat)

q_mat_all <- q_mat_all %>% 
  pivot_longer(cols = starts_with("K"), names_to = "K_cluster", values_to = "q") 
q_mat_all$linreg <- paste0(q_mat_all$lin, "_", q_mat_all$reg)

t <- q_mat_all %>% 
  dplyr::group_by(pop, lin, reg, K_cluster) %>% 
  reframe(q_sum = sum(q))

t1 <- t %>% 
  dplyr::group_by(pop, lin, reg) %>% 
  reframe(q_sum_all = sum(q_sum))

t2 <- merge(t,t1, by=c("pop", "lin", "reg"))
t2$q_frac <- t2$q_sum/t2$q_sum_all

# drop populations not interested in
t3 <- t2 %>% 
  filter(pop != "VA96" & pop != "NC90" & pop != "NC109B" & pop != "NC109C" & pop != "NC109D" & pop != "VA71")



### PLOTTING

# plot of individuals separated out
ggplot2::ggplot(data=q_mat_all)+
  geom_col(aes(x=id, y=q, fill=K_cluster))+
  theme(text = element_text(size = 12), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio=0.33)

# ancestry collated by population
tt <- t3
tt <- tt %>% mutate_at(vars(pop, lin, reg, K_cluster), as.factor)
# # by lineage and region: allo, nc, pa, va
tt$pop <- factor(tt$pop, levels=c("GA22", "KY51", "OH64", "VA86", "PA27", "PA94", "PA103", 
                                  "VA73", "MD5", "PA95", "PA101", "PA102", "PA104", "VA111", "VA131L",
                                  "VA93", "VA112", "VA85", 
                                  "NC105", "NC107", "NC108", "NC109A", "NC91", "NC109E", "NC110", "TN92"
                                  ))
# # color scheme for everything 80%
# # for K=6
col17 <- c("aquamarine", "skyblue2", "darkgreen", "purple3", "brown4",
           "deeppink", "antiquewhite3", "skyblue2", "azure4", "darkkhaki", "green3",
           "deeppink", "chocolate3", "coral1", "cyan4", "aquamarine", "darkorchid")
colors <- CreatePalette(color.vector=col17[1:best.K], palette.length = best.K)


# jpeg("sNMF_structure/plots/80perc/barplot_noVA71.jpeg", res=300, height=5, width=10, units="in")
ggplot2::ggplot(data=tt, aes(x=pop, y=q_frac, fill=K_cluster))+
  geom_col() +
  scale_fill_manual(values=col17)+
  ylab("Ancestry Proportion")+
  xlab("")+
  scale_x_discrete(position = "top")+
  # facet_grid(reg~lin)+
  # scale_y_continuous(labels = scales::percent)+
  theme(text = element_text(size = 15), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio=0.33, axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 90))
# dev.off()



# mapping ancestry

world <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))

locations <- read.csv("sNMF_structure/group_coords.csv", stringsAsFactors = T)
names(locations)[1] <- "pop"
# double check that q_frac adds appropriately
t3 %>% group_by(pop) %>% reframe(sum = sum(q_frac))

t4 <- merge(t3, locations[,c(1,3,4)]) # add in lat/long pop data

# pivot wide for scatterpie to work properly
t5 <- t4 %>% 
  pivot_wider(id_cols = c(pop, lin, reg, lat, long), names_from = K_cluster, values_from = q_frac)


# LOCATION DATA MERGE
# adding cordinate information
temp <- merge(locations %>% dplyr::select(pop, long, lat),
              pop_linreg %>% dplyr::select(id, pop), 
              by="pop")

# setting up coordinate matrix without NA's (filtering out pop's with NA's)
coord.mat <- temp %>%
  dplyr::filter(pop != "NC109B" & pop != "NC109C" & pop != "NC109D"  & pop != "NC90" & pop != "VA96" & pop != "VA71") %>% 
  dplyr::group_by(pop) %>% 
  dplyr::slice(1) %>%
  dplyr::ungroup() %>% 
  dplyr::select(long, lat) %>%
  as.matrix()
colnames(coord.mat) <- c("V1", "V2")
# usa <- rworldmap::getMap("United States", resolution="low")

# text label data
temp2 <- merge(locations %>% dplyr::select(pop, long, lat),
               pop_linreg %>% dplyr::select(id, pop, lin), 
               by="pop")
text.dat <- temp2 %>%
  dplyr::filter(pop != "NC109B" & pop != "NC109C" & pop != "NC109D"  & pop != "NC90" & pop != "VA96" & pop != "VA71") %>% 
  dplyr::group_by(pop) %>% 
  dplyr::slice(1) %>%
  as.data.frame()

# adding gene flow data from Dsuite
temp.unique <- temp2 %>% 
  group_by(pop) %>% 
  slice(1) %>% 
  dplyr::select(pop, long, lat)
names(temp.unique)[1] <- "P2"

gg <- merge(gene.bbaa, temp.unique, by="P2")
names(temp.unique)[1] <- "P3"
gg <- merge(gg, temp.unique, by="P3")

# gg nc
ggnc <- gg %>% mutate(lat.x = if_else(P2 == "NC109A", lat.x-0.06, lat.x)) %>%  # trying to bump the NC109A pop down a bit so you can see the gene flow line
  mutate(long.x = if_else(P2 == "NC109A", long.x+0.015, long.x)) %>% 
  mutate(lat.y = if_else(P3 == "NC109A", lat.y-0.06, lat.y)) %>% 
  mutate(long.y = if_else(P3 == "NC109A", long.y+0.015, long.y)) 

# gg va
ggva <- gg %>% mutate(lat.x = if_else(P2 == "VA131L", lat.x-0.06, lat.x)) %>%  # trying to bump the NC109A pop down a bit so you can see the gene flow line
  mutate(long.x = if_else(P2 == "VA131L", long.x+0.015, long.x)) %>% 
  mutate(lat.y = if_else(P3 == "VA131L", lat.y-0.06, lat.y)) %>% 
  mutate(long.y = if_else(P3 == "VA131L", long.y+0.015, long.y)) 

# plot everyone together

ggplot(data = world) +
  geom_sf() +
  geom_sf(data=lakes, fill = "gray50") +
  geom_sf(data = states, fill = NA) +
  coord_sf(xlim = c(-85, -76), ylim = c(34, 42), expand = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Longitude")+
  ylab("Latitude")+
  # ggnewscale::new_scale_fill() +
  geom_segment(data=subset(gg, lin_col=="between"), aes(x=long.x, y=lat.x, xend=long.y, yend=lat.y,
                                                        linewidth=no.sig, alpha=no.sig/max(gg$no.sig)))+
  scatterpie::geom_scatterpie(data=t5, aes(x=long, y=lat, group=lin, color=lin), 
                              cols=paste0("K", 1:best.K), pie_scale = 1)+
  scale_fill_manual(values=col17)+
  scale_color_manual(values=c("green4", "blue1", "purple1"))



# regions separated
size.t = 12
size.label = 4
alpha.fac = 10
width.fac = 2

gg$mod.alpha = gg$mean.Dstat*alpha.fac
ggnc$mod.alpha = ggnc$mean.Dstat*alpha.fac
ggva$mod.alpha = ggva$mean.Dstat*alpha.fac

gg$mod.sig = gg$no.sig*width.fac
ggnc$mod.sig = ggnc$no.sig*width.fac
ggva$mod.sig = ggva$no.sig*width.fac

# Pennsylvania
map_pa <-
  ggplot(data = world) +
  geom_sf() +
  geom_sf(data=lakes, fill = "gray50") +
  geom_sf(data = states, fill = "white", color = "white") +
  coord_sf(xlim = c(-80.7,-79.3), ylim = c(40.1, 41.7), expand = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(" ") + # Longitude
  ylab(" ") + # Latitude
  # ggnewscale::new_scale_fill() +
  geom_segment(data=subset(gg, lin_col=="between"), aes(x=long.x, y=lat.x, xend=long.y, yend=lat.y), 
               linewidth=gg$mod.sig, alpha=gg$mod.alpha)+
  scatterpie::geom_scatterpie(data=t5, aes(x=long, y=lat, group=lin), 
                              cols=paste0("K", 1:best.K), pie_scale = 0.5, color="black") +
  scale_fill_manual(values=col17)+
  ggrepel::geom_label_repel(data=subset(t5, reg=="PA"), aes(x=long, y=lat, label=pop, color=lin), size=size.label,
                            # arrow = arrow(length = unit(0.02, "npc")),
                            box.padding = 1.5)+
  scale_color_manual(values=c("green4", "purple1")) +
  # scale_color_manual(values=c("green4", "blue1", "purple1"))+
  theme(text = element_text(size = size.t), legend.position = "none", # right
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio=1, axis.text = element_blank()) # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# map_pa

# North Carolina
t6 <- t5 %>% mutate(lat = if_else(pop == "NC109A", lat-0.06, lat)) %>%  # trying to bump the NC109A pop down a bit so you can see the gene flow line
  mutate(long = if_else(pop == "NC109A", long+0.015, long)) 

map_nc <-
  ggplot(data = world) +
  geom_sf() +
  geom_sf(data=lakes, fill = "gray50") +
  geom_sf(data = states, fill = "white", color = "white") +
  coord_sf(xlim = c(-83.4, -82.7), ylim = c(35.4, 36.15), expand = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(" ") + # Longitude
  ylab(" ") + # Latitude
  # ggnewscale::new_scale_fill() +
  geom_segment(data=subset(ggnc, lin_col=="between"), aes(x=long.x, y=lat.x, xend=long.y, yend=lat.y), 
               linewidth=ggnc$mod.sig, alpha=ggnc$mod.alpha)+
  scatterpie::geom_scatterpie(data=t6, aes(x=long, y=lat, group=lin), 
                              cols=paste0("K", 1:best.K), pie_scale = 0.25, color="black") +
  scale_fill_manual(values=col17)+
  ggrepel::geom_label_repel(data=subset(t6, reg=="NC"), aes(x=long, y=lat, label=pop, color=lin), size=size.label,
                            # arrow = arrow(length = unit(0.02, "npc")),
                            max.overlaps = 20,
                            box.padding = 0.5)+
  scale_color_manual(values=c("green4", "purple1")) +
  # scale_color_manual(values=c("green4", "blue1", "purple1"))
  theme(text = element_text(size = size.t), legend.position = "none", # right
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio=1, axis.text = element_blank()) # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# map_nc

# Virginia
t7 <- t5 %>% mutate(lat = if_else(pop == "VA131L", lat-0.06, lat)) %>%  # trying to bump the NC109A pop down a bit so you can see the gene flow line
  mutate(long = if_else(pop == "VA131L", long+0.015, long)) %>% 
  filter(pop != "VA73")

map_va <- 
  ggplot(data = world) +
  geom_sf() +
  geom_sf(data=lakes, fill = "gray50") +
  geom_sf(data = states, fill = "white", color = "white") +
  coord_sf(xlim = c(-80.5,-78.75), ylim = c(36.675, 38.445), expand = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(" ") + # Longitude
  ylab(" ") + # Latitude
  # ggnewscale::new_scale_fill() +
  geom_segment(data=subset(ggva, lin_col=="between"), aes(x=long.x, y=lat.x, xend=long.y, yend=lat.y), 
               linewidth=ggva$mod.sig, alpha=ggva$mod.alpha)+
  scatterpie::geom_scatterpie(data=t7, aes(x=long, y=lat, group=lin), 
                              cols=paste0("K", 1:best.K), pie_scale = 1, color="black") +
  scale_fill_manual(values=col17)+
  ggrepel::geom_label_repel(data=subset(t7, reg=="VA" | pop == "VA73"), aes(x=long, y=lat, label=pop, color=lin), size=size.label,
                            # arrow = arrow(length = unit(0.02, "npc")),
                            box.padding = 1.5)+
  # scale_color_manual(values=c("green4", "skyblue3")) +
  # ggnewscale::new_scale_color()+
  scale_color_manual(values=c("green4", "skyblue4", "purple1"))+
  theme(text = element_text(size = size.t), legend.position = "none", # right
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio=1, axis.text = element_blank()) # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# map_va

# jpeg("sNMF_structure/plots/80perc/full_snmf_regions_GENEFLOW_k6_noVA71_white.jpeg", res=500, width=4, height=10, unit="in")
ggpubr::ggarrange(map_pa, map_va, map_nc, ncol=1, nrow=3, common.legend = F, legend = "none")
# dev.off()



### TESS3R mapping

# set up LEA q matrix as tess3 q matrix and getting rid of unwanted pops
q_prep <- t2 %>% 
  pivot_wider(id_cols="pop", names_from = "K_cluster", values_from = "q_frac") %>%
  dplyr::filter(pop != "NC109B" & pop != "NC109C" & pop != "NC109D" & pop != "VA71" & pop != "VA96" & pop != "NC90") %>% 
  dplyr::select(-pop) %>% 
  as.matrix() %>% 
  as.qmatrix()

# in base R plot
# jpeg("sNMF_structure/plots/full_snmf_qmat_map_k8.jpeg", res=200, width=10, height=10, unit="in")
# plot(x=q_prep, coord=coord.mat, method = "map.max",
#      resolution = c(400,400), col.palette=colors,
#      # raster.filename = usa,
#      interpolation.model = FieldsKrigModel(10), cex = 0,
#      xlab = "Longitude", ylab= "Latitude", main = "Ancestry coefficients")
# mark(x = text.dat$long, y = text.dat$lat, labels = text.dat$pop, cex = 0.8, col_bg = "gray90",
#      col = ifelse(text.dat$lin == "Eastern", "skyblue3", ifelse(text.dat$lin == "Western", "purple3", "green4")))
# dev.off()

# in ggplot
map.polygon <- rworldmap::getMap(resolution = "high")
# buff <- 0.75

# jpeg("sNMF_structure/plots/80perc/full_snmf_qmat_map_k6_noVA71_white.jpeg", res=500, width=10, height=10, unit="in")
ggtess3Q(q_prep, coord.mat, map.polygon = map.polygon, col.palette = colors, resolution = c(1e3,1e3))+
  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  geom_sf() +
  geom_sf(data=lakes, fill = "white") +
  geom_sf(data = states, fill = NA) +
  coord_sf(xlim = c(-85, -76.2), ylim = c(34, 42), expand = FALSE) +
  # xlim(c(min(coord.mat[,1]-buff), max(coord.mat[,1]+buff)))+
  # ylim(c(min(coord.mat[,2]-buff), max(coord.mat[,2]+buff)))+
  theme_bw()+
  xlab("Longitude")+
  ylab("Latitude")+
  ggnewscale::new_scale_fill()+
  scatterpie::geom_scatterpie(data=t5, aes(x=long, y=lat, group=lin), 
                              cols=paste0("K", 1:best.K), pie_scale = 1, color="black") +
  scale_fill_manual(values=col17)+
  # ggrepel::geom_label_repel(data=t5, aes(x=long, y=lat, label=pop, color=lin), size=4,
  #                           max.overlaps = 20,
  #                           # arrow = arrow(length = unit(0.02, "npc")),
  #                           box.padding = 0.25)+
  # scale_color_manual(values=c("green4", "skyblue3", "purple1"))+
  theme(text = element_text(size = 18), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"))
# dev.off()
