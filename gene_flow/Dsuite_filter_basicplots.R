
# gene flow using statistics calculated in Dsuite
library(ggplot2)

setwd("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/")

# outputs of Dsuite (no missing SNP filter)
bbaa <- read.csv("gene_flow/gene_flow_ALL_all_pops_ALLSNP_BBAA.txt", stringsAsFactors = T, sep="\t")
dmin <- read.csv("gene_flow/gene_flow_ALL_all_pops_ALLSNP_Dmin.txt", stringsAsFactors = T, sep="\t")

# population region information
pop_linreg <- read.csv("./VCF/noTperf_idpoplin.txt", sep="\t") %>% 
  tidyr::separate(lin, into=c("lin", "reg"), sep="_")


# plot(bbaa$p.value.adj, ylab="p value",xlab="trio number", ylim=c(0,0.05))
# plot(bbaa$f4.ratio, ylab="f4-ratio",xlab="trio number", ylim=c(0,1))

bbaa <- bbaa %>% 
  filter(P2 != "VA96" & P2 != "NC109B" & P2 != "NC109C" & P2 != "NC109D" & P2 != "NC90") %>% 
  filter(P3 != "VA96" & P3 != "NC109B" & P3 != "NC109C" & P3 != "NC109D" & P3 != "NC90")



# North Carolina subset
pop_nc <- pop_linreg %>% filter(reg=="NC")

bbaa_nc <- bbaa %>% 
  filter(P2 %in% pop_nc$pop & P3 %in% pop_nc$pop) %>% 
  mutate(lin_col = case_when((
    P2 %in% subset(pop_nc, lin == "Western")$pop & P3 %in% subset(pop_nc, lin == "Western")$pop) ~ "Western",
    (P2 %in% subset(pop_nc, lin == "Appalachian")$pop & P3 %in% subset(pop_nc, lin == "Appalachian")$pop) ~ "Appalachian",    (P2 %in% subset(pop_nc, lin == "Appalachian")$pop & P3 %in% subset(pop_nc, lin == "Appalachian")$pop) ~ "Appalachian",
    (P2 %in% subset(pop_nc, lin == "Eastern")$pop & P3 %in% subset(pop_nc, lin == "Eastern")$pop) ~ "Eastern",
    (P2 %in% subset(pop_nc, lin == "Appalachian")$pop & P3 %in% subset(pop_nc, lin == "Eastern")$pop) ~ "between",
    (P2 %in% subset(pop_nc, lin == "Eastern")$pop & P3 %in% subset(pop_nc, lin == "Appalachian")$pop) ~ "between",
    (P2 %in% subset(pop_nc, lin == "Appalachian")$pop & P3 %in% subset(pop_nc, lin == "Western")$pop) ~ "between",
    (P2 %in% subset(pop_nc, lin == "Western")$pop & P3 %in% subset(pop_nc, lin == "Eastern")$pop) ~ "between",
    (P2 %in% subset(pop_nc, lin == "Eastern")$pop & P3 %in% subset(pop_nc, lin == "Western")$pop) ~ "between",
    (P2 %in% subset(pop_nc, lin == "Appalachian")$pop & P3 %in% subset(pop_nc, lin == "Western")$pop) ~ "between",
    (P2 %in% subset(pop_nc, lin == "Western")$pop & P3 %in% subset(pop_nc, lin == "Appalachian")$pop) ~ "between"))

# Pennsylvania
pop_pa <- pop_linreg %>% filter(reg=="PA")

bbaa_pa <- bbaa %>% 
  filter(P2 %in% pop_pa$pop & P3 %in% pop_pa$pop) %>% 
  mutate(lin_col = case_when((
    P2 %in% subset(pop_pa, lin == "Western")$pop & P3 %in% subset(pop_pa, lin == "Western")$pop) ~ "Western",
    (P2 %in% subset(pop_pa, lin == "Appalachian")$pop & P3 %in% subset(pop_pa, lin == "Appalachian")$pop) ~ "Appalachian",    (P2 %in% subset(pop_pa, lin == "Appalachian")$pop & P3 %in% subset(pop_pa, lin == "Appalachian")$pop) ~ "Appalachian",
    (P2 %in% subset(pop_pa, lin == "Eastern")$pop & P3 %in% subset(pop_pa, lin == "Eastern")$pop) ~ "Eastern",
    (P2 %in% subset(pop_pa, lin == "Appalachian")$pop & P3 %in% subset(pop_pa, lin == "Eastern")$pop) ~ "between",
    (P2 %in% subset(pop_pa, lin == "Eastern")$pop & P3 %in% subset(pop_pa, lin == "Appalachian")$pop) ~ "between",
    (P2 %in% subset(pop_pa, lin == "Appalachian")$pop & P3 %in% subset(pop_pa, lin == "Western")$pop) ~ "between",
    (P2 %in% subset(pop_pa, lin == "Western")$pop & P3 %in% subset(pop_pa, lin == "Eastern")$pop) ~ "between",
    (P2 %in% subset(pop_pa, lin == "Eastern")$pop & P3 %in% subset(pop_pa, lin == "Western")$pop) ~ "between",
    (P2 %in% subset(pop_pa, lin == "Appalachian")$pop & P3 %in% subset(pop_pa, lin == "Western")$pop) ~ "between",
    (P2 %in% subset(pop_pa, lin == "Western")$pop & P3 %in% subset(pop_pa, lin == "Appalachian")$pop) ~ "between"))

# Virginia
pop_va <- pop_linreg %>% filter(reg=="VA")

bbaa_va <- bbaa %>% 
  filter(P2 %in% pop_va$pop & P3 %in% pop_va$pop) %>% 
  mutate(lin_col = case_when((
    P2 %in% subset(pop_va, lin == "Western")$pop & P3 %in% subset(pop_va, lin == "Western")$pop) ~ "Western",
    (P2 %in% subset(pop_va, lin == "Appalachian")$pop & P3 %in% subset(pop_va, lin == "Appalachian")$pop) ~ "Appalachian",    
    (P2 %in% subset(pop_va, lin == "Eastern")$pop & P3 %in% subset(pop_va, lin == "Eastern")$pop) ~ "Eastern",
    (P2 %in% subset(pop_va, lin == "Appalachian")$pop & P3 %in% subset(pop_va, lin == "Eastern")$pop) ~ "between",
    (P2 %in% subset(pop_va, lin == "Eastern")$pop & P3 %in% subset(pop_va, lin == "Appalachian")$pop) ~ "between",
    (P2 %in% subset(pop_va, lin == "Appalachian")$pop & P3 %in% subset(pop_va, lin == "Western")$pop) ~ "between",
    (P2 %in% subset(pop_va, lin == "Western")$pop & P3 %in% subset(pop_va, lin == "Eastern")$pop) ~ "between",
    (P2 %in% subset(pop_va, lin == "Eastern")$pop & P3 %in% subset(pop_va, lin == "Western")$pop) ~ "between",
    (P2 %in% subset(pop_va, lin == "Appalachian")$pop & P3 %in% subset(pop_va, lin == "Western")$pop) ~ "between",
    (P2 %in% subset(pop_va, lin == "Western")$pop & P3 %in% subset(pop_va, lin == "Appalachian")$pop) ~ "between"))


bbaa <- rbind(bbaa_nc, bbaa_pa, bbaa_va) %>% as.data.frame()

# restrict to comparisons of interest (W(WA)O - AAWO - E(EA)O - A(AE)O). O is implied. (XX) should be from same region. X() should be allopatric

pop_linreg <- read.csv("./VCF/noTperf_idpoplin.txt", sep="\t")
pop_linreg <- pop_linreg %>% 
  tidyr::separate(lin, into=c("lin", "reg"), sep="_") %>% 
  dplyr::select(-id) %>% 
  group_by(pop) %>% 
  slice(1) %>% 
  ungroup()
# merge pop_linreg into bbaa
names(pop_linreg)[1] <- "P1"
bbaa <- merge(bbaa, pop_linreg, by="P1") %>% group_by(P1,P2,P3) %>% slice(1) %>% ungroup()
names(pop_linreg)[1] <- "P2"
bbaa <- merge(bbaa, pop_linreg, by="P2") %>% group_by(P1,P2,P3) %>% slice(1) %>% ungroup()
names(pop_linreg)[1] <- "P3"
bbaa <- merge(bbaa, pop_linreg, by="P3") %>% group_by(P1,P2,P3) %>% slice(1) %>% ungroup()
names(bbaa)[12:ncol(bbaa)] <- c("P1.lin", "P1.reg", "P2.lin", "P2.reg", "P3.lin", "P3.reg")
# filter to interest
bbaa.temp <- bbaa %>% 
  filter(
    P1.reg == "Allopatry" & P2.reg != "Allopatry" & P3.reg != "Allopatry" & # make sure only the within lineage comparison is allopatric
      P2.reg == P3.reg & P2.lin != P3.lin # require population comparisons to be from same region but not same lineage
  )

bbaa.temp$p.value.adj <- p.adjust(bbaa.temp$p.value, method="BH")


# collate separate df and restrict by sig p-values
bbaa.sig.reg <- bbaa.temp %>% 
  as.data.frame() %>% 
  filter(p.value.adj <= 0.05)

write.csv(bbaa.sig.reg, "gene_flow/sig_bbaa_reg_CLDmethod.csv", row.names = F)

# # F4 ratio
# ggplot(data=(bbaa_nc %>% filter(p.value.adj <= 0.05)), aes(x=P2, y=P3))+
#   geom_tile(aes(alpha=f4.ratio, fill=lin_col))+
#   scale_fill_manual(values=c("Western" = "purple3", "Appalachian" = "green4", "between" = "blue3"))+
#   theme_bw()+
#   theme(text = element_text(size = 12), legend.position = "right",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.background = element_blank(), axis.line = element_line(colour = "black"),
#         aspect.ratio=1, axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# # D statistic
# ggplot(data=(bbaa_nc %>% filter(p.value.adj <= 0.05)), aes(x=P2, y=P3))+
#   geom_tile(aes(alpha=Dstatistic, fill=lin_col))+
#   scale_fill_manual(values=c("Western" = "purple3", "Appalachian" = "green4", "between" = "blue3"))+
#   theme_bw()+
#   theme(text = element_text(size = 12), legend.position = "right",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.background = element_blank(), axis.line = element_line(colour = "black"),
#         aspect.ratio=1, axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

