
library(tidyr)
library(readxl)
library(dplyr)
# make the metadata file for the Moments runs
# requires cols for: iterations=50, Pair, pop1, pop2, proj1, proj2

### POP MAP FOR COLLATED REGION/LINEAGE
# popmap_old <- read.csv("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/VCF/vcf_popmap.txt", header=F, stringsAsFactors = T, sep="\t")
# names(popmap_old) <- c("id", "pop")
# 
# pop_linreg <- read.csv("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/VCF/noTperf_idpoplin.txt", sep="\t") %>% tidyr::separate(lin, into=c("lin", "reg"), sep="_") # population information: id, pop, lineage, region
# 
# popmap <- merge(popmap_old, pop_linreg, by=c("id", "pop"))
# popmap$linreg <- paste0(popmap$lin, "_", popmap$reg)
# write.table(popmap[,c("id", "linreg")], "~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/VCF/vcf_popmap_linreg.txt", col.names = F, sep="\t", row.names = F, quote = F)

### easySFS data
setwd("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/")
x <- readxl::read_excel("./demo_infer/easySFS/easySFS_max_retain.xlsx", sheet = 5) %>% as.data.frame() # easySFS. sheet==1 is 65% SNP filter; sheet=2 is 80% SNP filter, sheet=4 is 65% filter & region/lineage collated; sheet=5 is 80% filter & region/lineage collated
pop_linreg <- read.csv("./VCF/noTperf_idpoplin.txt", sep="\t") %>% tidyr::separate(lin, into=c("lin", "reg"), sep="_") # population information: id, pop, lineage, region

pop_linreg <- pop_linreg %>% 
  group_by(pop, reg) %>% 
  slice(1) %>% 
  as.data.frame()

y <- expand.grid(pop1 = x$POP, pop2 = x$POP)
head(y)

# merge in projection data 
x_sub <- x[,c("POP", "BEST")]
names(x_sub)[1] <- "pop1"
y <- merge(x_sub, y, by="pop1")
names(y)[2] <- "proj1"
names(x_sub)[1] <- "pop2"
y <- merge(x_sub, y, by="pop2")
names(y)[2] <- "proj2"

# reorder
y <- y[,c(3,1,4,2)]
head(y)

# combine in region information to slim down total list of pairs
# # for metadata not using collated: 
# pops <- pop_linreg[,c("pop", "reg")] 

# for metadata using collated:
pop_linreg$linreg <- paste0(pop_linreg$lin, "_", pop_linreg$reg)
pops <- pop_linreg[,c("linreg", "reg")]

names(pops)[1] <- "pop1"
z <- merge(y, pops, by="pop1")

names(pops)[1] <- "pop2"
z <- merge(z, pops, by="pop2")

# keep only pairs with NC_NC, PA_PA, VA_VA
z$reg.xy <- paste0(z$reg.x, "_", z$reg.y)

temp <- z %>% 
  filter(reg.xy == "NC_NC" | reg.xy == "PA_PA" | reg.xy == "VA_VA") %>% 
  dplyr::select(pop1, pop2, proj1, proj2)

temp$Pair <- paste0(temp$pop1, ".", temp$pop2)
temp$iterations = 10

# reorder
metadata <- temp[,c("iterations", "Pair", "pop1", "pop2", "proj1", "proj2")]
head(metadata)

metadata <- metadata %>% 
  filter(pop1 != pop2)

# filter out reciprocals (e.g., NC105.NC108, NC108.NC105)
temp <- metadata

# function splits the pair name and sorts them alphanumerically before pasting them back together
extract_reciprocal <- function(x) {
  sorted_pair <- sort(unlist(strsplit(x, "\\.")))
  paste(sorted_pair, collapse = ".")
}

# remove pairs that match regardless of their order
metadata <- temp %>%
  mutate(reciprocal = sapply(Pair, extract_reciprocal)) %>%
  distinct(reciprocal, .keep_all = TRUE) %>%
  dplyr::select(-reciprocal)

metadata %>% head()

write.table(metadata, "demo_infer/All_metadat_80_linreg.txt", sep="\t", row.names=F, quote=F) 

