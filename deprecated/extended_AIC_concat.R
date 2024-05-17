### model concat for extended demography models ###

library(dplyr) # idk. presumably I'll use it though
library(tidyr) # for separate() model column into super models

# to do list:
# - read in all data
# - retain only best iteration per (pair & model)
# - create new super model column
# - convert param of interest to human-readable units
# - add column for region of interest (i.e. contact zone)
# - add column for bw vs wi

dir <- c("MIG", "SC", "NM") # list of final directories
dat <- data.frame(Pair_name = character(),
                  reversed = character(),
                  model = character(),
                  nu1 = numeric(),
                  nu2 = numeric(),
                  nu_ae = numeric(),
                  s = numeric(),
                  Ts = numeric(),
                  Tae = numeric(),
                  Tsc = numeric(),
                  m12 = numeric(),
                  m21 = numeric(),
                  theta = numeric(),
                  ll_model = numeric(),
                  aic = numeric())

# collects all diles in output directories, reads them in, then rbinds them into single file
for (i in dir) {
  setwd(sprintf("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/Moments/moments/extended/%s/", i)) # output directories
  files <- list.files(pattern = "\\.txt") # list all files in output directories
  tables <- lapply(files, read.csv, header = TRUE, sep="\t") # read them in as .txt's
  res <- do.call(rbind, tables) # rbind into single data frame
  dat <- rbind(dat, res) # rbinds different super models together
}

# take best model iteration per model and population pair, create new column for "super model" (i.e. SC, IM...)
t <- aggregate(aic~Pair_name+model, data=dat, FUN="min") # best model by AIC for each model and pair
t <- merge(t, dat, drop=F)
super_model <- separate(data = t, col = model, into = c("super_model", "right")) # split model col to get super order. ignore warning.
super_model <- super_model[,2] 
t <- cbind(t, super_model)

# bw vs wi addition

# first split pair into 2
# assign W or A or E lineage
pairs <- separate(data=t, col=Pair_name, into=c("pop1", "pop2"))
pairs <- pairs[,c(1,2)]
t <- cbind(t, pairs)

t$pop1.lin <- NA #doing pop 1 first in the duo
t <- t %>% 
  mutate(pop1.lin = case_when(
    pop1 == "GA22" | pop1 == "KY51" | pop1 == "NC105" | pop1 == "NC107" | pop1 == "NC108" | pop1 == "NC109A" |
    pop1 == "OH64" | pop1 == "PA103" | pop1 == "PA27" | pop1 == "PA94" | pop1 == "NC105" | pop1 == "VA86" ~ "W",
    pop1 == "MD5" | pop1 == "NC109E" | pop1 == "NC110" | pop1 == "NC91" | pop1 == "PA101" | pop1 == "PA102" |
    pop1 == "PA104" | pop1 == "PA95" | pop1 == "TN92" | pop1 == "VA111" | pop1 == "VA131L" | pop1 == "VA73" ~ "A",
    pop1 == "VA85" ~ "E"
  )
)

t$pop2.lin <- NA #doing pop 1 first in the duo
t <- t %>% 
  mutate(pop2.lin = case_when(
    pop2 == "KY51" | pop2 == "NC105" | pop2 == "NC107" | pop2 == "NC108" | pop2 == "NC109A" | pop2 == "OH64" | 
    pop2 == "PA103" | pop2 == "PA27" | pop2 == "PA94" | pop2 == "NC105" | pop2 == "VA86" ~ "W",
    pop2 == "MD5" | pop2 == "NC109E" | pop2 == "NC110" | pop2 == "NC91" | pop2 == "PA101" | pop2 == "PA102" |
    pop2 == "PA104" | pop2 == "PA95" | pop2 == "TN92" | pop2 == "VA111" | pop2 == "VA131L" | pop2 == "VA73" ~ "A",
    pop2 == "VA85" | pop2 == "VA93" ~ "E"
  )
)

# W&W==W, A&A==A, E&E==E, else B
t$lin.pair <- NA
t <- t %>%
  mutate(lin.pair = case_when(
    pop1.lin == "A" & pop2.lin == "A" ~ "A",
    pop1.lin == "W" & pop2.lin == "W" ~ "W",
    pop1.lin == "E" & pop2.lin == "E" ~ "E",
    TRUE ~ "B" # else case, where if pop1.lin != pop2.lin then lin.pair == "B"
  )
)

# region addition for CZ

# separate pop into state and id tags for easier sorting
temp <- separate(t, pop1, into=c("st1", "id1"), sep="(?<=[A-Za-z])(?=[0-9])")
temp <- separate(temp, pop2, into=c("st2", "id2"), sep="(?<=[A-Za-z])(?=[0-9])")
temp <- temp[,c("st1", "id1", "st2", "id2")]
t <- cbind(t, temp)

# identifying by state 
t$region <- NA
t <- t %>% 
  mutate(region = case_when(
    (st1 == "NC" & st2 == "NC") | (st1 == "TN" & st2 == "NC") | (st1 == "NC" & st2 == "TN") ~ "NC",
    st1 == "PA" & st2 == "PA" & pop1 != "PA95" & pop2 != "PA95" ~ "PA",
    (pop1 == "VA111" & pop2 == "VA131L") | (pop1 == "VA111" & pop2 == "VA85") | 
    (pop1 == "VA131L" & pop2 == "VA85") ~ "VA"
  )
)

write.csv(t, "~/Desktop/Documents/Research/Paper Code/Debban-Lamb/Moments/moments/extended/data/AIC_concat_extended.csv")  

