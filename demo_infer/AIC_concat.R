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

# perc: percentage SNP filter (80/65)
# meta: whether metapopulation data (T) is being aggregated or at the level of the population (F)
# cyclic: whether cyclic models are being aggregated (T) or not (F)
perc <- 80
meta <- F
cyclic <- F

if (meta == T) {
  # directory for determining metapopulation models
  dir <- c("migration", "secondary_migration", "no_migration", "ancient_migration") # list of final directories
}

if (meta == F) {
  # directory for determining migration and time of split from pop pairwise best model:
  dir <- c("migration", "secondary_migration", "no_migration", "ancient_migration") # list of final directories
}

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
# change directory for percentage missing SNP
# 2c folders require dummy no_migration folder
for (i in dir) {
  if (meta == T & cyclic == T) {
    setwd(sprintf("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/%s/data/80perc_linreg_2c", i)) # output directories. requires dummy folder in no_migration
  }
  if (meta == T & cyclic == F) {
    setwd(sprintf("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/%s/data/80perc_linreg", i)) 
  }
  if (meta == F & cyclic == T) {
    setwd(sprintf("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/%s/data/80perc_2c", i)) # requires dummy folder in no_migration
  }
  if (meta == F & cyclic == F) {
    setwd(sprintf("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/%s/data/80perc", i))
  }
  
  files <- list.files(pattern = "\\.txt") # list all files in output directories
  tables <- lapply(files, read.csv, header = TRUE, sep="\t") # read them in as .txt's
  res <- do.call(rbind, tables) # rbind into single data frame
  dat <- rbind(dat, res) # rbinds different super models together
}

# fixing issue in migration runs where model name did not save properly (fixed now, unnecessary for runs newer than 3/31/24)
dat <- dat %>% 
  separate(model, into=c("functions", "model", "at", "extraneous"), sep = " ") %>% 
  mutate(model = if_else(is.na(model) == T, functions, model)) %>% 
  dplyr::select(-functions, -at, -extraneous)

# different strategy of finding conserved model
t <- dat %>%
  group_by(Pair_name, model) %>% 
  arrange(aic) %>%
  slice(1)

# generate a 'super model' column that aggregates models by demographic theme (e.g., IM, SC, etc...)
super_model <- separate(data = t, col = model, into = c("super_model", "right")) # split model col to get super order. ignore warning.
super_model <- super_model[,"super_model"] 
t <- cbind(t, super_model)


# assign W or A or E lineage
if (meta == T) {
  pairs <- separate(data=t, col=Pair_name, into=c("pop1.lin", "st1", "pop2.lin", "st2"))
  pairs <- pairs[,c(1:4)]
}

if (meta == F) {
  pairs <- separate(data=t, col=Pair_name, into=c("pop1", "pop2"))
  pairs <- pairs[,c(1,2)]
}

t <- cbind(t, pairs)

if (meta == T) {
  # W&W==W, A&A==A, E&E==E, else B
  t$lin.pair <- NA
  t <- t %>%
    mutate(lin.pair = case_when(
      # # or lin and reg separated to pop comparisons: 
      pop1.lin == "A" & pop2.lin == "A" ~ "A",
      pop1.lin == "W" & pop2.lin == "W" ~ "W",
      pop1.lin == "E" & pop2.lin == "E" ~ "E",
      # # for lin and reg combined:
      pop1.lin == "Appalachian" & pop2.lin == "Appalachian" ~ "A",
      pop1.lin == "Western" & pop2.lin == "Western" ~ "W",
      pop1.lin == "Eastern" & pop2.lin == "Eastern" ~ "E",
      TRUE ~ "B" # else case, where if pop1.lin != pop2.lin then lin.pair == "B"
    )
  )
}

if (meta == F) {
  t$pop1.lin <- NA #doing pop 1 first in the duo
  t <- t %>%
    mutate(pop1.lin = case_when(
      pop1 == "GA22" | pop1 == "KY51" | pop1 == "NC105" | pop1 == "NC107" | pop1 == "NC108" | pop1 == "NC109A" |
        pop1 == "OH64" | pop1 == "PA103" | pop1 == "PA27" | pop1 == "PA94" | pop1 == "NC105" | pop1 == "VA86" ~ "W",
      pop1 == "MD5" | pop1 == "NC109E" | pop1 == "NC110" | pop1 == "NC91" | pop1 == "PA101" | pop1 == "PA102" |
        pop1 == "PA104" | pop1 == "PA95" | pop1 == "TN92" | pop1 == "VA111" | pop1 == "VA131L" | pop1 == "VA73" ~ "A",
      pop1 == "VA85" | pop1 == "VA93" | pop1 == "VA112" ~ "E"
    )
    )
  t$pop2.lin <- NA #doing pop 1 first in the duo
  t <- t %>%
    mutate(pop2.lin = case_when(
      pop2 == "KY51" | pop2 == "NC105" | pop2 == "NC107" | pop2 == "NC108" | pop2 == "NC109A" | pop2 == "OH64" |
        pop2 == "PA103" | pop2 == "PA27" | pop2 == "PA94" | pop2 == "NC105" | pop2 == "VA86" ~ "W",
      pop2 == "MD5" | pop2 == "NC109E" | pop2 == "NC110" | pop2 == "NC91" | pop2 == "PA101" | pop2 == "PA102" |
        pop2 == "PA104" | pop2 == "PA95" | pop2 == "TN92" | pop2 == "VA111" | pop2 == "VA131L" | pop2 == "VA73" ~ "A",
      pop2 == "VA85" | pop2 == "VA93" | pop2 == "VA112" ~ "E"
    )
  )
}

# include region information
if (meta == T) {
  # for linreg only:
  t$region <- NA
  t <- t %>%
    mutate(region = if_else(st1 == st2, st1, NA))
  
  # write data to file
  if (cyclic == T) {
    # metapopulation (all super models) data:
    write.csv(t, "~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/data_concat/AIC_concat_80perc_linreg_2c.csv", row.names = F)
  }
  
  if (cyclic == F) {
    # metapopulation (all super models) data:
    write.csv(t, "~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/data_concat/AIC_concat_80perc_linreg.csv", row.names = F)
    
  }
  
}

if (meta == F) {
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
      (st1 == "VA" & st2 == "VA" & pop1 != "VA71" & pop2 != "VA71" & pop1 != "VA73" & pop2 != "VA73") ~ "VA"

    )
  )
  
  # write data to file
  if (cyclic == T) {
    write.csv(t, "~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/data_concat/AIC_concat_80perc_2c.csv", row.names = F)
  }
  if (cyclic == F) {
    write.csv(t, "~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/data_concat/AIC_concat_80perc.csv", row.names = F)
  }
  
}
