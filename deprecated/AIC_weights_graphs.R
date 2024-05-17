# build boxplot of model outcomes

library(ggplot2)
library(ggpubr)
library(forcats)
library(ggpattern)
library(dplyr)

setwd("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/Moments/moments")

t <- read.csv("./data/AIC_weights_all.csv")

t.s <- subset(t, aic.ll < 1.25) # retain only models with 80% chance or greater of being best model
t.s <- subset(t.s, region == "NC" | region == "PA" | region == "VA")

t.s.b <- subset(t.s, type =="B")
t.s.a <- subset(t.s, type == "A")
t.s.w <- subset(t.s, type == "W")

supp.labs <- c("NC"="North Carolina", 
               "PA"="Pennsylvania", 
               "VA"="Virginia")
size.t <- 5

# by CZ best by AIC raw
p1 <-
  ggplot(t.s.b, aes(x=type, fill=fct_rev(model)))+
    geom_bar(position="fill")+
    theme_bw()+
    xlab("between lineage comparisons")+
    ylab("proportion")+
    scale_fill_manual("legend", values = c(no_mig = "#00B81F", sc = "#00A5FF", mig = "firebrick2", am = "orchid2"))+
    facet_wrap(~ region, labeller = labeller(region = supp.labs))+
    theme(strip.text.x = element_text(size=size.t), axis.text=element_text(size=size.t))

# by CZ best by AIC with contest
p2 <- 
  ggplot(t.s.a, aes(x=type, fill=fct_rev(model)))+
    geom_bar(position="fill")+
    theme_bw()+
    xlab("appalachian lineage comparisons")+
    ylab("proportion")+
    scale_fill_manual("legend", values = c(no_mig = "#00B81F", sc = "#00A5FF", mig = "firebrick2", am = "orchid2"))+
    facet_wrap(~ region, labeller = labeller(region = supp.labs))+
    theme(strip.text.x = element_text(size=size.t), axis.text=element_text(size=size.t))

# by CZ best by AIC raw
p3 <- 
  ggplot(t.s.w, aes(x=type, fill=fct_rev(model)))+
    geom_bar(position="fill")+
    theme_bw()+
    xlab("western lineage comparisons")+
    ylab("proportion")+
    scale_fill_manual("legend", values = c(no_mig = "#00B81F", sc = "#00A5FF", mig = "firebrick2", am = "orchid2"))+
    facet_wrap(~ region, labeller = labeller(region = supp.labs))+
    theme(strip.text.x = element_text(size=size.t), axis.text=element_text(size=size.t))
  
jpeg("./graphs/AIC_weights/temp_cz_best_waic_1dot2x.jpeg", width=10, height=4, units="in", res=250)
ggarrange(p1,p2,p3, nrow=1, ncol=3, common.legend=T)
dev.off()



# using model quality instead

# t.s <- subset(t, model_score > 0.8) # retain only models with 80% chance or greater of being best model
t.s <- subset(t.s, region == "VA" | region == "NC" | region == "PA")

t.s.b <- subset(t.s, type =="B")
t.s.a <- subset(t.s, type == "A")
t.s.w <- subset(t.s, type == "W")

# by CZ best by AIC raw
p1 <-
  ggplot(t.s.b, aes(x=type, fill=fct_rev(model)))+
    geom_bar(position="fill")+
    theme_bw()+
    xlab("between lineage comparisons")+
    ylab("proportion")+
    scale_fill_manual("legend", values = c(no_mig = "#00B81F", sc = "#00A5FF", mig = "firebrick2", am = "orchid2"))+
    facet_wrap(~ region, labeller = labeller(region = supp.labs))+
    theme(strip.text.x = element_text(size=size.t), axis.text=element_text(size=size.t))

# by CZ best by AIC with contest
p2 <- 
  ggplot(t.s.a, aes(x=type, fill=fct_rev(model)))+
    geom_bar(position="fill")+
    theme_bw()+
    xlab("appalachian lineage comparisons")+
    ylab("proportion")+
    scale_fill_manual("legend", values = c(no_mig = "#00B81F", sc = "#00A5FF", mig = "firebrick2", am = "orchid2"))+
    facet_wrap(~ region, labeller = labeller(region = supp.labs))+
    theme(strip.text.x = element_text(size=size.t), axis.text=element_text(size=size.t))

# by CZ best by AIC raw
p3 <- 
  ggplot(t.s.w, aes(x=type, fill=fct_rev(model)))+
    geom_bar(position="fill")+
    theme_bw()+
    xlab("western lineage comparisons")+
    ylab("proportion")+
    scale_fill_manual("legend", values = c(no_mig = "#00B81F", sc = "#00A5FF", mig = "firebrick2", am = "orchid2"))+
    facet_wrap(~ region, labeller = labeller(region = supp.labs))+
    theme(strip.text.x = element_text(size=size.t), axis.text=element_text(size=size.t))

jpeg("./graphs/AIC_weights/temp_cz_best_waic_modelscore_0dot80.jpeg", width=10, height=4, units="in", res=250)
ggarrange(p1,p2,p3, nrow=1, ncol=3, common.legend = T)
dev.off()





# alternate presentation that does a better job showing nuance of score


# calculate akaike weights
test <- read.csv("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/Moments/moments/data/all_data.csv")
aic.min <- aggregate(aic~pair, data=test, FUN="min") # find min aic for pop pair
names(aic.min)[2] <- "aic.min" # min score for delta(aic) calculation
t <- merge(test, aic.min, all=T) 
t$aic.delta <- t$aic - t$aic.min

t <- subset(t, aic.delta < 10) # remove highly improbable models

t$aic.d.exp <- exp((-t$aic.delta)/2) # numerator of aic.w calculation as per wagenmakers and farrell 2004
aic.sum <- aggregate(aic.d.exp~pair, data=t, FUN="sum") # denominator of aic.w equation
names(aic.sum)[2] <- "aic.d.exp.sum"
t <- merge(t, aic.sum, all=T)

t$aic.w <- t$aic.d.exp/t$aic.d.exp.sum # weighted aic

# ratio of weights to see likelihood
aic.w.best <- aggregate(aic.w~pair, data=t, FUN="max") # best weight per pair
names(aic.w.best)[2] <- "aic.w.best"
t <- merge(t, aic.w.best, all=T)
t$aic.ll <- t$aic.w/t$aic.w.best

t.s <- subset(t, region == "VA" | region == "NC" | region == "PA")
t.s <- subset(t.s, (region != "VA" & type != "A") | (region != "VA" & type != "W"))
t.s <- subset(t.s, model != "one_pop")
library(dplyr)
t.s %>% group_by(region, type) %>% count(model)

type.labs <- c("A" = "Appalachian", "B" = "Between Lineages", "W" = "Western")

# weighted AIC
vlines <- data.frame(group = unique(t.s$region))

jpeg("./graphs/AIC_weights/AICw_continuous.jpeg", width=8, height=4, units="in", res=250)
ggplot(data=t.s, aes(x=region, y=aic.w, fill=fct_rev(model)))+
  geom_boxplot(position = position_dodge2(preserve = "single"))+
  facet_wrap(~type, labeller = labeller(type = type.labs))+
  scale_fill_discrete(labels=c("Ancient Migration", 'Migration', 'No Migration', 'Secondary Migration'))+
  scale_fill_manual("Demographic Models", values = c(no_mig = "#00B81F", sc = "#00A5FF", mig = "firebrick2", am = "orchid2"),
                    labels=c('Secondary Migration', "No Migration", "Migration", "Ancient Migration"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Weighted AIC")+
  xlab("Contact Zone")+
  geom_vline(data = vlines, aes(xintercept = 1.5), color = "black", linetype = "solid", linewidth=0.25)
dev.off()

NC <- ggplot(data=subset(t.s, region=="NC"), aes(x=region, y=aic.w, fill=fct_rev(model)))+
        geom_boxplot(position = position_dodge2(preserve = "single"))+
        facet_wrap(~type, labeller = labeller(type = type.labs))+
        scale_fill_discrete(labels=c("Ancient Migration", 'Migration', 'No Migration', 'Secondary Migration'))+
        scale_fill_manual("Demographic Models", values = c(no_mig = "#00B81F", sc = "#00A5FF", mig = "firebrick2", am = "orchid2"),
                          labels=c('Secondary Migration', "No Migration", "Migration", "Ancient Migration"))+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.background = element_blank(), axis.line = element_line(colour = "black"))+
        ylab("Weighted AIC")+
        xlab("Contact Zone")

PA <- ggplot(data=subset(t.s, region=="PA"), aes(x=region, y=aic.w, fill=fct_rev(model)))+
        geom_boxplot(position = position_dodge2(preserve = "single"))+
        facet_wrap(~type, labeller = labeller(type = type.labs))+
        scale_fill_discrete(labels=c("Ancient Migration", 'Migration', 'No Migration', 'Secondary Migration'))+
        scale_fill_manual("Demographic Models", values = c(no_mig = "#00B81F", sc = "#00A5FF", mig = "firebrick2", am = "orchid2"),
                          labels=c('Secondary Migration', "No Migration", "Migration", "Ancient Migration"))+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.background = element_blank(), axis.line = element_line(colour = "black"))+
        ylab("Weighted AIC")+
        xlab("Contact Zone")

jpeg("./graphs/AIC_weights/aic_continuous.jpeg", width=4, height=8, units="in", res=250)
ggarrange(NC, PA, nrow=2, common.legend=T)
dev.off()

# model score 
# t.s.all <- subset(t, region == "NC" | region == "PA" | region == "VA")
# t.s.all <- subset(t.s.all, region != "VA")
# ggplot(data=t.s.all, aes(x=region, y=model_score, color=model))+
#   geom_boxplot(position = position_dodge2(width = 0.9, preserve = "single"))+
#   facet_wrap(~type, labeller = labeller(type = type.labs))+
#   scale_fill_discrete(labels=c('Migration', 'No Migration', 'Secondary Migration'))+
#   theme_bw()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.background = element_blank(), axis.line = element_line(colour = "black"))+
#   ylab("Model Score")+
#   xlab("Contact Zone")





