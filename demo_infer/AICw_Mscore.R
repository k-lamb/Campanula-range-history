library(ggplot2)
library(dplyr)
library(ggpubr)

# plotting AIC weights and Model scores
setwd("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/")
t <- read.csv("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/data_concat/AIC_weights_all_80perc_linreg_ext_2c_ALL.csv", stringsAsFactors = TRUE)

mprop <- seq(0.01,1.0, by=0.01)

# reassign the super model for cyclic models so they are independent of non-cyclic models
t <- t %>% 
  mutate(super_model = if_else((model == "IM_2C" | model == "IM_2C_ae"), "IM2C", 
                               if_else((model == "SC_2C" | model == "SC_2C_ae"), "SC2C", 
                                       if_else((model == "SC_AM" | model == "SC_AM_ae"), "A_SC", super_model))))

t.s.b <- subset(t, lin.pair =="B")
# t.s.a <- subset(t.s, lin.pair == "A")
# t.s.w <- subset(t.s, lin.pair == "W")
# t.s.e <- subset(t.s, lin.pair == "E")

# runs of linreg models do not require division of t.s.a, t.s.w, or t.s.e. commenting out these blocks
test <- data.frame(super_model=factor(), region=factor(), aicw_median=numeric(), mprop=numeric(), lin=factor())
for (i in 1:length(mprop)) {
  t.s <- subset(t, model_score >= mprop[i]) # retain only models where models within 10% of best
  t.s <- subset(t.s, region != "NA")
  # t.s <- subset(t.s, region != "VA")
  
  # t.s.sub <- t.s %>%
  #   filter(model == "SC" | model == "IM" | model == "NM")

  t.s.b <- subset(t.s, lin.pair =="B")
  # t.s.a <- subset(t.s, lin.pair == "A")
  # t.s.w <- subset(t.s, lin.pair == "W")
  # t.s.e <- subset(t.s, lin.pair == "E")
  
  t.s.b.mu <- t.s.b %>%
    group_by(super_model, region) %>%
    # filter(aic.w >= 0.001) %>% 
    dplyr::reframe(aicw_mean = mean(aic.w),
                   aicw_max = max(aic.w),
                   aic.w = aic.w,
                   model = model) %>%
    as.data.frame()

  # t.s.a.mu <- t.s.a %>%
  #   group_by(super_model, region) %>%
  #   dplyr::reframe(aicw_median = median(aic.w),
  #                  aicw_max = max(aic.w),
  #                  aic.w = aic.w,
  #                  model = model) %>%
  #   as.data.frame()
  # 
  # t.s.w.mu <- t.s.w %>%
  #   group_by(super_model, region) %>%
  #   dplyr::reframe(aicw_median = median(aic.w),
  #                  aicw_max = max(aic.w),
  #                  aic.w = aic.w,
  #                  model = model) %>%
  #   as.data.frame()
  # 
  # t.s.e.mu <- t.s.e %>%
  #   group_by(super_model, region) %>%
  #   dplyr::reframe(aicw_median = median(aic.w),
  #                  aicw_max = max(aic.w),
  #                  aic.w = aic.w,
  #                  model = model) %>%
  #   as.data.frame()

  t.s.b.mu$mprop <- mprop[i]
  t.s.b.mu$lin <- "b"
  # t.s.a.mu$mprop <- mprop[i]
  # t.s.a.mu$lin <- "a"
  # t.s.w.mu$mprop <- mprop[i]
  # t.s.w.mu$lin <- "w"
  # t.s.e.mu$mprop <- mprop[i]
  # t.s.e.mu$lin <- "e"

  # cb <- rbind(t.s.b.mu, t.s.a.mu, t.s.w.mu, t.s.e.mu)

  # test <- rbind(test, cb)
  test <- rbind(test, t.s.b.mu)

}

cz.labs <- c("NC"="Rear Edge", 
               "PA"="Leading Edge", 
               "VA"="Virginia")

lin.labs <- c("a"="Within Appalachian", 
               "b"="Between Lineages", 
               "w"="Within Western",
               "e"="Within Eastern")

size.t <- 18
axes.size.t <- 14

test <- test %>%
  mutate(lin2 = case_when(
    lin == "b" ~ "Between Lineages",
    lin == "a" ~ "Within Appalachian",
    lin == "w" ~ "Within Western",
    lin == "e" ~ "Within Eastern"
  ))


# to add 0.1 binned dots to graph below
# dot.panel <- test %>% 
#   subset(round(mprop,2) == 0.10 | mprop == 0.20 | mprop == 0.30 | mprop == 0.40 | mprop == 0.50 |
#          mprop == 0.60 | round(mprop,2) == 0.70 | mprop == 0.80 | mprop == 0.90 | mprop == 1.00)

# better display of the effect of model score on aicw
# aicw_median and aicw_max
ggplot(test, aes(y=aicw_mean, x=mprop, color=super_model))+
  # geom_histogram(position="identity", alpha=0.5)+
  geom_line(linewidth=0.5)+
  # geom_point(data=dot.panel, aes(y=aicw_mean, x=mprop, color=super_model, shape=super_model))+
  # geom_line(data=test, aes(x=mprop, y=aic.w, group=model, color=super_model),linewidth=0.25)+
  theme_bw()+
  xlab("Model score")+
  ylab("AIC weights")+
  ylim(c(-0.05,1.05))+
  xlim(c(-0.05,1.05))+
  # ggtitle("Region-wide SFS (80% missing)")+
  scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                          "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
  facet_grid(lin2 ~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
  theme(text = element_text(size=size.t-2), axis.text=element_text(size=axes.size.t),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), aspect.ratio=1,
        legend.position = "none")



# # model score
# t.s.b.mod <- subset(t, lin.pair =="B")
# t.s.b.mod <- t.s.b.mod %>% 
#   group_by(region) %>% 
#   arrange(model_score) %>% 
#   # slice(21:50) %>% 
#   mutate(model = factor(model, levels = unique(model))) %>% 
#   ungroup()
# 
# nc <- 
# ggplot(t.s.b.mod %>% filter(region=="NC")%>% 
#          group_by(region) %>% 
#          arrange(model_score) %>% 
#          mutate(model = factor(model, levels = unique(model))) %>% 
#          ungroup(), 
#        aes(y=model_score, x=model, color=super_model))+
#   # geom_histogram(position="identity", alpha=0.5)+
#   geom_point()+
#   stat_summary(fun = "median", geom = "crossbar")+
#   # geom_smooth(method="lm", se=F)+
#   theme_bw()+
#   xlab("")+
#   ylab("")+
#   ylim(c(0,1))+
#   scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
#                                           "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
#   # scale_fill_manual("legend", values = c(NM = "#00B81F", SC = "#00A5FF", IM = "firebrick2"))+
#   facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
#   theme(text = element_text(size=size.t-2), axis.text=element_text(size=axes.size.t),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# pa <- 
#   ggplot(t.s.b.mod %>% filter(region=="PA") %>% 
#            group_by(region) %>% 
#            arrange(model_score) %>%
#            mutate(model = factor(model, levels = unique(model))) %>% 
#            ungroup(), 
#          aes(y=model_score, x=model, color=super_model))+
#   # geom_histogram(position="identity", alpha=0.5)+
#   geom_point()+
#   stat_summary(fun = "median", geom = "crossbar")+
#   # geom_smooth(method="lm", se=F)+
#   theme_bw()+
#   xlab("")+
#   ylab("Model Score")+
#   ylim(c(0,1))+
#   scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
#                                           "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
#   # scale_fill_manual("legend", values = c(NM = "#00B81F", SC = "#00A5FF", IM = "firebrick2"))+
#   facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
#   theme(text = element_text(size=size.t-2), axis.text=element_text(size=axes.size.t),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# va <- 
#   ggplot(t.s.b.mod %>% filter(region=="VA") %>% 
#            group_by(region) %>% 
#            arrange(model_score) %>%
#            mutate(model = factor(model, levels = unique(model))) %>% 
#            ungroup(), 
#          aes(y=model_score, x=model, color=super_model))+
#   # geom_histogram(position="identity", alpha=0.5)+
#   geom_point()+
#   stat_summary(fun = "median", geom = "crossbar")+
#   # geom_smooth(method="lm", se=F)+
#   theme_bw()+
#   xlab("Model")+
#   ylab("")+
#   ylim(c(0,1))+
#   scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
#                                           "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
#   # scale_fill_manual("legend", values = c(NM = "#00B81F", SC = "#00A5FF", IM = "firebrick2"))+
#   facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
#   theme(text = element_text(size=size.t-2), axis.text=element_text(size=axes.size.t),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# jpeg("./plots/model_score_80perc_subset.jpeg", res=500, width=14, height=5, unit="in")
# ggarrange(pa,va,nc, ncol=3, common.legend = T)
# dev.off()



# # AIC
t.s.b.mod <- subset(t, lin.pair =="B")
t.s.b.mod <- t.s.b.mod %>% 
  group_by(region) %>% 
  arrange(aic) %>% 
  # slice(1:10) %>% # for main text subset of best models
  mutate(model = factor(model, levels = unique(model))) %>% 
  ungroup()
axes.size.t <- 6

nc <- 
  ggplot(t.s.b.mod %>% filter(region=="NC")%>% 
           group_by(region) %>% 
           arrange(desc(aic)) %>% 
           mutate(model = factor(model, levels = unique(model))) %>% 
           ungroup(), 
         aes(y=aic, x=model, color=super_model))+
  # geom_histogram(position="identity", alpha=0.5)+
  geom_point()+
  stat_summary(fun = "mean", geom = "crossbar")+
  # geom_smooth(method="lm", se=F)+
  theme_bw()+
  xlab("")+
  ylab("")+
  scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                          "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
  # ylim(c(0,1))+
  # scale_fill_manual("legend", values = c(NM = "#00B81F", SC = "#00A5FF", IM = "firebrick2"))+
  facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
  theme(text = element_text(size=size.t-2), axis.text=element_text(size=axes.size.t),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pa <- 
  ggplot(t.s.b.mod %>% filter(region=="PA") %>% 
           group_by(region) %>% 
           arrange(desc(aic)) %>% 
           mutate(model = factor(model, levels = unique(model))) %>% 
           ungroup(), 
         aes(y=aic, x=model, color=super_model))+
  # geom_histogram(position="identity", alpha=0.5)+
  geom_point()+
  stat_summary(fun = "mean", geom = "crossbar")+
  # geom_smooth(method="lm", se=F)+
  theme_bw()+
  scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                          "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"),
                     drop = FALSE)+
  xlab("")+
  ylab("AIC")+
  # ylim(c(0,1))+
  # scale_fill_manual("legend", values = c(NM = "#00B81F", SC = "#00A5FF", IM = "firebrick2"))+
  facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
  theme(text = element_text(size=size.t-2), axis.text=element_text(size=axes.size.t),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

va <- 
  ggplot(t.s.b.mod %>% filter(region=="VA") %>% 
           # group_by(region) %>% 
           arrange(desc(aic)) %>% 
           mutate(model = factor(model, levels = unique(model))) %>% 
           ungroup(), 
         aes(y=aic, x=model, color=super_model))+
  geom_point()+
  stat_summary(fun = "mean", geom = "crossbar")+
  # geom_bar(position="dodge", stat="identity")+
  # geom_smooth(method="lm", se=F)+
  theme_bw()+
  xlab("Model")+
  ylab("")+
  # ylim(c(0,1))+
  scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                          "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
  # scale_fill_manual("legend", values = c(NM = "#00B81F", SC = "#00A5FF", IM = "firebrick2"))+
  facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
  theme(text = element_text(size=size.t-2), axis.text=element_text(size=axes.size.t),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# jpeg("./plots/aic_80perc.jpeg", res=5e2, width=10, height=5, unit="in")
ggarrange(pa,va,nc, ncol=3, common.legend = T, legend = "right")
# dev.off()



# # # AIC weights
# t.s.b.mod <- subset(t, lin.pair =="B")
# t.s.b.mod <- t.s.b.mod %>% 
#   group_by(region) %>% 
#   arrange(aic.w) %>% 
#   mutate(model = factor(model, levels = unique(model))) %>% 
#   ungroup()
# 
# nc <- 
#   ggplot(t.s.b.mod %>% filter(region=="NC")%>% 
#            group_by(region) %>% 
#            arrange((aic.w)) %>% 
#            mutate(model = factor(model, levels = unique(model))) %>% 
#            ungroup(), 
#          aes(y=aic.w, x=model, color=super_model))+
#   # geom_histogram(position="identity", alpha=0.5)+
#   geom_point()+
#   stat_summary(fun = "median", geom = "crossbar")+
#   # geom_smooth(method="lm", se=F)+
#   theme_bw()+
#   xlab("")+
#   ylab("")+
#   # ylim(c(0,1))+
#   scale_fill_discrete(drop = FALSE) +  # Show all levels in the legend
#   scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
#                                           "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
#   # scale_fill_manual("legend", values = c(NM = "#00B81F", SC = "#00A5FF", IM = "firebrick2"))+
#   facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
#   theme(text = element_text(size=size.t-2), axis.text=element_text(size=axes.size.t),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# pa <- 
#   ggplot(t.s.b.mod %>% filter(region=="PA") %>% 
#            group_by(region) %>% 
#            arrange((aic.w)) %>% 
#            mutate(model = factor(model, levels = unique(model))) %>% 
#            ungroup(), 
#          aes(y=aic.w, x=model, color=super_model))+
#   # geom_histogram(position="identity", alpha=0.5)+
#   geom_point()+
#   stat_summary(fun = "median", geom = "crossbar")+
#   # geom_smooth(method="lm", se=F)+
#   theme_bw()+
#   xlab("")+
#   ylab("AIC weights")+
#   scale_fill_discrete(drop = FALSE) +  # Show all levels in the legend
#   # ylim(c(0,1))+
#   scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
#                                           "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
#   # scale_fill_manual("legend", values = c(NM = "#00B81F", SC = "#00A5FF", IM = "firebrick2"))+
#   facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
#   theme(text = element_text(size=size.t-2), axis.text=element_text(size=axes.size.t),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# va <- 
#   ggplot(t.s.b.mod %>% filter(region=="VA") %>% 
#            group_by(region) %>% 
#            arrange((aic.w)) %>% 
#            mutate(model = factor(model, levels = unique(model))) %>% 
#            ungroup(), 
#          aes(y=aic.w, x=model, color=super_model))+
#   geom_point()+
#   stat_summary(fun = "median", geom = "crossbar")+
#   # geom_bar(position="dodge", stat="identity")+
#   # geom_smooth(method="lm", se=F)+
#   theme_bw()+
#   xlab("Model")+
#   ylab("")+
#   # ylim(c(0,1))+
#   # scale_fill_manual("legend", values = c(NM = "#00B81F", SC = "#00A5FF", IM = "firebrick2"))+
#   facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
#   theme(text = element_text(size=size.t-2), axis.text=element_text(size=axes.size.t),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# va
# 
# jpeg("./plots/aic_weights_80perc_subset.jpeg", res=5e2, width=14, height=5, unit="in")
# ggarrange(pa,va,nc, ncol=3, common.legend = T)
# dev.off()
