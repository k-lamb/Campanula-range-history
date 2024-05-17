# test for differences in migration rates, t split, t mig

library(ggplot2)
library(ggpubr)
library(forcats)
library(ggpattern)
library(dplyr)
library(car) # Anova significance tests
library(effects) # standard error calculations from model means

setwd("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/")

meta = T # meta==F if working with populations aggregated by lineage and region; TRUE if working with populations separated
mu = 2.8e-9
L <- read.csv("./easySFS/pop_easySFS.csv") # SNPs retained from easySFS
perc = 80 # percentage SNP missing filter

# read in and grab best model where applicable
if (meta == T) {
  if (perc == 65) {
    t.s <- read.csv("./data_concat/AIC_weights_all_65perc_linreg_ext_2c_ALL.csv", stringsAsFactors = T)
  }
  if (perc == 80) {
    t.s <- read.csv("./data_concat/AIC_weights_all_80perc_linreg_ext_2c_ALL.csv", stringsAsFactors = T)
  }
}

# best model by region for METAPOPULATION MODELS
if (meta == F) {
  if (perc == 65) {
    t.s <- read.csv("./data_concat/AIC_weights_all_65perc.csv", stringsAsFactors = T)
  }
  if (perc == 80) {
    t.s <- read.csv("./data_concat/AIC_weights_all_80perc.csv", stringsAsFactors = T)
  }
  # for POPULATION-LEVEL MODELS
  # t.s <- t %>%
  #   subset(model == "IM_ae")
  
  # NEEDED FOR POP PAIRS BUT NOT LIN REG PAIRS
  # fix issue where mij/mji are flipped because populations in metadata did not respect lineage order
  t.s <- t.s %>% 
    mutate(m12T = if_else((pop1.lin != pop2.lin & ((pop1.lin == "W" & pop2.lin == "A") | (pop1.lin == "E" & pop2.lin == "A"))), m21, m12)) %>% 
    mutate(m21T = if_else((pop1.lin != pop2.lin & ((pop1.lin == "W" & pop2.lin == "A") | (pop1.lin == "E" & pop2.lin == "A"))), m12, m21)) %>% 
    mutate(p1 = if_else((pop1.lin != pop2.lin & ((pop1.lin == "W" & pop2.lin == "A") | (pop1.lin == "E" & pop2.lin == "A"))), pop2, pop1)) %>% 
    mutate(p2 = if_else((pop1.lin != pop2.lin & ((pop1.lin == "W" & pop2.lin == "A") | (pop1.lin == "E" & pop2.lin == "A"))), pop1, pop2)) %>% 
    mutate(Pair_name = paste0(p1,".",p2),
           pop1 = p1,
           pop2 = p2,
           m12=m12T,
           m21=m21T) %>% 
    dplyr::select(-p1, -p2, -m12T, -m21T)

  # fix issue where mij/mji are flipped for reverse models
  temp.rev <- subset(t.s, reversed == "TRUE")
  temp.for <- subset(t.s, reversed == "FALSE" | is.na(reversed) == T)
  temp.mij <- temp.rev$mji
  temp.mji <- temp.rev$mij
  temp.rev$mij <- temp.mij
  temp.rev$mji <- temp.mji
  t.s <- rbind(temp.rev, temp.for)
  
}

# fix variable names and make L appropriate (as possible)
if (meta == T) {
  L = 3.6e3
}

if (meta == F) {
  L <- L %>% dplyr::select(POP, BEST.SNP)
  names(L) <- c("pop1", "L1")
  t.s <- merge(t.s, L, by="pop1")
  names(L) <- c("pop2", "L2")
  t.s <- merge(t.s, L, by="pop2")
  t.s <- t.s %>% mutate(L = if_else(L1 <= L2, L1, L2))
  t.s$L <- t.s$L*(159+169)
}

t.s <- t.s %>% 
  mutate(super_model = if_else((model == "IM_2C" | model == "IM_2C_ae"), "IM2C", 
                               if_else((model == "SC_2C" | model == "SC_2C_ae"), "SC2C", 
                                       if_else((model == "SC_AM" | model == "SC_AM_ae"), "A_SC", super_model)))) %>% 
  mutate(model = if_else((model == "SC_AM"), "A_SC",
                          if_else(model == "SC_AM_ae", "A_SC_ae", model)))


if (meta == T) {
  # convert variables out of 
  t.s$Tsplit <- ((2*t.s$theta)/(4*mu*L))*t.s$Ts
  t.s$Tmig <-  ((2*t.s$theta)/(4*mu*L))*t.s$Tsc
  t.s$mij <- t.s$m12/(2*(t.s$theta/(4*mu*L)))
  t.s$mji <- t.s$m21/(2*(t.s$theta/(4*mu*L)))
  t.s$mig1 <- t.s$mij*(2*(t.s$theta/(4*mu*L)))*t.s$nu1
  t.s$mig2 <- t.s$mji*(2*(t.s$theta/(4*mu*L)))*t.s$nu2
  # t.s <- subset(t.s, region == "NC" | region == "PA" )
  
  t.s.b <- subset(t.s, lin.pair =="B")
}

if (meta == F) {
  # convert variables out of 
  t.s$Tsplit <- ((2*t.s$theta)/(4*mu*t.s$L))*t.s$Ts
  t.s$Tmig <-  ((2*t.s$theta)/(4*mu*t.s$L))*t.s$Tsc
  t.s$mij <- t.s$m12/(2*(t.s$theta/(4*mu*t.s$L)))
  t.s$mji <- t.s$m21/(2*(t.s$theta/(4*mu*t.s$L)))
  t.s$mig1 <- t.s$mij*(2*(t.s$theta/(4*mu*t.s$L)))*t.s$nu1
  t.s$mig2 <- t.s$mji*(2*(t.s$theta/(4*mu*t.s$L)))*t.s$nu2
  # t.s <- subset(t.s, region == "NC" | region == "PA" )
  
  t.s.b <- t.s %>% 
    mutate(lin.pair = if_else(pop1.lin != pop2.lin, "B", pop1.lin)) %>% 
    group_by(Pair_name) %>% 
    mutate(aic.w = MuMIn::Weights(aic)) %>% # overwriting existing aic.w to try to figure out why order of aic.w != order of aic
    mutate_at(vars(lin.pair), as.factor) %>% 
    mutate_at(vars(aic.w), as.numeric)
}

# t.s.a <- subset(t.s, lin.pair == "A")
# t.s.w <- subset(t.s, lin.pair == "W")
# t.s.w.e <- subset(t.s, lin.pair == "E")

# tables by mode best model
if (meta == F) {
  table(t.s.b %>% filter(aic < aic.min+2 & region=="NC" & lin.pair=="B") %>% dplyr::select(model)) %>% colSums()
  table(t.s.b %>% filter(aic < aic.min+2 & region=="PA" & lin.pair=="B") %>% dplyr::select(model)) %>% colSums()
  table(t.s.b %>% filter(aic < aic.min+2 & region=="VA" & lin.pair=="B") %>% dplyr::select(model)) %>% colSums()
}



# whether to get best models by mode AICmin (mode == T) or by mean AICmin (mode == F)
mode = F

if (meta == T) {
  # NC: AM, IM_ae, NM_ae, SC -- PA: AM_ae_br, IM_ae, NM_ae, SC_ae -- VA: AM, IM_ae, NM_ae, SC_ae
  if (perc == 65) {
    t.s.b.1 <- t.s.b %>%
      filter((region == "NC" & model == "AM") | (region == "NC" & model == "IM_ae") | (region == "NC" & model == "SC") | (region == "NC" & model == "NM") | (region == "NC" & model == "SC_2C") | (region == "NC" & model == "IM_2C_ae") |
               (region == "VA" & model == "AM") | (region == "VA" & model == "IM_ae") |  (region == "VA" & model == "SC_ae") | (region == "VA" & model == "NM_ae") |  (region == "VA" & model == "SC_2C") | (region == "VA" & model == "IM_2C_ae") |
               (region == "PA" & model == "AM_ae_br") | (region == "PA" & model == "IM_ae") | (region == "PA" & model == "SC_ae") | (region == "PA" & model == "NM_ae_br")  | (region == "PA" & model == "SC_2C") | (region == "PA" & model == "IM_2C_ae"))
  }
  if (perc == 80) {
    t.s.b.1 <- t.s.b %>%
      filter((region == "NC" & model == "AM") | (region == "NC" & model == "SC") | (region == "NC" & model == "IM_ae") | (region == "NC" & model == "NM_ae") | (region == "NC" & model == "SC_2C") | (region == "NC" & model == "IM_2C_ae") |
               (region == "VA" & model == "AM") | (region == "VA" & model == "SC_ae") |  (region == "VA" & model == "IM_2C_ae") | (region == "VA" & model == "SC_2C") |  (region == "VA" & model == "IM_ae") | (region == "VA" & model == "NM_ae") |
               (region == "PA" & model == "AM_ae_br") | (region == "PA" & model == "IM_ae") | (region == "PA" & model == "SC_ae") | (region == "PA" & model == "NM_ae")  | (region == "PA" & model == "IM_2C_ae") | (region == "PA" & model == "SC_2C"))
  }
}

if (meta == F) {
  if (perc == 65) {
    # mean AIC only: 
    # NC: NM_ae_br, IM, SC -- PA: IM_ae, NM_ae, SC -- VA: SC, IM, NM_ae
      t.s.b.1 <- t.s.b %>%
        filter((region == "NC" & model == "SC") | (region == "NC" & model == "IM") | (region == "NC" & model == "NM_ae_br") |
               (region == "VA" & model == "SC") | (region == "VA" & model == "IM") | (region == "VA" & model == "NM_ae") |
               (region == "PA" & model == "IM_ae") | (region == "PA" & model == "SC") | (region == "PA" & model == "NM_ae"))
  }
  if (perc == 80) {
    if (mode == T) {
      # mode AIC
      t.s.b.1 <- t.s.b %>%
        filter((region == "NC" & model == "AM_ae") | (region == "NC" & model == "SC") | (region == "NC" & model == "IM_ae") | (region == "NC" & model == "NM_ae_b") |
                 (region == "VA" & model == "AM_br") | (region == "VA" & model == "SC") | (region == "VA" & model == "IM") | (region == "VA" & model == "NM") |
                 (region == "PA" & model == "AM_br") | (region == "PA" & model == "IM_ae") | (region == "PA" & model == "SC") | (region == "PA" & model == "NM_ae"))
    }
    if (mode == F) {
      # mean AIC
      t.s.b.1 <- t.s.b %>%
        filter((region == "NC" & model == "AM_br") | (region == "NC" & model == "SC") | (region == "NC" & model == "IM_ae") | (region == "NC" & model == "NM_ae_b") | (region == "NC" & model == "SC_2C") | (region == "NC" & model == "IM_2C") |
                 (region == "VA" & model == "AM_ae_br") | (region == "VA" & model == "SC") | (region == "VA" & model == "IM") | (region == "VA" & model == "NM") | (region == "VA" & model == "IM_2C") | (region == "VA" & model == "SC_2C") |
                 (region == "PA" & model == "AM_br") | (region == "PA" & model == "IM_ae") | (region == "PA" & model == "SC") | (region == "PA" & model == "NM_ae") | (region == "PA" & model == "IM_2C_ae") | (region == "PA" & model == "SC_2C"))
    }
  }
}



# PLOTTING AIC 
aic.w = T # plot AIC weights
aic = T # plot AIC (raw)
model.score = T

cz.labs <- c("NC"="Rear Edge", 
               "PA"="Leading Edge", 
               "VA"="Virginia")

if (meta == meta) { # originally had this plotting for only meta == F but now I want it to be able to plot either
    axes.size.t = 17
    size.t = 20
    aspect = 0.5
    # AIC Weights plotting
    if (aic.w == T) {
      t.s.b.mod <- subset(t.s.b.1, lin.pair =="B")
      t.s.b.mod <- t.s.b.mod %>%
        group_by(Pair_name, super_model) %>% 
        arrange(aic) %>% 
        slice(1) %>% # select only one the best model within the super_model
        ungroup() %>% 
        group_by(Pair_name) %>% 
        mutate(aic.w = MuMIn::Weights(aic)) %>% # recalculate weights without additional models interfering with selection
        mutate_at(vars(aic.w), as.numeric) %>% 
        ungroup()
      t.s.b.mod <- t.s.b.mod %>% 
        group_by(region, lin.pair, model, super_model) %>% 
        reframe(mean.aic = mean(aic.w),
                median.aic = median(aic.w)) %>% 
        ungroup()

      
      nc.w <- 
        ggplot(t.s.b.mod %>% filter(region=="NC")%>% 
                 arrange((mean.aic)) %>% 
                 mutate(model = factor(model, levels = unique(model))) %>%
                 ungroup(), 
               aes(y=mean.aic, x=model, color=super_model))+
        geom_point()+
        stat_summary(fun = "mean", geom = "crossbar")+
        # geom_point(data=t.s.b %>% 
        #              filter(lin.pair == "B" & region == "NC"),
        #            aes(x=model, y=aic.w), alpha=0.25)+
        theme_bw()+
        xlab("")+
        ylab("")+
        # geom_hline(data=t.s.b.mod %>% filter(region=="NC"), aes(yintercept = min(mean.aic)+2), color="red", linetype=2)+
        scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                                "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
        ylim(c(0,1))+
        facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
        theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=axes.size.t), aspect.ratio = aspect)
      
      pa.w <- 
        ggplot(t.s.b.mod %>% filter(region=="PA") %>% 
                 arrange((mean.aic)) %>% 
                 mutate(model = factor(model, levels = unique(model))) %>%
                 ungroup(), 
               aes(y=mean.aic, x=model, color=super_model))+
        geom_point()+
        stat_summary(fun = "mean", geom = "crossbar")+
        # geom_point(data=t.s.b %>%
        #              filter(lin.pair == "B" & region == "PA"),
        #            aes(x=model, y=aic.w), alpha=0.25)+
        theme_bw()+
        scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                                "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"),
                           drop = FALSE)+
        xlab("")+
        ylab("AIC Weights")+
        # geom_hline(data=t.s.b.mod %>% filter(region=="PA"), aes(yintercept = min(mean.aic)+2), color="red", linetype=2)+
        ylim(c(0,1))+
        facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
        theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=axes.size.t), aspect.ratio = aspect)
      
      va.w <- 
        ggplot(t.s.b.mod %>% filter(region=="VA") %>% 
                 arrange((mean.aic)) %>% 
                 mutate(model = factor(model, levels = unique(model))) %>%
                 ungroup(), 
               aes(y=mean.aic, x=model, color=super_model))+
        geom_point()+
        stat_summary(fun = "mean", geom = "crossbar")+
        # geom_point(data=t.s.b %>% 
        #              filter(lin.pair == "B" & region == "VA"),
        #            aes(x=model, y=aic.w), alpha=0.25)+
        theme_bw()+
        scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                                "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"),
                           drop = FALSE)+
        xlab("")+
        ylab("")+
        # geom_hline(data=t.s.b.mod %>% filter(region=="PA"), aes(yintercept = min(mean.aic)+2), color="red", linetype=2)+
        ylim(c(0,1))+
        facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
        theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=axes.size.t), aspect.ratio = aspect)
      
      # jpeg("./plots/aic_pop_80perc.jpeg", res=5e2, width=10, height=5, unit="in")
      # ggarrange(pa.w,va.w,nc.w, ncol=3, common.legend = T, legend = "right")
      # dev.off()
  }
  
  if (aic == T) {
      t.s.b.mod <- subset(t.s.b, lin.pair =="B")
      t.s.b.mod <- t.s.b.mod %>% 
        group_by(region, lin.pair, model, super_model) %>% 
        reframe(mean.aic = mean(aic),
                mean.aicw = mean(aic.w),
                median.aic = median(aic),
                median.aicw = median(aic.w),
                min.aic = min(aic),
                max.aic = max(aic),
                sd.aic = sd(aic)) %>% 
        ungroup()
      
      nc <- 
        ggplot(t.s.b.mod %>% filter(region=="NC")%>% 
                 arrange(desc(mean.aic)) %>% 
                 mutate(model = factor(model, levels = unique(model))) %>%
                 ungroup(), 
               aes(y=mean.aic, x=model, color=super_model))+
        geom_point()+
        stat_summary(fun = "mean", geom = "crossbar")+
        # geom_point(data=t.s.b %>%
        #              filter(region=="NC"),
        #            aes(x=model, y=aic, color=super_model), alpha=0.25)+
        theme_bw()+
        xlab("")+
        ylab("")+
        geom_hline(data=t.s.b.mod %>% filter(region=="NC"), aes(yintercept = min(mean.aic)+2), color="red", linetype=2)+
        scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                                "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
        facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
        theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=axes.size.t), aspect.ratio = aspect)
      
      pa <- 
        ggplot(t.s.b.mod %>% filter(region=="PA") %>% 
                 arrange(desc(mean.aic)) %>% 
                 mutate(model = factor(model, levels = unique(model))) %>%
                 ungroup(), 
               aes(y=mean.aic, x=model, color=super_model))+
        geom_point()+
        stat_summary(fun = "mean", geom = "crossbar")+
        # geom_point(data=t.s.b %>%
        #              filter(region=="PA"),
        #            aes(x=model, y=aic, color=super_model), alpha=0.25)+
        theme_bw()+
        scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                                "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"),
                           drop = FALSE)+
        xlab("")+
        ylab("AIC")+
        geom_hline(data=t.s.b.mod %>% filter(region=="PA"), aes(yintercept = min(mean.aic)+2), color="red", linetype=2)+
        facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
        theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=axes.size.t), aspect.ratio = aspect)
      
      va <- 
        ggplot(t.s.b.mod %>% filter(region=="VA") %>% 
                 arrange(desc(mean.aic)) %>% 
                 mutate(model = factor(model, levels = unique(model))) %>%
                 ungroup(), 
               aes(y=mean.aic, x=model, color=super_model))+
        geom_point()+
        stat_summary(fun = "mean", geom = "crossbar")+
        geom_hline(data=t.s.b.mod %>% filter(region=="VA"), aes(yintercept = min(mean.aic)+2), color="red", linetype=2)+
        # geom_point(data=t.s.b %>%
        #              filter(region=="VA"),
        #            aes(x=model, y=aic, color=super_model), alpha=0.25)+
        theme_bw()+
        xlab("Model")+
        ylab("")+
        scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                                "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
        facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
        theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=axes.size.t), aspect.ratio = aspect)
      
      
      # jpeg("./plots/aic_pop_80perc.jpeg", res=5e2, width=10, height=5, unit="in")
      # ggarrange(pa,va,nc, ncol=3, common.legend = T, legend = "right")
      # dev.off()
  }
    
    if (model.score == T) {
      t.s.b.mod <- subset(t.s.b, lin.pair =="B")
      t.s.b.mod <- t.s.b.mod %>% 
        group_by(region, lin.pair, model, super_model) %>% 
        reframe(mean.ms = mean(model_score),
                median.ms = median(model_score)) %>% 
        ungroup()
      
      nc.ms <- 
        ggplot(t.s.b.mod %>% filter(region=="NC")%>% 
                 arrange((mean.ms)) %>% 
                 mutate(model = factor(model, levels = unique(model))) %>%
                 ungroup(), 
               aes(y=mean.ms, x=model, color=super_model))+
        geom_point()+
        stat_summary(fun = "mean", geom = "crossbar")+
        # geom_point(data=t.s.b %>%
        #              filter(region=="NC"),
        #            aes(x=model, y=aic, color=super_model), alpha=0.25)+
        theme_bw()+
        xlab("")+
        ylab("")+
        ylim(c(0,1))+
        # geom_hline(data=t.s.b.mod %>% filter(region=="NC"), aes(yintercept = min(mean.ms)), color="red", linetype=2)+
        scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                                "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
        facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
        theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=axes.size.t), aspect.ratio = aspect)
      
      pa.ms <- 
        ggplot(t.s.b.mod %>% filter(region=="PA") %>% 
                 arrange((mean.ms)) %>% 
                 mutate(model = factor(model, levels = unique(model))) %>%
                 ungroup(), 
               aes(y=mean.ms, x=model, color=super_model))+
        geom_point()+
        stat_summary(fun = "mean", geom = "crossbar")+
        # geom_point(data=t.s.b %>%
        #              filter(region=="PA"),
        #            aes(x=model, y=aic, color=super_model), alpha=0.25)+
        theme_bw()+
        scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                                "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"),
                           drop = FALSE)+
        xlab("")+
        ylab("Model Score")+
        ylim(c(0,1))+
        # geom_hline(data=t.s.b.mod %>% filter(region=="PA"), aes(yintercept = minX(mean.aic)+2), color="red", linetype=2)+
        facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
        theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=axes.size.t), aspect.ratio = aspect)
      
      va.ms <- 
        ggplot(t.s.b.mod %>% filter(region=="VA") %>% 
                 arrange((mean.ms)) %>% 
                 mutate(model = factor(model, levels = unique(model))) %>%
                 ungroup(), 
               aes(y=mean.ms, x=model, color=super_model))+
        geom_point()+
        stat_summary(fun = "mean", geom = "crossbar")+
        # geom_hline(data=t.s.b.mod %>% filter(region=="VA"), aes(yintercept = min(mean.aic)+2), color="red", linetype=2)+
        # geom_point(data=t.s.b %>%
        #              filter(region=="VA"),
        #            aes(x=model, y=aic, color=super_model), alpha=0.25)+
        theme_bw()+
        xlab("Model")+
        ylab("")+
        ylim(c(0,1))+
        scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                                "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
        facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
        theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=axes.size.t), aspect.ratio = aspect)
      
      
      # jpeg("./plots/aic_pop_80perc.jpeg", res=5e2, width=10, height=5, unit="in")
      # ggarrange(pa,va,nc, ncol=3, common.legend = T, legend = "right")
      # dev.off()
    }
    
}


# if (meta == T) {
#   jpeg("./plots/comparison/meta65.jpeg", res=5e2, width=10, height=10, unit="in")
  # ggarrange(pa.w,va.w,nc.w, # AIC Weights models
  #           pa,va,nc, # AIC models
  #           pa.ms, va.ms, nc.ms,
  #           ncol=3, nrow=3, common.legend = T, legend = "none")
#   dev.off()
# }
# 
# if (meta == F) {
#   jpeg("./plots/comparison/pop65.jpeg", res=5e2, width=10, height=10, unit="in")
#   ggarrange(pa.w,va.w,nc.w, # AIC Weights models
#             pa,va,nc, # AIC models
#             pa.ms, va.ms, nc.ms,
#             ncol=3, nrow=3, common.legend = T, legend = "none")
#   dev.off()
# }



# MODEL DEMOGRAPHIC PARAMETERS

supp.labs <- c("NC"="Rear Edge",
               "PA"="Leading Edge", 
               "VA"="Virginia")

# get the best model in each region
# for AIC: PA=NM_ae; NC=IM_ae; VA=AM_ae_br
# for AIC weights: PA=IM_ae; NC=AM_ae_b; VA=AM_b
aicw=T
if (meta == T) {
  if (perc == 65) {
    # AIC W == AIC for meta models
    pa.sub <- t.s.b %>% filter(region == "PA" & model == "IM_ae" & lin.pair == "B")
    nc.sub <- t.s.b %>% filter(region == "NC" & model == "IM_ae" & lin.pair == "B")
    va.sub <- t.s.b %>% filter(region == "VA" & model == "IM_ae" & lin.pair == "B")
  }
  if (perc == 80) {
    # AIC W == AIC for meta models
    pa.sub <- t.s.b %>% filter(region == "PA" & model == "IM_ae" & lin.pair == "B")
    nc.sub <- t.s.b %>% filter(region == "NC" & model == "IM_ae" & lin.pair == "B")
    va.sub <- t.s.b %>% filter(region == "VA" & model == "IM_ae" & lin.pair == "B")
  }
 }

if (meta == F) {
  if (perc == 65) {
    pa.sub <- t.s.b %>% filter(region == "PA" & model == "IM_ae" & lin.pair == "B")
    nc.sub <- t.s.b %>% filter(region == "NC" & model == "NM_ae_br" & lin.pair == "B")
    va.sub <- t.s.b %>% filter(region == "VA" & model == "SC" & lin.pair == "B")
  }
  
  if (perc == 80) {
    if (aicw==T) {
      pa.sub <- t.s.b %>% filter(region == "PA" & model == "IM_ae" & lin.pair == "B")
      nc.sub <- t.s.b %>% filter(region == "NC" & model == "IM_ae" & lin.pair == "B")
      va.sub <- t.s.b %>% filter(region == "VA" & model == "AM_ae_br" & lin.pair == "B")
    }
    if (aicw==F) {
      pa.sub <- t.s.b %>% filter(region == "PA" & model == "NM_ae" & lin.pair == "B")
      nc.sub <- t.s.b %>% filter(region == "NC" & model == "IM_ae" & lin.pair == "B")
      va.sub <- t.s.b %>% filter(region == "VA" & model == "AM_ae_br" & lin.pair == "B")
    }
  }
}

t.s.b.reg <- rbind(pa.sub, nc.sub, va.sub)

# Time of split between lineages
t.s.b.reg <- t.s.b.reg %>% mutate(Tsplit2 = if_else(is.na(Tmig)==T, Tsplit, Tmig+Tsplit)) # Tmig + Tsplit = the time of lineage divergence
mu.Ts <- mean(t.s.b.reg$Tsplit2)
t.s.b.reg <- t.s.b.reg %>% 
  mutate(Tsplit.rel = Tsplit2/mu.Ts)
t.s.b.reg$placeholder <- "Relative time of lineage divergence"

t.s.b.reg <- t.s.b.reg %>% 
  mutate(name = case_when(region == "NC" ~ "Rear Edge",
                          region == "PA" ~ "Leading Edge",
                          region == "VA" ~ "Virginia")) %>% 
  mutate(name = factor(name, levels=c("Leading Edge", "Virginia", "Rear Edge")))

if (meta == T) {
  time <- ggplot()+
    geom_hline(yintercept = 1, linetype=2, color="gray50")+
    geom_boxplot(data=t.s.b.reg, aes(x=name, y=Tsplit.rel, group=region), alpha=0.4, width=0.25)+
    geom_point(data=t.s.b.reg, aes(x=name, y=Tsplit.rel, group=region), size=2)+
    stat_summary(fun = "median", geom = "crossbar")+
    theme_bw()+
    ylab("Relative time of divergence")+
    xlab(" ")+
    facet_wrap(~ placeholder)+
    theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 0, size=axes.size.t), aspect.ratio = aspect)
}

if (meta == F) {
  tsplit.lm <- lm(Tsplit.rel ~ region, data=t.s.b.reg)
  car::Anova(tsplit.lm, type=3)
  
  lm.se <- effect("region", tsplit.lm) %>% as.data.frame()
  lm.se <- lm.se %>% 
    mutate(region.name = case_when(region == "NC" ~ "Rear Edge",
                                   region == "VA" ~ "Virginia",
                                   region == "PA" ~ "Leading Edge"))
  
  lm.se$region.name <- factor(lm.se$region.name, levels = c("Leading Edge", "Virginia", "Rear Edge"))
  
  time <- ggplot()+
    geom_hline(yintercept = 1, linetype=2, color="gray50")+
    geom_point(data=lm.se, aes(x=region.name, y=fit, color=region.name), size=4)+
    geom_errorbar(data=lm.se, aes(x=region.name, ymin=(fit-se), ymax=(fit+se), color=region.name), width=0.2)+
    theme_bw()+
    # ylim(c(0,1))+
    # ylab("Time of Split")+
    # xlab("Contact Zone")+
    facet_wrap(~ placeholder)+
    theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 0, size=axes.size.t), aspect.ratio = aspect)
}

# ggplot(t.s.b, aes(x=super_model, y=Tsplit, fill=super_model))+
#   # geom_violin(width=0.5, alpha=0.4)+
#   geom_point(position = position_dodge2(width = 0.9, preserve = "single"))+
#   stat_summary(fun = "median", geom = "crossbar")+
#   theme_bw()+
#   xlab("between lineage comparisons")+
#   ylab("Tsplit")+
#   # ylim(c(0,5e6))+
#   scale_fill_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
#                                           "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
#   facet_wrap(~ region, labeller = labeller(region = supp.labs))+
#   theme(strip.text.x = element_text(size=size.t), axis.text=element_text(size=size.t),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.background = element_blank(), axis.line = element_line(colour = "black"),
#         legend.position = "none")

### MIGRATION
if (meta == F) {
  nc.dat <- t.s.b.reg %>% filter(region=="NC")
  t.test(nc.dat$mij, nc.dat$mji, paired=T, alternative = "greater")
  
  va.dat <- t.s.b.reg %>% filter(region=="VA")
  t.test(va.dat$mij, va.dat$mji, paired=T, alternative = "greater")
  
  pa.dat <- t.s.b.reg %>% filter(region=="PA")
  t.test(pa.dat$mij, pa.dat$mji, paired=T, alternative = "greater")
}

# ggplot(t.s.b.mod %>% filter(region=="NC")%>% 
#          arrange((mean.aic)) %>% 
#          mutate(model = factor(model, levels = unique(model))) %>%
#          ungroup(), 
#        aes(y=mean.aic, x=model, color=super_model))+
#   geom_point()+
#   stat_summary(fun = "mean", geom = "crossbar")+
#   # geom_point(data=t.s.b %>% 
#   #              filter(lin.pair == "B" & region == "NC"),
#   #            aes(x=model, y=aic.w), alpha=0.25)+
#   theme_bw()+
#   xlab("")+
#   ylab("")+
#   # geom_hline(data=t.s.b.mod %>% filter(region=="NC"), aes(yintercept = min(mean.aic)+2), color="red", linetype=2)+
#   scale_color_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
#                                           "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
#   ylim(c(0,1))+
#   facet_wrap(~ region, labeller = labeller(region = cz.labs, lin2=label_wrap_gen(width=10)))+
#   theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=axes.size.t), aspect.ratio = aspect)

#mij between lineages

size.t.mig = 12
t.s.b.reg$placeholder <- "Ratio of migration between lineages"

a <-
ggplot(data=t.s.b.reg, aes(x=name, y=mij, fill=super_model))+
  # geom_violin(width=0.5, alpha=0.4)+
  geom_point(position = position_dodge2(width = 0.9, preserve = "single"))+
  stat_summary(fun = "median", geom = "crossbar")+
  theme_bw()+
  xlab("between lineage comparisons")+
  ylab("Migration i <- j")+
  ylim(c(0,2e-6))+ # meta=T & perc=80
  # ylim(c(0,6e-7))+ # meta=T & perc=65
  scale_fill_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                         "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
  # facet_wrap(~ placeholder)+
  theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, size=axes.size.t), aspect.ratio = aspect)

# mji between lineages
b <-
ggplot(data=t.s.b.reg, aes(x=name, y=mji, fill=super_model))+
  # geom_violin(width=0.5, alpha=0.4)+
  geom_point(position = position_dodge2(width = 0.9))+
  stat_summary(fun = "median", geom = "crossbar")+
  theme_bw()+
  xlab("between lineage comparisons")+
  ylab("Migration j <- i")+
  ylim(c(0,2e-6))+ # meta=T & perc=80
  # ylim(c(0,6e-7))+ # meta=T & perc=65
  scale_fill_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                         "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+  
  # facet_wrap(~ placeholder)+
  theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, size=axes.size.t), aspect.ratio = aspect)

# mig1 between lineages
c <-
  ggplot(data=t.s.b.reg, aes(x=name, y=mig1, fill=super_model))+
  # geom_violin(width=0.5, alpha=0.4)+
  geom_point(position = position_dodge2(width = 0.9, preserve = "single"))+
  stat_summary(fun = "median", geom = "crossbar")+
  theme_bw()+
  xlab("between lineage comparisons")+
  ylab("Number of Migrants i <- j")+
  ylim(c(0,5))+ # meta=T & perc=80
  # ylim(c(0,11))+ # meta=T & perc=65
  scale_fill_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                         "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
  # facet_wrap(~ placeholder)+
  theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0,  size=axes.size.t), aspect.ratio = aspect)

# mji between lineages
d <-
  ggplot(data=t.s.b.reg, aes(x=name, y=mig2, fill=super_model))+
  # geom_violin(width=0.5, alpha=0.4)+
  geom_point(position = position_dodge2(width = 0.9))+
  stat_summary(fun = "median", geom = "crossbar")+
  theme_bw()+
  xlab("between lineage comparisons")+
  ylab("Number of Migrants j <- i")+
  ylim(c(0,5))+  # meta=T & perc=80
  # ylim(c(0,11))+ # meta=T & perc=65
  scale_fill_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                         "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+
  # facet_wrap(~ placeholder)+
  theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, size=axes.size.t), aspect.ratio = aspect)

# ggpubr::ggarrange(a,b,c,d, ncol=2, nrow=2, legend = "none")



### RATIO APPROACH
ratio <- t.s.b.reg
ratio$mig.ratio <- ratio$mig1/ratio$mig2
ratio$mij.ratio <- ratio$mij/ratio$mji

ratio$placeholder <- "Ratio of migration between lineages"

e <- 
ggplot(ratio, aes(x=name, y=mij.ratio, fill=super_model))+
  geom_hline(yintercept = 1, linetype=2)+
  geom_point(data=ratio, aes(x=name, y=mij.ratio, group=region), size=2)+
  geom_boxplot(alpha=0.4, width=0.25)+
  # geom_point(position=position_dodge2(width=0.25))+
  # stat_summary(fun = "mean", geom = "crossbar")+
  # ggrepel::geom_label_repel(data=ratio %>% filter(region != "PA"), aes(x=region, y=mij.ratio, label=Pair_name), alpha=0.75, size=0.5)+
  theme_bw()+
  xlab(" ")+
  ylab("Ratio of migration rates")+
  # ylim(c(0,1.25e-6))+
  # ylim(c(0,10))+
  # scale_fill_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
  #                                        "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+  
  facet_wrap(~ placeholder)+
  theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, size=axes.size.t), aspect.ratio = aspect)


f <- 
ggplot(ratio, aes(x=name, y=mig.ratio, fill=super_model))+
  geom_hline(yintercept = 1, linetype=2)+
  geom_boxplot(alpha=0.4, width=0.25)+
  # geom_point(position=position_dodge2(width=0.25))+
  # stat_summary(fun = "mean", geom = "crossbar")+
  # ggrepel::geom_label_repel(data=ratio %>% filter(region != "PA"), aes(x=region, y=mij.ratio, label=Pair_name), alpha=0.75, size=0.5)+
  theme_bw()+
  # xlab("Between lineage comparisons")+
  # ylab("Ratio of Migration")+
  # ylim(c(0,1.25e-6))+
  # ylim(c(0,10))+
  scale_fill_manual("legend", values = c("NM" = "#00B81F", "SC" = "#00A5FF", "IM" = "firebrick2",
                                         "SC2C" = "blue4", "IM2C" = "red4", "AM" = "orchid2", "A_SC" = "purple3"))+  
  # facet_wrap(~ region, labeller = labeller(region = supp.labs))+
  theme(text = element_text(size=size.t), axis.text=element_text(size=axes.size.t),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, size=axes.size.t), aspect.ratio = aspect)

# ggarrange(e,f, ncol=2)

alpha <- 
ggarrange(# pa.w,va.w,nc.w, # AIC Weights models
          pa, va, nc,
          pa.ms, va.ms, nc.ms,
          ncol=3, nrow=2, common.legend = T, legend = "none")

beta <-
ggarrange(time,e, ncol=2, nrow=2, common.legend = T, legend = "none") # c,d are throw away

gamma <- 
ggarrange(a,b, ncol=2,nrow=2, legend="none")


# model selection
jpeg("./plots/comparison/meta80_AIC_MS.jpeg", res=5e2, width=25, height=16, unit="in")
alpha
dev.off()

# model parameters
jpeg("./plots/comparison/meta80_param_axes_text.jpeg", res=5e2, width=15, height=10, unit="in")
beta
dev.off()

# SI model parameters
jpeg("./plots/comparison/meta80_param_SI_axes_text.jpeg", res=5e2, width=15, height=10, unit="in")
gamma
dev.off()
