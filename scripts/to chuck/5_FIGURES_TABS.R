



# OUTPUT RESULTS ----

load("./output/model_output/vf_2304_VICTOR_d118_CHAIN.RData")

# Model summaries
summary(meta.a.secretome)
summary(meta.a.biofilm)
summary(meta.a.sid)
summary(meta.a.ab)
summary(meta.a.qs)
summary(meta.a.ssyst)
summary(meta.a.coop) # all classes together

plot(meta.a.secretome)
plot(meta.a.biofilm)
plot(meta.a.coop)


# Getting rounded coeffs, back to OR scale
round(exp(summary(meta.a.secretome)$solutions), 2)
round(exp(summary(meta.a.biofilm)$solutions), 2)
round(exp(summary(meta.a.sid)$solutions), 2)
round(exp(summary(meta.a.ab)$solutions), 2)
round(exp(summary(meta.a.qs)$solutions), 2)
round(exp(summary(meta.a.ssyst)$solutions), 2)
round(exp(summary(meta.a.coop)$solutions), 2)  # all classes together



# FIGURES ----

# ... OR meta-analysis ----

meta.results<- rbind(
  data.frame(summary(meta.a.secretome)$solutions, trait = 'Secretome'),
  data.frame(summary(meta.a.biofilm)$solutions, trait = 'Biofilm'),
  data.frame(summary(meta.a.sid)$solutions, trait = 'Siderophores'),
  data.frame(summary(meta.a.ab)$solutions, trait = 'Antib. degr.'),
  data.frame(summary(meta.a.qs)$solutions, trait = 'Quorum-sensing'),
  data.frame(summary(meta.a.ssyst)$solutions, trait = 'Secr. syst.'),
  data.frame(summary(meta.a.coop)$solutions, trait = 'Cooperation (any)'))


meta.results$trait<- factor(meta.results$trait,
                            levels = c('Quorum-sensing',
                                       'Biofilm',
                                       'Siderophores',
                                       'Antib. degr.',
                                       'Secr. syst.',
                                       'Secretome',
                                       'Cooperation (any)'))

# turn it back to normal scale
meta.results$post.mean = exp(meta.results$post.mean)
meta.results$l.95..CI = exp(meta.results$l.95..CI)
meta.results$u.95..CI = exp(meta.results$u.95..CI)



p.vf.meta<- ggplot(meta.results, aes(x = `post.mean`, y = trait, col = 'black'))+
  xlab('Odds ratio (posterior distribution)')+ylab('Cooperation category')+
  geom_point()+
  #xlim(-1,5)+
  geom_errorbarh(aes(xmin = `l.95..CI`, xmax = `u.95..CI`), height = 0.1)+
  geom_vline(xintercept = 1, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l.95..CI`, xmax = `u.95..CI`), height = 0.1, col = 'black', size = 0.3)+
  scale_color_manual(values = 'white')+
  theme_bw()+
  theme(axis.title = element_text(face = 'bold', size = 5),
        axis.text = element_text(size = 5),
        panel.grid = element_line(size = 0.1),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 0),
        legend.title = element_blank())


# ... OR overall  ----

getor<- function(trait){
  data.frame(
    trait = trait, 
    or = ctab.all(trait)$or$es,
    lci = ctab.all(trait)$or$ci.lo,
    uci = ctab.all(trait)$or$ci.hi,
    n = ctab.all(trait)$or$totaln,
    chi = ctab.all(trait)$chi$statistic,
    p = ctab.all(trait)$chi$p.value)
}

or.all<- rbind(
  getor('secretome'),
  getor('biofilm'),
  getor('siderophores'),
  getor('quorum_sensing'),
  getor('antibiotic_degradation'),
  getor('secretion_system'),
  getor('is_coop'))

or.all<- mutate(or.all,
                trait2 = c('Secretome', 'Biofilm', 'Siderophores', 'Quorum-sensing', 'Antib. degr.', 'Secr. syst.', 'Cooperation (any)'))


or.all$trait2<- factor(or.all$trait2,
                       levels = c('Quorum-sensing',
                                  'Biofilm',
                                  'Siderophores',
                                  'Antib. degr.',
                                  'Secr. syst.',
                                  'Secretome',
                                  'Cooperation (any)'))


p.vf.or<- ggplot(or.all, aes(x = or, y = trait2, col = 'black'))+
  xlab('Odds ratio')+ylab('Cooperation category')+
  geom_point()+
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.1)+
  geom_vline(xintercept = 1, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.1, col = 'black', size = 0.3)+
  scale_color_manual(values = 'white')+
  theme_bw()+
  theme(axis.title = element_text(face = 'bold', size = 5),
        axis.text = element_text(size = 5),
        panel.grid = element_line(size = 0.1),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 0),
        legend.title = element_blank())


# ... Coop (any) funnel  ----


head(df.all)


p.vf.funnel.is_coop<- 
  df.all %>%
  filter(cooperative_trait == 'is_coop') %>%
  ggplot(., aes(x = logOR, y = 1/se_logOR, col = 'black'))+
  ggtitle('Cooperation (any)')+
  xlab('Species-specific logOR')+
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  geom_rect(ymin = 0, ymax = 5,
            xmin = summary(meta.a.coop)$solutions[,2],
            xmax = summary(meta.a.coop)$solutions[,3], fill = 'gray89')+
  geom_point(size = 0.5)+
  geom_point(col = 'black', size = 0.5)+
  scale_color_manual(values = 'white')+
  theme_bw()+
  theme(axis.title = element_text(face = 'bold', size = 5),
        axis.text = element_text(size = 5),
        panel.grid = element_line(size = 0.1),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 0),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 6))



# ... paper figure ----

library(patchwork)

pdf('./output/figures/Fig2.main.pdf', width = 4.8, height = 0.55*4.33)

{(p.vf.funnel.is_coop | p.vf.meta)} + plot_layout(ncol = 2, widths = c(0.5, 0.5))+
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 7, face = 'bold'))

dev.off()



# suppelemt figs ----


plot.funnel<- function(df, trait, title, model){
  
  df %>%
    filter(cooperative_trait == trait) %>%
    ggplot(., aes(x = logOR, y = 1/se_logOR, col = 'black'))+
    xlab('Species-specific logOR')+
    ggtitle(title)+
    geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
    geom_rect(ymin = 0, ymax = 5,
              xmin = summary(model)$solutions[,2],
              xmax = summary(model)$solutions[,3], fill = 'gray89')+
    geom_point(size = 0.5)+
    geom_point(col = 'black', size = 0.5)+
    scale_color_manual(values = 'white')+
    theme_bw()+
    theme(axis.title = element_text(face = 'bold', size = 5),
          axis.text = element_text(size = 5),
          panel.grid = element_line(size = 0.1),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", size = 0.2),
          legend.position = 'none', #c(0.9, 0.1),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 0),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = 'bold', size = 6))
} 


pdf('./output/figures/Fig2.supp1.pdf', width = 4.33, height = 4.33)

grid.arrange(ncol = 3,
             plot.funnel(df.all, 'is_coop', 'Cooperative (any)', meta.a.coop),
             plot.funnel(df.all, 'secretome', 'Secretome', meta.a.secretome),
             plot.funnel(df.all, 'qs', 'Quorum-sensing', meta.a.qs),
             plot.funnel(df.all, 'sid', 'Siderophores', meta.a.sid),
             plot.funnel(df.all, 'ssyst', 'Secretion-systems', meta.a.ssyst),
             plot.funnel(df.all, 'biofilm', 'Biofilm', meta.a.biofilm),
             plot.funnel(df.all, 'ab', 'Antibiotic degradation', meta.a.ab))

dev.off()



# SUPP. TABLE (model summary) ----

get.effects.M2<- function(model.output, model = ''){
  extract.d<- rbind(
    rbind(summary(model.output)$Gcovariances, summary(model.output)$Rcovariances) %>% as.data.frame() %>% mutate(effect = rownames(.)) %>% mutate(pMCMC = NA),
    summary(model.output)$solutions %>% as.data.frame() %>% mutate(effect = rownames(.)))
  extract.d<- extract.d %>%
    select(effect, post.mean, `l-95% CI`, `u-95% CI`, eff.samp, pMCMC) %>%
    mutate(model_name = model)
  
  
  extract.d$post.mean<- ifelse(round(extract.d$post.mean, 3) > 0,
                               round(extract.d$post.mean, 3),
                               format(extract.d$post.mean, scientific = TRUE, big.mark = ",", digit = 2))
  
  
  extract.d$`l-95% CI`<- ifelse(round(extract.d$`l-95% CI`, 3) > 0,
                                round(extract.d$`l-95% CI`, 3),
                                format(extract.d$`l-95% CI`, scientific = TRUE, big.mark = ",", digit = 2))
  
  extract.d$`u-95% CI`<- ifelse(round(extract.d$`u-95% CI`, 3) > 0,
                                round(extract.d$`u-95% CI`, 3),
                                format(extract.d$`u-95% CI`, scientific = TRUE, big.mark = ",", digit = 2))
  
  extract.d$eff.samp<- round(extract.d$eff.samp, 3)
  
  
  extract.d$pMCMC<- ifelse(round(extract.d$pMCMC, 3) > 0,
                           round(extract.d$pMCMC, 3),
                           format(extract.d$pMCMC, scientific = TRUE, big.mark = ",", digit = 2))
  
  
  return(extract.d)
}


summary.tab.vf<- rbind(
  get.effects.M2(meta.a.secretome, 'secretome'),
  get.effects.M2(meta.a.biofilm, 'biofilm'),
  get.effects.M2(meta.a.sid, 'siderophores'),
  get.effects.M2(meta.a.ab, 'antibiotic degradation'),
  get.effects.M2(meta.a.qs, 'quorum-sensing'),
  get.effects.M2(meta.a.ssyst, 'secretion systems'),
  get.effects.M2(meta.a.coop, 'cooperation (any)')
)


write.table(summary.tab.vf, './output/summary_tables_raw/model_virulenceFactors.csv', sep = '\t', col.names = TRUE,quote = FALSE,  row.names = FALSE)


summary.tab.vf %>%
  filter(effect == '(Intercept)')


meta.results




