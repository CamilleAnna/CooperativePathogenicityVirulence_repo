local_project_dir="/Users/s1687811/Documents/PhD/Research/"
setwd(local_project_dir)
setwd('./CooperativePathogenicityVirulence/')
source('./scripts/sourced_packages.R')
#source('./scripts/sourced_ggthemes.R')

library(ape)
library(MCMCglmm)


# CFR DATASET ----
d.cfr<- read.table('./data/processed_tables/assembled_CFR_SPECIES_annotation.txt', header=TRUE,
                   sep = '\t',
                   colClasses = c('character', 'character', 'character',
                                  'numeric', 'numeric',
                                  'character', 'character', 'character',
                                  'numeric', 'numeric',
                                  rep('character', 4), rep('numeric', 8)))



d.cfr<- rename(d.cfr,
               supfam_id = matching_supfam_id,
               nb_cds = total_cds)
# filter out:
# missing CFR (3)
#d.cfr<- d.cfr[!is.na(d.cfr$case_fatality_rate),] 

# not bacteria (5) or not found in supfam (3) 
#not.b<- c('Cryptosporidium parvum', 'Entamoeba histolytica', 'Giardia lamblia', 'Histoplasma capsulatum', 'Plasmodium falciparum')
#not.in.SUPFAM<- c('Escherichia coli EIEC', 'Escherichia coli EPEC', 'Plesiomonas shigelloides')
#d.cfr<- d.cfr[!d.cfr$species_leggett %in% c(not.b, not.in.SUPFAM),]

# not p or n gram profile (2)
#d.cfr<- d.cfr[!is.na(d.cfr$gram_profile),]


# This reproduces Figure 3a of Leggett et al
# So the CFR in the data must indeed be a PERCENTAGE
# Also confirms the data I have scrapped for infective dose is correct
# ggplot(d.cfr, aes(x = log(infective_dose), y = case_fatality_rate, col = infection_route))+
#   geom_point()+
#   scale_color_manual(values = c('lightgrey', 'darkgrey', 'black'))
# 
# d.cfr<- d.cfr %>%
#   filter(!is.na(case_fatality_rate), # some species do not have CFR data in their SI table
#          !is.na(supfam_id)) %>%   # to keep only legget species that I did find in SUPFAM / processed with psortb and panzzer
#   mutate(generations_per_week = (24*7)/generation_time_h) # include growth  rate in mode, which is in generation per weeks in their original model ... ? Actually I'll keep it as a generation tim in hours, as they give it
# 
# d.cfr$species<- gsub(' ', '_', d.cfr$species)



# PHYLOGENY
tree<- read.tree('./data/phylogeny_files/supfam_cfr_tree.newick')
tree$tip.label
tree<- drop.tip(tree, tree$tip.label[!tree$tip.label %in% d.cfr$supfam_id])
tree$tip.label<- d.cfr$species_id[match(tree$tip.label, d.cfr$supfam_id)]
tree<- chronopl(tree, lambda = 1)
plot(tree) # ultrametricised tree
tree<-makeNodeLabel(tree)
Ainv<-inverseA(tree, scale=FALSE)$Ainv


# TRANSFORM RESPONSE
d.cfr$log_cfr_plus1<- log(d.cfr$case_fatality_rate + 1)


# PRIOR AND RUN SPECS
prior <- list(R = list(R1 = list(V = diag(1), nu = 0.002)), # non-pylo / residual
              G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = 100))) # phylo

div = 2
nitt = 10500000/div
burnin = 500000/div
thin = ceiling(5000/div)
nitt
(nitt-burnin)/thin



# REPLICATE LEGGET ----
m3.legget <- MCMCglmm(
  log_cfr_plus1  ~ -1 + infection_route + generation_time_h,
  random = ~species_id,
  ginverse = list(species_id= Ainv), 
  data = d.cfr,
  prior = prior,
  family = c("gaussian"),
  start = list(QUASI = FALSE), 
  DIC = TRUE, nitt = nitt, thin = thin, burnin = burnin)

summary(m3.legget)


# Posteriors contrasts:

summary(m3.legget$Sol[,c('infection_routeInhalation')] - m3.legget$Sol[,c('infection_routeIngestion')])
inh_ing<- m3.legget$Sol[,c('infection_routeInhalation')] - m3.legget$Sol[,c('infection_routeIngestion')]
summary(inh_ing)
(sum(inh_ing < 0)/nrow(m3.legget$Sol))

skin_ing<- m3.legget$Sol[,c('infection_routeSkin')] - m3.legget$Sol[,c('infection_routeIngestion')]
summary(skin_ing)
(sum(skin_ing < 0)/nrow(m3.legget$Sol))


inh_skin<- m3.legget$Sol[,c('infection_routeInhalation')] - m3.legget$Sol[,c('infection_routeSkin')]
summary(inh_skin)
(sum(inh_skin < 0)/nrow(m3.legget$Sol))



# UNIVARIATE ----

# code wrapper for univariate models with transmission route and growth rate included
m3 <-function(df, focal_trait){
  
  d.tmp<- df
  colnames(d.tmp)[which(colnames(d.tmp) == focal_trait)]<- 'focal_trait'

  m<- MCMCglmm(
  log_cfr_plus1  ~ -1 + infection_route + generation_time_h + nb_cds + focal_trait,
  random = ~species_id,
  ginverse = list(species_id= Ainv), 
  data = d.tmp,
  prior = prior,
  family = c("gaussian"),
  start = list(QUASI = FALSE), 
  DIC = TRUE, nitt = nitt, thin = thin, burnin = burnin,
  verbose = FALSE)
  
  return(m)
  
}

# code wrapper for univariate models with cooperation only
m3.coop_only <-function(df, focal_trait){
  
  d.tmp<- df
  colnames(d.tmp)[which(colnames(d.tmp) == focal_trait)]<- 'focal_trait'
  
  m<- MCMCglmm(
    log_cfr_plus1  ~ -1 + nb_cds + focal_trait,
    random = ~species_id,
    ginverse = list(species_id= Ainv), 
    data = d.tmp,
    prior = prior,
    family = c("gaussian"),
    start = list(QUASI = FALSE), 
    DIC = TRUE, nitt = nitt, thin = thin, burnin = burnin,
    verbose = FALSE)
  
  return(m)
  
}


m3.bio<- m3(df = d.cfr, focal_trait = 'biofilm')
m3.qs<- m3(df = d.cfr, focal_trait = 'quorum_sensing')
m3.ab<- m3(df = d.cfr, focal_trait = 'antibiotic_degradation')
m3.sid<- m3(df = d.cfr, focal_trait = 'siderophores')
m3.ssy<- m3(df = d.cfr, focal_trait = 'secretion_system')
m3.sec<- m3(df = d.cfr, focal_trait = 'secretome')
m3.vf<- m3(df = d.cfr, focal_trait = 'is_victor_vf')

save.image("./output/model_output/cfr_2604.RData")


m3.bio.coop_onl<- m3.coop_only(df = d.cfr, focal_trait = 'biofilm')
m3.qs.coop_onl<- m3.coop_only(df = d.cfr, focal_trait = 'quorum_sensing')
m3.ab.coop_onl<- m3.coop_only(df = d.cfr, focal_trait = 'antibiotic_degradation')
m3.sid.coop_onl<- m3.coop_only(df = d.cfr, focal_trait = 'siderophores')
m3.ssy.coop_onl<- m3.coop_only(df = d.cfr, focal_trait = 'secretion_system')
m3.sec.coop_onl<- m3.coop_only(df = d.cfr, focal_trait = 'secretome')
m3.vf.coop_onl<- m3.coop_only(df = d.cfr, focal_trait = 'is_victor_vf')

save.image("./output/model_output/cfr_2604.RData")


# MULTIVARIATE MODEL ----

# z-transform predictors
z.transform<- function(x){(x - mean(x))/sd(x)}

d.cfr$z_nb_cds<- z.transform(d.cfr$nb_cds)
d.cfr$z_biofilm<- z.transform(d.cfr$biofilm)
d.cfr$z_quorum_sensing<- z.transform(d.cfr$quorum_sensing)
d.cfr$z_antibiotic_degradation<- z.transform(d.cfr$antibiotic_degradation)
d.cfr$z_secretion_system<- z.transform(d.cfr$secretion_system)
d.cfr$z_siderophores<- z.transform(d.cfr$siderophores)
d.cfr$z_secretome<- z.transform(d.cfr$secretome)
d.cfr$z_is_victor_vf<- z.transform(d.cfr$is_victor_vf)



m3.multi<- MCMCglmm(
  log_cfr_plus1  ~ -1 + infection_route + generation_time_h + z_nb_cds + z_biofilm + z_quorum_sensing + z_antibiotic_degradation + z_secretion_system + z_siderophores + z_secretome,
  random = ~species_id,
  ginverse = list(species_id= Ainv), 
  data = d.cfr,
  prior = prior,
  family = c("gaussian"),
  start = list(QUASI = FALSE), 
  DIC = TRUE, nitt = nitt, thin = thin, burnin = burnin,
  verbose = FALSE)

summary(m3.multi)



m3.multi.with_vf<- MCMCglmm(
  log_cfr_plus1  ~ -1 + infection_route + generation_time_h + z_nb_cds + z_biofilm + z_quorum_sensing + z_antibiotic_degradation + z_secretion_system + z_siderophores + z_secretome + z_is_victor_vf,
  random = ~species_id,
  ginverse = list(species_id= Ainv), 
  data = d.cfr,
  prior = prior,
  family = c("gaussian"),
  start = list(QUASI = FALSE), 
  DIC = TRUE, nitt = nitt, thin = thin, burnin = burnin,
  verbose = FALSE)



# SAVE
save.image("./output/model_output/cfr_2604.RData")


# MODEL OUTPUT ----

#load("./output/model_output/CFR_models_withOUT_interactions.RData")
load("./output/model_output/cfr_2604.RData")

summary(m3.legget)

round(summary(m3.bio)$solutions, 3)
round(summary(m3.qs)$solutions, 3)
round(summary(m3.sec)$solutions, 3)
round(summary(m3.ab)$solutions, 3)
round(summary(m3.sid)$solutions, 3)
round(summary(m3.ssy)$solutions, 3)
round(summary(m3.vf)$solutions, 3)


summary(m3.bio.coop_onl) # .
summary(m3.qs.coop_onl)
summary(m3.ab.coop_onl)
summary(m3.sid.coop_onl)
summary(m3.ssy.coop_onl)
summary(m3.sec.coop_onl)
summary(m3.vf.coop_onl)


summary(m3.multi)
summary(m3.multi.with_vf)









# CHECK: per week generation time ----

m3.legget_week <- MCMCglmm(
  log_cfr_plus1  ~ -1 + infection_route + generations_per_week,
  random = ~species,
  ginverse = list(species= Ainv), 
  data = d.cfr,
  prior = prior,
  family = c("gaussian"),
  start = list(QUASI = FALSE), 
  DIC = TRUE, nitt = nitt, thin = thin, burnin = burnin)


m3_week <-function(df, focal_trait){
  
  d.tmp<- df
  colnames(d.tmp)[which(colnames(d.tmp) == focal_trait)]<- 'focal_trait'
  
  m<- MCMCglmm(
    log_cfr_plus1  ~ -1 + infection_route +  generations_per_week + nb_cds + focal_trait,
    random = ~species,
    ginverse = list(species= Ainv), 
    data = d.tmp,
    prior = prior,
    family = c("gaussian"),
    start = list(QUASI = FALSE), 
    DIC = TRUE, nitt = nitt, thin = thin, burnin = burnin)
  
  return(m)
  
}
m3.bio_week<- m3_week(df = d.cfr, focal_trait = 'biofilm')

summary(m3.legget)
summary(m3.legget_week)
plot(m3.legget_week$VCV)
plot(m3.legget$VCV)


#save.image("./output/model_output/CFR_models_withOUT_interactions.RData")

# FIGURE ----

#load("./output/model_output/CFR_models_withOUT_interactions.RData")


cfr.theme<-
  theme_bw()+
  theme(axis.title = element_text(face = 'bold', size = 5),
        axis.text = element_text(size = 4),
        panel.grid = element_line(size = 0.1),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        legend.position = 'none',
        #legend.position = 'top', #c(0.9, 0.1),
        #legend.key.size = unit(0.3, "cm"),
        #legend.text = element_text(size = 0),
        #legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 5, face = 'bold'))


pt.size = 0.5
pt.lwd = 0.12

p.cfr.bio<- ggplot(d.cfr, aes(x = biofilm/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = c('black', 'grey', 'white'))+
  xlab(' ')+ylab(' ')+ggtitle('Biofilm')+ cfr.theme
  
  
p.cfr.sec<- ggplot(d.cfr, aes(x = secretome/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = c('black', 'grey', 'white'))+xlab(' ')+ylab(' ')+ggtitle('Secretome')+ cfr.theme

p.cfr.qs<- ggplot(d.cfr, aes(x = quorum_sensing/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = c('black', 'grey', 'white'))+xlab(' ')+ylab(' ')+ggtitle('Quorum-sensing')+ cfr.theme

p.cfr.ab<- ggplot(d.cfr, aes(x = antibiotic_degradation/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = c('black', 'grey', 'white'))+xlab(' ')+ylab(' ')+ggtitle('Antib. degr.')+ cfr.theme

p.cfr.sid<- ggplot(d.cfr, aes(x = siderophores/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = c('black', 'grey', 'white'))+xlab(' ')+ylab(' ')+ggtitle('Siderophores')+ cfr.theme

p.cfr.ssy<- ggplot(d.cfr, aes(x = secretion_system/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = c('black', 'grey', 'white'))+xlab(' ')+ylab(' ')+ggtitle('Secr. Syst.')+ cfr.theme

p.cfr.vf<- ggplot(d.cfr, aes(x = is_victor_vf/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = c('black', 'grey', 'white'))+xlab(' ')+ylab(' ')+ggtitle('Virulence factors')+ cfr.theme



library(grid)
library(gridExtra)
library(scales)

pdf('./output/figures/Figure4_CFRscatterplots.pdf', width = 4.4, height = 4.4+0.2)

grid.arrange(grobs = c(list(p.cfr.bio, p.cfr.qs,p.cfr.ab),
                       list(p.cfr.sid, p.cfr.ssy, p.cfr.sec),
                       list(p.cfr.vf)
),
ncol = 3,
as.table = FALSE,
bottom = textGrob("Proportion of genes coding for trait",gp=gpar(fontsize=6,font=1)),
left = textGrob("Case fatality rate (%)",gp=gpar(fontsize=6,font=1), rot = 90))

dev.off()






p.cfr.legend<- ggplot(d.cfr, aes(x = biofilm/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(name = 'Infection route', values = c('black', 'grey', 'white'))+
  theme_bw()+
  theme(axis.title = element_text(face = 'bold', size = 5),
        axis.text = element_text(size = 4),
        panel.grid = element_line(size = 0.1),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        #legend.position = 'none',
        #legend.position = 'top', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.title = element_text(hjust = 0.5, size = 5, face = 'bold'))


library(grid)
library(gridExtra) 
legend <- cowplot::get_legend(p.cfr.legend)
pdf('./output/figures/Fig4.main.A.legend.pdf', width = 4.33, height = (0.48)*4.33)
grid.draw(legend)
dev.off()


# params of multivariate models
m3.outs.multi<- as.data.frame(summary(m3.multi.with_vf)$solutions)[-c(1:5),] %>%
  mutate(trait = c('Biofilm', 'Quorum-sensing', 'Antib. degr.', 'Secr. Syst.', 'Siderophores', 'Secretome', 'Virulence factors')) %>%
  select(trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)

m3.outs.multi$trait<- factor(m3.outs.multi$trait,
  levels = c('Virulence factors', 'Quorum-sensing', 'Biofilm', 'Siderophores', 'Antib. degr.', 'Secr. Syst.', 'Secretome'))


p.models.multi<- ggplot(m3.outs.multi, aes(x = post.mean, y = trait, col = 'black'))+
  xlab('Posterior mean (multivariate model)')+ylab('Cooperation category')+
  geom_point()+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1, col = 'black', size = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
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



# params of univariate models
m3.outs.uni<- rbind(
as.data.frame(summary(m3.sec)$solutions)['focal_trait',] %>% mutate(trait = 'Secretome'),
as.data.frame(summary(m3.ssy)$solutions)['focal_trait',] %>% mutate(trait = 'Secr. Syst.'),
as.data.frame(summary(m3.ab)$solutions)['focal_trait',] %>% mutate(trait = 'Antib. degr.'),
as.data.frame(summary(m3.sid)$solutions)['focal_trait',] %>% mutate(trait = 'Siderophores'),
as.data.frame(summary(m3.bio)$solutions)['focal_trait',] %>% mutate(trait = 'Biofilm'),
as.data.frame(summary(m3.qs)$solutions)['focal_trait',] %>% mutate(trait = 'Quorum-sensing'),
as.data.frame(summary(m3.vf)$solutions)['focal_trait',] %>% mutate(trait = 'Virulence factors')
) %>%
  select(trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)

m3.outs.uni$trait<- factor(m3.outs.uni$trait, levels = rev(m3.outs.uni$trait))


p.models.uni<- ggplot(m3.outs.uni, aes(x = post.mean, y = trait, col = 'black'))+
  xlab('Posterior mean (univariate models)')+ylab('Cooperation category')+
  geom_point()+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1, col = 'black', size = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
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



pdf('./output/figures/Fig4.main.B.pdf', width = 3, height = (0.8)*3)
p.models.multi
dev.off()


pdf('./output/figures/Fig4.supp.pdf', width = 3, height = (0.8)*3)
p.models.uni
dev.off()


grid.arrange(p.models.multi, p.models.uni, ncol = 2)


# SUPP. TABLE (model summary) ----

get.effects.M3<- function(model.output, model = ''){
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

summary.tab.cfr<- 
  rbind(
    get.effects.M3(m3.legget, 'Leggett'),
    get.effects.M3(m3.bio, 'biofilm'),
      get.effects.M3(m3.qs, 'quorum-sensing'),
      get.effects.M3(m3.sec, 'secretome'),
      get.effects.M3(m3.ab, 'antibiotic-resistance'),
      get.effects.M3(m3.sid, 'siderophores'),
      get.effects.M3(m3.ssy, 'secretion-systems'),
      get.effects.M3(m3.vf, 'virulence-factors'),
      get.effects.M3(m3.multi.with_vf, "multivariate-model"))


write.table(summary.tab.cfr, './output/summary_tables_raw/model_cfr.csv', sep = '\t', col.names = TRUE,quote = FALSE,  row.names = FALSE)




