local_project_dir="/Users/s1687811/Documents/PhD/Research/"
setwd(local_project_dir)
setwd('./CooperativePathogenicityVirulence//')
source('./scripts/sourced_packages.R')
#source('./scripts/sourced_ggthemes.R')

library(ape)
library(MCMCglmm)



# LOAD DATA ----
d<- read.table('./data/processed_tables/assembled_SPECIES_annotation.txt', header=TRUE, sep = '\t',
               colClasses = c(rep('character', 4),
                              rep('numeric', 9)))

d<- rename(d,
           nb_extracellular = secretome,
           species.patric = species,
           vf = is_victor_vf,
           #mp3_pathogenic = is_mp3_hybrid_pathogenic,
           nb_cds = total_cds,
           ab_degradation = antibiotic_degradation) %>%
  rename(species = species_id)

d$pathogen<- as.factor(d$pathogen)

# Sanity check
nrow(d) # 118
length(unique(d$species))
table(d$pathogen) # we have: 59 non-pathogen, 9 pathogen (species that were both already reclassified as pathogen)


#d<- d[d$gram_profile %in% c('p', 'n'),] # after this filtering, we have 90 species, 39 non-pathogen, 51 pathogen


ggplot(d, aes(x = pathogen, y = nb_extracellular, col = gram_profile))+
  geom_boxplot() # bear in mind there are just 4 OM+, three of which are Mycobacterium. So likely OM+ just is part of same distribution as p

ggplot(d, aes(x = pathogen, y = nb_extracellular, col = gram_profile))+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2))

ggplot(d, aes(x = as.factor(pathogen), y = biofilm))+
  geom_boxplot()


# PHYLOGENY ----
#tree<- read.tree('./data/species_info_files/hm_phylogeny.newick')
tree<- read.tree('./data/phylogeny_files/midas_tree_renamed.newick')
tree$tip.label
tree<- drop.tip(tree, tree$tip.label[!tree$tip.label %in% d$species])

phylogeny<- chronopl(tree, lambda = 1)
plot(phylogeny) # ultrametricised tree
phylogeny<-makeNodeLabel(phylogeny)
Ainv<-inverseA(phylogeny, scale=FALSE)$Ainv


# M1 (UNI): pathogen ~ n.o. cooperative genes ----

# PRIOR:
# R: fixed residual variance because binary trait
# G1: random effect  for  phylogeny,
# One observation per species, so  non-phylogenetic species variance is  residual, which is fixed

prior <- list(R = list(V = 1, fix = 1),
              G = list(G1 = list(V = diag(1), nu = 1000, alpha.mu = rep(0,1), alpha.V = diag(1))))


div = 2
nitt = 10500000/div
burnin = 500000/div
thin = ceiling(5000/div)
nitt
(nitt-burnin)/thin

d<- as.data.frame(d)

# Code wrapper to run models on ech trait (uses env variables for data and nitt/burnin/thin)
run_m1.go<- function(df, focal_trait){
  
  
  colnames(df)[which(colnames(df) == focal_trait)]<- 'focal_trait'
  
  m1.go<- MCMCglmm(pathogen  ~ 1 + nb_cds + focal_trait ,
                   random = ~species,
                   ginverse = list(species=Ainv),
                   data = df,
                   prior = prior,
                   family=c("categorical"), trunc=T,
                   start=list(QUASI=FALSE),
                   #pl = TRUE, pr = TRUE, nodes = 'ALL',
                   DIC = TRUE,
                   verbose = FALSE,
                   nitt=nitt, thin=thin, burnin=burnin)
  
  return(m1.go)
}
run_m1.go.gram<- function(df, focal_trait){
  
  
  colnames(df)[which(colnames(df) == focal_trait)]<- 'focal_trait'
  
  m1.go<- MCMCglmm(pathogen  ~ 1 + nb_cds + gram_profile + focal_trait ,
                   random = ~species,
                   ginverse = list(species=Ainv),
                   data = df,
                   prior = prior,
                   family=c("categorical"), trunc=T,
                   start=list(QUASI=FALSE),
                   #pl = TRUE, pr = TRUE, nodes = 'ALL',
                   DIC = TRUE,
                   nitt=nitt, thin=thin, burnin=burnin)
  
  return(m1.go)
}


m1.ss<- run_m1.go(d, focal_trait = 'nb_extracellular')
m1.biofilm<- run_m1.go(d, 'biofilm')
m1.siderophores<- run_m1.go(d,'siderophores')
m1.ab_degradation<- run_m1.go(d,'ab_degradation')
m1.secretion_system<- run_m1.go(d,'secretion_system')
m1.quorum_sensing<- run_m1.go(d,'quorum_sensing')
m1.vf<- run_m1.go(d, 'vf')
#m1.mp3<- run_m1.go(d[!is.na(d$mp3_pathogenic),], 'mp3_pathogenic')


summary(m1.biofilm)          # 
summary(m1.siderophores)     # . in multi
summary(m1.ab_degradation)   # **
summary(m1.quorum_sensing)   # 
summary(m1.secretion_system) # ***
summary(m1.ss) # **
summary(m1.vf) # *** 
#summary(m1.mp3) #



# MULTIVARIATE MODEL ----

# z-transform predictors
z.transform<- function(x){(x - mean(x))/sd(x)}

d$z_nb_cds<- z.transform(d$nb_cds)
d$z_biofilm<- z.transform(d$biofilm)
d$z_quorum_sensing<- z.transform(d$quorum_sensing)
d$z_ab_degradation<- z.transform(d$ab_degradation)
d$z_secretion_system<- z.transform(d$secretion_system)
d$z_siderophores<- z.transform(d$siderophores)
d$z_nb_extracellular<- z.transform(d$nb_extracellular)
d$z_vf<- z.transform(d$vf)
#d$z_mp3<- z.transform(d$mp3_pathogenic)


# DEPRECATED - Actually, standardize the nb_extracellular respective to the gram_profile, so that we can include it in the multivariate model without having gram profile along, which would not be releavant for the other predictors
#d$z_nb_extracellular_gram_spec<- NA

#d$z_nb_extracellular_gram_spec[which(d$gram_profile == 'p')]<- 
#(d$nb_extracellular[d$gram_profile == 'p'] - mean(d$nb_extracellular[d$gram_profile == 'p']))/sd(d$nb_extracellular[d$gram_profile == 'p'])

#d$z_nb_extracellular_gram_spec[which(d$gram_profile == 'n')]<- 
#  (d$nb_extracellular[d$gram_profile == 'n'] - mean(d$nb_extracellular[d$gram_profile == 'n']))/sd(d$nb_extracellular[d$gram_profile == 'n'])


# m1.multi.with_vf_WITHgram<- MCMCglmm(pathogen  ~ 1 + gram_profile + z_nb_cds +  z_biofilm + z_ab_degradation + z_quorum_sensing + z_secretion_system + z_siderophores + z_nb_extracellular + z_vf,
#                                    random = ~species,
#                                    ginverse = list(species=Ainv),
#                                    data = d,
#                                    prior = prior,
#                                    family=c("categorical"),trunc=T,
#                                    start=list(QUASI=FALSE),
#                                    #pl = TRUE, pr = TRUE, nodes = 'ALL',
#                                    DIC = TRUE,
#                                    nitt=nitt, thin=thin, burnin=burnin)



m1.multi.with_vf_NOgram<- MCMCglmm(pathogen  ~ 1 + z_nb_cds +  z_biofilm + z_ab_degradation + z_quorum_sensing + z_secretion_system + z_siderophores + z_nb_extracellular + z_vf,
                    random = ~species,
                    ginverse = list(species=Ainv),
                    data = d,
                    prior = prior,
                    family=c("categorical"),trunc=T,
                    start=list(QUASI=FALSE),
                    #pl = TRUE, pr = TRUE, nodes = 'ALL',
                    DIC = TRUE,
                    nitt=nitt, thin=thin, burnin=burnin)




# m1.multi.with_MP3_NOgram<- MCMCglmm(pathogen  ~ 1 + z_nb_cds +  z_biofilm + z_ab_degradation + z_quorum_sensing + z_secretion_system + z_siderophores + z_nb_extracellular + z_mp3,
#                                    random = ~species,
#                                    ginverse = list(species=Ainv),
#                                    data = d,
#                                    prior = prior,
#                                    family=c("categorical"),trunc=T,
#                                    start=list(QUASI=FALSE),
#                                    #pl = TRUE, pr = TRUE, nodes = 'ALL',
#                                    DIC = TRUE,
#                                    nitt=nitt, thin=thin, burnin=burnin)


#save.image('./output/model_output/model1_reps_genomes_v2.RData')
#save.image('./output/model_output/model_output_pathogens.RData')

#save.image('./output/model_output/pathogens_2304_VICTOR_d118.RData')


#summary(m1.multi.with_vf_WITHgram)
summary(m1.multi.with_vf_NOgram) # results are the same, except marginal significance for AB if gram not included
#summary(m1.multi.with_MP3_NOgram) # results are the same, except marginal significance for AB if gram not included


# CHECK MODEL OUTPUT ----

#load('./output/model_output/model1_reps_genomes_v2.RData')
load('./output/model_output/pathogens_2304_VICTOR_d118.RData')


summary(m1.biofilm)          # 
summary(m1.siderophores)     # (*, neg, in multivariate model)
summary(m1.ab_degradation)   # 
summary(m1.quorum_sensing)   # 
summary(m1.secretion_system) # *
summary(m1.ss) #  *
summary(m1.vf) #  ***
#summary(m1.mp3) #


summary(m1.multi.with_vf_NOgram)
#summary(m1.multi.with_MP3_NOgram)
plot(m1.multi.with_vf_NOgram)
plot(m1.ss$Sol)
plot(m1.ss$VCV)
summary(m1.ss)


get.effects.M1<- function(model.output, model = ''){
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






# CHECK EMBLEMATIC ----
d$cfr_cat<- NA
d$cfr_cat[which(d$species %in% c('Yersinia_pestis_56773', 'Rickettsia_prowazekii_51996', 'Burkholderia_pseudomallei_54177', 'Listeria_grayi_55743'))]<- 'top'
d$cfr_cat[which(d$species %in% c('Francisella_tularensis_57093', 'Vibrio_parahaemolyticus_57034', 'Escherichia_coli_58110', 'Campylobacter_jejuni_62119'))]<- 'mid'

d$cfr_cat[which(d$species %in% c('Ruminococcus_gnavus_57638', 'Prevotella_melaninogenica_58075', 'Veillonella_parvula_57794', 'Bacteroides_vulgatus_57955'))]<- 'a'

d.test<- d[!is.na(d$cfr_cat),]

d.test.long<- d.test %>%
  gather('trait', 'z_value', 19:24)

p.test.z<- ggplot(d.test.long, aes(x = z_value, y = trait, fill = cfr_cat))+
  geom_boxplot(width = 0.5)#+
#facet_wrap(~trait)




# FIGURE 3 ----


# params of univariate models
m1.outs.uni<- rbind(
as.data.frame(summary(m1.ss)$solutions)['focal_trait',] %>% mutate(trait = 'Secretome'),
as.data.frame(summary(m1.secretion_system)$solutions)['focal_trait',] %>% mutate(trait = 'Secr. Syst.'),
as.data.frame(summary(m1.ab_degradation)$solutions)['focal_trait',] %>% mutate(trait = 'Antib. degr.'),
as.data.frame(summary(m1.siderophores)$solutions)['focal_trait',] %>% mutate(trait = 'Siderophores'),
as.data.frame(summary(m1.biofilm)$solutions)['focal_trait',] %>% mutate(trait = 'Biofilm'),
as.data.frame(summary(m1.quorum_sensing)$solutions)['focal_trait',] %>% mutate(trait = 'Quorum-sensing'),
as.data.frame(summary(m1.vf)$solutions)['focal_trait',] %>% mutate(trait = 'Virulence factors')
#as.data.frame(summary(m1.mp3)$solutions)['focal_trait',] %>% mutate(trait = 'MP3')
) %>%
  select(trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)

m1.outs.uni$trait<- factor(m1.outs.uni$trait, levels = rev(m1.outs.uni$trait))

p.models.uni<- ggplot(m1.outs.uni, aes(x = post.mean, y = trait, col = 'black'))+
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




# params of multivariate models
m1.outs.multi<- as.data.frame(summary(m1.multi.with_vf_NOgram)$solutions)[-c(1,2),] %>%
  mutate(trait = c('Biofilm', 'Antib. degr.', 'Quorum-sensing', 'Secr. Syst.', 'Siderophores', 'Secretome', 'Virulence factors')) %>%
  #select(model, trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)
   select(trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)


m1.outs.multi$trait<- factor(m1.outs.multi$trait, levels = rev(m1.outs.uni$trait)) # same order as **UNI** table


p.models.multi<- ggplot(m1.outs.multi, aes(x = post.mean, y = trait, col = 'black'))+
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
        legend.position = 'bottom', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 0),
        legend.title = element_blank())


# boxplots
d.long<- d %>% gather('trait', 'z_value', grep('z_', colnames(d)))
d.long.raw<- d %>% gather('trait', 'value', gsub('z_', '', colnames(d)[grep('z_', colnames(d))]))

d.long$trait<- gsub('z_biofilm', 'Biofilm', d.long$trait)
d.long$trait<- gsub('z_quorum_sensing', 'Quorum-sensing', d.long$trait)
d.long$trait<- gsub('z_ab_degradation', 'Antib. degr.', d.long$trait)
d.long$trait<- gsub('z_secretion_system', 'Secr. Syst.', d.long$trait)
d.long$trait<- gsub('z_siderophores', 'Siderophores', d.long$trait)
d.long$trait<- gsub('z_nb_extracellular', 'Secretome', d.long$trait)
d.long$trait<- gsub('z_vf', 'Virulence factors', d.long$trait)
#d.long$trait<- gsub('z_mp3', 'MP3', d.long$trait)

d.long<- d.long[d.long$trait != 'z_nb_cds',]
d.long$trait<- factor(d.long$trait, levels = rev(m1.outs.uni$trait)) # order as models *UNI* table


d.long.raw$trait<- gsub('biofilm', 'Biofilm', d.long.raw$trait)
d.long.raw$trait<- gsub('quorum_sensing', 'Quorum-sensing', d.long.raw$trait)
d.long.raw$trait<- gsub('ab_degradation', 'Antib. degr.', d.long.raw$trait)
d.long.raw$trait<- gsub('secretion_system', 'Secr. Syst.', d.long.raw$trait)
d.long.raw$trait<- gsub('siderophores', 'Siderophores', d.long.raw$trait)
d.long.raw$trait<- gsub('nb_extracellular', 'Secretome', d.long.raw$trait)
d.long.raw$trait<- gsub('vf', 'Virulence factors', d.long.raw$trait)
#d.long.raw$trait<- gsub('mp3', 'MP3', d.long.raw$trait)

d.long.raw<- d.long.raw[d.long.raw$trait != 'nb_cds',]
d.long.raw$trait<- factor(d.long.raw$trait, levels = rev(m1.outs.uni$trait)) # order as models *UNI* table


p.z<- ggplot(d.long[d.long$trait != 'z_nb_cds',], aes(x = z_value, y = trait, fill = pathogen))+
  xlab('Z-number of genes')+ylab('Cooperation category')+
  geom_boxplot(width = 0.5, outlier.size = 0.05, lwd = 0.2)+
  theme_bw()+
  scale_fill_manual(values = c('dodgerblue', 'darkorange'), labels = c('Non-pathogen', 'Pathogen'))+
  theme(axis.title = element_text(face = 'bold', size = 5),
        axis.text = element_text(size = 5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        panel.grid = element_line(size = 0.1),
        legend.position = 'bottom', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank())


p.raw<- ggplot(d.long.raw[d.long.raw$trait != 'nb_cds',], aes(x = value, y = trait, fill = pathogen))+
  xlab('Number of genes')+ylab('Cooperation category')+
  geom_boxplot(width = 0.5, outlier.size = 0.05, lwd = 0.2)+
  theme_bw()+
  scale_fill_manual(values = c('dodgerblue', 'darkorange'), labels = c('Non-pathogen', 'Pathogen'))+
  theme(axis.title = element_text(face = 'bold', size = 5),
        axis.text = element_text(size = 5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        panel.grid = element_line(size = 0.1),
        legend.position = 'right', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank())

#library(gridExtra)


#ft<- 'Secretome'

make.boxplot.pathogen<- function(focal.trait){
d.long.raw %>%
  filter(trait == focal.trait) %>%
  mutate(pathogen = ifelse(pathogen == "1", 'P', 'NP')) %>%
  ggplot(., aes(x = value, y = pathogen, fill = pathogen))+
  ggtitle(focal.trait)+
  xlab('Number of genes')+ylab(' ')+
  geom_boxplot(width = 0.5, outlier.size = 0.05, lwd = 0.2)+
  theme_bw()+
  scale_fill_manual(values = c('dodgerblue', 'darkorange'), labels = c('Non-pathogen', 'Pathogen'))+
  theme(axis.title = element_text(size = 5),
        axis.text = element_text(size = 5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        panel.grid = element_line(size = 0.1),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 5, face = 'bold'))
}



p1<- make.boxplot.pathogen('Secretome')
p2<- make.boxplot.pathogen('Secr. Syst.')
p3<- make.boxplot.pathogen('Antib. degr.')
p4<- make.boxplot.pathogen('Siderophores')
p5<- make.boxplot.pathogen('Biofilm')
p6<- make.boxplot.pathogen('Quorum-sensing')
p7<- make.boxplot.pathogen('Virulence factors')
#p8<- make.boxplot.pathogen('MP3')


#library(devtools)
#install_github("thomasp85/patchwork")
#library(patchwork)



pdf('./output/figures/Fig3.supp1.pdf', width = 4.33, height = 4.33)
#grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 3)
grid.arrange(p1,p2,p3,p4,p5,p6,p7, ncol = 3)
dev.off()

pdf('./output/figures/Fig3.supp2.pdf', width = 3, height = (0.8)*3)
p.models.uni
dev.off()


library(grid)
library(gridExtra) 
legend <- cowplot::get_legend(p.z)
pdf('./output/figures/Fig3.legend.pdf', width = 4.33, height = (0.48)*4.33)
grid.draw(legend)
dev.off()


library(patchwork)
pdf('./output/figures/Fig3.main.pdf', width = 4.33, height = (0.48)*4.33)
(p.z + theme(legend.position = 'none') | p.models.multi + theme(legend.position = 'none', axis.text.y = element_blank(), axis.title.y = element_blank())) + plot_annotation(tag_levels = 'A') & 
   theme(plot.tag = element_text(size = 7, face = 'bold'))
dev.off()





# SUPP. TABLE (model summary) ----



summary.tab.pathogenicity<- rbind(
  get.effects.M1(m1.ss, 'secretome'),
  get.effects.M1(m1.biofilm, 'biofilm'),
  get.effects.M1(m1.siderophores, 'siderophores'),
  get.effects.M1(m1.ab_degradation, 'antibiotic degradation'),
  get.effects.M1(m1.quorum_sensing, 'quorum-sensing'),
  get.effects.M1(m1.secretion_system, 'secretion systems'),
  get.effects.M1(m1.vf, 'virulence factors'),
  get.effects.M1(m1.multi.with_vf_NOgram, 'Multivariate model')
)


write.table(summary.tab.pathogenicity, './output/summary_tables_raw/model_pathogenicity.csv', sep = '\t', col.names = TRUE,quote = FALSE,  row.names = FALSE)


