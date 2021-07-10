# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```~~~~~~~~~~~ #
#                   Script for comparative analysis of pathogenicity
#
#   Takes as input:
#         assembled_SPECIES_annotation.txt
#         midas_tree_renamed.newick
#   Runs analysis:
#         Takes the species Annotation
#         Run a binary response phylogenetic analysis to test whether each cooperative trait predicts pathogenicty
#         Run both univariate models (one for each trait)
#         And a multivariate model with all traits, z-transformed, included as predictors
#   Output:
#         3_model_output/CompAnalysis_pathogens_d118_CHAIN1.RData
#         3_model_output/CompAnalysis_pathogens_d118_CHAIN2.RData
#         3_model_output/CompAnalysis_pathogens_d118_CHAIN3.RData
#
# Code for extracting tables and producing figures in separate script file, using RData object
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#local_project_dir=path/to/cloned/repo #'/Users/s1687811/Documents/PhD/Research/CooperativePathogenicityVirulence/'

setwd(local_project_dir)
setwd('./CooperativePathogenicityVirulence_repo/')
source('./scripts/sourced_packages.R')
library(ape)
library(MCMCglmm)
library(dmetar)
library(esc)


# 0 - LOAD DATA  ----


d<- read.table('./output/1_processed_tables/2.1_assembled_SPECIES_annotation.txt', header=TRUE, sep = '\t',
               colClasses = c(rep('character', 4),
                              rep('numeric', 9)))

d<- rename(d,
           nb_extracellular = secretome,
           species.patric = species,
           vf = is_victor_vf,
           nb_cds = total_cds,
           ab_degradation = antibiotic_degradation) %>%
  rename(species = species_id)

d$pathogen<- as.factor(d$pathogen)

# SANITY CHECK
nrow(d) # 118
length(unique(d$species))
table(d$pathogen) # we have: 59 non-pathogen, 9 pathogen (species that were both already reclassified as pathogen)
unique(d$gram_profile) # Three gram profiles types



# QUICK OVERVIEW
ggplot(d, aes(x = pathogen, y = nb_extracellular, col = gram_profile))+
  geom_boxplot() # Only  4 OM+, three of which are Mycobacterium. Likely OM+ just is part of same distribution as p ...?

ggplot(d, aes(x = pathogen, y = nb_extracellular, col = gram_profile))+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2))

ggplot(d, aes(x = as.factor(pathogen), y = biofilm))+
  geom_boxplot()


ggplot(d, aes(x = pathogen, y = nb_cds, col = gram_profile))+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2))


# Thre gram p clearly have more nb_extracellular
# But they do not have more genes
# Not sure if they are more likely to be pathogens (as in, do I have more gram P pathogens than gram N pathogens?)
# If I have more gram P than gram N pathogens, and gram profile associated with more nb_extracellular
# could create a spurious correlation between pathogenicity and secretions, just as a result of gram profile covariate
# But the gram profile being unrelated to pathogenicity

# --> CHECK THAT: in fact I have more grame negative pathogens
# And equal number of gram P and N for non-pathogens
table(d$gram_profile, d$pathogen)

# --> run models with and without grame profile included as covariate
# found that results are the same, the same effects come out significant
# significance tends to be higher when including grame profile
# although DIC does not improve
# left code in script for records, but not in manuscript to focus on main message


# PHYLOGENY ----
tree<- read.tree('./data/5_phylogeny_files/midas_tree_renamed.newick')
tree$tip.label
tree<- drop.tip(tree, tree$tip.label[!tree$tip.label %in% d$species])

phylogeny<- chronopl(tree, lambda = 1)
#plot(phylogeny) # ultrametricised tree
phylogeny<-makeNodeLabel(phylogeny)
Ainv<-inverseA(phylogeny, scale=FALSE)$Ainv


# M1 (UNI): pathogen ~ n.o. cooperative genes ----

# PRIOR:
# R: fixed residual variance (not identifiable with binary response)
# G1: random effect  for  phylogeny
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
                   verbose = FALSE,
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


summary(m1.biofilm)          # DIC = 83
summary(m1.siderophores)     # DIC = 83
summary(m1.ab_degradation)   # DIC = 83
summary(m1.quorum_sensing)   # DIC = 83
summary(m1.secretion_system) # *, DIC = 85
summary(m1.ss) #  *, DIC = 79
summary(m1.vf) #  ***, DIC = 76



# DEPRECATED ~~~~~~~~  keep for records but remove from final analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CHECK: what does it look like if I include gram profile?

# m1.ss.gram<- run_m1.go.gram(d, focal_trait = 'nb_extracellular')
# m1.biofilm.gram<- run_m1.go.gram(d, 'biofilm')
# m1.siderophores.gram<- run_m1.go.gram(d,'siderophores')
# m1.ab_degradation.gram<- run_m1.go.gram(d,'ab_degradation')
# m1.secretion_system.gram<- run_m1.go.gram(d,'secretion_system')
# m1.quorum_sensing.gram<- run_m1.go.gram(d,'quorum_sensing')
# m1.vf.gram<- run_m1.go.gram(d, 'vf')

# summary(m1.biofilm.gram)          # DIC = 82
# summary(m1.siderophores.gram)     # DIC = 82
# summary(m1.ab_degradation.gram)   # DIC = 82
# summary(m1.quorum_sensing.gram)   # DIC = 82
# summary(m1.secretion_system.gram) # *, DIC = 84
# summary(m1.ss.gram) # **, DIC = 79
# summary(m1.vf.gram) # *** , DIC = 75

# DIC are the same, results are qualitatively the same
# Significance levels are the same, except for secretome which has higher significance when including gram
# convergence is fine in both cases
# keep model without gram to focus on main emssage in figures/text
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




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


# DEPRECATED ~~~~~~~~  keep for records but remove from final analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# standardize the nb_extracellular respective to the gram_profile, so that we can include it in the multivariate model without having gram profile along, which would not be releavant for the other predictors
# d$z_nb_extracellular_gram_spec<- NA
# 
# d$z_nb_extracellular_gram_spec[which(d$gram_profile == 'p')]<- 
# (d$nb_extracellular[d$gram_profile == 'p'] - mean(d$nb_extracellular[d$gram_profile == 'p']))/sd(d$nb_extracellular[d$gram_profile == 'p'])
# 
# d$z_nb_extracellular_gram_spec[which(d$gram_profile == 'n')]<- 
#   (d$nb_extracellular[d$gram_profile == 'n'] - mean(d$nb_extracellular[d$gram_profile == 'n']))/sd(d$nb_extracellular[d$gram_profile == 'n'])
# 
# 
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


m1.multi.with_vf_NOgram<- MCMCglmm(pathogen  ~ 1 + z_nb_cds +  z_biofilm + z_ab_degradation + z_quorum_sensing + z_secretion_system + z_siderophores + z_nb_extracellular + z_vf,
                    random = ~species,
                    ginverse = list(species=Ainv),
                    data = d,
                    prior = prior,
                    family=c("categorical"),trunc=T,
                    start=list(QUASI=FALSE),
                    #pl = TRUE, pr = TRUE, nodes = 'ALL',
                    DIC = TRUE,
                    verbose = TRUE,
                    nitt=nitt, thin=thin, burnin=burnin)



#summary(m1.multi.with_vf_WITHgram)
summary(m1.multi.with_vf_NOgram)
plot(m1.multi.with_vf_WITHgram)
plot(m1.multi.with_vf_NOgram)
# results w/wo gram profile are the same
# in fact significance levels even higher if include gram profile
# DIC are similar
# variance associated with species are similar
# chains converge well in both cases
# So let's keep model without gram to avoid unnecessary details and simpflify plotting


# SAVE & CHECK OUTPUT ----


#save.image("./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN1.RData")
#save.image("./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN2.RData")
#save.image("./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN3.RData")


# GELMAN RUBIN TESTS ----

# in figures script











