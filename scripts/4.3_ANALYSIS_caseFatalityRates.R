# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```~~~~~~~~~~~ #
#                   Script for comparative analysis of Case Fatality Rates
#
#   Takes as input:
#         assembled_CFR_SPECIES_annotation.txt
#         supfam_cfr_tree.newick
#   Runs analysis:
#         get cfr dataset assembled
#         compute log+1 measure to run gaussian response model
#         replicate legget's 2017 analysis, althought here very different modelling approach
#         but recover main result that transmission route has a strong effect
#         include transmission route as covariate in all subsequent models
#         run 6 univariate models to test effect of each form of cooperation on CFR
#         run on multivariate model with all traits, z-transformed, included as predictors
#   Output:
#         3_model_output/CompAnalysis_cfr_CHAIN1.RData
#         3_model_output/CompAnalysis_cfr_CHAIN2.RData
#         3_model_output/CompAnalysis_cfr_CHAIN3.RData
#
# Code for extracting tables and producing figures in separate script file, using RData object
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#local_project_dir=PATH/TO/CLONE/REPO # '/Users/s1687811/Documents/PhD/Research/CooperativePathogenicityVirulence/'

setwd(local_project_dir)
setwd('./CooperativePathogenicityVirulence_repo/')
source('./scripts/sourced_packages.R')
library(ape)
library(MCMCglmm)
library(dmetar)
library(esc)


# CFR DATASET ----
d.cfr<- read.table('./output/1_processed_tables/4.1_assembled_CFR_SPECIES_annotation.txt', header=TRUE,
                   sep = '\t',
                   colClasses = c('character', 'character', 'character',
                                  'numeric', 'numeric',
                                  'character', 'character', 'character',
                                  'numeric', 'numeric',
                                  rep('character', 4), rep('numeric', 8)))

d.cfr<- rename(d.cfr,
               supfam_id = matching_supfam_id,
               nb_cds = total_cds)


# This reproduces Figure 3a of Leggett et al
# So the CFR in the data must indeed be a PERCENTAGE
# Also confirms the data I have scrapped for infective dose is correct
 ggplot(d.cfr, aes(x = log(infective_dose), y = case_fatality_rate, col = infection_route))+
   geom_point()+
   scale_color_manual(values = c('lightgrey', 'darkgrey', 'black'))


 
# PHYLOGENY
tree<- read.tree('./data/5_phylogeny_files/supfam_cfr_tree.newick')
tree$tip.label
tree<- drop.tip(tree, tree$tip.label[!tree$tip.label %in% d.cfr$supfam_id])
tree$tip.label<- d.cfr$species_id[match(tree$tip.label, d.cfr$supfam_id)]
tree<- chronopl(tree, lambda = 1)
#plot(tree) # ultrametricised tree
tree<-makeNodeLabel(tree)
Ainv<-inverseA(tree, scale=FALSE)$Ainv


# TRANSFORM RESPONSE
# Legget's data are %death, but no actual number is given.
# Can't analyse with binomial regression
# So went for a log+1 gaussian response instead

d.cfr$log_cfr_plus1<- log(d.cfr$case_fatality_rate + 1)

ggplot(d.cfr, aes(x = infection_route, y = log_cfr_plus1))+
  geom_boxplot()

ggplot(d.cfr, aes(x = infection_route, y = log_cfr_plus1, fill = gram_profile))+
  geom_boxplot()


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


# as reported in Leggett's, Skin & Inhalation have significantly greater CFR
# their lower CI bound is higher than upper bound of CI of Ingestion
# However we do not recover result that generation time hase an effect
# However model strategy is different, their modelled a multivariate response
# and I think they transformed data for make up numbers corresponding to the %, because they say they used a binomial response model
# Anyway, main effect they report is the strong effect of transmission route
# Will include this in my model



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


m3.bio<- m3(df = d.cfr, focal_trait = 'biofilm')
m3.qs<- m3(df = d.cfr, focal_trait = 'quorum_sensing')
m3.ab<- m3(df = d.cfr, focal_trait = 'antibiotic_degradation')
m3.sid<- m3(df = d.cfr, focal_trait = 'siderophores')
m3.ssy<- m3(df = d.cfr, focal_trait = 'secretion_system')
m3.sec<- m3(df = d.cfr, focal_trait = 'secretome')
m3.vf<- m3(df = d.cfr, focal_trait = 'is_victor_vf')


# DIC clearly better when including the transmission route, down to ~170
# while DIC ~190 for all models when transmission route is not included
# no significance for any of the traits
summary(m3.bio) # .
summary(m3.qs)
summary(m3.ab)
summary(m3.sid)
summary(m3.ssy)
summary(m3.sec)
summary(m3.vf)


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

summary(m3.multi.with_vf)
plot(m3.multi.with_vf)



# SAVE ----

#save.image("./output/3_model_output/CompAnalysis_cfr_CHAIN1.RData") #(former cfr_2604.RData)
#save.image("./output/3_model_output/CompAnalysis_cfr_CHAIN2.RData")
#save.image("./output/3_model_output/CompAnalysis_cfr_CHAIN3.RData")

# GELMAN-RUBIN tests ----

# in figures script








