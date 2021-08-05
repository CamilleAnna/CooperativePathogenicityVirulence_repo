# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```~~~~~~~~~~~ #
#                   Script for comparative analysis of virulence factors
#
#   Takes as input:
#         2.1_assembled_SPECIES_annotation.txt
#         2.2_assembled_GENES_annotation.txt
#         midas_tree_renamed.newick
#   Runs analysis:
#         Use genes annotations to make contingency table VFxCoop
#         Compute Odds ratio and SE on the odds ratio
#         Use MCMCglmm to run phylogenetic meta-analysis
#   Output:
#         3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN3.RData
#         3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN3.RData
#         3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN3.RData
#         additional_figs/VF_ForestFunnel_VICTOR_dataset118.pdf
#
# Code for extracting tables and producing figures in separate script file, using RData object
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#local_project_dir=PATH/TO/CLONE/REPO # local_project_dir='/Users/s1687811/Documents/PhD/Research/CooperativePathogenicityVirulence/'

setwd(local_project_dir)
setwd('./CooperativePathogenicityVirulence_repo/')
source('./scripts/sourced_packages.R')
library(ape)
library(MCMCglmm)
library(dmetar)
library(esc)


# 0 - LOAD DATA  ----

# Annotations summarised at species level
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



# Detail annotation at gene level
d.annot<- read.csv('./output/1_processed_tables/2.2_assembled_GENES_annotation.txt', header=TRUE, sep = '\t',
               colClasses = c(rep('character', 7),
                              rep('numeric', 8)))

head(d.annot)

sum(is.na(d.annot$secretome)) # all genes are annotated

d.annot2<- d.annot %>%
  filter(!is.na(secretome)) %>% # in fact no case of this
  mutate(is_coop = biofilm+antibiotic_degradation+quorum_sensing+siderophores+secretion_system+secretome)%>%
  mutate(is_coop = ifelse(is_coop == 0, is_coop, 1)) %>%
  select(-genus, -species, -product_patric, -total_cds) %>%
  rename(is_vf = is_victor_vf,
         species = species_id)

head(d.annot2)




foo<- data.frame(species = rep(LETTERS[1:4], each = 5),
           biofilm = sample(c(0,1), replace = TRUE, 20),
           secretome = sample(c(0,1), replace = TRUE, 20),
           is_vf = sample(c(0,1), replace = TRUE, 20))




foo2<- foo %>%
  filter(is_vf == 1) %>%
  mutate(biofilm = ifelse(biofilm == 1, 'biofilm', '')) %>%
  mutate(secretome = ifelse(secretome == 1, 'secretome', '')) %>%
  mutate(traits = paste(biofilm, secretome, sep = '_'))



as.data.frame(table(foo2$traits)) %>%
ggplot(., aes(fill=Var1, y=Freq, x=1)) + 
  geom_bar(position="fill", stat="identity")



head(foo2)

get_prop_vf_coop<- function(df, focal_sp){
test<- df %>% filter(species == focal_sp)
dfoo<- data.frame(species = focal_sp,
           as.data.frame(table(test$traits)),
           nb_vf = nrow(test))
 
return(dfoo)          
}



yo<- data.frame(species = character(),
                Var1 = character(),
                Freq = numeric(),
                nb_vf = numeric())


for(i in 1:length(unique(foo$species))){
  
  yo<- rbind(yo, get_prop_vf_coop(foo2, unique(foo$species)[i]))
}


yo$prop = yo$Freq/yo$nb_vf
yo

ggplot(yo, aes(x = Var1, y = prop, fill = Var1))+
  geom_boxplot()





foo3<- d.annot2 %>%
  filter(is_vf == 1) %>%
  mutate(biofilm = ifelse(biofilm == 1, 'biofilm', '')) %>%
  mutate(antibiotic_degradation = ifelse(antibiotic_degradation == 1, 'antibDegr', '')) %>%
  mutate(quorum_sensing = ifelse(quorum_sensing == 1, 'quorumSensing', '')) %>%
  mutate(siderophores = ifelse(siderophores == 1, 'siderophores', '')) %>%
  mutate(secretion_system = ifelse(secretion_system == 1, 'secrSyst', '')) %>%
  mutate(secretome = ifelse(secretome == 1, 'secretome', '')) %>%
  mutate(traits = paste(biofilm, antibiotic_degradation, quorum_sensing, siderophores, secretion_system, secretome, sep = '_'))




yo<- data.frame(species = character(),
                Var1 = character(),
                Freq = numeric(),
                nb_vf = numeric())


for(i in 1:length(unique(foo3$species))){
  
  yo<- rbind(yo, get_prop_vf_coop(foo3, unique(foo3$species)[i]))
}


yo$prop = yo$Freq/yo$nb_vf
yo

ggplot(yo, aes(x = Var1, y = prop, fill = Var1))+
  geom_boxplot()








# 1 - QUICK OVERAL CHI-SQUARE ----

# Code wrapper:
# make a contingency table for a focus cooperative trait
# run a chi-squre test on it
# compute OR and CI

ctab.all<- function(focus_trait){
  
  tab<- d.annot2 %>% 
    select(focus_trait, is_vf)%>%
    rename(cooperative = 1) %>%
    #group_by(cooperative, is_vf) %>%
    group_by(is_vf,cooperative) %>%
    table()
  
  
  or<- esc_2x2(grp1yes = tab[2,2],
               grp1no = tab[2,1],
               grp2yes = tab[1,2],
               grp2no = tab[1,1],
               es.type = 'or')
  
  
  chi<- chisq.test(tab)
  
  return(list(tab = tab, or = or, chi = chi))
}

ctab.all.secretome<- ctab.all('secretome')
ctab.all.biofilm<- ctab.all('biofilm')
ctab.all.sid<- ctab.all('siderophores')  # marginally 
ctab.all.qs<- ctab.all('quorum_sensing')
ctab.all.ab<- ctab.all('antibiotic_degradation') # approx may be incorrect # ns
ctab.all.ssyst<- ctab.all('secretion_system')    # approx may be incorrect
ctab.all.any<- ctab.all('is_coop') 




d.annot2$biofilm2<- ifelse(d.annot2$biofilm == 0, '', 'biofilm')
d.annot2$quorum_sensing2<- ifelse(d.annot2$quorum_sensing == 0, '', 'qs')
d.annot2$siderophores2<- ifelse(d.annot2$siderophores == 0, '', 'sid')
d.annot2$antibiotic_degradation2<- ifelse(d.annot2$antibiotic_degradation == 0, '', 'ab')
d.annot2$secretion_system2<- ifelse(d.annot2$secretion_system == 0, '', 'SecSyst')
d.annot2$secretome2<- ifelse(d.annot2$secretome == 0, '', 'secretome')
d.annot2$VF2<- ifelse(d.annot2$is_vf == 0, '', 'VF')




foo<- d.annot2[d.annot2$is_vf == 1,]

foo$profile<- 
paste(foo$VF2, foo$secretome2, foo$secretion_system2, foo$antibiotic_degradation2, foo$siderophores2, foo$quorum_sensing2, foo$biofilm2)


unique(foo$profile)


foo.tab<- as.data.frame(table(foo$profile))
foo.tab$foo = 1

ggplot(foo.tab[foo.tab$Var1!='VF      ',], aes(x = foo, y = Freq, fill = Var1))+
  geom_bar(stat = "identity", position = "fill")




# 2 - ACCOUNTING FOR PHYLOGENY ----
# 2.1 - CONTINGENCY TABLES ----

# WORKS WITH ENVIRONMENT DATAFRAME 'd.annot2'

# Code wrappers:
# get_conttable:
#   makes a contingency table for a given trait for a given species
# get_effects:  
#   run metabin() on each 2x2 table
#   extarcts individual species effects, SE and CI --> then used in PMM
#   also runs a meta-analysis, although not phylogenetic


# NOTE: I use basic 'MH' method. Using  GLMM method gives different results for the meta-analysis but the individual effects exacted are the same
get_conttable<- function(focal_trait = ''){
  
  d.annot2[,c('species', focal_trait, 'is_vf')] %>%
    rename(renamed_focal_trait = 2) %>%
    mutate(treatment_outcome = paste(renamed_focal_trait, is_vf, sep = '_')) %>%
    select(species, treatment_outcome) %>%
    group_by(species, treatment_outcome) %>%
    table() %>%
    as.data.frame() %>%
    spread(treatment_outcome, Freq) %>%
    rename(coop0_vf0 = `0_0`,
           coop0_vf1 = `0_1`,
           coop1_vf0 = `1_0`,
           coop1_vf1 = `1_1`)
}
get_effects<- function(ctab=''){
  
  m<- metabin(data = ctab,
              event.e = coop1_vf1,
              n.e = coop1_vf1+coop1_vf0,
              event.c = coop0_vf1,
              n.c = coop0_vf1+coop0_vf0,
              sm = "or",
              incr = 0.5, MH.exact = FALSE,
              method = 'MH', 
              allstudies = FALSE)
  
  ctab$logOR = m$TE
  ctab$se_logOR = m$seTE
  ctab$ci_l = m$lower
  ctab$ci_u = m$upper
  ctab$z = m$zval
  ctab$pval = m$pval
  
  return(list(df = ctab, meta.out = m))
  
}

# Make contingency tables
ctab.iscoop<- get_conttable(focal_trait = 'is_coop')
ctab.secretome<- get_conttable(focal_trait = 'secretome')
ctab.ab<- get_conttable(focal_trait = 'antibiotic_degradation')
ctab.qs<- get_conttable(focal_trait = 'quorum_sensing')
ctab.sid<- get_conttable(focal_trait = 'siderophores')
ctab.ssyst<- get_conttable(focal_trait = 'secretion_system')
ctab.biofilm<- get_conttable(focal_trait = 'biofilm')

library(meta)

# Run metabin to extract effects
effects.iscoop<- get_effects(ctab.iscoop)
effects.secretome<- get_effects(ctab.secretome)
effects.ab<- get_effects(ctab.ab)
effects.qs<- get_effects(ctab.qs)
effects.sid<- get_effects(ctab.sid)
effects.ssyst<- get_effects(ctab.ssyst)
effects.biofilm<- get_effects(ctab.biofilm)


# ASSEMBLE ALL OUTPUT
df.all<- rbind(
  cbind(effects.iscoop$df,cooperative_trait = 'is_coop'),
  cbind(effects.secretome$df,cooperative_trait = 'secretome'),
  cbind(effects.ab$df,cooperative_trait = 'ab'),
  cbind(effects.qs$df,cooperative_trait = 'qs'),
  cbind(effects.sid$df,cooperative_trait = 'sid'),
  cbind(effects.ssyst$df,cooperative_trait = 'ssyst'),
  cbind(effects.biofilm$df,cooperative_trait = 'biofilm'))


df.all$species.ide = df.all$species
unique(df.all$species)


# SANITY CHECK ON THE PROCESS USED TO COMPUTE OR AND CI
t<- df.all[df.all$cooperative_trait == 'sid',] # Try on siderophore 
t2<- t[-which(rowSums(t[,2:5] == 0)>1),]  # filter those with > 1 zero cell
t2[which(rowSums(t2[,2:5] == 0) > 0), 2:5]<- t2[which(rowSums(t2[,2:5] == 0) > 0), 2:5]+0.5 # add 0.5 to studies with a zero cell

t2$myOR<- (t2$coop1_vf1/t2$coop1_vf0)/(t2$coop0_vf1/t2$coop0_vf0) # compute OR
t2$mylogOR<- log(t2$myOR) # transform OR into log(or)
t2$mySE<- sqrt((1/t2$coop0_vf0)+ (1/t2$coop0_vf1) + (1/t2$coop1_vf0) + (1/t2$coop1_vf1)) # compute SE

round(t2$mylogOR, 5) == round(t2$logOR, 5) # this matches
round(t2$mySE, 5) == round(t2$se_logOR, 5) # this matches 




# 2.2 - QUICK PLOT CHECK ----

sorted_forest<- function(df, trait){
  
  df$species<- factor(df$species, df$species[order(df$logOR)])
  
  #df$species<- factor(df$species, df$species[rev(order(df$ci_u - df$ci_l))])
  
  
  p<- ggplot(df, aes(x = logOR, y = species))+
    #geom_errorbarh(aes(xmin = ci_l, xmax = ci_u))+
    geom_errorbarh(aes(xmin = logOR-se_logOR, xmax = logOR+se_logOR))+
    geom_point()+#xlim(0,100)+
    ggtitle(trait)+
    #facet_wrap(~cooperative_trait)+
    geom_vline(xintercept = 0, col = 'darkgrey', linetype = 'dashed')+
    theme(axis.text.y = element_text(size = 3))
  
  return(p)
  
}

p.sf.iscoop<- sorted_forest(effects.iscoop$df, trait = 'is_coop')
p.sf.secretome<- sorted_forest(effects.secretome$df, trait = 'secretome')
p.sf.ab<- sorted_forest(effects.ab$df, trait = 'ab')
p.sf.qs<- sorted_forest(effects.qs$df, trait = 'qs')
p.sf.sid<- sorted_forest(effects.sid$df, trait = 'sid')
p.sf.ssyst<- sorted_forest(effects.ssyst$df, trait = 'ssyst')
p.sf.biofilm<- sorted_forest(effects.biofilm$df, trait = 'biofilm')



pdf('./output/2_figures/additional_figs/VF_ForestFunnel_VICTOR_dataset118.pdf', width = 8.75, height = 11.5)

grid.arrange(p.sf.iscoop, p.sf.secretome, p.sf.ab, p.sf.qs,
             p.sf.sid, p.sf.ssyst, p.sf.biofilm, ncol = 3)


ggplot(df.all, aes(x = logOR, y = 1/se_logOR))+
  geom_point(alpha = 0.5)+
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey')+
  facet_wrap(~cooperative_trait, scales = 'free_y')

dev.off()


# 2.3 - META-ANALYSIS ON ODD RATIO ----

library(MCMCglmm)
library(ape)

midas.tree<- read.tree('./data/5_phylogeny_files/midas_tree_renamed.newick')


div = 2000
nitt = 10500000/div
burnin = 500000/div
thin = ceiling(5000/div)
nitt
(nitt-burnin)/thin


a <- 1000
prior.a <- list(R = list(V = diag(1), nu = 0.002), # residual variance = non-phylo variance
                G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),  # phylogeny
                         G2 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a) # measurement error, estimate it rather than fix it to 1
                         )) 


# Code wrapper to run meta-analysis on a given trait
run_meta<- function(df = '', focal_trait, prior = ''){
  
  df.tmp = df[df$cooperative_trait %in% focal_trait,]
  
  phylogeny<- drop.tip(midas.tree, midas.tree$tip.label[which(!midas.tree$tip.label %in% df.tmp$species)])
  phylogeny<- chronopl(phylogeny, lambda = 0)
  phylogeny<-makeNodeLabel(phylogeny)
  Ainv<-inverseA(phylogeny, scale=FALSE)$Ainv
  
  randomtest <- MCMCglmm(logOR ~ 1,
                         random = ~species + idh(se_logOR):units,   # non-phylo species effect here is residual because one observation per species, and that's also equivalent to "study effect" in a usual meta-analysis
                         #random = ~species + se_logOR,             # non-phylo species effect here is residual because one observation per species, and that's also equivalent to "study effect" in a usual meta-analysis
                         ginverse = list(species=Ainv),
                         family = 'gaussian',
                         prior = prior,
                         data = df.tmp,
                         start = list(QUASI = FALSE),
                         verbose = FALSE,
                         nitt=nitt, thin=thin, burnin=burnin)
  
  return(randomtest)
  
}

# Run meta-analysis on each cooperation category
meta.a.secretome<- run_meta(df = df.all[!is.na(df.all$logOR),], focal_trait = 'secretome', prior = prior.a)
meta.a.biofilm<- run_meta(df = df.all[!is.na(df.all$logOR),], focal_trait = 'biofilm', prior = prior.a)
meta.a.sid<- run_meta(df = df.all[!is.na(df.all$logOR),], focal_trait = 'sid', prior = prior.a)
meta.a.ab<- run_meta(df = df.all[!is.na(df.all$logOR),], focal_trait = 'ab', prior = prior.a)
meta.a.qs<- run_meta(df = df.all[!is.na(df.all$logOR),], focal_trait = 'qs', prior = prior.a)
meta.a.ssyst<- run_meta(df = df.all[!is.na(df.all$logOR),], focal_trait = 'ssyst', prior = prior.a)

# Run meta-analysis on grouped categories under a single 'cooperative' category
meta.a.coop<- run_meta(df = df.all[!is.na(df.all$logOR),], focal_trait = 'is_coop', prior = prior.a)


meta.a.secretome.1<- run_meta(df = df.all[!is.na(df.all$logOR),], focal_trait = 'secretome', prior = prior.a)
meta.a.secretome.2<- run_meta(df = df.all[!is.na(df.all$logOR),], focal_trait = 'secretome', prior = prior.a)

summary(meta.a.secretome.1)$solutions
summary(meta.a.secretome.2)$solutions



# Check how many datapoint there were for in each model
# (if a species has two entries of the 2x2 table VF/COOP that are equal to zero, the odds ratio cannot be computed)
df = df.all[!is.na(df.all$logOR),]
nrow(df[df$cooperative_trait %in% 'is_coop',])
nrow(df[df$cooperative_trait %in% 'secretome',])
nrow(df[df$cooperative_trait %in% 'biofilm',])
nrow(df[df$cooperative_trait %in% 'sid',])
nrow(df[df$cooperative_trait %in% 'ab',])
nrow(df[df$cooperative_trait %in% 'qs',])
nrow(df[df$cooperative_trait %in% 'ssyst',])


# SAVE ----

# We run the analysis chains to then use Gelman-Rubin test to assess convergence
#save.image("./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN1.RData")
#save.image("./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN2.RData")
#save.image("./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN3.RData")


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



# GELMAN-RUBIN TESTS ----

# --> in figures script













