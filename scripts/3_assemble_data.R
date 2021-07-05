# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# This script assembles the annotation data t produce dataframe used in downstream statistical analyses codes
# Specifically:
# 1) Uses PANNZER output, to see for each CDS which is top hit assigned GO, and whether this GO falls within a social trait and if yes which one
# 2) Uses PATRIC features table to record total number of CDS for each species
# 3) Uses VICTOR tables to record number of genes that are virulence factors (+ gene-level records as well)
# 4) Uses PSORTB output to records number of secreted proteins (+ gene-level records as well)


#local_project_dir=PATH/TO/CLONE/REPO # '/Users/s1687811/Documents/PhD/Research/CooperativePathogenicityVirulence/'
setwd(local_project_dir)
library(dplyr)
library(tidyr)
library(readxl)

# (I) ASSEMBLE DATA FOR VF AND PATHOGENS ANALYSIS ----

# PANNZER OUTPUT (GO) ----
 
argot_filter = 1

path.to.pannzer.output<- paste0(local_project_dir, '/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/pannzer_output/')
files<- paste0(path.to.pannzer.output, list.files(path.to.pannzer.output))
sps<- list.files(path.to.pannzer.output)

social_go_list<- as.data.frame(read_excel('./CooperativePathogenicityVirulence_repo/data/0_species_info_files/social_go_list_final.xls'))
social_gos<- social_go_list
cooperative_behaviours<- c('biofilm', 'antibiotic_degradation', 'quorum_sensing', 'siderophores', 'secretion_system')

per_peg_annotations<- vector('list')


for(s in 1:length(sps)){
  
  print(s)
  species<- sps[s]
  pegs_assignation<- vector('list')
  
  
  for(k in 1:length(cooperative_behaviours)){
    
    focal_behaviour<- cooperative_behaviours[k]
    
    sp<- read.csv(files[s], header=TRUE, sep = '\t', colClasses = c(rep('character', 4), rep('numeric', 3)))
    sp<- sp[sp$ARGOT_rank <= argot_filter,]
    sp<- sp[,1:4]
    colnames(sp)<- c('peg', 'Ontology', 'GO_id', 'Description')
    sp$GO_id<- paste0('GO:', sp$GO_id)
    # Fix peg name to match that of PATRIC feature tables
    sp$peg<- paste0(do.call('rbind', strsplit(sp$peg, '\\|'))[,1], '|', do.call('rbind', strsplit(sp$peg, '\\|'))[,2])
    sp$peg<- gsub(' ', '', sp$peg)
    
    
    # open PATRIC features table of that species
    sp_cds<-
      read.csv(paste0(local_project_dir, '/CooperativePathogenicityVirulence_repo/data/1_patric/features/', gsub('GO.out', 'features', fixed = TRUE, sps[s])),
               header = TRUE, sep = '\t', colClasses = 'character') %>%
      filter(feature_type == 'CDS') %>%
      mutate(species = gsub('.GO.out', '', fixed = TRUE, sps[s])) %>%
      select(species, patric_id, product) %>%
      rename(peg = patric_id,
             product_patric = product)
    
    
    # Assemble both tables
    sp_cds_annot<- full_join(sp_cds, sp, 'peg')
    
    
    # Flag pegs that have a GO id part of focal_behaviour list
    sp_cds_annot$is_focal_behaviour <- sp_cds_annot$GO_id %in% social_gos[social_gos$behaviour == focal_behaviour,1]
    sp_cds_annot2<- sp_cds_annot[,c('peg', 'Description', 'is_focal_behaviour')]
    #test<- sp_cds_annot2 %>% group_by(peg) %>% summarise(focal_behaviour_counts = sum(is_focal_behaviour))
    test<- sp_cds_annot2 %>% group_by(peg) %>% mutate(focal_behaviour_counts = sum(is_focal_behaviour)) %>% select(peg, focal_behaviour_counts) %>% unique()
    test2<- data.frame(peg = test$peg, ifelse(test[,c('focal_behaviour_counts')] > 0, 1, 0))
    colnames(test2)<- c('peg', focal_behaviour)
    pegs_assignation[[k]]<- test2
    rm(test2)
    rm(test)
    rm(sp_cds_annot)
    rm(sp_cds_annot2)
    rm(sp)
    rm(sp_cds)
  }
  
  
  all<- left_join(pegs_assignation[[1]], pegs_assignation[[2]], by = 'peg') %>%
    left_join(pegs_assignation[[3]], 'peg') %>%
    left_join(pegs_assignation[[4]], 'peg') %>%
    left_join(pegs_assignation[[5]], 'peg')
  
  all$species<- gsub('.GO.out', '', sps[s], fixed = TRUE)
  
  all<- select(all,
               species, peg, biofilm, antibiotic_degradation, quorum_sensing, siderophores, secretion_system)
  
  per_peg_annotations[[s]]<- all
  
}

per_peg_annotations_flat<- do.call('rbind', per_peg_annotations)


# PATRIC FEATURES (all CDS) ----

path.to.patric.features<- paste0(local_project_dir, '/CooperativePathogenicityVirulence_repo/data/1_patric/features/')
files<- paste0(path.to.patric.features, list.files(path.to.patric.features))
sps<- list.files(path.to.patric.features)
patric.features<- list()
for(j in 1:length(sps)){
  print(j)
  
  patric.features[[j]]<- read.csv(files[j], header = TRUE, sep = '\t', colClasses = 'character') %>%
    filter(feature_type == 'CDS') %>%
    mutate(species = gsub('.features', '', fixed = TRUE, sps[j])) %>%
    select(species, patric_id, product) %>%
    rename(peg = patric_id,
           product_patric = product)
}
patric.features.df<- do.call('rbind', patric.features)
nrow(patric.features.df)



# VICTOR (virulence factors) ----
focus.species<- read.table('./CooperativePathogenicityVirulence_repo/output/1_processed_tables/1.1_pathogen_commensal_genomes_118.txt', header = TRUE, sep = '\t', colClasses = 'character')
focus.species<- focus.species %>%
  select(genome_id, species_id, gram_profile)
# Load the virulence factor tabs
vf.list<- vector('list')
for(i in 1:nrow(focus.species)){
  print(i) 
  sp<- focus.species[i,'species_id']
  
  d<- read.csv(paste0('./CooperativePathogenicityVirulence_repo/data/1_patric/virulenceFactors/', sp,'.vf'), header=TRUE, sep = '\t', colClasses = 'character') %>%
    select(genome_id, patric_id, gene, product, source) %>%
    rename(qpid = patric_id) %>%
    filter(source == 'Victors')
  
  if(nrow(d) == 0){
    d[1,]<- NA
  }
  
  
  vf.list[[i]]<- d %>%
    mutate(species = sp) %>%
    rename(peg = qpid) %>%
    select(species, peg, gene, product, source)
  
  rm(d)
  
}
vf.tab<- do.call('rbind', vf.list)
vf.tab<- select(vf.tab, species, peg, source)



# PSORTB (secretome) ---- 

# Gathering PSORTb output
path.to.psotb.output<- paste0(local_project_dir, '/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb_output/')
files<- paste0(path.to.psotb.output, list.files(path.to.psotb.output))
sps<- list.files(path.to.psotb.output)
psortb.outs<- list()
for(s in 1:length(sps)){
  
  coop<- read.csv(files[s], sep='\t', row.names = NULL)
  print(s)
  colnames(coop)<- c(colnames(coop)[-1], 'foo')
  coop.keep<- coop[,c('SeqID', 'Final_Localization')]
  colnames(coop.keep)<- c('peg', 'secretome')
  
  coop.keep$peg<- paste0(do.call('rbind', strsplit(coop.keep$peg, '\\|'))[,1], '|', do.call('rbind', strsplit(coop.keep$peg, '\\|'))[,2])
  coop.keep$peg<- gsub(' ', '', coop.keep$peg)
  coop.keep$secretome<- ifelse(coop.keep$secretome == 'Extracellular', 1, 0)
  
  coop.keep$species<- gsub('.psortb.out', '', sps[s], fixed = TRUE)
  
  psortb.outs[[s]]<- coop.keep
  
}
psortb.outs.df<- do.call('rbind', psortb.outs) %>%
  select(species, peg, secretome)



# ASSEMBLE ALL ----

genes_annotations<- left_join(patric.features.df, per_peg_annotations_flat) %>%
  left_join(psortb.outs.df) %>%
  #left_join(mp3.out.df) %>%
  left_join(vf.tab)


genes_annot.final<-
  rename(genes_annotations, is_victor_vf = source) %>%
  mutate(is_victor_vf = ifelse(is.na(is_victor_vf), 0, 1),
         total_cds = 1)


genes_annot.final<- genes_annot.final %>%
  #select(species, peg, product_patric, total_cds, biofilm, antibiotic_degradation, quorum_sensing, siderophores, secretion_system, secretome, is_mp3_hybrid_pathogenic, is_victor_vf) %>%
  select(species, peg, product_patric, total_cds, biofilm, antibiotic_degradation, quorum_sensing, siderophores, secretion_system, secretome, is_victor_vf) %>%
  rename(species_id = species)
  


# Summarise at Species level annotation
summary.annot<- genes_annot.final %>%
  select(-peg, -product_patric) %>%
  gather('trait', 'flag', 2:ncol(.)) %>%
  #rename(species_id = species)%>%
  group_by(species_id, trait) %>%
  summarise(n = sum(flag, na.rm = TRUE)) %>% # just two species that have one NA each on their entire genome. Just ignore it.
  spread(trait, n) %>%
  #select(species_id, total_cds, biofilm, antibiotic_degradation, quorum_sensing, siderophores, secretion_system, secretome, is_mp3_hybrid_pathogenic, is_victor_vf)
  select(species_id, total_cds, biofilm, antibiotic_degradation, quorum_sensing, siderophores, secretion_system, secretome, is_victor_vf)



# Add gram profile to both tables
genomes_selected<- read.table('./CooperativePathogenicityVirulence_repo/output/1_processed_tables/1.1_pathogen_commensal_genomes_118.txt', header=TRUE, sep = '\t')
genomes_selected2<- genomes_selected %>%
  select(species_id, genus, species, gram_profile, pathogen)

summary.annot.final<- left_join(genomes_selected2, summary.annot)

genes_annot.final2<- left_join(genomes_selected2, genes_annot.final)



write.table(summary.annot.final, './CooperativePathogenicityVirulence_repo/output/1_processed_tables/2.1_assembled_SPECIES_annotation.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
write.table(genes_annot.final2, './CooperativePathogenicityVirulence_repo/output/1_processed_tables/2.2_assembled_GENES_annotation.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')



# (II) ASSEMBLE DATA FOR CFR ANALYSIS ----

# We repeated the process for downloading PATRIC files + annotating them with PANNZZER and PSORTB for a second sets of genomes
# corresponding to genomes previously analyzed by Leggett et al 2017
# OUtput for these are in 4_data_cfr_analysis
# Then assemble data exactly as above

#local_project_dir=PATH/TO/CLONE/REPO # '/Users/s1687811/Documents/PhD/Research/CooperativePathogenicityVirulence/'
setwd(local_project_dir)
library(dplyr)
library(tidyr)
library(readxl)



# PANNZER OUTPUT (GO) ----

argot_filter = 1

path.to.pannzer.output<- paste0(local_project_dir, '/CooperativePathogenicityVirulence_repo/data/4_data_cfr_analysis/pannzer/')
files<- paste0(path.to.pannzer.output, list.files(path.to.pannzer.output))
sps<- list.files(path.to.pannzer.output)

social_go_list<- as.data.frame(read_excel('./CooperativePathogenicityVirulence_repo/data/0_species_info_files/social_go_list_final.xls'))
social_gos<- social_go_list
cooperative_behaviours<- c('biofilm', 'antibiotic_degradation', 'quorum_sensing', 'siderophores', 'secretion_system')

per_peg_annotations<- vector('list')



for(s in 1:length(sps)){
  
  print(s)
  species<- sps[s]
  pegs_assignation<- vector('list')
  
  
  for(k in 1:length(cooperative_behaviours)){
    
    focal_behaviour<- cooperative_behaviours[k]
    
    sp<- read.csv(files[s], header=TRUE, sep = '\t', colClasses = c(rep('character', 4), rep('numeric', 3)))
    sp<- sp[sp$ARGOT_rank <= argot_filter,]
    sp<- sp[,1:4]
    colnames(sp)<- c('peg', 'Ontology', 'GO_id', 'Description')
    sp$GO_id<- paste0('GO:', sp$GO_id)
    # Fix peg name to match that of PATRIC feature tables
    sp$peg<- paste0(do.call('rbind', strsplit(sp$peg, '\\|'))[,1], '|', do.call('rbind', strsplit(sp$peg, '\\|'))[,2])
    sp$peg<- gsub(' ', '', sp$peg)
    
    
    # open PATRIC features table of that species
    sp_cds<-
      read.csv(paste0(local_project_dir, '/CooperativePathogenicityVirulence_repo/data/4_data_cfr_analysis/features/', gsub('GO.out', 'features', fixed = TRUE, sps[s])),
               header = TRUE, sep = '\t', colClasses = 'character') %>%
      filter(feature_type == 'CDS') %>%
      mutate(species = gsub('.GO.out', '', fixed = TRUE, sps[s])) %>%
      select(species, patric_id, product) %>%
      rename(peg = patric_id,
             product_patric = product)
    
    
    # Assemble both tables
    sp_cds_annot<- full_join(sp_cds, sp, 'peg')
    
    
    # Flag pegs that have a GO id part of focal_behaviour list
    sp_cds_annot$is_focal_behaviour <- sp_cds_annot$GO_id %in% social_gos[social_gos$behaviour == focal_behaviour,1]
    sp_cds_annot2<- sp_cds_annot[,c('peg', 'Description', 'is_focal_behaviour')]
    #test<- sp_cds_annot2 %>% group_by(peg) %>% summarise(focal_behaviour_counts = sum(is_focal_behaviour))
    test<- sp_cds_annot2 %>% group_by(peg) %>% mutate(focal_behaviour_counts = sum(is_focal_behaviour)) %>% select(peg, focal_behaviour_counts) %>% unique()
    test2<- data.frame(peg = test$peg, ifelse(test[,c('focal_behaviour_counts')] > 0, 1, 0))
    colnames(test2)<- c('peg', focal_behaviour)
    pegs_assignation[[k]]<- test2
    rm(test2)
    rm(test)
    rm(sp_cds_annot)
    rm(sp_cds_annot2)
    rm(sp)
    rm(sp_cds)
  }
  
  
  all<- left_join(pegs_assignation[[1]], pegs_assignation[[2]], by = 'peg') %>%
    left_join(pegs_assignation[[3]], 'peg') %>%
    left_join(pegs_assignation[[4]], 'peg') %>%
    left_join(pegs_assignation[[5]], 'peg')
  
  all$species<- gsub('.GO.out', '', sps[s], fixed = TRUE)
  
  all<- select(all,
               species, peg, biofilm, antibiotic_degradation, quorum_sensing, siderophores, secretion_system)
  
  per_peg_annotations[[s]]<- all
  
}

per_peg_annotations_flat<- do.call('rbind', per_peg_annotations)


# PATRIC FEATURES (all CDS) ----

path.to.patric.features<- paste0(local_project_dir, '/CooperativePathogenicityVirulence_repo/data/4_data_cfr_analysis/features/')
files<- paste0(path.to.patric.features, list.files(path.to.patric.features))
sps<- list.files(path.to.patric.features)
patric.features<- list()
for(j in 1:length(sps)){
  print(j)
  
  patric.features[[j]]<- read.csv(files[j], header = TRUE, sep = '\t', colClasses = 'character') %>%
    filter(feature_type == 'CDS') %>%
    mutate(species = gsub('.features', '', fixed = TRUE, sps[j])) %>%
    select(species, patric_id, product) %>%
    rename(peg = patric_id,
           product_patric = product)
}
patric.features.df<- do.call('rbind', patric.features)
nrow(patric.features.df)


# VICTOR (virulence factors) ----
focus.species<- read.table('./CooperativePathogenicityVirulence_repo/output/1_processed_tables/3.2_cfr_SUPFAM_match.txt', header=TRUE, sep = '\t', colClasses = 'character')

focus.species<- focus.species %>%
  select(genome_id, Legget_species, file_name_record, matching_supfam_id, gram_profile)  %>%
  filter(!is.na(file_name_record))
# Load the virulence factor tabs
vf.list<- vector('list')
for(i in 1:nrow(focus.species)){
  print(i) 
  sp<- focus.species[i,'file_name_record']
  
  d<- read.csv(paste0('./CooperativePathogenicityVirulence_repo/data/4_data_cfr_analysis/virulenceFactors/', sp,'.vf'), header=TRUE, sep = '\t', colClasses = 'character') %>%
    select(genome_id, patric_id, gene, product, source) %>%
    rename(qpid = patric_id) %>%
    filter(source == 'Victors')
  
  if(nrow(d) == 0){
    d[1,]<- NA
  }
  
  
  vf.list[[i]]<- d %>%
    mutate(species = sp) %>%
    rename(peg = qpid) %>%
    select(species, peg, gene, product, source)
  
  rm(d)
  
}
vf.tab<- do.call('rbind', vf.list)
vf.tab<- select(vf.tab, species, peg, source)



# PSORTB (secretome) ---- 

# Gathering PSORTb output
path.to.psotb.output<- paste0(local_project_dir, '/CooperativePathogenicityVirulence_repo/data/4_data_cfr_analysis/psortb/')
files<- paste0(path.to.psotb.output, list.files(path.to.psotb.output))
sps<- list.files(path.to.psotb.output)
psortb.outs<- list()
for(s in 1:length(sps)){
  
  coop<- read.csv(files[s], sep='\t', row.names = NULL)
  print(s)
  colnames(coop)<- c(colnames(coop)[-1], 'foo')
  coop.keep<- coop[,c('SeqID', 'Final_Localization')]
  colnames(coop.keep)<- c('peg', 'secretome')
  
  coop.keep$peg<- paste0(do.call('rbind', strsplit(coop.keep$peg, '\\|'))[,1], '|', do.call('rbind', strsplit(coop.keep$peg, '\\|'))[,2])
  coop.keep$peg<- gsub(' ', '', coop.keep$peg)
  coop.keep$secretome<- ifelse(coop.keep$secretome == 'Extracellular', 1, 0)
  
  coop.keep$species<- gsub('.psortb.out', '', sps[s], fixed = TRUE)
  
  psortb.outs[[s]]<- coop.keep
  
}
psortb.outs.df<- do.call('rbind', psortb.outs) %>%
  select(species, peg, secretome)


# ASSEMBLE ALL ----

genes_annotations<- left_join(patric.features.df, per_peg_annotations_flat) %>%
  left_join(psortb.outs.df) %>%
  #left_join(mp3.out.df) %>%
  left_join(vf.tab)


genes_annot.final<-
  rename(genes_annotations, is_victor_vf = source) %>%
  mutate(is_victor_vf = ifelse(is.na(is_victor_vf), 0, 1),
         total_cds = 1)


genes_annot.final<- genes_annot.final %>%
  #select(species, peg, product_patric, total_cds, biofilm, antibiotic_degradation, quorum_sensing, siderophores, secretion_system, secretome, is_mp3_hybrid_pathogenic, is_victor_vf) %>%
  select(species, peg, product_patric, total_cds, biofilm, antibiotic_degradation, quorum_sensing, siderophores, secretion_system, secretome, is_victor_vf) %>%
  rename(species_id = species)



# Summarise at Species level annotation
summary.annot<- genes_annot.final %>%
  select(-peg, -product_patric) %>%
  gather('trait', 'flag', 2:ncol(.)) %>%
  #rename(species_id = species)%>%
  group_by(species_id, trait) %>%
  summarise(n = sum(flag, na.rm = TRUE)) %>% # just two species that have one NA each on their entire genome. Just ignore it.
  spread(trait, n) %>%
  #select(species_id, total_cds, biofilm, antibiotic_degradation, quorum_sensing, siderophores, secretion_system, secretome, is_mp3_hybrid_pathogenic, is_victor_vf)
  select(species_id, total_cds, biofilm, antibiotic_degradation, quorum_sensing, siderophores, secretion_system, secretome, is_victor_vf)


# Add genomes and species details info + Leggett data to these

spInfoAppend<- focus.species %>%
  select(file_name_record, Legget_species, genome_id, matching_supfam_id, gram_profile) %>%
  rename(species_id = file_name_record)

legget_data<- read.table('./CooperativePathogenicityVirulence_repo/data/4_data_cfr_analysis/Leggett_data_bacteria_only_matching_genomes_only.txt', header=TRUE, sep = '\t')

legget_data.edit<- legget_data %>% # 3 species dropped because no CFR data for them --> 50 species left
  filter(!is.na(case_fatality_rate),
         !is.na(infection_route)) %>%
  rename(Legget_species = species2)

cfr.d<- left_join(legget_data.edit, spInfoAppend, by = 'Legget_species')



summary.annot.final<- left_join(cfr.d, summary.annot, by = 'species_id')
nrow(summary.annot) # 53
nrow(summary.annot.final) # 50, because three species dropped because missing CFR data


head(genes_annot.final)
nrow(genes_annot.final)
length(unique(genes_annot.final$species_id))
cfr.d$species_id %in% unique(genes_annot.final$species_id)
genes_annot.final2<- left_join(cfr.d, genes_annot.final, by = 'species_id')
nrow(genes_annot.final2)
length(unique(genes_annot.final2$species_id)) # 50 species left


write.table(summary.annot.final, './CooperativePathogenicityVirulence_repo/output/1_processed_tables/4.1_assembled_CFR_SPECIES_annotation.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
write.table(genes_annot.final2, './CooperativePathogenicityVirulence_repo/output/1_processed_tables/4.2_assembled_CFR_GENES_annotation.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


