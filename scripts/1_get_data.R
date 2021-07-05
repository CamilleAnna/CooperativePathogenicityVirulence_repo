# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# This scripts identifies a list of human associated pathogens and non-pathogens
# and outputs a table with listed species and identified representative genome
# from the PATRIC database to proceed with analysis
#
# identify PATHOGENS: For each "bacteria pathogen" listed in PATRIC database, we retreived the list of representative genomes, isolate from Human hosts = list of Human pathogen species
# identify NON-PATHOGENS: Use HMP dataset in curatedMetagenomicData() package to identify species present in >80% samples = set of core commensal Human species
#
# For each of these two lists we:
# - use MIDAS database to identify the most common representative strain for each species
# - For a few species where a candidate genome could not be found directly by grep matching MIDAS database, we manually browsed the MIDAS database. Details of genomes that were manualy identified are given in excel file "pathogen_commensals_MIDASmatch"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Path to cloned repository
#local_project_dir=PATH/TO/CLONE/REPO # '/Users/s1687811/Documents/PhD/Research/CooperativePathogenicityVirulence/'
setwd(paste0(local_project_dir, '/CooperativePathogenicityVirulence_repo'))
source('./scripts/sourced_packages.R')


# PATHOGENS ---- 
# ~~~~~~~~~~~~~~~#

setwd(local_project_dir)
source('./scripts/sourced_packages.R')

# PATRIC release files (from: https://docs.patricbrc.org/user_guides/ftp.html, accessed January 2021)
patric.lineages<- read.csv('./data/1_patric/1.1_patric_release_notes/genome_lineage', sep = '\t', colClasses = 'character')
patric.metadata<- read.csv('./data/1_patric/1.1_patric_release_notes/genome_metadata', colClasses = 'character', sep = '\t', quote='"')
patricdb<- left_join(patric.lineages, patric.metadata) %>%
  select(genus,species,genome_id,genome_name,taxon_id,genome_status,host_name,gram_stain,disease,strain,serovar,biovar,pathovar,type_strain,refseq_accessions,chromosomes,plasmids,contigs,patric_cds) %>%
  filter(genome_status %in% c('WGS', 'Complete'))


# Load tables of representative genomes of pathogens taken from patric, under 'pathogen' section
files<- list.files('./data/1_patric/1.2_genomes_human_disease/')
files # 24 pathogen genera
patric.tabs<- vector('list')
for(i in 1:length(files)){
  patric.tabs[[i]]<- read.csv(paste0('./data/1_patric/1.2_genomes_human_disease/', files[i]), colClasses=c("Genome.ID"="character")) %>%
    select(Genome.ID,Genome.Name,Genome.Status,Strain,Serovar,Biovar,Pathovar,Type.Strain,Chromosomes,Plasmids,Contigs,PATRIC.CDS,Host.Name,Gram.Stain,Disease, RefSeq.Accessions) %>%
    filter(Genome.Status != 'Plasmid')%>%
    rename(genome_id = Genome.ID) %>%
    left_join(patric.lineages[,c('genome_id', 'species')])
}
patric.tabs.df<- do.call('rbind', patric.tabs)
length(unique(patric.tabs.df$species)) # identifies 112 pathogen species
unique(do.call('rbind', strsplit(unique(patric.tabs.df$species), ' '))[,1]) # 25 genera (Because both Borrelia table gives Borrelia and Borreliella genus)


# That's 2845 genomes, for 112 pathogen species, 25 genera
# We will keep one representative genome per species
# Yet for some genera, many closely related species. Results in large dataset but few information.
# Must do some trimming
# For genera with many closely related species, keep most common pathogens to avoid
strep.remove<- c("Streptococcus infantis","Streptococcus pseudopneumoniae","Streptococcus urinalis","Streptococcus iniae","Streptococcus intermedius","Streptococcus anginosus","Streptococcus dysgalactiae subsp. equisimilis","Streptococcus equi subsp. zooepidemicus","Streptococcus gallolyticus subsp. gallolyticus","Streptococcus suis","Streptococcus sobrinus")
campylo.remove<- c("Campylobacter jejuni subsp. jejuni", "Campylobacter rectus", "Campylobacter upsaliensis")
staph.remove<-c("Staphylococcus arlettae","Staphylococcus aureus","Staphylococcus massiliensis","Staphylococcus saprophyticus subsp. saprophyticus","Staphylococcus capitis","Staphylococcus warneri","Staphylococcus lugdunensis")
burkohl.remove<- c("Burkholderia ambifaria", "Burkholderia glumae")

# Idea with this approach is having to most objective way as possible to select a list of pathogens
# e.g. not bias analysis towards nastiest pathogens
# Yet here thre Streptococcus + one Staphylococcus species show up that are known as typical commensal ...
# I'll remove these manually. Leave everything else in.
rm.commensals<- c('Staphylococcus epidermidis', 'Streptococcus mitis', 'Streptococcus mutans', 'Streptococcus sanguinis')
filter.out<- c(strep.remove, campylo.remove, staph.remove, burkohl.remove, rm.commensals)
patric.tabs.df<- patric.tabs.df[!patric.tabs.df$species %in% filter.out,]


length(unique(patric.tabs.df$species)) # identifies 85 pathogen species
unique(do.call('rbind', strsplit(unique(patric.tabs.df$species), ' '))[,1]) # from 25 genera
nrow(patric.tabs.df) # 1981 genomes

# Now identify a single representative genome for each species
# in MIDAS database, ref_genome corresponds to the type strain of a cluster of closely related genomes
# This can be used as the type strain for each species
genome_info<- read.csv('./data/2_midas_files/genome_info.txt', sep = '\t', colClasses = 'character')
genome_tax<- read.csv('./data/2_midas_files/genome_taxonomy.txt', sep = '\t', colClasses = 'character')
species_info<- read.csv('./data/2_midas_files/species_info.txt', sep = '\t', colClasses = 'character') %>% rename(genome_id = rep_genome_id)

midas_genomes<- genome_info %>%
  select(species_id, genome_name, genome_id, is_rep_genome) %>%
  left_join(species_info) %>%
  left_join(genome_tax[,c('genome_id', 'genus', 'species')], 'genome_id')


# List of pathogen species to look for in MIDAS database
pathogens<- unique(patric.tabs.df$species) # That's 85 pathogen species
pathogens<- gsub('Borreliella', 'Borrelia', pathogens) # Borreliella named Borrelia in MIDAS. Substitute it to ensure grep matching pattern matches MIDAS database


# Look for those in MIDAS database
pathogensMIDAS<- midas_genomes[midas_genomes$species %in% pathogens,] %>%
  arrange(species_id, species) # That's all matching genomes

# Now get the type strain
pathogensMIDAS.reps<- pathogensMIDAS %>%
  #select(species_id, is_rep_genome) %>%
  filter(is_rep_genome == 1)

nrow(pathogensMIDAS.reps)
length(unique(pathogensMIDAS.reps$species))
length(unique(pathogensMIDAS.reps$genus))
head(pathogensMIDAS.reps)


# It's still a lot of genomes, but in fact a lot come from the same species,
# i.e. some species have a lot of genomes sequences and this forms a lot of different strains
# Some species have a lot of distinct strains.
# Yet most of the time those are very rare strains: sometimes it's just a single genome, indicated it's just one rarely isolate strain
as.data.frame(table(pathogensMIDAS.reps$species)) %>% arrange(Freq)
# H. pylori in particular has many strains

# Focus on the most common strain by taking the type strain (representative genome) that has the most individual genomes clustering into it
midas_genomes.pathogens_typeStrains.MostCommon<-
  pathogensMIDAS.reps %>%
  mutate(count_genomes = as.numeric(count_genomes)) %>%
  arrange(species, -count_genomes) %>%
  mutate(first = !duplicated(species)) %>%
  filter(first == TRUE)


midas_genomes.pathogens_typeStrains.MostCommon[,c('species', 'species_id')] %>%
  arrange(species) %>%
  full_join(data.frame(species = pathogens))


# Found  a match for about half of the species
# Some species are trickier to find a match by simple grep pattern but there is a representative genome for them
# manually browsed MIDAS database to found those
# EXCEL FILE "pathogen_commensals_MIDASmatch" gives details about which species a representative genome was manually added, and which were dropped because none could be found.



# NON-PATHOGENS HUMAN COMMENSALS ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(curatedMetagenomicData)
library(stringi)
library(stringr)
library(phytools)
library(phyloseq)


# Code wrapper to process datasets
identify_core_species_HMP<- function(body.site, cutoff, data  = 'HMP_2012.metaphlan_bugs_list.'){
  
  dataset<- paste0('HMP_2012.metaphlan_bugs_list.', body.site)
  
  # select dataset
  suppressPackageStartupMessages(library(curatedMetagenomicData))
  focus <- curatedMetagenomicData(dataset, dryrun = FALSE)
  focus.eset <- focus[[1]] # object contains otu_table(),  sample_metadata(), tax_table()
  
  
  # keep only taxa identified at species level
  df<- exprs( focus.eset )[grep("|s__", fixed = TRUE, rownames( focus.eset)), ] %>%
    as.data.frame()
  
  # add prevalence as number of samples  with non-zero abundance
  df$prev<- rowSums((df > 0))
  
  # extract core set of species: those with prevalence > cutoff, extract species name from the rowname 
  core<- sort(do.call('rbind', strsplit(rownames(df[df$prev >= round(ncol(df)*cutoff),]), fixed = TRUE, '|'))[,7])
  
  
  # Get full taxonomy info from corresponding tax_table()
  tax.table<- data.frame(as(tax_table(ExpressionSet2phyloseq(focus.eset)), "matrix"))
  
  tax.table.filtered<- tax.table[rownames(tax.table) %in% core,] %>%
    rownames_to_column() %>%
    select(Family, Genus, Species)  %>%
    mutate(body_site = body.site) %>%
    arrange(Species)
  
  return(tax.table.filtered)
  
}

# Available datasets:
curatedMetagenomicData("HMP_*metaphlan_bugs_list.*", dryrun = TRUE)

# Extract species present in >80% samples
hmp.stool<- identify_core_species_HMP(body.site = 'stool', cutoff = 0.80)
hmp.nasalcavity<- identify_core_species_HMP(body.site = 'nasalcavity', cutoff = 0.80)
hmp.oralcavity<- identify_core_species_HMP(body.site = 'oralcavity', cutoff = 0.80)
hmp.skin<- identify_core_species_HMP(body.site = 'skin', cutoff = 0.80)
hmp.vagina<- identify_core_species_HMP(body.site = 'vagina', cutoff = 0.80)

hmp.core.species<- rbind(hmp.stool, hmp.nasalcavity, hmp.skin, hmp.vagina, hmp.oralcavity)
table(hmp.core.species$body_site)
length(unique(hmp.core.species$Species)) # 67 unique species
length(unique(hmp.core.species$Genus)) # 33 unique genus

commensals<- gsub('_', ' ', unique(hmp.core.species$Species))

# SOME CLEANUP
# Edit some of the species names to make sure it matches MIDAS by grepping patterns
commensals<- gsub('Clostridium leptum', '[Clostridium] leptum', commensals)
commensals<- gsub('Eubacterium eligens', '[Eubacterium] eligens', commensals)
commensals<- gsub('Eubacterium hallii', '[Eubacterium] hallii', commensals)
commensals<- gsub('Lachnospiraceae bacterium 3 1 46FAA', 'Lachnospiraceae bacterium 3_1_46FAA', commensals)
commensals<- gsub('Lachnospiraceae bacterium 5 1 63FAA', 'Lachnospiraceae bacterium 5_1_63FAA', commensals)
commensals<- gsub('Lachnospiraceae bacterium 7 1 58FAA', 'Lachnospiraceae bacterium 7_1_58FAA', commensals)
commensals<- gsub('Lachnospiraceae bacterium ICM7', 'Lachnoanaerobaculum sp. ICM7', commensals)
commensals<- gsub('Porphyromonas sp oral taxon 279', 'Porphyromonas sp. oral taxon 279', commensals)
commensals<- gsub('Ruminococcus obeum', '[Ruminococcus] obeum', commensals)
commensals<- gsub('Ruminococcus torques', '[Ruminococcus] torques', commensals)
commensals<- gsub('Staphylococcus caprae capitis', 'Staphylococcus capitis', commensals)
commensals<- gsub('Streptococcus mitis oralis pneumoniae', 'Streptococcus mitis', commensals); commensals<- c(commensals, 'Streptococcus oralis') # S. pneumoniae is well known pathogen, leave it out
commensals<- gsub('Subdoligranulum unclassified', 'Subdoligranulum variabile', commensals)

# Streptococcus has A LOT of species.
# Do not want too much data to come from a single genus.
# As for the pathogens, focus on the main well known commensal ones, remove the other
commensals<- commensals[!commensals %in% c('Streptococcus australis', 'Streptococcus gordonii', 'Streptococcus infantis', 'Streptococcus mitis', 'Streptococcus oralis')]

# Remove remaining "unclassified" species + Malassezia globosa which is a fungi
commensals<- commensals[!commensals %in% c('Actinobacillus unclassified','Bilophila unclassified','Capnocytophaga unclassified','Malassezia globosa','Oscillibacter unclassified','Riemerella unclassified','Veillonella unclassified')]

# Agin purpose here is to have a list of commensal species obtained in the most objective way as possible
# Yet a few well known commensal do not show up here which I know were present and prevalent in the HMP healthy cohort dataset from Kin slection paper (Simonet & McNally, 2021)
# Adding those manually
commensals<- c(commensals, c('Akkermansia muciniphila', 'Bifidobacterium longum subsp. longum', 'Bacteroides thetaiotaomicron', 'Acidaminococcus intestini'))



# We have a list of commensals.
# Again, use MIDAS database to identify representative genome for each, and if several strains for a single species, take the most common strain
commensalsMIDAS<- midas_genomes[midas_genomes$species %in% commensals,] %>%
  arrange(species_id, species)

# It's a lot of genomes, so let's focus on representative genomes, i.e. the type strain
commensalsMIDAS.reps<- commensalsMIDAS %>%
  #select(species_id, is_rep_genome) %>%
  filter(is_rep_genome == 1)

# As before:
# It's still a lot of genomes, but in fact a lot come from the same species,
# i.e. some species have a lot of genomes sequences and this forms a lot of different strains
# Some species have a lot of distinct strains.
# Yet most of the time those are very rare strains: sometimes it's just a single genome, indicated it's just one rarely isolate strain
as.data.frame(table(commensalsMIDAS.reps$species)) %>% arrange(Freq)


# Focus on the most common strain by taking the type strain (representative genome) that has the most individual genomes clustering into it
midas_genomes.commensals_typeStrains.MostCommon<-
  commensalsMIDAS.reps %>%
  mutate(count_genomes = as.numeric(count_genomes)) %>%
  arrange(species, -count_genomes) %>%
  mutate(first = !duplicated(species)) %>%
  filter(first == TRUE)



midas_genomes.commensals_typeStrains.MostCommon[,c('species', 'species_id')] %>%
  arrange(species) %>%
  full_join(data.frame(species = commensals))


# Found  a match for  almost all the species
# Again, completed by manually searching MIDAS to check if there were a relevant match.
# details of added genome in excel file "pathogen_commensals_MIDASmatch.xlsx"



# WRITE TABLE FOR GENOMES PROCESSING ----

library(readxl)

midasmatch.pathogens<- read_excel('./output/1_processed_tables/1.0_pathogen_commensals_MIDASmatch.xlsx', sheet = 'pathogen') %>%
  select(species, species_id, pathogen) %>%
  filter(!is.na(species_id))


midasmatch.commensals<- read_excel('./output/1_processed_tables/1.0_pathogen_commensals_MIDASmatch.xlsx', sheet = 'commensal') %>%
  select(species, species_id, pathogen) %>%
  filter(!is.na(species_id))


nrow(midasmatch.pathogens)
length(unique(midasmatch.pathogens$species))
length(unique(midasmatch.pathogens$species_id))


nrow(midasmatch.commensals)
length(unique(midasmatch.commensals$species))
length(unique(midasmatch.commensals$species_id))

# two species are both in pathogen and non-pathogen list
sum(midasmatch.pathogens$species_id %in% midasmatch.commensals$species_id)

intersect(midasmatch.pathogens$species_id, midasmatch.commensals$species_id)
# These are two species of Campylobacter.
# However keeping those as non-pathogen, because can rarely cause infection but most of the time part of normal mouth microbiome
# In contrast the other three Campylobacter
# Campylobacter_coli_58237
# Campylobacter_fetus_53261
# Campylobacter_jejuni_57666
# Are very common pathogens

midasmatch.pathogens<- midasmatch.pathogens[-which(midasmatch.pathogens$species_id %in% intersect(midasmatch.pathogens$species_id, midasmatch.commensals$species_id)
),]

nrow(midasmatch.pathogens)
nrow(midasmatch.commensals)

# we end up with 59  species of each, so 118 species total


genome_info<- read.csv('./data/2_midas_files/genome_info.txt', sep = '\t', colClasses = 'character')
genome_tax<- read.csv('./data/2_midas_files/genome_taxonomy.txt', sep = '\t', colClasses = 'character')
species_info<- read.csv('./data/2_midas_files/species_info.txt', sep = '\t', colClasses = 'character') %>% rename(genome_id = rep_genome_id)

midas_genomes<- genome_info %>%
  select(species_id, genome_name, genome_id, is_rep_genome) %>%
  left_join(species_info) %>%
  left_join(genome_tax[,c('genome_id', 'genus', 'species')], 'genome_id')

midas.selection<- c(midasmatch.pathogens$species_id, midasmatch.commensals$species_id)

midas.selection.sp.info<- species_info[species_info$species_id %in% midas.selection,] %>% arrange(genome_id)
midas.selection.genome.info<- genome_info[genome_info$genome_id %in% midas.selection.sp.info$genome_id,]%>% arrange(genome_id)
#midas.selection.genome.tax<- genome_tax[genome_tax$genome_id %in% midas.selection.sp.info$genome_id,]

assembled<- cbind(midas.selection.sp.info, midas.selection.genome.info[,c('genome_name', 'is_rep_genome')])

assembled2<- assembled %>%
  left_join(genome_tax[,c('genome_id', 'genus', 'species')]) %>%
  left_join(rbind(midasmatch.pathogens, midasmatch.commensals)[,c('species_id', 'pathogen')]) %>%
  select(species_id, genome_name, genome_id, is_rep_genome, count_genomes, genus, species, pathogen)


# Add gram profiles
grams<- read.table('./data/0_species_info_files/grams_stains.txt', header=TRUE, sep = '\t')
pathogen_commensal_genomes<- left_join(assembled2, grams) %>%
  arrange(pathogen, species_id, species)

# Some gram profiles missing in the gram_stain table. Add manually
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Acidaminococcus_intestini_54097",'gram_profile']<- 'n'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Akkermansia_muciniphila_55290",'gram_profile']<- 'n'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Anaerostipes_hadrus_55206",'gram_profile']<- 'p'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Bacteroides_thetaiotaomicron_56941",'gram_profile']<- 'n'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Bifidobacterium_longum_57796",'gram_profile']<- 'p'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Clostridiales_bacterium_56470",'gram_profile']<- 'p'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Lachnospiraceae_bacterium_51870",'gram_profile']<- 'p'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Porphyromonas_sp_57899",'gram_profile']<- 'n'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Prevotella_nanceiensis_44721",'gram_profile']<- 'n'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Bacillus_cereus_57918",'gram_profile']<- 'p'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Campylobacter_fetus_53261",'gram_profile']<- 'n'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Helicobacter_bilis_58212",'gram_profile']<- 'n'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Mycobacterium_avium_56142",'gram_profile']<- 'OM+'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Mycobacterium_marinum_56420",'gram_profile']<- 'OM+'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Rickettsia_rickettsii_57256",'gram_profile']<- 'n'
pathogen_commensal_genomes[pathogen_commensal_genomes$species_id == "Vibrio_owensii_57591",'gram_profile']<- 'n'


# some genomesin PATRICdb have changed genome  id since MIDAS was published (or it is errors in MIDAS db)
# fix those manually

pathogen_commensal_genomes[which(pathogen_commensal_genomes$genome_id == '1313.2455'),'genome_id']<- '1313.10903'

# write it
write.table(pathogen_commensal_genomes, file = './output/1_processed_tables/1.1_pathogen_commensal_genomes_118.txt', quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')


