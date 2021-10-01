# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# SCRIPT FOR FIGURES AND TABLES #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#local_project_dir=PATH/TO/CLONE/REPO
#local_project_dir='/Users/s1687811/Documents/PhD/Research/CooperativePathogenicityVirulence/'

setwd(local_project_dir)
setwd('./CooperativePathogenicityVirulence_repo/')
source('./scripts/sourced_packages.R')
source('./scripts/sourced_ggthemes.R')
library(ape)
library(MCMCglmm)
library(dmetar)
library(esc)
library(ggtree)
library(grid)
library(gridExtra)
library(scales)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(gtable)
library(ggpubr)
library(aplot)

# Figure 1 (web of science) ----

# table with assembled tables of the three Web of Science searches
d<- read_excel('./data/6_webofscience_search/WebOfScience_microbial_cooperation.xlsx', sheet = 1)
d<- d[!is.na(d$`Publication Years`),]

colnames(d)<- c('year', 'microbe_general', 'cooperation', 'coop_and_health', 'prop_cooperation', 'MA5Y_prop_cooperation')

d$coop_and_health[which(d$coop_and_health == 0)]<- NA # for bubbles to not appear on plot if n = 0

d$prop_coop_perc<- d$prop_cooperation*100 # get % instead of proprotions 
d$MA5Y_prop_cooperation_perc<- d$MA5Y_prop_cooperation*100 # get % instead of proprotions 

# rename to get neat legend title
d<- d %>% rename(`Number of papers \nmentionning \nhealth impact` = coop_and_health)


# main plot
p<- ggplot(d, aes(x = year, y = prop_coop_perc, col))+
  ylab('Percentage of research on microbial cooperation')+
  xlab('Year')+
  xlim(2000, 2021)+
  geom_line(col = 'grey', size = 0.3)+
  geom_line(aes(x = year, y = MA5Y_prop_cooperation_perc), size = 0.5)+
  geom_point(aes(x = year, y = c(d$MA5Y_prop_cooperation_perc[1:20],d$MA5Y_prop_cooperation_perc[c(20,20)]), size = `Number of papers \nmentionning \nhealth impact`), shape = 21, fill = 'dodgerblue', col = 'black', alpha = .5)+
  theme_bw()+
  theme_bw()+
  scale_size_continuous(breaks = c(1, 2, 4, 8))+
  scale_color_manual(values = c('black', 'grey'))+
  theme(axis.title = element_text(face = 'bold', size = 8),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.3),
        axis.ticks = element_line(colour = "black", size = 0.3),
        axis.text = element_text(colour = "black", size = 7),
        #legend.position = 'none', #c(0.9, 0.1),
        legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5)
  )

p

# Make plot to extract lines legend
p.legend.lines<- d %>%
  select(year, prop_coop_perc, MA5Y_prop_cooperation_perc) %>%
  rename(`per year` = prop_coop_perc,
         `5-year moving average` = MA5Y_prop_cooperation_perc) %>%
  gather('measure', 'value', 2:3) %>%
  ggplot(., aes(x = year, y = value, col = measure, size = measure))+
  geom_line()+
  scale_size_manual(values=c(`per year` = 0.3,`5-year moving average` = 0.4))+
  scale_color_manual(values = c('black', 'grey'))+
  theme_bw()+
  theme(legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_blank())

legend.lines <- cowplot::get_legend(p.legend.lines)


# p+ annotate("text", x = 2014, y = 0.046,
#             label = "First principles of\nHamiltonian medicine\n(Crespi, 2005)",
#             size = 2, face = 'bold', hjust = 0, fontface = 'italic')+
#   annotate("segment", x = 2014, xend = 2014, y = 0.05, yend = 0.05+0.013, size=0.4, alpha=0.6, arrow=arrow(length = unit(0.15,"cm")))+
#   annotate("text", x = 2005, y = 0.014,
#            label = "Hamiltonian Medicine:\nWhy the Social Lives\nof Pathogens Matter\n(Foster, 2005)",
#            size = 2, face = 'bold', hjust = 0, fontface = 'italic')+
#   annotate("segment", x = 2005, xend = 2005, y = 0.015, yend = 0.015+0.008, size=0.4, alpha=0.6, arrow=arrow(length = unit(0.15,"cm")))



# output figure
pdf('./output/2_figures/Figure1_WebOfScience.pdf', width = 8.2/2.55, height = 8.6/2.55)
p+
  theme(legend.position = 'top',
        legend.justification = c("left"))+
  annotation_custom(grob = legend.lines, xmax = 2007, ymin = 0.08)

dev.off()



# Figure 2 ----

library(ggtree)
library(viridis)

# 2A (Phylogeny) ----

df<- read.table('./output/1_processed_tables/2.1_assembled_SPECIES_annotation.txt', header=TRUE, sep = '\t',
                colClasses = c(rep('character', 4),
                               rep('numeric', 9)))

# Make shorter names for plot
df$species_id_short<- gsub('_', ' ', substr(df$species_id,1,nchar(df$species_id)-6))
df$species_id_short<- gsub('Streptococcus','Streptoco.', df$species_id_short)
df$species_id_short<- gsub('bacterium ','bact. ', df$species_id_short)
df$species_id_short<- gsub('Staphylococcus ','Staphyloco.', df$species_id_short)
df$species_id_short<- gsub('Actinomyces ','Actinom.', df$species_id_short)
df$species_id_short<- gsub('thetaiotaomicron','thetaiota.', df$species_id_short)
df$species_id_short<- gsub('Chlamydophila pneumoniae','Chlamydophila pneumo.', df$species_id_short)
df$species_id_short<- gsub('Burkholderia','Burkhold.', df$species_id_short)
df$species_id_short<- gsub('Pseudomonas','Pseudom.', df$species_id_short)
df$species_id_short<- gsub('Acinetobacter','Acineto.', df$species_id_short)
df$species_id_short<- gsub('Haemophilus ','Haemoph. ', df$species_id_short)


# Fix some species names to avoid confusion in phylogeny plot due to strains
df[which(df$species_id == 'Acinetobacter_baumannii_56211'),'species_id_short']<- 'Acineto. baumannii (56211)'
df[which(df$species_id == 'Acinetobacter_baumannii_57014'),'species_id_short']<- 'Acineto. baumannii (57014)'
df[which(df$species_id == 'Bacillus_cereus_57918'),'species_id_short']<- 'Bacillus cereus (57918)'
df[which(df$species_id == 'Bacillus_cereus_57982'),'species_id_short']<- 'Bacillus cereus (57982)'
df[which(df$species_id == 'Mycobacterium_tuberculosis_57144'),'species_id_short']<-'Mycobact. colombiense'#df[which(df$species_id == 'Mycobacterium_tuberculosis_57144'),'species']
df[which(df$species_id == 'Porphyromonas_sp_57899'),'species_id_short']<- "Porphyromonas KLE1280" #df[which(df$species_id == 'Porphyromonas_sp_57899'),'species']
df[which(df$species_id == 'Lachnospiraceae_bacterium_62587'),'species_id_short']<- 'Lachno. bact. ICM7' # df[which(df$species_id == 'Lachnospiraceae_bacterium_62587'),]
df[which(df$species_id == 'Lachnospiraceae_bacterium_51870'),'species_id_short']<- 'Lachno. bact. 8.1.57FAA' # df[which(df$species_id == 'Lachnospiraceae_bacterium_51870'),]


df<- df %>%
  mutate(species_id2 = species_id,
         species_id = species_id_short) %>%
  rename(id = species_id,
         ab_degradation = antibiotic_degradation,
         vf = is_victor_vf,
         nb_extracellular = secretome) %>%
  mutate(#species_id_short  = id,
    pathogen = as.factor(pathogen)) %>%
  as.data.frame()


midas.tree<- read.tree('./data/5_phylogeny_files/midas_tree_renamed.newick')
phylogeny<- drop.tip(midas.tree, midas.tree$tip.label[which(!midas.tree$tip.label %in% df$species_id2)])
#phylogeny<- chronopl(phylogeny, lambda = 0)
#plot(phylogeny)
new.tipslabs<- df$id[match(phylogeny$tip.label, df$species_id2)]
phylogeny$tip.label<- new.tipslabs

# MAKE THE TREE
p<- ggtree(phylogeny, size = 0.3) %<+% df + 
  geom_tiplab(aes(label=species_id_short, col = pathogen),
              align=T,
              offset=0, hjust=0, linesize = 0.2, size = 1.5) +
  scale_color_manual(values = c('dodgerblue','darkorange'))+
  theme(legend.position = 'none')+
  xlim(0,8)
p




# 2E (boxplot variation) ----

coop.data<- 
  df %>%
  select(-genus, -species, -gram_profile, -pathogen, -total_cds, -species_id_short, -species_id2) %>%
  rename(label = id) %>%
  gather('category', 'value', 2:8)
head(coop.data)


coop.data.props<- 
  df %>%
  select(id, total_cds) %>%
  rename(label = id) %>%
  left_join(coop.data, by = 'label') %>%
  select(label, category, value, total_cds) %>%
  mutate(prop = value/total_cds) %>%
  mutate(perc = prop*100)


# Having a look at how each form of cooperation varies
hist(coop.data.props[which(coop.data.props$category == 'nb_extracellular'),'value'])
summary(coop.data.props[which(coop.data.props$category == 'nb_extracellular'),'value'])
var(coop.data.props[which(coop.data.props$category == 'siderophores'),'value'])
range(coop.data.props[which(coop.data.props$category == 'siderophores'),'perc'])


# renaming the categories for figure
unique(coop.data.props$category)
coop.data.props$category[which(coop.data.props$category == 'ab_degradation')]<- 'Antib. degr.'
coop.data.props$category[which(coop.data.props$category == 'biofilm')]<- 'Biofilm'
coop.data.props$category[which(coop.data.props$category == 'quorum_sensing')]<- 'QS'
coop.data.props$category[which(coop.data.props$category == 'siderophores')]<- 'Siderophores'
coop.data.props$category[which(coop.data.props$category == 'secretion_system')]<- 'Secr. Syst.'
coop.data.props$category[which(coop.data.props$category == 'nb_extracellular')]<- 'Secretome'
coop.data.props$category[which(coop.data.props$category == 'vf')]<- 'VF'

unique(coop.data.props$category)
coop.data.props$category<- factor(coop.data.props$category,
                                  levels = c('VF', 'QS', 'Biofilm', 'Siderophores', 'Antib. degr.', 'Secr. Syst.', 'Secretome'))



pdf('./output/2_figures/Figure2_DataOverview_E.pdf', width = 4.04/2.55, height = 3.38/2.55)

coop.data.props %>%
  filter(category != 'VF') %>%
  ggplot(., aes(y = category, x = perc))+
  geom_boxplot(width = 0.5, outlier.size = 0.05, lwd = 0.15)+
  ylab('Cooperation category')+xlab('Percentage of CDS coding for\ncooperation category')+
  theme_bw()+
  theme(axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 5, colour = 'black'),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.grid = element_blank(),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(colour = 'black', size = 0.3),
        plot.margin=unit(c(0,0,0,0), "null")
  )

dev.off()



# 2B (Heatmap) ----

plotHeatmap<- function(focal.trait, pal.type = 'viridis', leg.pos = 'none', margins = c(0,-0.09,0,-0.09)){
  
  plot<- coop.data.props %>%
    filter(category %in% focal.trait) %>%
    ggplot(., aes(x = category, y = label, fill = perc))+
    geom_tile()+
    scale_fill_viridis(option=pal.type)+
    scale_x_discrete(position = "top")+
    theme(axis.text.x = element_text(angle = 45, size = 5, hjust = 0))+
    xlab('')+ylab('')+
    #theme_bw()+
    # theme(axis.ticks.y = element_blank(),
    #       axis.text.y = element_blank(),
    #       panel.grid = element_blank(),
    #       panel.border = element_blank(),
    #       panel.background = element_blank())+
    theme_void()+
    theme(legend.position = leg.pos)+
    theme(plot.margin=unit(margins, "null"))
  
  return(plot)
  
}
plotHeatmapLegends<- function(focal.trait, pal.type = 'viridis', leg.pos = 'right', margins = c(0.1,0.1,0.1,0.1), breaks){
  
  plot<- coop.data.props %>%
    filter(category %in% focal.trait) %>%
    ggplot(., aes(x = category, y = label, fill = perc))+
    geom_tile()+
    scale_fill_viridis(option=pal.type, breaks = breaks, name = '')+
    theme_void()+
    theme(legend.position = leg.pos)+
    theme(plot.margin=unit(margins, "null"),
          legend.key.height = unit(0.2, "cm"),
          legend.key.width = unit(0.2, "cm"),
          legend.spacing.x = unit(0.04, 'cm'),
          #legend.spacing = unit(0, "cm"),
          legend.text = element_text(size = 4))
  
  return(plot)
  
}


l1<- get_legend(plotHeatmapLegends('Secr. Syst.', breaks = c(0, 0.5, 1.03)))
l2<- get_legend(plotHeatmapLegends('Secretome', breaks = c(0.2, 1.8, 3.5)))
l3<- get_legend(plotHeatmapLegends('Siderophores', breaks = c(0, 0.3, 0.58)))
l4<- get_legend(plotHeatmapLegends('Biofilm', breaks = c(0, 0.6, 1.15)))
l5<- get_legend(plotHeatmapLegends('Antib. degr.', breaks = c(0, 0.25, 0.47)))
l6<- get_legend(plotHeatmapLegends('QS', breaks = c(0, 0.13, 0.25)))
l7<- get_legend(plotHeatmapLegends('VF', breaks = c(0, 5, 10, 15)))


#pdf('~/Desktop/HM1Figs/legend.bars.pdf', width = 3.9/2.55, height = 1.88/2.55)
#grid.arrange(l1, l2, l3, l4, l5, l6, l7, ncol = 7)
#dev.off()


# (heatmap plot output in next section with %VF & composition)


# 2C (VF vs nVF Coop, per species) ----

# Get and format the annotation data

# Annotations summarised at species level
d.annot<- read.csv('./output/1_processed_tables/2.2_assembled_GENES_annotation.txt', header=TRUE, sep = '\t',
                   colClasses = c(rep('character', 7),
                                  rep('numeric', 8)))

sum(is.na(d.annot$secretome)) # all genes are annotated

d.annot2<- d.annot %>%
  filter(!is.na(secretome)) %>% # in fact no case of this
  mutate(is_coop = biofilm+antibiotic_degradation+quorum_sensing+siderophores+secretion_system+secretome)%>%
  mutate(is_coop = ifelse(is_coop == 0, is_coop, 1)) %>%
  select(-genus, -species, -product_patric, -total_cds) %>%
  rename(is_vf = is_victor_vf,
         species = species_id)

# Code wrapper to make the 2x2 table for each species for a given focal_trait

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

data.overall.cooperative.proportions <- 
  get_conttable(focal_trait = 'is_coop') %>%
  mutate(propCoopVF = coop1_vf1/(coop1_vf1+coop0_vf1),
         propCoopNonVF = coop1_vf0/(coop1_vf0+coop0_vf0)) %>%
  rename(species_id2 = species) %>%
  left_join(df[,c('species_id2', 'id')], by = 'species_id2') %>%
  select(id, propCoopVF, propCoopNonVF, species_id2, coop0_vf0, coop0_vf1, coop1_vf0, coop1_vf1) %>%
  rename(label = id)


# Replace proportion of 1 for Campylobacter showae, by 0.55
# Plot the bar up to 0.55 but replace 0.55 axis graduation by 1
# and will make the axis break display manually when assembling the figure
max(data.overall.cooperative.proportions$propCoopVF, na.rm = TRUE)
data.overall.cooperative.proportions[which(data.overall.cooperative.proportions$propCoopVF == 1),]
data.overall.cooperative.proportions$propCoopVF[which(data.overall.cooperative.proportions$propCoopVF == 1)]<- 0.55

# Separate barplots, on for VF one for non-VF
p.barplot.cooperative.prop.VF<-
  ggplot(data.overall.cooperative.proportions, aes(x = propCoopVF, y = label))+
  geom_bar(stat="identity", fill = 'grey')+
  theme_tree2()+
  xlim(0, 0.55)+
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.55), labels = c(0, 0.1, 0.2, 0.3, 0.4, 1))+
  theme(axis.text.x = element_text(size = 5, colour = 'black', angle = 45, hjust = 1),
        axis.ticks.x = element_line(colour = 'black', size = 0.3))


p.barplot.cooperative.prop.NonVF<-
  ggplot(data.overall.cooperative.proportions, aes(x = propCoopNonVF, y = label))+
  geom_bar(stat="identity", fill = 'black')+
  theme_tree2()+
  theme(axis.text.x = element_text(size = 5, colour = 'black', angle = 45, hjust = 1),
        axis.ticks.x = element_line(colour = 'black', size = 0.3))

# Single barplot with both on same scale
p.barplot.single<- 
  data.overall.cooperative.proportions %>%
  gather('category', 'prop', 2:3) %>%
  #rename(label = id) %>%
  ggplot(., aes(x = prop, y = label, fill = category))+
  #xlim(0, 0.5)+
  geom_bar(stat = 'identity', position = 'dodge')+
  scale_fill_manual(values = c('black', 'grey'))+
  theme(#legend.justification = 'top',
    axis.text.y = element_blank())+
  theme_tree2()+
  theme(legend.position = 'none')+
  xlim(0, 0.55)+
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.55), labels = c(0, 0.1, 0.2, 0.3, 0.4, 1))+
  #annotate('rect', xmin = 0.52, xmax = 0.53, ymin = -Inf, ymax = 50, col = 'red')
  theme(axis.text.x = element_text(size = 7, colour = 'black', angle = 0),
        axis.ticks.x = element_line(colour = 'black', size = 0.3))



# 2D (VF composition, per species) ----

# BARPLOT of VF, all genes (cooperative and non-cooperative)

dat.barplot.per.species.vf<- 
  d.annot2 %>%
  mutate(code = paste0(biofilm, antibiotic_degradation,quorum_sensing, siderophores, secretion_system, secretome, is_vf)) %>%
  filter(is_vf == 1) %>%
  select(species, gram_profile, pathogen, code) %>%
  arrange(species) %>%
  group_by(species,code) %>%
  summarise(number = n()) %>%
  rename(species_id2 = species) %>%
  left_join(df[,c('species_id2', 'id')], by = 'species_id2') %>%
  ungroup() %>%
  select(id, code, number)

pal.coop<- c('lightgrey', '#E69F00', '#F0E442', '#D55E00', '#CC79A7', '#009E73', '#56B4E9', 'firebrick')


barplot.species.vf<- ggplot(dat.barplot.per.species.vf, aes(fill=code, y=number, x=id)) + 
  #geom_bar(stat="identity")+
  geom_bar(position="fill", stat="identity")+
  coord_flip()+
  theme_bw()+
  scale_fill_manual(
    'Virulence factor\nannotation',
    values = pal.coop,
    labels = c('non-cooperative', 'secretome', 'secretion syst.', 'siderophore',
               'quorum-sensing', 'antibiotic degr.', 'biofilm', 'biofilm & secretome'))+
  theme(legend.justification = 'top',
        axis.text.y = element_blank())+
  theme_tree2()


# 2F (BARPLOT of NON-VF, all genes) ----

dat.barplot.per.species.NONvf<- 
  d.annot2 %>%
  mutate(code = paste0(biofilm, antibiotic_degradation,quorum_sensing, siderophores, secretion_system, secretome, is_vf)) %>%
  filter(is_vf == 0) %>%
  select(species, gram_profile, pathogen, code) %>%
  arrange(species) %>%
  group_by(species,code) %>%
  summarise(number = n()) %>%
  rename(species_id2 = species) %>%
  left_join(df[,c('species_id2', 'id')], by = 'species_id2') %>%
  ungroup() %>%
  select(id, code, number)


dat.barplot.per.species.NONvf$code<- factor(dat.barplot.per.species.NONvf$code,
                                            levels = c("0000000","0000010","0000100","0001000","0010000","0100000","1000000","0000110","0100010","1000010","1000100"))

pal.coop.2<- c('lightgrey','#E69F00', '#F0E442', '#D55E00', '#CC79A7', '#009E73', '#56B4E9', 'firebrick', 'black', 'green', 'purple')

barplot.species.NONvf<- ggplot(dat.barplot.per.species.NONvf, aes(fill=code, y=number, x=id)) + 
  #geom_bar(stat="identity")+
  geom_bar(position="fill", stat="identity")+
  coord_flip()+
  theme_bw()+
  scale_fill_manual(
    'Virulence factor\nannotation',
    values = pal.coop.2,
    labels = c('non-cooperative', 'secretome', 'secretion syst.', 'siderophore',
               'quorum-sensing', 'antibiotic degr.', 'biofilm',
               "secretome & secretion syst.",
               "secretome & antibiotic degr.",
               'biofilm & secretome',
               'biofilm & secretion syst.'))+
  theme(legend.justification = 'top',
        axis.text.y = element_blank())+
  theme_tree2()


# 2G (BARPLOTS OF VF and non-VF, cooperative genes only) ----

barplot.species.vf.CoopOnly<-
  dat.barplot.per.species.vf %>% filter(code != '0000001') %>%
  ggplot(., aes(fill=code, y=number, x=id)) + 
  #geom_bar(stat="identity")+
  geom_bar(position="fill", stat="identity")+
  coord_flip()+
  theme_bw()+
  scale_fill_manual(
    'Virulence factor\nannotation',
    values = pal.coop[-1],
    labels = c(#'non-cooperative',
      'secretome', 'secretion syst.', 'siderophore',
      'quorum-sensing', 'antibiotic degr.', 'biofilm', 'biofilm & secretome'))+
  theme_tree2()+
  theme(legend.justification = 'top',
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7, colour = 'black', angle = 0),
        axis.ticks.x = element_line(colour = 'black', size = 0.3))+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, '', 0.5, '', 1))


barplot.species.NONvf.CoopOnly<- 
  dat.barplot.per.species.NONvf %>% filter(code != '0000000') %>%
  ggplot(., aes(fill=code, y=number, x=id)) + 
  #geom_bar(stat="identity")+
  geom_bar(position="fill", stat="identity")+
  coord_flip()+
  theme_bw()+
  scale_fill_manual(
    'Virulence factor\nannotation',
    values = pal.coop.2[-1],
    labels = c(#'non-cooperative',
      'secretome', 'secretion syst.', 'siderophore',
      'quorum-sensing', 'antibiotic degr.', 'biofilm',
      "secretome & secretion syst.",
      "secretome & antibiotic degr.",
      'biofilm & secretome',
      'biofilm & secretion syst.'))+
  theme_tree2()+
  theme(legend.justification = 'top',
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7, colour = 'black', angle = 0),
        axis.ticks.x = element_line(colour = 'black', size = 0.3))+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, '', 0.5, '', 1))


# Assemble ----

# PART 1, sizing tailored to get phylogeny and heatmap
pdf('./output/2_figures/Figure2_DataOverview_A_B.pdf', width = 20/2.55, height = 17/2.55) # 20, 17

plotHeatmap('Secretome') %>%
  insert_left(p, width = 10) %>% # 10
  insert_right(plotHeatmap('Secr. Syst.')) %>%
  insert_right(plotHeatmap('Antib. degr.')) %>%
  insert_right(plotHeatmap('Siderophores')) %>%
  insert_right(plotHeatmap('Biofilm')) %>%
  insert_right(plotHeatmap('QS')) %>%
  insert_right(plotHeatmap('VF', pal.type = 'viridis')) %>%
  insert_right(p.barplot.single, width = 3) %>% # 3
  insert_right(barplot.species.vf.CoopOnly+theme(legend.position = 'none'), width = 2) %>% # 2
  insert_right(barplot.species.NONvf.CoopOnly+theme(legend.position = 'none'), width = 2) # 2

dev.off()

# PART 2, sizing tailored to get VF composition data
pdf('./output/2_figures/Figure2_DataOverview_C_D.pdf', width = 20/2.55, height = 17/2.55)

  plotHeatmap('Secretome') %>%
  insert_left(p) %>%
  insert_right(p.barplot.single) %>%
  insert_right(barplot.species.vf.CoopOnly+theme(legend.position = 'none'), width = 0.5) %>%
  insert_right(barplot.species.NONvf.CoopOnly+theme(legend.position = 'none'), width = 0.5)

dev.off() 


# 2F (VF vs non-VF Coop, overall) ----

tab.props.overall<- 
  d.annot2 %>% 
  select(is_coop, is_vf)%>%
  rename(cooperative = 1) %>%
  group_by(cooperative, is_vf) %>%
  table()

p.vf<- tab.props.overall[2,2]/sum(tab.props.overall[,2]) # Proportion of VF that are cooperative
p.nonVf<- tab.props.overall[2,1]/sum(tab.props.overall[,1]) # Proportion of non-VF that are cooperative

data.propCoop.Overall<- data.frame(category = c('Non-VF', 'VF'),
                                   PropsCooperative = c(p.nonVf, p.vf))



pdf('./output/2_figures/Figure2_DataOverview_F.pdf', width = 3.35/2.55, height = 3.38/2.55)


ggplot(data.propCoop.Overall, aes(y = category, x = PropsCooperative, fill = category))+
  geom_bar(stat = 'identity', width = 0.5)+
  xlab('Proportion of CDS which \nare cooperative')+ylab('')+
  scale_fill_manual(values = c('black', 'grey'))+
  scale_x_continuous(limits = c(0, 0.055), expand = c(0, 0))+
  theme_bw()+
  theme(axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 5, colour = 'black', angle = 90,  hjust=0.5),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.grid = element_blank(),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(colour = 'black', size = 0.3),
        plot.margin=unit(c(0,0,0,0), "null")
        )
  
dev.off()


# 2G (VF composition, overall) ----

d.annot2$biofilm2<- ifelse(d.annot2$biofilm == 0, '', 'biofilm')
d.annot2$quorum_sensing2<- ifelse(d.annot2$quorum_sensing == 0, '', 'qs')
d.annot2$siderophores2<- ifelse(d.annot2$siderophores == 0, '', 'sid')
d.annot2$antibiotic_degradation2<- ifelse(d.annot2$antibiotic_degradation == 0, '', 'ab')
d.annot2$secretion_system2<- ifelse(d.annot2$secretion_system == 0, '', 'SecSyst')
d.annot2$secretome2<- ifelse(d.annot2$secretome == 0, '', 'secretome')
d.annot2$VF2<- ifelse(d.annot2$is_vf == 0, '', 'VF')

#foo<- d.annot2[d.annot2$is_vf == 1,]
foo<- d.annot2[d.annot2$is_vf %in% c(0,1),]

foo$profile<- 
  paste(foo$VF2, foo$secretome2, foo$secretion_system2, foo$antibiotic_degradation2, foo$siderophores2, foo$quorum_sensing2, foo$biofilm2)
unique(foo$profile)

foo.tab<- as.data.frame(table(foo$profile))


# remove the first and 12th row which is non-cooperative and non-VF
foo.tab<- foo.tab[!foo.tab$Var1 %in% c('VF      ', "      "),]
foo.tab$category<- c(rep('non-VF', 10), rep('VF', 7))
foo.tab$cooperative_trait<- c("biofilm",
                              "quorum-sensing",
                              "siderophore",
                              "antibiotic-degration",
                              "secretion system",
                              "secretion system & biofilm",
                              "secretome",
                              "secretome & biofilm",
                              "secretome & antibiotic-degration",
                              "secretome & secretion system",
                              "biofilm",
                              "quorum-sensing",
                              "siderophore",
                              "antibiotic-degration",
                              "secretion system",
                              "secretome",
                              "secretome & biofilm")

foo.tab$cooperative_trait<- factor(foo.tab$cooperative_trait,
                                   levels = c("secretome", "secretion system", "siderophore", "quorum-sensing", "antibiotic-degration", "biofilm",
                                              "secretome & secretion system", "secretome & antibiotic-degration", "secretome & biofilm", "secretion system & biofilm"))

library(viridis)
#foo.tab$foo = 1

pal.coop<- c('#E69F00', '#F0E442', '#D55E00', '#CC79A7', '#009E73', '#56B4E9', 'firebrick', 'black', 'green', 'purple')


pdf('./output/2_figures/Figure2_DataOverview_G.pdf', width = 3.35/2.55, height = 3.38/2.55)

ggplot(foo.tab, aes(y = category, x = Freq, fill = cooperative_trait))+
  geom_bar(stat = "identity", position = "fill", width = 0.5)+
  xlab('Proportion of cooperative CDS\nof each cooperation type')+ylab('')+
  scale_x_continuous(expand = c(0, 0))+
  scale_fill_manual(
    'Annotation',
    values = pal.coop,
    labels = c(
      'secretome',
      'secretion syst.',
      'siderophore',
      'quorum-sensing',
      'antibiotic degr.',
      'biofilm',
      'secretome & secretion. syst',
      'secretome & antibiotic degr.',
      'biofilm & secretome',
      'biofilm & secretion syst'))+
  theme(legend.justification = 'top',
        legend.position = 'none')+
  theme_bw()+
  theme(axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 5, colour = 'black', angle = 90,  hjust=0.5),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.grid = element_blank(),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(colour = 'black', size = 0.3),
        plot.margin=unit(c(0,0,0,0), "null")
  )


dev.off()



pdf('./output/2_figures/Figure2_DataOverview_legend.pdf', width = 6/2.55, height = 6/2.55)

ggplot(foo.tab, aes(x = category, y = Freq, fill = cooperative_trait))+
  geom_bar(stat = "identity", position = "fill")+
  ylab('Proportion of cooperative genes\nof each type')+xlab('')+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(
    'Annotation',
    values = pal.coop,
    labels = c(
      'secretome',
      'secr. syst.',
      'siderophore',
      'quorum-sensing',
      'antib. degr.',
      'biofilm',
      'secretome & secr. syst.',
      'secretome & antib. degr.',
      'biofilm & secretome',
      'biofilm & secr. syst.'))+
  theme(legend.justification = 'top',
        legend.position = 'none')+
theme_bw()+
  theme(axis.title = element_text(size = 5),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.grid = element_line(size = 0),
        legend.position = 'bottom', #c(0.9, 0.1),
        legend.justification = 'left',
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(colour = 'black', size = 0.3),
        plot.margin=unit(c(0.025,0,0,0), "null"))+
  guides(fill=guide_legend(nrow=5))

dev.off()



# Figure 3 ----

load("./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN1.RData")

# ... (A) observed OR ----

df.all2<- left_join(df.all, unique(d.annot2[,c('species', 'pathogen')]), by = 'species')

unique(df.all2$cooperative_trait)
head(df.all2)


# renaming the categories for figure
unique(df.all2$cooperative_trait)
df.all2$cooperative_trait[which(df.all2$cooperative_trait == 'ab')]<- 'Antib. degr.'
df.all2$cooperative_trait[which(df.all2$cooperative_trait == 'biofilm')]<- 'Biofilm'
df.all2$cooperative_trait[which(df.all2$cooperative_trait == 'qs')]<- 'QS'
df.all2$cooperative_trait[which(df.all2$cooperative_trait == 'sid')]<- 'Sideroph.'
df.all2$cooperative_trait[which(df.all2$cooperative_trait == 'ssyst')]<- 'Secr. Syst.'
df.all2$cooperative_trait[which(df.all2$cooperative_trait == 'secretome')]<- 'Secretome'

unique(df.all2$cooperative_trait)
df.all2$cooperative_trait<- factor(df.all2$cooperative_trait,
                                   levels = c('QS', 'Biofilm', 'Sideroph.', 'Antib. degr.', 'Secr. Syst.', 'Secretome', 'is_coop'))


p.ORobs<- df.all2 %>%
  filter(cooperative_trait != 'is_coop') %>%
  ggplot(., aes(x = exp(logOR), y = cooperative_trait))+
  geom_boxplot(width = 0.5, outlier.size = 0.05, lwd = 0.2)+
  xlab('Observed odds ratios')+ylab('')+
  xlim(0,1000)+
  theme_bw()+
  theme(axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 5, colour = 'black'),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.grid = element_blank(),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(colour = 'black', size = 0.3),
        plot.margin=unit(c(0,0,0,0), "null"))

# ... (B) estimated OR ----


#load("./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN1.RData")

meta.results<- rbind(
  data.frame(summary(meta.a.secretome)$solutions, trait = 'Secretome'),
  data.frame(summary(meta.a.biofilm)$solutions, trait = 'Biofilm'),
  data.frame(summary(meta.a.sid)$solutions, trait = 'Siderophores'),
  data.frame(summary(meta.a.ab)$solutions, trait = 'Antib. degr.'),
  data.frame(summary(meta.a.qs)$solutions, trait = 'QS'),
  data.frame(summary(meta.a.ssyst)$solutions, trait = 'Secr. syst.'),
  data.frame(summary(meta.a.coop)$solutions, trait = 'Cooperation (any)'))


meta.results$trait<- factor(meta.results$trait,
                            levels = c('QS',
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

meta.results<- meta.results %>% filter(trait != 'Cooperation (any)')

p.vf.meta<- ggplot(meta.results, aes(x = `post.mean`, y = trait, col = 'black'))+
  xlab('Estimated odds ratio\n  ')+ylab('')+
  #ggtitle('Are cooperative products\nvirlence factors?')+
  geom_point()+
  #xlim(-1,5)+
  geom_errorbarh(aes(xmin = `l.95..CI`, xmax = `u.95..CI`), height = 0.1, size = 0.25)+
  geom_vline(xintercept = 1, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l.95..CI`, xmax = `u.95..CI`), height = 0.1, col = 'black', size = 0.25)+
  scale_color_manual(values = 'white')+
  theme_bw()+
  theme(axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 5, colour = 'black'),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(colour = 'black', size = 0.3),
        plot.margin=unit(c(0,0,0,0), "null"))


pdf('./output/2_figures/Figure3_analysis_VF.pdf', width = 12/2.55, height = 5/2.55)
grid.arrange(ncol = 2, p.ORobs, p.vf.meta)
dev.off()



# Figure 4 ----
# ... (A) P vs NP boxplots  ----

# Load model run output which has the dataframe with the z-transformed number of genes
load('./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN1.RData')

# boxplots
d.long<- d %>% gather('trait', 'z_value', grep('z_', colnames(d)))

d.long$trait<- gsub('z_biofilm', 'Biofilm', d.long$trait)
d.long$trait<- gsub('z_quorum_sensing', 'QS', d.long$trait)
d.long$trait<- gsub('z_ab_degradation', 'Antib. degr.', d.long$trait)
d.long$trait<- gsub('z_secretion_system', 'Secr. Syst.', d.long$trait)
d.long$trait<- gsub('z_siderophores', 'Sideroph.', d.long$trait)
d.long$trait<- gsub('z_nb_extracellular', 'Secretome', d.long$trait)
d.long$trait<- gsub('z_vf', 'VF', d.long$trait)
#d.long$trait<- gsub('z_mp3', 'MP3', d.long$trait)

d.long<- d.long[d.long$trait != 'z_nb_cds',]
d.long$trait<- factor(d.long$trait,
                      levels = c('VF', 'QS', 'Biofilm', 'Sideroph.', 'Antib. degr.', 'Secr. Syst.', 'Secretome'))

# REPLACING VALUE 7.78 by 5.78, modify x-axis labels to display 8 and will manually display axis break
d.long[which(d.long$z_value > 6), 'z_value']<- 5.78

#pdf('~/Desktop/HM1Figs/PvsNP_boxplot_2.pdf', width = 6/2.55, height = 5/2.55)
PNP_boxplot<- ggplot(d.long[d.long$trait != 'z_nb_cds',], aes(x = z_value, y = trait, fill = pathogen))+
  xlab('Z-number of genes')+ylab('')+
  xlim(-1,6)+
  geom_boxplot(width = 0.5, outlier.size = 0.05, lwd = 0.2)+
  scale_fill_manual(values = c('dodgerblue', 'darkorange'), labels = c('Non-pathogen', 'Pathogen'))+
  theme_bw()+
  theme(axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 5, colour = 'black'),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.grid = element_blank(),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(colour = 'black', size = 0.3),
        plot.margin=unit(c(0,0,0,0), "null"))+
  scale_x_continuous(breaks = c(0, 2, 4, 6), labels = c(0, 2, 4, 8))
#dev.off()


# ... (B) P vs NP estimates ----

m1.outs.multi<- as.data.frame(summary(m1.multi.with_vf_NOgram)$solutions)[-c(1,2),] %>%
  mutate(trait = c('Biofilm', 'Antib. degr.', 'QS', 'Secr. Syst.', 'Sideroph.', 'Secretome', 'VF')) %>%
  #select(model, trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)
  select(trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)

m1.outs.multi$trait<- factor(m1.outs.multi$trait,
                             levels = c('VF', 'QS', 'Biofilm', 'Sideroph.', 'Antib. degr.', 'Secr. Syst.', 'Secretome'))

p.patho.meta<- ggplot(m1.outs.multi, aes(x = post.mean, y = trait, col = 'black'))+
  xlab('Estimated regression coefficients\n(multivariate model)')+ylab('')+
  #ggtitle('Is cooperation predictive\nof pathogenicity?')+
  geom_point()+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1, size = 0.25)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1, col = 'black', size = 0.25)+
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  scale_color_manual(values = 'white')+
  theme_bw()+
  theme(axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 5, colour = 'black'),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(colour = 'black', size = 0.3),
        plot.margin=unit(c(0,0,0,0), "null"))



pdf('./output/2_figures/Figure4_analysis_Pathogenicity.pdf', width = 12/2.55, height = 5/2.55)
grid.arrange(ncol = 2, PNP_boxplot, p.patho.meta)
dev.off()


# Getting out plot with legend as well
pdf('./output/2_figures/Figure4_analysis_Pathogenicity_legend.pdf', width = 6/2.55, height = 6/2.55)
PNP_boxplot+
  theme(legend.position = 'top', #c(0.9, 0.1),
        legend.justification = 'left', 
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 4),
        legend.title = element_blank())
dev.off()



# Figure 5 ----
# ... (A) CFR phylogeny and heatmap ----

load("./output/3_model_output/CompAnalysis_cfr_CHAIN1.RData")

df<- d.cfr %>%
  mutate(species_id2 = species_id,
         Legget_species2 = Legget_species
  ) %>%
  rename(id = Legget_species,
         ab_degradation = antibiotic_degradation,
         vf = is_victor_vf,
         nb_extracellular = secretome) %>%
  as.data.frame()

head(df)


# Make shorter names for plot
df$Legget_species2<- gsub(' subsp. enterica serovar', '',df$Legget_species2)
df$id<- gsub(' subsp. enterica serovar', '',df$id)


coop.data.cfr<- 
  df %>%
  select(id, biofilm, ab_degradation, quorum_sensing, siderophores, secretion_system, nb_extracellular, vf) %>%
  rename(label = id) %>%
  gather('category', 'value', 2:8)
head(coop.data.cfr)

coop.data.cfr.props<- 
  df %>%
  select(id, nb_cds) %>%
  rename(label = id) %>%
  left_join(coop.data.cfr, by = 'label') %>%
  select(label, category, value, nb_cds) %>%
  mutate(prop = value/nb_cds) %>%
  mutate(perc = prop*100)

head(coop.data.cfr.props)

coop.data.cfr.props<-
  rbind(coop.data.cfr.props,
        df[,c('id', 'case_fatality_rate')] %>%
          gather('category', 'perc', 2) %>%
          mutate(value = NA, nb_cds = NA, prop = NA) %>%
          rename(label = id) %>%
          select(label, category, value, nb_cds, prop, perc))



library(viridis)
plotHeatmap<- function(focal.trait, pal.type = 'viridis', leg.pos = 'none', margins = c(0,-0.09,0,-0.09)){
  
  plot<- coop.data.cfr.props %>%
    filter(category %in% focal.trait) %>%
    ggplot(., aes(x = category, y = label, fill = perc))+
    geom_tile()+
    scale_fill_viridis(option=pal.type)+
    scale_x_discrete(position = "top")+
    theme(axis.text.x = element_text(angle = 45, size = 5, hjust = 0))+
    xlab('')+ylab('')+
    #theme_bw()+
    # theme(axis.ticks.y = element_blank(),
    #       axis.text.y = element_blank(),
    #       panel.grid = element_blank(),
    #       panel.border = element_blank(),
    #       panel.background = element_blank())+
    theme_void()+
    theme(legend.position = leg.pos)+
    theme(plot.margin=unit(margins, "null"))
  
  return(plot)
  
}


# Get the phylogeny
tree<- read.tree('./data/5_phylogeny_files/supfam_cfr_tree.newick')
tree$tip.label
tree<- drop.tip(tree, tree$tip.label[!tree$tip.label %in% df$supfam_id])
tree$tip.label<- df$Legget_species2[match(tree$tip.label, df$supfam_id)]



# MAKE THE TREE
p<- ggtree(tree, size = 0.3) %<+% df + 
  geom_tiplab(aes(label=Legget_species2, col = infection_route),
              align=T,
              offset=0, hjust=0, linesize = 0.2, size = 1.5) +
  scale_color_manual(values = c('darkgoldenrod1','forestgreen', 'plum'))+
  theme(legend.position = 'none')+
  xlim(0,1)
p


pdf('./output/2_figures/Figure5_CFR_A1_phylogeny.pdf', width = 8/2.55, height = 8/2.55)
p
dev.off()


pdf('./output/2_figures/Figure5_CFR_A2_heatmap.pdf', width = 8/2.55, height = 8/2.55)
plotHeatmap('nb_extracellular') %>%
  insert_left(p, width = 9) %>%
  insert_right(plotHeatmap('secretion_system')) %>%
  insert_right(plotHeatmap('ab_degradation')) %>%
  insert_right(plotHeatmap('siderophores')) %>%
  insert_right(plotHeatmap('biofilm')) %>%
  insert_right(plotHeatmap('quorum_sensing')) %>%
  insert_right(plotHeatmap('vf')) %>%
  insert_right(plotHeatmap('case_fatality_rate', pal.type = 'magma'))
dev.off()


# ... (B) CFR scatterplots ----

#load("./output/3_model_output/CompAnalysis_cfr_CHAIN1.RData")

cfr.theme<-
  theme_bw()+
  theme(axis.title = element_text(face = 'bold', size = 5),
        axis.text = element_text(size = 4, colour = 'black'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.25),
        axis.ticks = element_line(colour = "black", size = 0.25),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 5, face = 'bold'))

pt.size = 0.9
pt.lwd = 0.16

cfr.cols<- c('darkgoldenrod1','forestgreen', 'plum')


p.cfr.bio<- ggplot(d.cfr, aes(x = biofilm/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = cfr.cols)+
  xlab(' ')+ylab(' ')+ggtitle('Biofilm')+ cfr.theme


p.cfr.sec<- ggplot(d.cfr, aes(x = secretome/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = cfr.cols)+xlab(' ')+ylab(' ')+ggtitle('Secretome')+ cfr.theme

p.cfr.qs<- ggplot(d.cfr, aes(x = quorum_sensing/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = cfr.cols)+xlab(' ')+ylab(' ')+ggtitle('Quorum-sensing')+ cfr.theme

p.cfr.ab<- ggplot(d.cfr, aes(x = antibiotic_degradation/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = cfr.cols)+xlab(' ')+ylab(' ')+ggtitle('Antib. degr.')+ cfr.theme

p.cfr.sid<- ggplot(d.cfr, aes(x = siderophores/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = cfr.cols)+xlab(' ')+ylab(' ')+ggtitle('Siderophores')+ cfr.theme

p.cfr.ssy<- ggplot(d.cfr, aes(x = secretion_system/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = cfr.cols)+xlab(' ')+ylab(' ')+ggtitle('Secr. Syst.')+ cfr.theme

p.cfr.vf<- ggplot(d.cfr, aes(x = is_victor_vf/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(values = cfr.cols)+xlab(' ')+ylab(' ')+ggtitle('Virulence factors')+ cfr.theme



p.cfr.legend<- ggplot(d.cfr, aes(x = biofilm/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(name = 'Infection route', values = cfr.cols)+
  theme_bw()+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))


cfr.legend <- cowplot::get_legend(p.cfr.legend)



pdf('./output/2_figures/Figure5_CFR_B_scatterplot.pdf', width = 4.4, height = 4.4+0.2)

grid.arrange(grobs = c(#list(p.cfr.bio, p.cfr.qs,p.cfr.ab),
  #list(p.cfr.sid, p.cfr.ssy, p.cfr.sec),
  #list(p.cfr.vf, cfr.legend)
  list(p.cfr.sec, p.cfr.sid, p.cfr.vf),
  list(p.cfr.ssy, p.cfr.bio,cfr.legend),
  list(p.cfr.ab, p.cfr.qs)
),
ncol = 3,
as.table = FALSE,
bottom = textGrob("Proportion of genes coding for trait",gp=gpar(fontsize=6,font=1)),
left = textGrob("Case fatality rate (%)",gp=gpar(fontsize=6,font=1), rot = 90))

dev.off()




# ... (C)Plot of estimates ----

#load("./output/3_model_output/CompAnalysis_cfr_CHAIN1.RData")

# params of multivariate models
m3.outs.multi<- as.data.frame(summary(m3.multi.with_vf)$solutions)[-c(1:5),] %>%
  mutate(trait = c('Biofilm', 'QS', 'Antib. degr.', 'Secr. Syst.', 'Siderophores', 'Secretome', 'VF')) %>%
  select(trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)


m3.outs.multi$trait<- factor(m3.outs.multi$trait,
                             levels = c('VF', 'QS', 'Biofilm', 'Siderophores', 'Antib. degr.', 'Secr. Syst.', 'Secretome'))


p.cfr.meta<- ggplot(m3.outs.multi, aes(x = post.mean, y = trait, col = 'black'))+
  xlab('Estimated regression coefficients\n(multivariate model)')+ylab('')+
  #ggtitle('Are more cooperative pathogens\nmore virulent?')+
  geom_point()+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1, size = 0.25)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1, col = 'black', size = 0.25)+
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  scale_color_manual(values = 'white')+
  theme_bw()+
  theme(axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 5, colour = 'black'),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(colour = 'black', size = 0.3),
        plot.margin=unit(c(0,0,0,0), "null"))


pdf('./output/2_figures/Figure5_CFR_C_effectEstimates.pdf', width = 6/2.55, height = 5/2.55)
p.cfr.meta
dev.off()



# Figure 6 (literature) ----


# cheking out overall what literature on microbial cooperation is talking about

library(ggplot2)
library(readxl)
library(dplyr)
library(plyr)


d<- read_excel('./data/6_webofscience_search/WebOfScience_microbial_cooperation.xlsx', sheet = 3, skip = 2) %>%
  filter(!is.na(Title))

length(d$Title) # 405 rows
length(unique(d$Title)) # but 398 titles. This is because we duplicated a few records where different studies were looking independently at different social traits


d.filtered<- d %>%
  filter(!trait %in% c('zNA', 'zNA(duplicate record)'))

length(d.filtered$Title) # 283 rows
length(unique(d.filtered$Title)) # but 277 titles. This is because we duplicated a few records where different studies were looking independently at different social traits

# Among experimental papers, 277/398 = 70% that were about a specific microbial social trait
sort(unique(d.filtered$trait))



dat.freq<- d.filtered %>%
  select(study_type, trait) %>%
  #filter(!study_type %in% c('review', 'computational/theoretical')) %>%
  #filter(study_type %in% c('experiment/data analysis')) %>%
  group_by(trait) %>%
  table() %>%
  as.data.frame() %>%
  arrange(Freq)

head(dat.freq)
unique(dat.freq$trait)

dat.freq$trait<- factor(dat.freq$trait,
                        levels = unique(as.character(dat.freq$trait)))

p.barplot.literature<- ggplot(dat.freq, aes(y = trait, x = Freq))+
  geom_bar(stat="identity")+
  ylab('Social trait studied\n ')+xlab('Number of studies')+
  theme_bw()+
  theme(axis.title = element_text(size = 5),
        axis.text.y = element_text(size = 5, colour = 'black'),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.grid = element_blank(),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(colour = 'black', size = 0.3),
        plot.margin=unit(c(0,0,0,0), "null"))


pdf('./output/2_figures/Figure6_LiteratureReview_experiments.pdf', height = 9/2.55, width = 9/2.55)
p.barplot.literature
dev.off()



# SUPPLEMENTARY MATERIAL ----


# Figure S1 (funnel funnel plot) ----


load("./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN1.RData")


plot.funnel<- function(df, trait, title, model){
  
  df %>%
    filter(cooperative_trait == trait) %>%
    ggplot(., aes(x = logOR, y = 1/se_logOR))+
    xlab(' ')+ylab(' ')+
    ggtitle(title)+
    #scale_color_manual(values = 'white')+
    geom_rect(ymin = 0, ymax = 5,
              xmin = summary(model)$solutions[,2],
              xmax = summary(model)$solutions[,3],
              fill = 'gray89', color = NA)+
    geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.3)+
    geom_point(size = 0.5, alpha = .5)+
    #geom_point(col = 'black', size = 0.5)+
    #scale_color_manual(values = c('dodgerblue', 'darkorange'))+
    geom_vline(xintercept = summary(model)$solutions[,1], linetype = 'dashed', col = 'firebrick', size = 0.5)+
    
    theme_bw()+
    theme(axis.title = element_text(face = 'bold', size = 5),
          axis.text = element_text(size = 5.5, colour = "black"),
          #panel.grid = element_line(size = 0.1),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", size = 0.3),
          axis.ticks = element_line(colour = "black", size = 0.3),
          legend.position = 'none', #c(0.9, 0.1),
          #legend.key.size = unit(0.3, "cm"),
          #legend.text = element_text(size = 0),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = 'bold', size = 6))
} 


df.all2<- left_join(df.all, unique(d.annot2[,c('species', 'pathogen')]), by = 'species')

library(gridExtra)
library(grid)
library(gridtext)

pdf('./output/2_figures/FigureS1_funnelPlots.pdf', width = 3.43, height = 5.2)

grid.arrange(ncol = 2,
             #plot.funnel(df.all, 'is_coop', 'Cooperative (any)', meta.a.coop),
             plot.funnel(df.all2, 'secretome', 'Secretome', meta.a.secretome),
             plot.funnel(df.all2, 'ssyst', 'Secretion-systems', meta.a.ssyst),
             plot.funnel(df.all2, 'ab', 'Antibiotic degradation', meta.a.ab),
             plot.funnel(df.all2, 'sid', 'Siderophores', meta.a.sid),
             plot.funnel(df.all2, 'biofilm', 'Biofilm', meta.a.biofilm),
             plot.funnel(df.all2, 'qs', 'Quorum-sensing', meta.a.qs),
             
             bottom = textGrob("Log odds ratio estimated for each species,\n and estimate of the mean odds ratio accross species (meta-analysis)",gp=gpar(fontsize=6,font=2)),
             left = textGrob("1/SE(logOR)",gp=gpar(fontsize=6,font=2), rot = 90))

dev.off()



# Figure S2 (boxplots raw) ----

load('./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN1.RData')

# boxplots
d.long.raw<- d %>% gather('trait', 'value', gsub('z_', '', colnames(d)[grep('z_', colnames(d))]))

d.long.raw$trait<- gsub('biofilm', 'Biofilm', d.long.raw$trait)
d.long.raw$trait<- gsub('quorum_sensing', 'QS', d.long.raw$trait)
d.long.raw$trait<- gsub('ab_degradation', 'Antib. degr.', d.long.raw$trait)
d.long.raw$trait<- gsub('secretion_system', 'Secr. Syst.', d.long.raw$trait)
d.long.raw$trait<- gsub('siderophores', 'Siderophores', d.long.raw$trait)
d.long.raw$trait<- gsub('nb_extracellular', 'Secretome', d.long.raw$trait)
d.long.raw$trait<- gsub('vf', 'VF', d.long.raw$trait)

d.long.raw<- d.long.raw[d.long.raw$trait != 'nb_cds',]
d.long.raw$trait<- factor(d.long.raw$trait,
                          levels = c('VF', 'QS', 'Biofilm', 'Siderophores', 'Antib. degr.', 'Secr. Syst.', 'Secretome')) # order as models *UNI* table


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
          axis.text = element_text(size = 5, colour = 'black'),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", size = 0.25),
          axis.ticks = element_line(colour = "black", size = 0.25),
          panel.grid = element_blank(),
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
p6<- make.boxplot.pathogen('QS')
p7<- make.boxplot.pathogen('VF')


pdf('./output/2_figures/FigureS2_boxplotRaw.pdf', width = 4.4, height = 4.4+0.2)
#grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 3)
grid.arrange(p1,p2,p3,p4,p5,p6,p7, ncol = 3)
dev.off()



# Figure S3 (univariate model estimates) ----

theme_results<-   theme_bw()+
  theme(axis.title = element_text(size = 5),
        axis.text = element_text(size = 5, colour = 'black'),
        #panel.grid = element_line(size = 0.1),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks = element_line(colour = "black", size = 0.2),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 0),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 6, face = 'bold'))


# PATHOGENICITY MODEL
load("./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN1.RData")

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

p.patho.uni<- ggplot(m1.outs.uni, aes(x = post.mean, y = trait, col = 'black'))+
  xlab('Estimated regression coefficients\n(univariate models)')+ylab('')+
  geom_point()+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1, col = 'black', size = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  scale_color_manual(values = 'white')+
  theme_bw()+
  theme_results


# CFR model 

load("./output/3_model_output/CompAnalysis_cfr_CHAIN1.RData")


# Model version with growth rate included

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

p.cfr.uni<- ggplot(m3.outs.uni, aes(x = post.mean, y = trait, col = 'black'))+
  xlab('Estimated regression coefficients\n(univariate models)')+ylab('')+
  geom_point()+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1, col = 'black', size = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  scale_color_manual(values = 'white')+
  theme_results



# Model version WITHOUT growth rate included

# params of univariate models
m3.outs.uni.noGrowth<- rbind(
  as.data.frame(summary(m3.sec.noGrowth)$solutions)['focal_trait',] %>% mutate(trait = 'Secretome'),
  as.data.frame(summary(m3.ssy.noGrowth)$solutions)['focal_trait',] %>% mutate(trait = 'Secr. Syst.'),
  as.data.frame(summary(m3.ab.noGrowth)$solutions)['focal_trait',] %>% mutate(trait = 'Antib. degr.'),
  as.data.frame(summary(m3.sid.noGrowth)$solutions)['focal_trait',] %>% mutate(trait = 'Siderophores'),
  as.data.frame(summary(m3.bio.noGrowth)$solutions)['focal_trait',] %>% mutate(trait = 'Biofilm'),
  as.data.frame(summary(m3.qs.noGrowth)$solutions)['focal_trait',] %>% mutate(trait = 'Quorum-sensing'),
  as.data.frame(summary(m3.vf.noGrowth)$solutions)['focal_trait',] %>% mutate(trait = 'Virulence factors')
) %>%
  select(trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)

m3.outs.uni.noGrowth$trait<- factor(m3.outs.uni.noGrowth$trait, levels = rev(m3.outs.uni.noGrowth$trait))

p.cfr.uni.noGrowth<- ggplot(m3.outs.uni.noGrowth, aes(x = post.mean, y = trait, col = 'black'))+
  xlab('Estimated regression coefficients\n(univariate models)')+ylab('')+
  geom_point()+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1, col = 'black', size = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  scale_color_manual(values = 'white')+
  theme_results



library(patchwork)


pdf('./output/2_figures/FigureS3_modelResults_univariate.pdf', width = (5.8)/2.55, height = (4.5*3)/2.55)

(p.patho.uni/p.cfr.uni/p.cfr.uni.noGrowth)


dev.off()



# CODE WRAPPERS FOR SUPP. TABLES ----

format_effects.HM1<- function(string){
  
  string_edited<- string %>%
    gsub(pattern = "(Intercept)", replacement = "Intercept", ., fixed = TRUE) %>%
    gsub(pattern = "infection_routeIngestion", replacement = "Infection route: Ingestion", ., fixed = TRUE) %>%
    gsub(pattern = "infection_routeInhalation", replacement = "Infection route: Inhalation", ., fixed = TRUE) %>%
    gsub(pattern = "infection_routeSkin", replacement = "Infection route: Skin", ., fixed = TRUE) %>%
    gsub(pattern = "generation_time_h", replacement = "Generation time (h)", ., fixed = TRUE) %>%
    gsub(pattern = "z_nb_cds", replacement = "z(Proteome size)", ., fixed = TRUE) %>%
    gsub(pattern = "z_biofilm", replacement = "z(Biofilm)", ., fixed = TRUE) %>%
    gsub(pattern = "z_ab_degradation", replacement = "z(Antib. degr.)", ., fixed = TRUE) %>%
    gsub(pattern = "z_antibiotic_degradation", replacement = "z(Antib. degr.)", ., fixed = TRUE) %>%
    gsub(pattern = "z_quorum_sensing", replacement = "z(Quorum-sensing)", ., fixed = TRUE) %>%
    gsub(pattern = "z_siderophores", replacement = "z(Siderophores)", ., fixed = TRUE) %>%
    gsub(pattern = "z_nb_extracellular", replacement = "z(Secretome)", ., fixed = TRUE) %>%
    gsub(pattern = "z_secretome", replacement = "z(Secretome)", ., fixed = TRUE) %>%
    gsub(pattern = "z_vf", replacement = "z(Virulence Factors)", ., fixed = TRUE) %>%
    gsub(pattern = "z_is_victor_vf", replacement = "z(Virulence Factors)", ., fixed = TRUE) %>%
    gsub(pattern = "z_secretion_system", replacement = "z(Secretion systems)", ., fixed = TRUE) %>%
    gsub(pattern = "nb_cds", replacement = "Proteome size", ., fixed = TRUE) %>%
    gsub(pattern = "species_id", replacement = "Species, phylogenetic", ., fixed = TRUE) %>%
    gsub(pattern = "species", replacement = "Species, phylogenetic", ., fixed = TRUE) %>%
    gsub(pattern = "se_logOR.units", replacement = "Measurement error", ., fixed = TRUE) %>%
    gsub(pattern = "focal_trait", replacement = "Focal trait", ., fixed = TRUE) %>%
    gsub(pattern = "units", replacement = "Residual", ., fixed = TRUE)
  
  return(string_edited)  
}
format_full_summary.HM1<- function(model, trait){
  
  fixed<- summary(model)$solutions %>%
    as.data.frame()%>%
    rownames_to_column('Effect') %>%
    mutate(Structure = 'Fixed effect',
           pMCMC = ifelse(pMCMC<0.01,
                          formatC(pMCMC, digit = 2, format = 'e'),
                          formatC(pMCMC, digit = 3, format = 'f')))
  
  random<- summary(model)$Gcovariances %>%
    as.data.frame()%>%
    rownames_to_column('Effect') %>%
    mutate(pMCMC = '') %>%
    mutate(Structure = '(Co)variance')
  
  
  unit<- summary(model)$Rcovariances %>%
    as.data.frame()%>%
    rownames_to_column('Effect') %>%
    mutate(pMCMC = '') %>%
    mutate(Structure = '(Co)variance')
  
  
  tab<- rbind(fixed, random, unit) %>%
    mutate(Effect = format_effects.HM1(Effect),
           Model = trait) %>%
    mutate(Structure = ifelse(duplicated(Structure) == TRUE, '', Structure),
           Model = ifelse(duplicated(Model) == TRUE, '', Model)) %>%
    select(Model, Structure, Effect, post.mean,`l-95% CI`, `u-95% CI`, eff.samp, pMCMC) %>%
    rename(`Posterior\n mean` = post.mean,
           `CI95% lower` = `l-95% CI`,
           `CI95% upper` = `u-95% CI`,
           `Effective\n sampling` = eff.samp) %>%
    as.data.frame() %>%
    mutate(`Effective\n sampling` = formatC(`Effective\n sampling`, digit = 0, format = 'f')) %>%
    mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'f')))
  
  
  return(tab)
  
}




# SUPP. TABLE 1 (global chi squares) ----


make.ctab.all.Table<- function(focal.ctab, focal.trait){
  
  supptab.ctab<- data.frame(Trait = focal.trait, 
                            `VF\n Coop` = focal.ctab$tab[2,2],
                            `VF_non\n Coop` = focal.ctab$tab[2,1],
                            `non-VF\n Coop` = focal.ctab$tab[1,2],
                            `non-VF\n non-Coop` = focal.ctab$tab[1,1],
                            OR = as.data.frame(focal.ctab$or)[,'es'],
                            `l-95% CI` = as.data.frame(focal.ctab$or)[,'ci.lo'],
                            `u-95% CI` = as.data.frame(focal.ctab$or)[,'ci.hi'],
                            Chisq = focal.ctab$chi$statistic,
                            df = focal.ctab$chi$parameter,
                            pval = focal.ctab$chi$p.value)
  
  rownames(supptab.ctab)<- NULL
  
  return(supptab.ctab)
}

# Need to run section 1 of script 4.1 before this ("QUICK OVERAL CHI-SQUARE")

ctabSuppsTable<- rbind(
  make.ctab.all.Table(ctab.all.secretome, 'Secretome'),
  make.ctab.all.Table(ctab.all.biofilm, 'Biofilm'),
  make.ctab.all.Table(ctab.all.sid, 'Siderophore'),
  make.ctab.all.Table(ctab.all.qs, 'QS'),
  make.ctab.all.Table(ctab.all.ab, 'Antib. degr.'),
  make.ctab.all.Table(ctab.all.ssyst, 'Secr. Syst.'),
  make.ctab.all.Table(ctab.all.any, 'All'))


roundit<- function(vec){
  
  ifelse(round(vec, 2) > 0,
         round(vec, 2),
         format(vec, scientific = TRUE, big.mark = ",", digit = 2))
}


ctabSuppsTable$OR<- roundit(ctabSuppsTable$OR)
ctabSuppsTable$`l.95..CI`<- roundit(ctabSuppsTable$`l.95..CI`)
ctabSuppsTable$`u.95..CI`<- roundit(ctabSuppsTable$`u.95..CI`)
ctabSuppsTable$Chisq<- roundit(ctabSuppsTable$Chisq)
ctabSuppsTable$pval<- roundit(ctabSuppsTable$pval)

colnames(ctabSuppsTable)<- c('Trait', '(VF,Coop)', '(VF,non-Coop)', '(non-VF,Coop)', '(non-VF,non-Coop)',
                             'OR', 'l-95% CI', 'u-95% CI', 'Chisq', 'df', 'pvalue')

write.table(ctabSuppsTable, './output/4_summary_tables/TABLE_S1.hm1.csv', sep = '\t', col.names = TRUE,quote = FALSE,  row.names = FALSE)


library(kableExtra)

tab.s1<- kable(ctabSuppsTable, "latex", booktabs = T,
               caption = 'Cooperative and non-cooperative virulence factors in the global pool of coding sequences.
               Corresponding odds ratio (OR) and chi-square test testing whether the OR is significantly greater than 1.
               The odds ratio is as defined in the main text methods: ratio of the odds of a virulence factor being cooperative to the odds of a non-virulence factor being cooperative.
               Abbreviations: VF = Virulence Factor, non-VF = not virulence factor, Coop = cooperative, non-Coop = not cooperative, QS = quorum-sensing, All = when a coding sequence is considered cooperative if it is annotated as any of the 6 forms of cooperation.')

fileConn<-file("./output/4_summary_tables/TABLE_S1.hm1.tex")
writeLines(tab.s1, fileConn)
close(fileConn)



# SUPP. TABLE 2 (summary VF analysis) ----

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


load('./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN1.RData')

summary.tab.vf<- rbind(
  get.effects.M2(meta.a.secretome, 'secretome'),
  get.effects.M2(meta.a.biofilm, 'biofilm'),
  get.effects.M2(meta.a.sid, 'siderophores'),
  get.effects.M2(meta.a.ab, 'antibiotic degradation'),
  get.effects.M2(meta.a.qs, 'quorum-sensing'),
  get.effects.M2(meta.a.ssyst, 'secretion systems'),
  get.effects.M2(meta.a.coop, 'cooperation (any)')
)


write.table(summary.tab.vf, './output/4_summary_tables/TABLE_S2.hm1.csv', sep = '\t', col.names = TRUE,quote = FALSE,  row.names = FALSE)


library(kableExtra)


tab<-rbind(
  format_full_summary.HM1(meta.a.secretome, ''),
  format_full_summary.HM1(meta.a.biofilm, ''),
  format_full_summary.HM1(meta.a.sid, ''),
  format_full_summary.HM1(meta.a.ab, ''),
  format_full_summary.HM1(meta.a.qs, ''),
  format_full_summary.HM1(meta.a.ssyst, '')
)


tab.s2<- kable(tab, "latex", booktabs = T, caption = 'Model summaries of virulence factors phylogenetic meta-analyses (estimates in model summary are on log scale)') %>%
  footnote(c("CI95%: 95% credible interval of the posterior distribution",
             "pMCMC: taken as twice the posterior probability that the estimate is negative"),
           fixed_small_size = TRUE, general_title = "") %>%
  kable_styling() %>%
  pack_rows("Secretome", 1, 4) %>%
  pack_rows("Biofilm", 5, 8) %>%
  pack_rows("Siderophores", 9, 12) %>%
  pack_rows("Antibiotic degradation", 13, 16) %>%
  pack_rows("Quorum sensing", 17, 20) %>%
  pack_rows("Secretion systems", 21, 24)


fileConn<-file("./output/4_summary_tables/TABLE_S2.hm1.tex")
writeLines(tab.s2, fileConn)
close(fileConn)



# SUPP. TABLE 3 (summary PATHOGENICITY analysis) ----


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


load('./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN1.RData')


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


write.table(summary.tab.pathogenicity, './output/4_summary_tables/TABLE_s3.hm1.csv', sep = '\t', col.names = TRUE,quote = FALSE,  row.names = FALSE)


tab3<-rbind(
  format_full_summary.HM1(m1.ss, ''),
  format_full_summary.HM1(m1.biofilm, ''),
  format_full_summary.HM1(m1.siderophores, ''),
  format_full_summary.HM1(m1.ab_degradation, ''),
  format_full_summary.HM1(m1.quorum_sensing, ''),
  format_full_summary.HM1(m1.secretion_system, ''),
  format_full_summary.HM1(m1.vf, ''),
  format_full_summary.HM1(m1.multi.with_vf_NOgram, '')
)



tab.s3<- kable(tab3, "latex", booktabs = T, caption = 'Model summaries of pathogenicity comparative analyses. This includes: 6 univariate models with each cooperative traits as predictor, one univariate model with number of virulence factors as predictor, and one multivariate model with all predictors included)') %>%
  footnote(c("CI95%: 95% credible interval of the posterior distribution",
             "pMCMC: taken as twice the posterior probability that the estimate is negative",
             "z(): For multivariate model, predictors were z-transformed"),
           fixed_small_size = TRUE, general_title = "") %>%
  kable_styling() %>%
  pack_rows("Secretome", 1, 5) %>%
  pack_rows("Biofilm", 6, 10) %>%
  pack_rows("Siderophores", 11, 15) %>%
  pack_rows("Antibiotic degradation", 16, 20) %>%
  pack_rows("Quorum sensing", 21, 25) %>%
  pack_rows("Secretion systems", 26, 30) %>%
  pack_rows("Virulence factors", 31, 35) %>%
  pack_rows("Multivariate model", 36, 46)


fileConn<-file("./output/4_summary_tables/TABLE_S3.hm1.tex")
writeLines(tab.s3, fileConn)
close(fileConn)



# SUPP. TABLE 4 (summary CFR analysis) ----

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


load('./output/3_model_output/CompAnalysis_cfr_CHAIN1.RData')

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


write.table(summary.tab.cfr, './output/4_summary_tables/TABLE_s4.hm1.csv', sep = '\t', col.names = TRUE,quote = FALSE,  row.names = FALSE)


tab4<-rbind(
  format_full_summary.HM1(m3.legget, ''),
  format_full_summary.HM1(m3.bio, ''),
  format_full_summary.HM1(m3.qs, ''),
  format_full_summary.HM1(m3.sec, ''),
  format_full_summary.HM1(m3.ab, ''),
  format_full_summary.HM1(m3.sid, ''),
  format_full_summary.HM1(m3.ssy, ''),
  format_full_summary.HM1(m3.vf, ''),
  format_full_summary.HM1(m3.multi.with_vf, "")
)



tab.s4<- kable(tab4, "latex", booktabs = T, caption = 'Model summaries of virulence comparative analyses (case fatality rate). This includes: one model replicating Leggett\'s et al (2017) main results (non-cooperative traits), 6 univariate models with each cooperative traits as predictor, one univariate model with number of virulence factors as predictor, and one multivariate model with all predictors included)') %>%
  footnote(c("CI95%: 95% credible interval of the posterior distribution",
             "pMCMC: taken as twice the posterior probability that the estimate is negative",
             "z(): For multivariate model, predictors were z-transformed"),
           fixed_small_size = TRUE, general_title = "") %>%
  kable_styling() %>%
  pack_rows("Non-cooperative traits", 1, 6) %>%
  pack_rows("Secretome", 7, 14) %>%
  pack_rows("Biofilm", 15, 22) %>%
  pack_rows("Siderophores", 23, 30) %>%
  pack_rows("Antibiotic degradation", 31, 38) %>%
  pack_rows("Quorum sensing", 39, 46) %>%
  pack_rows("Secretion systems", 47, 54) %>%
  pack_rows("Virulence factors", 55, 62) %>%
  pack_rows("Multivariate model", 63, 76)



fileConn<-file("./output/4_summary_tables/TABLE_S4.hm1.tex")
writeLines(tab.s4, fileConn)
close(fileConn)



# SUPP. TABLE 5 (summary CFR analysis, without growth rate) ----

summary.tab.cfr.noGrowth<- 
  rbind(
    get.effects.M3(m3.bio.noGrowth, 'biofilm'),
    get.effects.M3(m3.qs.noGrowth, 'quorum-sensing'),
    get.effects.M3(m3.sec.noGrowth, 'secretome'),
    get.effects.M3(m3.ab.noGrowth, 'antibiotic-resistance'),
    get.effects.M3(m3.sid.noGrowth, 'siderophores'),
    get.effects.M3(m3.ssy.noGrowth, 'secretion-systems'),
    get.effects.M3(m3.vf.noGrowth, 'virulence-factors'),
    get.effects.M3(m3.multi.with_vf.NOgrowthRate, "multivariate-model"))


write.table(summary.tab.cfr.noGrowth, './output/4_summary_tables/TABLE_s5.hm1.csv', sep = '\t', col.names = TRUE,quote = FALSE,  row.names = FALSE)


tab5<-rbind(
  format_full_summary.HM1(m3.bio.noGrowth, ''),
  format_full_summary.HM1(m3.qs.noGrowth, ''),
  format_full_summary.HM1(m3.sec.noGrowth, ''),
  format_full_summary.HM1(m3.ab.noGrowth, ''),
  format_full_summary.HM1(m3.sid.noGrowth, ''),
  format_full_summary.HM1(m3.ssy.noGrowth, ''),
  format_full_summary.HM1(m3.vf.noGrowth, ''),
  format_full_summary.HM1(m3.multi.with_vf.NOgrowthRate, "")
)



tab.s5<- kable(tab5, "latex", booktabs = T,
               caption = 'Summaries of CFR models similar to those described in table S4, but ommiting the growth rate predictor.') %>%
  footnote(c("CI95%: 95% credible interval of the posterior distribution",
             "pMCMC: taken as twice the posterior probability that the estimate is negative",
             "z(): For multivariate model, predictors were z-transformed"),
           fixed_small_size = TRUE, general_title = "") %>%
  kable_styling() %>%
  #pack_rows("Non-cooperative traits", 1, 6) %>%
  pack_rows("Secretome", 1, 7) %>%
  pack_rows("Biofilm", 8, 14) %>%
  pack_rows("Siderophores", 15, 21) %>%
  pack_rows("Antibiotic degradation", 22, 28) %>%
  pack_rows("Quorum sensing", 29, 35) %>%
  pack_rows("Secretion systems", 36, 42) %>%
  pack_rows("Virulence factors", 43, 49) %>%
  pack_rows("Multivariate model", 50, 62)



fileConn<-file("./output/4_summary_tables/TABLE_S5.hm1.tex")
writeLines(tab.s5, fileConn)
close(fileConn)



# Gelman rubin tests ----


do.gelman.vf<- function(focus.model, rm.units){
  
  
  load("./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN1.RData")
  chain1<- get(focus.model)
  
  load("./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN2.RData")
  chain2<- get(focus.model)
  
  load("./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN3.RData")
  chain3<- get(focus.model)
  
  
  # Sol
  combinedchains.Sol = mcmc.list(chain1$Sol, chain2$Sol, chain3$Sol)
  psrf.sol<- gelman.diag(combinedchains.Sol)$psrf[,2] # extracts upper CI of the estimate of potential scale reduction factor
  
  psrf.sol<- gelman.diag(combinedchains.Sol)$psrf %>%
    as.data.frame()%>%
    rownames_to_column(var = 'effect') %>%
    mutate(model = focus.model,
           structure = 'Fixed effects') %>%
    select(model, structure, effect, `Point est.`, `Upper C.I.`)
  
  
  # VCV
  if(rm.units == TRUE){
    # VCV - in this case, unit variance (residual) was fixed. Remove from the evaluated chain
    
    combinedchains.vcv = mcmc.list(chain1$VCV[,1], chain2$VCV[,1], chain3$VCV[,1])
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf[,2]
    
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf %>%
      as.data.frame()%>%
      rownames_to_column(var = 'effect') %>%
      mutate(model = focus.model,
             structure = 'Variances') %>%
      select(model, structure, effect, `Point est.`, `Upper C.I.`)
    
    psrf.vcv[,'effect']<- colnames(chain1$VCV)[1]
    
  }else{
    combinedchains.vcv = mcmc.list(chain1$VCV, chain2$VCV, chain3$VCV)
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf[,2]
    
    
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf %>%
      as.data.frame()%>%
      rownames_to_column(var = 'effect') %>%
      mutate(model = focus.model,
             structure = 'Variances') %>%
      select(model, structure, effect, `Point est.`, `Upper C.I.`)
    
  }
  
  gel.tab<- rbind(psrf.sol, psrf.vcv) %>%
    mutate(structure = ifelse(duplicated(structure) == TRUE, '', structure),
           model = ifelse(duplicated(model) == TRUE, '', model)) 
  
  return(gel.tab)
  
}



do.gelman.patho<- function(focus.model, rm.units){
  
  
  load("./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN1.RData")
  chain1<- get(focus.model)
  
  load("./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN2.RData")
  chain2<- get(focus.model)
  
  load("./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN3.RData")
  chain3<- get(focus.model)
  
  
  # Sol
  combinedchains.Sol = mcmc.list(chain1$Sol, chain2$Sol, chain3$Sol)
  psrf.sol<- gelman.diag(combinedchains.Sol)$psrf[,2] # extracts upper CI of the estimate of potential scale reduction factor
  
  psrf.sol<- gelman.diag(combinedchains.Sol)$psrf %>%
    as.data.frame()%>%
    rownames_to_column(var = 'effect') %>%
    mutate(model = focus.model,
           structure = 'Fixed effects') %>%
    select(model, structure, effect, `Point est.`, `Upper C.I.`)
  
  
  # VCV
  if(rm.units == TRUE){
    # VCV - in this case, unit variance (residual) was fixed. Remove from the evaluated chain
    
    combinedchains.vcv = mcmc.list(chain1$VCV[,1], chain2$VCV[,1], chain3$VCV[,1])
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf[,2]
    
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf %>%
      as.data.frame()%>%
      rownames_to_column(var = 'effect') %>%
      mutate(model = focus.model,
             structure = 'Variances') %>%
      select(model, structure, effect, `Point est.`, `Upper C.I.`)
    
    psrf.vcv[,'effect']<- colnames(chain1$VCV)[1]
    
  }else{
    combinedchains.vcv = mcmc.list(chain1$VCV, chain2$VCV, chain3$VCV)
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf[,2]
    
    
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf %>%
      as.data.frame()%>%
      rownames_to_column(var = 'effect') %>%
      mutate(model = focus.model,
             structure = 'Variances') %>%
      select(model, structure, effect, `Point est.`, `Upper C.I.`)
    
  }
  
  gel.tab<- rbind(psrf.sol, psrf.vcv) %>%
    mutate(structure = ifelse(duplicated(structure) == TRUE, '', structure),
           model = ifelse(duplicated(model) == TRUE, '', model)) 
  
  return(gel.tab)
  
}


do.gelman.cfr<- function(focus.model, rm.units){
  
  
  load("./output/3_model_output/CompAnalysis_cfr_CHAIN1.RData")
  chain1<- get(focus.model)
  
  load("./output/3_model_output/CompAnalysis_cfr_CHAIN2.RData")
  chain2<- get(focus.model)
  
  load("./output/3_model_output/CompAnalysis_cfr_CHAIN3.RData")
  chain3<- get(focus.model)
  
  
  # Sol
  combinedchains.Sol = mcmc.list(chain1$Sol, chain2$Sol, chain3$Sol)
  psrf.sol<- gelman.diag(combinedchains.Sol)$psrf[,2] # extracts upper CI of the estimate of potential scale reduction factor
  
  psrf.sol<- gelman.diag(combinedchains.Sol)$psrf %>%
    as.data.frame()%>%
    rownames_to_column(var = 'effect') %>%
    mutate(model = focus.model,
           structure = 'Fixed effects') %>%
    select(model, structure, effect, `Point est.`, `Upper C.I.`)
  
  
  # VCV
  if(rm.units == TRUE){
    # VCV - in this case, unit variance (residual) was fixed. Remove from the evaluated chain
    
    combinedchains.vcv = mcmc.list(chain1$VCV[,1], chain2$VCV[,1], chain3$VCV[,1])
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf[,2]
    
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf %>%
      as.data.frame()%>%
      rownames_to_column(var = 'effect') %>%
      mutate(model = focus.model,
             structure = 'Variances') %>%
      select(model, structure, effect, `Point est.`, `Upper C.I.`)
    
    psrf.vcv[,'effect']<- colnames(chain1$VCV)[1]
    
  }else{
    combinedchains.vcv = mcmc.list(chain1$VCV, chain2$VCV, chain3$VCV)
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf[,2]
    
    
    psrf.vcv<- gelman.diag(combinedchains.vcv)$psrf %>%
      as.data.frame()%>%
      rownames_to_column(var = 'effect') %>%
      mutate(model = focus.model,
             structure = 'Variances') %>%
      select(model, structure, effect, `Point est.`, `Upper C.I.`)
    
  }
  
  gel.tab<- rbind(psrf.sol, psrf.vcv) %>%
    mutate(structure = ifelse(duplicated(structure) == TRUE, '', structure),
           model = ifelse(duplicated(model) == TRUE, '', model)) 
  
  return(gel.tab)
  
}




# Run the gelman tests for VF analysis
gelman.tab.vf<- rbind(
  do.gelman.vf('meta.a.secretome', rm.units = TRUE),
  do.gelman.vf('meta.a.biofilm', rm.units = TRUE),
  do.gelman.vf('meta.a.sid', rm.units = TRUE),
  do.gelman.vf('meta.a.ab', rm.units = TRUE),
  do.gelman.vf('meta.a.qs', rm.units = TRUE),
  do.gelman.vf('meta.a.ssyst', rm.units = TRUE)
)


# Run the gelman tests for PATHOGENICITY analysis
gelman.tab.pathogenicity<- rbind(
  do.gelman.patho('m1.ss', rm.units = TRUE),
  do.gelman.patho('m1.biofilm', rm.units = TRUE),
  do.gelman.patho('m1.siderophores', rm.units = TRUE),
  do.gelman.patho('m1.ab_degradation', rm.units = TRUE),
  do.gelman.patho('m1.quorum_sensing', rm.units = TRUE),
  do.gelman.patho('m1.secretion_system', rm.units = TRUE),
  do.gelman.patho('m1.vf', rm.units = TRUE),
  do.gelman.patho('m1.multi.with_vf_NOgram', rm.units = TRUE)
)


# Run the gelman tests for CFR analysis --> Something wrong to fix here!
# gelman.tab.cfr<- rbind(
#   do.gelman.cfr('m3.bio', rm.units = FALSE),
#   do.gelman.cfr('m3.qs', rm.units = FALSE),
#   do.gelman.cfr('m3.ab', rm.units = FALSE),
#   do.gelman.cfr('m3.sid', rm.units = FALSE),
#   do.gelman.cfr('m3.ssy', rm.units = FALSE),
#   do.gelman.cfr('m3.sec', rm.units = FALSE),
#   do.gelman.cfr('m3.vf', rm.units = FALSE),
#   do.gelman.cfr('m3.multi.with_vf', rm.units = FALSE)
# )
# 
# 
# gelman.tab.cfr.NoGrowth<- rbind(
#   do.gelman.cfr('m3.bio.noGrowth', rm.units = FALSE),
#   do.gelman.cfr('m3.qs.noGrowth', rm.units = FALSE),
#   do.gelman.cfr('m3.ab.noGrowth', rm.units = FALSE),
#   do.gelman.cfr('m3.sid.noGrowth', rm.units = FALSE),
#   do.gelman.cfr('m3.ssy.noGrowth', rm.units = FALSE),
#   do.gelman.cfr('m3.sec.noGrowth', rm.units = FALSE),
#   do.gelman.cfr('m3.vf.noGrowth', rm.units = FALSE),
#   do.gelman.cfr('m3.multi.with_vf.NOgrowthRate', rm.units = FALSE)
# )




# Format the gelman rubin tests table for manuscript
dgl<- rbind(gelman.tab.vf,gelman.tab.pathogenicity)#,gelman.tab.cfr)

range(dgl$`Upper C.I.`) # max scale reduction factor is 1.017





