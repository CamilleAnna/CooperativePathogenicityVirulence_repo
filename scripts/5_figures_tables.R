# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# SCRIPT FOR FIGURES AND TABLES #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#local_project_dir=PATH/TO/CLONE/REPO # '/Users/s1687811/Documents/PhD/Research/CooperativePathogenicityVirulence/'

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


# Figure 2 ----
# ... A (phylogeny) ----
df<- read.table('./output/1_processed_tables/2.1_assembled_SPECIES_annotation.txt', header=TRUE, sep = '\t',
                colClasses = c(rep('character', 4),
                               rep('numeric', 9)))


# Make shorter names for plot
df$species_id_short<- gsub('_', ' ', substr(df$species_id,1,nchar(df$species_id)-6))

# Fix some species names to avoid confusion in phylogeny plot due to strains
df[which(df$species_id == 'Acinetobacter_baumannii_56211'),'species_id_short']<- 'Acinetobacter baumannii (56211)'
df[which(df$species_id == 'Acinetobacter_baumannii_57014'),'species_id_short']<- 'Acinetobacter baumannii (57014)'
df[which(df$species_id == 'Bacillus_cereus_57918'),'species_id_short']<- 'Bacillus cereus (57918)'
df[which(df$species_id == 'Bacillus_cereus_57982'),'species_id_short']<- 'Bacillus cereus (57982)'
df[which(df$species_id == 'Mycobacterium_tuberculosis_57144'),'species_id_short']<- df[which(df$species_id == 'Mycobacterium_tuberculosis_57144'),'species']
df[which(df$species_id == 'Porphyromonas_sp_57899'),'species_id_short']<- df[which(df$species_id == 'Porphyromonas_sp_57899'),'species']
df[which(df$species_id == 'Lachnospiraceae_bacterium_62587'),'species_id_short']<- df[which(df$species_id == 'Lachnospiraceae_bacterium_62587'),'species']
df[which(df$species_id == 'Lachnospiraceae_bacterium_51870'),'species_id_short']<- 'Lachnospiraceae bacterium'


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


# PREP HEATMAP DATA
# Set color scales for each cooperative trait for the heatmap
# For each trait, scales goes from light grey to black according to min and max number of genes coding for that trait across the 118 species

heatmapData<-  df[,c(7:13)] # GOS + secretome
#heatmapData[which(df$vf > 100),'vf']<- '100p'
heatmapData <- as.data.frame(sapply(heatmapData, as.character))
rownames(heatmapData) <- df$id


light_dark<- c("lightgrey", "black")
cols.gos.common<- colorRampPalette(light_dark)(length(seq(0, max(df[,7:13], na.rm = TRUE))))
cols.gos.common<- setNames(cols.gos.common, seq(0, max(df[,7:13], na.rm = TRUE)))


cols.biofilm<- colorRampPalette(light_dark)(length(seq(0, max(df[,'biofilm'], na.rm = TRUE))))
cols.biofilm<- setNames(cols.biofilm, paste0('biofilm.',seq(0, max(df[,'biofilm'], na.rm = TRUE))))
heatmapData[,'biofilm']<- paste0('biofilm.',heatmapData[,'biofilm'])

cols.ab_degradation<- colorRampPalette(light_dark)(length(seq(0, max(df[,'ab_degradation'], na.rm = TRUE))))
cols.ab_degradation<- setNames(cols.ab_degradation, paste0('ab_degradation.',seq(0, max(df[,'ab_degradation'], na.rm = TRUE))))
heatmapData[,'ab_degradation']<- paste0('ab_degradation.',heatmapData[,'ab_degradation'])

cols.quorum_sensing<- colorRampPalette(light_dark)(length(seq(0, max(df[,'quorum_sensing'], na.rm = TRUE))))
cols.quorum_sensing<- setNames(cols.quorum_sensing, paste0('quorum_sensing.',seq(0, max(df[,'quorum_sensing'], na.rm = TRUE))))
heatmapData[,'quorum_sensing']<- paste0('quorum_sensing.',heatmapData[,'quorum_sensing'])

cols.secretion_system<- colorRampPalette(light_dark)(length(seq(0, max(df[,'secretion_system'], na.rm = TRUE))))
cols.secretion_system<- setNames(cols.secretion_system, paste0('secretion_system.',seq(0, max(df[,'secretion_system'], na.rm = TRUE))))
heatmapData[,'secretion_system']<- paste0('secretion_system.',heatmapData[,'secretion_system'])

cols.siderophores<- colorRampPalette(light_dark)(length(seq(0, max(df[,'siderophores'], na.rm = TRUE))))
cols.siderophores<- setNames(cols.siderophores, paste0('siderophores.',seq(0, max(df[,'siderophores'], na.rm = TRUE))))
heatmapData[,'siderophores']<- paste0('siderophores.',heatmapData[,'siderophores'])

cols.secretome<- colorRampPalette(light_dark)(length(seq(0, max(df[,'nb_extracellular'], na.rm = TRUE))))
cols.secretome<- setNames(cols.secretome, paste0('secretome.',seq(0, max(df[,'nb_extracellular'], na.rm = TRUE))))
heatmapData[,'nb_extracellular']<- paste0('secretome.',heatmapData[,'nb_extracellular'])

#cols.gram_profile<- c(n = 'gray33', p =  'gray85', gram0 = 'snow')
#cols.vf<- colorRampPalette(c("mistyrose", "firebrick"))(101+1)
#cols.vf<- setNames(cols.vf, c(paste0('vf.',seq(0, 100)), 'vf.100p'))
#heatmapData[,'vf']<- paste0('vf.',heatmapData[,'vf'])

cols.vf<- colorRampPalette(c("mistyrose", "firebrick"))(length(seq(0, max(df[,'vf'], na.rm = TRUE))))
cols.vf<- setNames(cols.vf, paste0('vf.',seq(0, max(df[,'vf'], na.rm = TRUE))))
heatmapData[,'vf']<- paste0('vf.',heatmapData[,'vf'])


#cols.heatmap<- c(cols.gos.common,cols.secretome)

cols.heatmap<- c(cols.biofilm, cols.ab_degradation, cols.siderophores, cols.quorum_sensing,
                 cols.secretion_system, cols.secretome,
                 #cols.gram_profile,
                 cols.vf)


colnames(heatmapData)<- c('Biofilm', 'Antib. degr.', 'QS', 'Siderophores', 'Secr. syst.', 'Secretome', 'VF')

heatmapData<- heatmapData[,c('Secretome', 'Secr. syst.', 'Antib. degr.',
                             'Siderophores', 'Biofilm', 'QS', 
                             'VF')]


# Get the phylogeny
midas.tree<- read.tree('./data/5_phylogeny_files/midas_tree_renamed.newick')
phylogeny<- drop.tip(midas.tree, midas.tree$tip.label[which(!midas.tree$tip.label %in% df$species_id2)])
#phylogeny<- chronopl(phylogeny, lambda = 0)

new.tipslabs<- df$id[match(phylogeny$tip.label, df$species_id2)]
phylogeny$tip.label<- new.tipslabs


# MAKE THE TREE
p<- ggtree(phylogeny, size = 0.3) %<+% df + 
  geom_tiplab(aes(label=species_id_short, col = pathogen),
              align=T,
              offset=0, hjust=0, linesize = 0.2, size = 1.9) +
  scale_color_manual(values = c('dodgerblue','darkorange'))+
  #theme(legend.position = 'none')+
  xlim(0,30)


# PLOT THE HEATMAP
p2<- gheatmap(p, heatmapData, offset = 3, color=NULL, 
              colnames_position="top", 
              colnames_angle=45,
              colnames_offset_y = 0,
              hjust=0, font.size=2, width = 1) +
  scale_fill_manual(values=cols.heatmap)+ #+
  theme(legend.position = 'none')+
  annotate(geom="text", x=0.55, y=120, label="Pathogen", color="darkorange", size = 2, hjust = 0)+
  annotate(geom="text", x=0.55, y=118, label="Commensal", color="dodgerblue", size = 2, hjust = 0)+
  annotate(geom='segment', x = 0, y = 120, xend = 0.5, yend = 120, linetype = 2, color = 'darkorange', size = 0.3)+
  annotate(geom='segment', x = 0, y = 118, xend = 0.5, yend = 118, linetype = 2, color = 'dodgerblue', size = 0.3)

  
pdf('./output/2_figures/Figure2A_phylogeny.pdf', width = (17.3/2.54), height = 1.4*(17.3/2.54))
print(p2+ ylim(0,120)+ xlim(0,15))
dev.off()


p.2A<- p2+ ylim(0,120)+ xlim(0,15)


# ... B (funnels) ----

load("./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN1.RData")


plot.funnel<- function(df, trait, title, model){
  
  df %>%
    filter(cooperative_trait == trait) %>%
    ggplot(., aes(x = logOR, y = 1/se_logOR, col = as.factor(pathogen)))+
    xlab(' ')+ylab(' ')+
    ggtitle(title)+
    #scale_color_manual(values = 'white')+
    geom_rect(ymin = 0, ymax = 5,
              xmin = summary(model)$solutions[,2],
              xmax = summary(model)$solutions[,3],
              fill = 'gray89', color = NA)+
    geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.4)+
    geom_point(size = 0.5, alpha = .5)+
    #geom_point(col = 'black', size = 0.5)+
    scale_color_manual(values = c('dodgerblue', 'darkorange'))+
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

pdf('./output/2_figures/Figure2B_funnelPlots.pdf', width = 3.43, height = 5.2)

grid.arrange(ncol = 2,
             #plot.funnel(df.all, 'is_coop', 'Cooperative (any)', meta.a.coop),
             plot.funnel(df.all2, 'secretome', 'Secretome', meta.a.secretome),
             plot.funnel(df.all2, 'ssyst', 'Secretion-systems', meta.a.ssyst),
             plot.funnel(df.all2, 'ab', 'Antibiotic degradation', meta.a.ab),
             plot.funnel(df.all2, 'sid', 'Siderophores', meta.a.sid),
             plot.funnel(df.all2, 'biofilm', 'Biofilm', meta.a.biofilm),
             plot.funnel(df.all2, 'qs', 'Quorum-sensing', meta.a.qs),

             bottom = textGrob("Species logOR",gp=gpar(fontsize=6,font=2)),
             left = textGrob("1/SE(logOR)",gp=gpar(fontsize=6,font=2), rot = 90))

dev.off()


# ... C (boxplots) ----
load('./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN1.RData')

# boxplots
d.long<- d %>% gather('trait', 'z_value', grep('z_', colnames(d)))

d.long$trait<- gsub('z_biofilm', 'Biofilm', d.long$trait)
d.long$trait<- gsub('z_quorum_sensing', 'QS', d.long$trait)
d.long$trait<- gsub('z_ab_degradation', 'Antib. degr.', d.long$trait)
d.long$trait<- gsub('z_secretion_system', 'Secr. Syst.', d.long$trait)
d.long$trait<- gsub('z_siderophores', 'Siderophores', d.long$trait)
d.long$trait<- gsub('z_nb_extracellular', 'Secretome', d.long$trait)
d.long$trait<- gsub('z_vf', 'VF', d.long$trait)
#d.long$trait<- gsub('z_mp3', 'MP3', d.long$trait)

d.long<- d.long[d.long$trait != 'z_nb_cds',]
d.long$trait<- factor(d.long$trait,
                      levels = c('VF', 'QS', 'Biofilm', 'Siderophores', 'Antib. degr.', 'Secr. Syst.', 'Secretome'))


p.z<- ggplot(d.long[d.long$trait != 'z_nb_cds',], aes(x = z_value, y = trait, fill = pathogen))+
  xlab('Z-number of genes')+ylab('')+
  geom_boxplot(width = 0.5, outlier.size = 0.05, lwd = 0.25)+
  theme_bw()+
  scale_fill_manual(values = c('dodgerblue', 'darkorange'), labels = c('Non-pathogen', 'Pathogen'))+
  theme(axis.title = element_text(face = 'bold', size = 5),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", size = 0.25),
        axis.ticks = element_line(colour = "black", size = 0.25),
        legend.position = c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_blank())


pdf('./output/2_figures/Figure2C_pathogensBoxplot.pdf', width = 3.43, height = 3)
print(p.z)
dev.off()


# Figure 3 ----

theme_results<-   theme_bw()+
  theme(axis.title = element_text(face = 'bold', size = 5),
        axis.text = element_text(size = 5, colour = 'black'),
        #panel.grid = element_line(size = 0.1),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 0.25),
        axis.ticks = element_line(colour = "black", size = 0.25),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 0),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 5, face = 'bold'))


# ... A (virulence factor) ----

load("./output/3_model_output/CompAnalysis_VirulenceFactors_d118_CHAIN1.RData")

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
  ggtitle('Are cooperative products virlence factors?')+
  geom_point()+
  #xlim(-1,5)+
  geom_errorbarh(aes(xmin = `l.95..CI`, xmax = `u.95..CI`), height = 0.1)+
  geom_vline(xintercept = 1, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l.95..CI`, xmax = `u.95..CI`), height = 0.1, col = 'black', size = 0.3)+
  scale_color_manual(values = 'white')+
  theme_results


# ... B (pathogenicity) ----

load('./output/3_model_output/CompAnalysis_pathogens_d118_CHAIN1.RData')

m1.outs.multi<- as.data.frame(summary(m1.multi.with_vf_NOgram)$solutions)[-c(1,2),] %>%
  mutate(trait = c('Biofilm', 'Antib. degr.', 'QS', 'Secr. Syst.', 'Siderophores', 'Secretome', 'VF')) %>%
  #select(model, trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)
  select(trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)


m1.outs.multi$trait<- factor(m1.outs.multi$trait,
                             levels = c('VF', 'QS', 'Biofilm', 'Siderophores', 'Antib. degr.', 'Secr. Syst.', 'Secretome'))


p.patho.meta<- ggplot(m1.outs.multi, aes(x = post.mean, y = trait, col = 'black'))+
  xlab('Estimated regression coefficients\n(multivariate model)')+ylab('')+
  ggtitle('Is cooperation predictive of pathogenicity?')+
  geom_point()+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1, col = 'black', size = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  scale_color_manual(values = 'white')+
  theme_results



# ... C (cfr) ----

load("./output/3_model_output/CompAnalysis_cfr_CHAIN1.RData")

# params of multivariate models
m3.outs.multi<- as.data.frame(summary(m3.multi.with_vf)$solutions)[-c(1:5),] %>%
  mutate(trait = c('Biofilm', 'QS', 'Antib. degr.', 'Secr. Syst.', 'Siderophores', 'Secretome', 'VF')) %>%
  select(trait, post.mean, `l-95% CI`, `u-95% CI`, `eff.samp`, pMCMC)


m3.outs.multi$trait<- factor(m3.outs.multi$trait,
                             levels = c('VF', 'QS', 'Biofilm', 'Siderophores', 'Antib. degr.', 'Secr. Syst.', 'Secretome'))


p.cfr.meta<- ggplot(m3.outs.multi, aes(x = post.mean, y = trait, col = 'black'))+
  xlab('Estimated regression coefficients\n(multivariate model)')+ylab('')+
  ggtitle('Are more cooperative pathogens more virulent?')+
  geom_point()+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1)+
  geom_point(col = 'black', size = 0.5)+
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.1, col = 'black', size = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'darkgrey', size = 0.2)+
  scale_color_manual(values = 'white')+
  theme_results


pdf('./output/2_figures/Figure3_modelResults_main2.pdf', width = (5.8*3)/2.55, height = 5.8/2.55)

(p.vf.meta|p.patho.meta|p.cfr.meta)+
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 7, face = 'bold'))

dev.off()






# Figure 4 (phylogeny CFR) ----

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




# Get the phylogeny
tree<- read.tree('./data/5_phylogeny_files/supfam_cfr_tree.newick')
tree$tip.label
tree<- drop.tip(tree, tree$tip.label[!tree$tip.label %in% d.cfr$supfam_id])
tree$tip.label<- d.cfr$Legget_species[match(tree$tip.label, d.cfr$supfam_id)]


# PREP HEATMAP DATA
df$case_fatality_rate<- round(df$case_fatality_rate) # round CFR for color scale to work

heatmapData<-  df[,c(5,7,16:22)]
heatmapData <- as.data.frame(sapply(heatmapData, as.character))
rownames(heatmapData) <- df$id

light_dark<- c("lightgrey", "black")

cols.biofilm<- colorRampPalette(light_dark)(length(seq(0, max(df[,'biofilm'], na.rm = TRUE))))
cols.biofilm<- setNames(cols.biofilm, paste0('biofilm.',seq(0, max(df[,'biofilm'], na.rm = TRUE))))
heatmapData[,'biofilm']<- paste0('biofilm.',heatmapData[,'biofilm'])

cols.ab_degradation<- colorRampPalette(light_dark)(length(seq(0, max(df[,'ab_degradation'], na.rm = TRUE))))
cols.ab_degradation<- setNames(cols.ab_degradation, paste0('ab_degradation.',seq(0, max(df[,'ab_degradation'], na.rm = TRUE))))
heatmapData[,'ab_degradation']<- paste0('ab_degradation.',heatmapData[,'ab_degradation'])

cols.quorum_sensing<- colorRampPalette(light_dark)(length(seq(0, max(df[,'quorum_sensing'], na.rm = TRUE))))
cols.quorum_sensing<- setNames(cols.quorum_sensing, paste0('quorum_sensing.',seq(0, max(df[,'quorum_sensing'], na.rm = TRUE))))
heatmapData[,'quorum_sensing']<- paste0('quorum_sensing.',heatmapData[,'quorum_sensing'])

cols.secretion_system<- colorRampPalette(light_dark)(length(seq(0, max(df[,'secretion_system'], na.rm = TRUE))))
cols.secretion_system<- setNames(cols.secretion_system, paste0('secretion_system.',seq(0, max(df[,'secretion_system'], na.rm = TRUE))))
heatmapData[,'secretion_system']<- paste0('secretion_system.',heatmapData[,'secretion_system'])

cols.siderophores<- colorRampPalette(light_dark)(length(seq(0, max(df[,'siderophores'], na.rm = TRUE))))
cols.siderophores<- setNames(cols.siderophores, paste0('siderophores.',seq(0, max(df[,'siderophores'], na.rm = TRUE))))
heatmapData[,'siderophores']<- paste0('siderophores.',heatmapData[,'siderophores'])

cols.secretome<- colorRampPalette(light_dark)(length(seq(0, max(df[,'nb_extracellular'], na.rm = TRUE))))
cols.secretome<- setNames(cols.secretome, paste0('secretome.',seq(0, max(df[,'nb_extracellular'], na.rm = TRUE))))
heatmapData[,'nb_extracellular']<- paste0('secretome.',heatmapData[,'nb_extracellular'])


cols.vf<- colorRampPalette(c("mistyrose", "firebrick"))(length(seq(0, max(df[,'vf'], na.rm = TRUE))))
cols.vf<- setNames(cols.vf, paste0('vf.',seq(0, max(df[,'vf'], na.rm = TRUE))))
heatmapData[,'vf']<- paste0('vf.',heatmapData[,'vf'])


cols.case_fatality_rate<- colorRampPalette(c("cadetblue2", "navy"))(length(seq(0, max(df[,'case_fatality_rate'], na.rm = TRUE))))
cols.case_fatality_rate<- setNames(cols.case_fatality_rate, paste0('case_fatality_rate.',seq(0, max(df[,'case_fatality_rate'], na.rm = TRUE))))
heatmapData[,'case_fatality_rate']<- paste0('case_fatality_rate.',heatmapData[,'case_fatality_rate'])



cols.infection_route<- colorRampPalette(c('forestgreen', 'gold', 'plum'))(3)
cols.infection_route<- setNames(cols.infection_route, paste0('infection_route.',c('Skin', 'Ingestion', 'Inhalation')))
heatmapData[,'infection_route']<- paste0('infection_route.',heatmapData[,'infection_route'])




cols.heatmap<- c(cols.case_fatality_rate,cols.infection_route,
                 cols.secretion_system, cols.secretome, cols.siderophores,
                 cols.biofilm, cols.ab_degradation,cols.quorum_sensing,
                 cols.vf)



heatmapData<- heatmapData[,c('nb_extracellular', 'secretion_system', 'ab_degradation',
                             'siderophores', 'biofilm', 'quorum_sensing', 
                             'vf',
                             'infection_route', 'case_fatality_rate')]


colnames(heatmapData)<- c('Secretome', 'Secr. syst.', 'Antib. degr.',
                          'Siderophores', 'Biofilm', 'QS', 
                          'VF',
                          'Infection route', 'CFR')


heatmapData<- heatmapData %>%
  select(-`Infection route`)

library(RColorBrewer)


# MAKE THE TREE
p<- ggtree(tree, size = 0.3) %<+% df + 
  geom_tiplab(aes(label=Legget_species2, col = infection_route),
              align=T,
              offset=0, hjust=0, linesize = 0.2, size = 1.9) +
  scale_color_manual(values = c('darkgoldenrod1','forestgreen', 'plum'))+
  #scale_color_brewer(paletter = 'Set1')
  #theme(legend.position = 'none')+
  xlim(0,1)+
  annotate(geom="text", x=0.04, y=53, label="Ingestion", color="darkgoldenrod1", size = 2, hjust = 0)+
  annotate(geom="text", x=0.04, y=52, label="Inhalation", color="forestgreen", size = 2, hjust = 0)+
  annotate(geom="text", x=0.04, y=51, label="Skin", color="plum", size = 2, hjust = 0)+
  annotate(geom='segment', x = 0, y = 53, xend = 0.03, yend = 53, linetype = 2, color = 'darkgoldenrod1', size = 0.3)+
  annotate(geom='segment', x = 0, y = 52, xend = 0.03, yend = 52, linetype = 2, color = 'forestgreen', size = 0.3)+
  annotate(geom='segment', x = 0, y = 51, xend = 0.03, yend = 51, linetype = 2, color = 'plum', size = 0.3)

p


p2<- gheatmap(p, heatmapData, offset = 0.3, color=NULL, 
              colnames_position="top", 
              colnames_angle=45,
              colnames_offset_y = 0,
              hjust=0, font.size=2, width = 3) +
  scale_fill_manual(values=cols.heatmap)+ #+
  theme(legend.position = 'none')


pdf('./output/2_figures/Figure4_phylogenyCFR.pdf', width = (17.3/2.54), height = 1*(17.3/2.54))
print(p2 + ylim(0, 53)+ xlim(0,1))
dev.off()


# SUPPLEMENTARY MATERIAL ----

# Figure S1 (boxplots raw) ----

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


pdf('./output/2_figures/FigureS1_boxplotRaw2.pdf', width = 4.4, height = 4.4+0.2)
#grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 3)
grid.arrange(p1,p2,p3,p4,p5,p6,p7, ncol = 3)
dev.off()



# Figure S2 (cfr scatterplots) ----


load("./output/3_model_output/CompAnalysis_cfr_CHAIN1.RData")

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



p.cfr.legend<- ggplot(d.cfr, aes(x = biofilm/nb_cds, y = case_fatality_rate, fill = infection_route))+
  geom_point(shape=21, size=pt.size, stroke = pt.lwd)+
  scale_fill_manual(name = 'Infection route', values = c('black', 'grey', 'white'))+
  theme_bw()+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))


cfr.legend <- cowplot::get_legend(p.cfr.legend)



pdf('./output/2_figures/FigureS2_CFRscatterplots.pdf', width = 4.4, height = 4.4+0.2)

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


# Figure S3 (univariate model estimates) ----

# PATHOGENICITY MODEL
load("./output/3_model_output/CompAnalysis_cfr_CHAIN1.RData")
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



pdf('./output/2_figures/FigureS3_modelResults_univariate.pdf', width = (5.8*2)/2.55, height = 4.8/2.55)
grid.arrange(p.patho.uni, p.cfr.uni, 
             ncol = 2)
dev.off()


# SUPP. TABLE 1 (summary VF analysis) ----

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


write.table(summary.tab.vf, './output/4_summary_tables/model_virulenceFactors.csv', sep = '\t', col.names = TRUE,quote = FALSE,  row.names = FALSE)


# SUPP. TABLE 2 (summary PATHOGENICITY analysis) ----


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


write.table(summary.tab.pathogenicity, './output/4_summary_tables/model_pathogenicity.csv', sep = '\t', col.names = TRUE,quote = FALSE,  row.names = FALSE)



# SUPP. TABLE 3 (summary CFR analysis) ----

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


write.table(summary.tab.cfr, './output/4_summary_tables/model_cfr.csv', sep = '\t', col.names = TRUE,quote = FALSE,  row.names = FALSE)





