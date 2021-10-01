# CooperativePathogenicityVirulence_repo
 Comparative analysis of cooperative pathogenicity and virulence
 
 # data
 
 Mostly, files used for genomes processing pipeline and cooperative traits annotations.
 
- **0_species_info_files**: collated list of gram stains and list of social GO used to perform cooperative GO annotation. Methods to establish social GO comes from Simonet & McNally 2021 (PNAS)

- **1_patric**: fasta, features and virulence factor tables downloaded from PATRIC, for the 118 microbial organisms included in the analysis. 1.1_patric_release_notes are the complete PATRIC databased tab delimited format list of genomes and their annotations. Downloaded from PATRIC April 2021. 1.2_genomes_human_disease are the list of representative genomes listed by PATRIC under pathogen when filtering for human hosts, as of April 2021.

- **2_midas_files**: files from MIDAS database, used in pipeline

- **3_genomes_processing**: output of genomes processing with PANNZER (GO annotation) and PSORTb (secretome annotation)

- **4_data_cfr_analysis**: all data for the CFR dataset. List fo genomes used  + PATRIC data (fasta/features/virulence factors) + PANNZER + PSORTB output for the CFR dataset

- **5_phylogeny_files**: phylogenies used for comparative analyses

- **6_webofscience_search**: data records from Web of Science search to see patterns of research in microbial cooperation (figure 1 and figure 6). Searched perform September 2021.



 # output
 
All output from scripts in the "script" directory (data assembled, model output RData, summary tables in csv and tex formats, figures in pdf format)

# scripts

- **Script 1**: establish a list of pathogen and non-pathogen Human associated microbes and identify most common strain representative genome in MIDAS database

- **Script 2**: shell script to process genomes and annotate with PANNZER and MOCAT. Likely would need to be adapted for each system. All output of these goes in ./data/3_genomes_processing

- **Script 3**: assembling all output in 3_genomes_processing to make dataset to feed in analysis scripts. Output are tables 2.1, 2.2, 4.1 and 4.2 in ./output/1_processed_tables

- **Script 4**: Statistical analyses for *(4.1)* virulence factor phylogenetic meta-analysis, *(4.2)* phylogenetic binomial response model for comparative analysis of pathogenicity and *(4.3)* phylogenetic gaussian response model comparative analysis of case fatality rate. All output save as RData object in ./output/3_model_output. Output explored and plotted in script 5

- **Script 5**: loads RData object and look at model summaries, produces figures and tables, runs Gelman-Rubin tests

- other: miscellaneous sourced script to load packages and ggplot themes

