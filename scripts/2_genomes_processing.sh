# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# This scripts processes the genomes for analysis 
#
# 1) download PATRIC files:
# 		fasta (protein sequences in fasta format) [+ rename sequences for compatibility with PSORTb downstream]
#		features (table of CDS)
# 		vf (virulence factors listed from Victor database)
#
# 2) Run psortB on FASTA
#
# 3) run Pannzer on fasta [in fact ended up running this manually using the web interface]
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Download files  PATRIC      #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # For each genome, get:
  #   - Genomes amino-acids sequences
  #   - Corresponding feature tabs
  #   - Corresponding virulence factors annotation

  # define local_project_dir='path/to/cloned/repository'

  cd $local_project_dir/CooperativePathogenicityVirulence_repo/data/1_patric

  cat $local_project_dir/CooperativePathogenicityVirulence_repo/output/1_processed_tables/1.1_pathogen_commensal_genomes_118.txt | sed '1d' | while read line
  do
  SPECIES=$(echo $line | cut -f 1 -d' ')
  GENOME=$(grep $SPECIES $local_project_dir/CooperativePathogenicityVirulence_repo/output/1_processed_tables/1.1_pathogen_commensal_genomes_118.txt | cut -f 3)

  #FILE=$(echo $line | cut -f 2 -d'-')

  wget -N "ftp://ftp.patricbrc.org/genomes/$GENOME/$GENOME.PATRIC.faa";
  wget -N "ftp://ftp.patricbrc.org/genomes/$GENOME/$GENOME.PATRIC.features.tab";
  wget -N "ftp://ftp.patricbrc.org/genomes/$GENOME/$GENOME.PATRIC.spgene.tab"
  grep 'Virulence Factor' $GENOME.PATRIC.spgene.tab > tmp.2 && head -n 1 $GENOME.PATRIC.spgene.tab > tmp.1 && cat tmp.* > $GENOME.PATRIC.spgene.tab && rm tmp.* # Keep virulence factors only


  mv $GENOME.PATRIC.faa $SPECIES\.fasta
  mv $GENOME.PATRIC.features.tab $SPECIES\.features
  mv $GENOME.PATRIC.spgene.tab $SPECIES\.vf

  done

  # SOME HOUSEKEEPING
  #mkdir fasta
  #mkdir features
  #mkdir virulenceFactors
  mv *.fasta ./fasta
  mv *.features ./features
  mv *.vf ./virulenceFactors


  # SEQUENCES RENAMING
  # cutting sequences names to keep only first part of it, to avoid crash downstream with PSORTb
  cd ./fasta
  ls | while read line
  do
  sed 's/   .*//' $line > temp.txt && mv temp.txt $line
  done


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Run PSORTb - (perso note: code desktop Mac)  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # GENERATE PSORTB COMMANDS
  # defined $installation_dir pointing to installatioon directory of psortb

  installation_dir='/Users/s1687811/Applications/psortb'

  cd $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb
  # mkdir psortb_output

  cat $local_project_dir/CooperativePathogenicityVirulence_repo/output/1_processed_tables/1.1_pathogen_commensal_genomes_118.txt | sed '1d' | while read line
  do


  FILE=$(echo $line | cut -f 1 -d ' ')
  GRAM=$(grep $FILE $local_project_dir/CooperativePathogenicityVirulence_repo/output/1_processed_tables/1.1_pathogen_commensal_genomes_118.txt | cut -f 9)

  if [[ $GRAM == *n* ]]; then echo "$installation_dir/psortb -i $local_project_dir/CooperativePathogenicityVirulence_repo/data/1_patric/fasta/$FILE.fasta -r $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb_output/$FILE -n -o long && rm $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb_output/$FILE/*\.fasta && mv $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb_output/$FILE/*.txt $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb_output/$FILE/$FILE\.psortb.out && echo 'done with $FILE' >> $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb.log" >> psortb_commands.sh;
  elif [[ $GRAM == *p* ]]; then echo "$installation_dir/psortb -i $local_project_dir/CooperativePathogenicityVirulence_repo/data/1_patric/fasta/$FILE.fasta -r $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb_output/$FILE -p -o long && rm $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb_output/$FILE/*\.fasta && mv $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb_output/$FILE/*.txt $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb_output/$FILE/$FILE\.psortb.out && echo 'done with $FILE' >> $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb.log" >> psortb_commands.sh;
  else echo "For $FILE, special gram, use precomputed genome if available" >> SpecialGram.txt;
  fi
  done


  # RUN PSORTB COMMANDS
  
  cd $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb_output

  cat $local_project_dir/CooperativePathogenicityVirulence_repo/output/1_processed_tables/1.1_pathogen_commensal_genomes_118.txt | sed '1d' | while read line
  do
  FILE=$(echo $line | cut -f 1 -d ' ')
  mkdir $FILE
  done

  # you may run into a 'Permission denied' error when running the following, related to github tracking
  # either run in a separate local directory, or unlink directory from github by running
  # rm -rf .git*
  cd $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb
  sudo su # required to run docker image of  psortb if you don't have administrative rights
  sh psortb_commands.sh
  exit


  # housekeeping
  cd $local_project_dir/CooperativePathogenicityVirulence_repo/data/3_genomes_processing/psortb/psortb_output
  mv ./*/*.psortb.out .
  rmdir * # removes NOW EMPTY directories


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   PANNZER                                      #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
	# --> Manual run of Pannzer at: http://ekhidna2.biocenter.helsinki.fi/sanspanz/
	# input with cds FASTA files
	# output stored in ./data/3_genomes_processing/pannzer


