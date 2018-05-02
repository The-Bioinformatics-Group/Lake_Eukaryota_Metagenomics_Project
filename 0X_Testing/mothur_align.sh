## mothur "#align.seqs(candidate=seq.fasta, template=../extractSilva/silva.nr_v132.align)" #command line
set.dir(input=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/, output=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/)
align.seqs(candidate=./clean_and_matched_2/seq.fasta, template=./extractSilva/silva.nr_v132.align)

filter.seqs(fasta=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.align, trump=.) #Removing leading and trailing dots
dist.seqs(fasta=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.filter.fasta) #creating distance matrix