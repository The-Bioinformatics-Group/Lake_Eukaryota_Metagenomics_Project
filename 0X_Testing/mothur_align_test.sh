## mothur "#align.seqs(candidate=seq.fasta, template=../extractSilva/silva.nr_v132.align)" #command line
set.dir(input=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/, output=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/)
align.seqs(candidate=./clean_and_matched_2/seq.fasta, template=./extractSilva/silva.nr_v132.align)

#filter.seqs(fasta=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.align, trump=.) #Removing leading and trailing dots
dist.seqs(fasta=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.align, output=lt) #creating distance matrix

list.seqs(fasta=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.align)
make.group(fasta=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.align, groups=A)
#count.seqs(name=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.accnos, group=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/groups)
count.seqs(name=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.accnos)
cluster(phylip=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.phylip.dist, count=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.accnos, cutoff=0.10)
make.shared(list=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.phylip.opti_mcc.list, group=/home/emmah/Lake_Eukaryota_Metagenomics_Project/00_Data/mothur_out_2/seq.count_table, label=0.03)
tree.shared(shared=abrecovery.fn.shared)
