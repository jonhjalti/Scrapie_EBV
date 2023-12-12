# Scrapie_EBV
Repository for codes used in simulation for testing method for calculating breeding values for relative scrapie resistance.

simgt_5.f95 Reads a pedigree file and control.txt and similates genotypes for a single locus for the animals in the pedigree.
ruglgt.f95 introduces errors in genotype information. gtrugl_control.txt has control parameters for ruglgt
ruglped.f95 introduces errors to pedigree file. rugl_control.txt has control parameters for ruglped
fromsol_8.f95 reads the solutions files from dmu and creates files with predicted and true values. 

rktest90.DIR is an input file for DMU5 for assumed heritability=0.90
rktest95.DIR is an input file for DMU5 for assumed heritability=0.95
rktest99.DIR is an input file for DMU5 for assumed heritability=0.99

meta1.R contains R codes used in the meta-analysis
Text file ahq_ahq.txt and similar names refering to PRNP genotypes have numbers from case-control studies. 
The columns are: paper, scrapie cases AHQ/AHQ, healthy AHQ/AHQ, scrapie cases ARQ/ARQ, healthy ARQ/ARQ 
