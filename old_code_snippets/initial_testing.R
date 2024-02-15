





################################################################################
# How to read stuff from the human reference genome

library(seqinr) # install.packages('seqinr')

# E.g. the first entry is a mutation in PRDM16, in chromosome 1.
3411765 - 3411766
fa_file = '/Users/m.wehrens/Data/__resources/ClinVar/NC_000001.11.fasta'
seq = read.fasta(fa_file, seqtype="DNA", as.string = F)

# Now obtain a sequence from the clinvar
seq$NC_000001.11[3411765:3411766]
# plusminus 10
theSequencePluMin10 = paste0(seq$NC_000001.11[(3411765-10):(3411766+10)], collapse = '')
paste0(seq$NC_000001.11[(3411765-3):(3411766+3)], collapse = '')
theSequencePluMin10
toupper(theSequencePluMin10)



################################################################################
#####
# Let's see how convenient it is to read in the entire genome:
path_to_genome = '/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/GRCh38.p14/'
genomefile     = 'GCF_000001405.40_GRCh38.p14_genomic.fna' 
seq = read.fasta(paste0(path_to_genome, genomefile), seqtype="DNA", as.string = F)




################################################################################
# Just playing aroud with reading/checkgin ref seq data

chromosome_data_current = read.fasta(paste0(path_to_genome, mychromo_files[1], '.fna'))
    # chromosome_data_current$NC_000001.11[1:10]

# So let's now see whether we can get information about the 
SNV_IDX = 1
current_SNV_chr = ClinVar_Table_CM_SNV$GRCh38Chromosome[SNV_IDX]
current_SNV_loc = as.double( ClinVar_Table_CM_SNV$GRCh38Location[SNV_IDX] )
current_SNV_nt_WT   = strsplit( ClinVar_Table_CM_SNV$Canonical.SPDI[SNV_IDX] , split=':')[[1]][3]
current_SNV_nt_var  = strsplit( ClinVar_Table_CM_SNV$Canonical.SPDI[SNV_IDX] , split=':')[[1]][4]

# SNV_IDX = 1 regards the PRDM16 gene
# This is the mutation c.1627C>T
# paste0(chromosome_data_current[[mychromo_files[current_SNV_chr]]][(current_SNV_loc-5):(current_SNV_loc+5)], collapse="")

chromosome_data_current[[mychromo_files[current_SNV_chr]]][current_SNV_loc]

# so let's see whether it has a PAM, for a specific editor



################################################################################

# Look for correct PAM region

# grep gives a hit -->     

names(the_ClinVar_Table)
View(the_ClinVar_Table[, c('Name','Canonical.SPDI','SAKKH-Abe8_candidate_forward','SAKKH-Abe8_PAM_region_fwd','SAKKH-Abe8_PAM_locations_fwd','SAKKH-Abe8_PAM_bystander_count_fwd','MW_rowidx')])

View(the_ClinVar_Table[the_ClinVar_Table$Name=='NM_001103.4(ACTN2):c.355G>A (p.Ala119Thr)',])

################################################################################
