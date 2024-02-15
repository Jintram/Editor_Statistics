
# See also document 
# /Users/m.wehrens/Documents/Labwork/Digital_journal_ELAB_BIOINF/2024-02_Thomas-Editors-ClinVar.docx
# for notes.

library(seqinr) # install.packages('seqinr')

revcom_lookup = c(A='T',T='A',C='G',G='C')
revco = revcom_lookup

get_revcom_charvector = function(charvector) {
    charvector_revcom = 
        rev(revcom_lookup[sapply(charvector,toupper)])
    return(charvector_revcom)
}
get_revcom_sequence = function(sequence) {
    sequence_revcom = paste0(
        revcom_lookup[rev(strsplit(sequence, split="")[[1]])]
        , collapse="")
    return(sequence_revcom)
}
get_complement_sequence = function(sequence) {
    sequence_revcom = paste0(
        revcom_lookup[strsplit(sequence, split="")[[1]]]
        , collapse="")
    return(sequence_revcom)
}

ClinVar_Table_CM =
    read.table('/Users/m.wehrens/Data/__resources/ClinVar/clinvar_result_cardiomyopathy_pathogenic.txt', sep='\t', header=1)
View(ClinVar_Table_CM)

ClinVar_Table_CM_SNV = ClinVar_Table_CM[ClinVar_Table_CM$Variant.type=='single nucleotide variant',]

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


#####
# Let's see how convenient it is to read in the entire genome:
path_to_genome = '/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/GRCh38.p14/'
genomefile     = 'GCF_000001405.40_GRCh38.p14_genomic.fna' 
seq = read.fasta(paste0(path_to_genome, genomefile), seqtype="DNA", as.string = F)


mychromos = 
    read.table('/Users/m.wehrens/Data/__resources/ClinVar/chromosome_files_overview.tsv', sep = '\t', header=1)

mychromo_files = mychromos$RefSeq
names(mychromo_files) = mychromos$Chromosome
mychromo_files

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

# Editor TO DO --> create a list with editor_list[['SAKKH-Abe8']]$FROM='A' ect?
EDITOR_NAME = 'SAKKH-Abe8'
EDITOR_FROM = 'A'
EDITOR_TO = 'G'
MY_PAM = 'NNNRRT'
MY_PAM_GREPL = '...[G|A][G|A]T'
X = 7; Y = 8; L_g = 22; L_pam=nchar(MY_PAM)

# Check whether the mutation is editable
if ((EDITOR_FROM == current_SNV_nt_var & current_SNV_nt_WT == EDITOR_TO) |     # template strand
    (EDITOR_FROM == revco[current_SNV_nt_var] & revco[current_SNV_nt_WT] == EDITOR_TO)){      # rev com strand

    nucleotide_applicable = T 
} else { nucleotide_applicable = F }
    
# Update table with which locations can be edited
# First define the from and tos
ClinVar_Table_CM_SNV$MW_nt_WT   = sapply( ClinVar_Table_CM_SNV$Canonical.SPDI, function(S) {strsplit( S , split=':')[[1]][3] })
ClinVar_Table_CM_SNV$MW_nt_var  = sapply( ClinVar_Table_CM_SNV$Canonical.SPDI, function(S) {strsplit( S , split=':')[[1]][4] })

# Now check where the editor can be applied
ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_candidate_forward')]]=
    (EDITOR_FROM == ClinVar_Table_CM_SNV$MW_nt_var & ClinVar_Table_CM_SNV$MW_nt_WT == EDITOR_TO)
ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_candidate_revcom')]]=
    (EDITOR_FROM == revco[ClinVar_Table_CM_SNV$MW_nt_var] & revco[ClinVar_Table_CM_SNV$MW_nt_WT] == EDITOR_TO)
    # Note that the strandedness is important for genes in the table.
    # For example, entry 105, in the cDNA has a G>A change, 
    # NM_001276345.2(TNNT2):c.891G>A (p.Trp297Ter)
    # whilst the primary strand of the NCBI annotated genome has C>T.
    # NC_000001.11:201359215:C:T
    # However, when looking at which positions are editable, this doesn't matter
    # as long as its consistent, and the "NC_000001.11:201359215:C:T" is consistent
    # with the refseq genome.
    #
    # Sanity check, both should never be true, so no 2 values expected when summing
    # table(ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_candidate_forward')]]+ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_candidate_revcom')]])

if (F) {
    # Convenient to identify the nucleotide of interest
    toupper(paste0( chromosome_data_current[[mychromo_files[current_SNV_chr]]][
            ( current_SNV_loc - 10 ):(current_SNV_loc + 10)]    , collapse=""))
    toupper(paste0( chromosome_data_current[[mychromo_files[current_SNV_chr]]][
            ( current_SNV_loc - 3 ):(current_SNV_loc + 3)]    , collapse=""))    
}


# ===
# Now retrieve the accompanying PAM sequences
# Let's do this with a simple loop for clarity of the code
seq_skipped=0
ClinVar_Table_CM_SNV$MW_rowidx=1:nrow(ClinVar_Table_CM_SNV)
ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_region_fwd')]]=NA
ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_locations_fwd')]]=NA
ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_bystander_count_fwd')]]=NA
ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_region_rev')]]=NA
ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_locations_rev')]]=NA
ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_bystander_count_rev')]]=NA
last_chromosome='none' # last_chromosome = '1'; current_chromosome='1'
for (row_idx in 1:nrow(ClinVar_Table_CM_SNV)) {
    # row_idx=1
    # row_idx=105
    # row_idx=37
    # row_idx=40

    if (ClinVar_Table_CM_SNV$Canonical.SPDI[row_idx]=="") {seq_skipped=seq_skipped+1; next}
    
    current_chromosome = ClinVar_Table_CM_SNV$GRCh38Chromosome[row_idx]
    
    # DEBUGGING FEATURE, run the thing only on chromo 1
    # if (current_chromosome!="1") {print('STOPPING NOW'); break} # just testing stuff REMOVE
    
    # Check whether we should load a new chromosome
    # This assumes the table is sorted by chromosomal locations, as then it
    # will work efficiently
    if (current_chromosome != last_chromosome) {
        new_chromo_file=mychromo_files[current_chromosome]
        print(paste0('Loading new chromosome (',new_chromo_file,')'))
        chromosome_data_current = read.fasta(paste0(path_to_genome, new_chromo_file, '.fna'))
        print('Loading done, getting to work..')
    }
    
    if (ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_candidate_forward')]][row_idx]) {
        # print('testing forward')
        # FORWARD CASE
        
        current_SNV_loc = as.double(ClinVar_Table_CM_SNV$GRCh38Location[row_idx])
        
        # First obtain the region where potential PAM is located
        thePAMregion_coordinates = c(( current_SNV_loc + L_g - Y + 1 ), (current_SNV_loc + L_g + L_pam - X))
        thePAMregion = toupper(paste0( chromosome_data_current[[mychromo_files[current_chromosome]]][
            thePAMregion_coordinates[1]:thePAMregion_coordinates[2]]    , collapse=""))
        # Save that sequence
        ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_region_fwd')]][row_idx] = thePAMregion

        # now annotate whether it's editable given the PAM sequence
        PAM_locations = gregexpr(pattern = MY_PAM_GREPL, text = thePAMregion)[[1]][1]
        
        # Save hits; these are now relative to the PAM region
        ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_locations_fwd')]][row_idx] = 
            if (PAM_locations==-1) {NA} else {toString(PAM_locations)}
        
        # now count the number of bystander edits
        # This needs to be done for the region corresponding to the PAM hit
        if (PAM_locations != -1) {
            #for (current_PAM_pos in PAM_locations) {
            bystander_counts=
                unlist(sapply(PAM_locations, function(current_PAM_pos) {
                    #current_PAM_pos=PAM_locations[1] # testing purposes
                    PAM_loc_start = thePAMregion_coordinates[1]+current_PAM_pos-1
                    guide_locations = (PAM_loc_start-1-L_g+1):(PAM_loc_start-1)
                    guide_nucleotides = chromosome_data_current[[mychromo_files[current_chromosome]]][guide_locations]
                    # paste0(guide_nucleotides, collapse="")
                    edit_target_sequence = guide_nucleotides[X:Y]
                    bystander_count = sum(sapply(edit_target_sequence, toupper)==EDITOR_FROM)
                        # Note that the ref sequence contains the wild type sequence.
                    return(bystander_count)
                }))
            #}
            ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_bystander_count_fwd')]][row_idx] = toString(bystander_counts)
        }
    }
        
    
    if (ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_candidate_revcom')]][row_idx]) {
        
        # row_idx=36
        
        current_SNV_loc = as.double(ClinVar_Table_CM_SNV$GRCh38Location[row_idx])
        
        # Get the region with potential PAM sequences
        thePAMregion_coordinates = c((current_SNV_loc - L_g - L_pam + X),( current_SNV_loc - L_g + Y - 1 ))

        thePAMregion = get_revcom_sequence(
                toupper(paste0( chromosome_data_current[[mychromo_files[current_chromosome]]][
                    (thePAMregion_coordinates[1]):(thePAMregion_coordinates[2])]    , collapse=""))
            )
        # Save that sequence
        ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_region_rev')]][row_idx] = thePAMregion

        
        # now annotate whether it's editable given the PAM sequence
        PAM_locations = gregexpr(pattern = MY_PAM_GREPL, text = thePAMregion)[[1]][1]
        
        # Save hits; these are now relative to the PAM region
        ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_locations_rev')]][row_idx] = 
            if (PAM_locations==-1) {NA} else {toString(PAM_locations)}

        
        # now count the number of bystander edits
        # This needs to be done for the region corresponding to the PAM hit
        # THIS NEEDS TO BE UPDATED FOR REVERSE COORDINATES ... blergh
        if (PAM_locations != -1) {
            #for (current_PAM_pos in PAM_locations) {
            bystander_counts=
                unlist(sapply(PAM_locations, function(current_PAM_pos) {
                    #current_PAM_pos=PAM_locations[1] # testing purposes
                    PAM_loc_start = thePAMregion_coordinates[2]-current_PAM_pos
                    guide_locations = (PAM_loc_start+1):(PAM_loc_start+1+L_g-1)
                    guide_nucleotides = get_revcom_charvector(
                        chromosome_data_current[[mychromo_files[current_chromosome]]][guide_locations])
                    # paste0(guide_nucleotides, collapse="")
                    edit_target_sequence = guide_nucleotides[X:Y]
                    bystander_count = sum(sapply(edit_target_sequence, toupper)==EDITOR_FROM)
                        # Note that the ref sequence contains the wild type sequence.
                    return(bystander_count)
                }))
            #}
            ClinVar_Table_CM_SNV[[paste0(EDITOR_NAME, '_PAM_bystander_count_rev')]][row_idx] = toString(bystander_counts)
        }
        
    } 
 
    if (row_idx %% 100 == 0) {
        #print('.')
        print(paste0(round(100*row_idx / nrow(ClinVar_Table_CM_SNV),1), '% done'))
    }
       
    last_chromosome = current_chromosome
    
}
if (seq_skipped>0) { warning(paste0('Sequences had lacking SPDI, and where skipped! -- ',seq_skipped, ' sequences.')) }


# Now create an output table


# Look for correct PAM region

# grep gives a hit -->     

names(ClinVar_Table_CM_SNV)
View(ClinVar_Table_CM_SNV[, c('Name','Canonical.SPDI','SAKKH-Abe8_candidate_forward','SAKKH-Abe8_PAM_region_fwd','SAKKH-Abe8_PAM_locations_fwd','SAKKH-Abe8_PAM_bystander_count_fwd','MW_rowidx')])

View(ClinVar_Table_CM_SNV[ClinVar_Table_CM_SNV$Name=='NM_001103.4(ACTN2):c.355G>A (p.Ala119Thr)',])

