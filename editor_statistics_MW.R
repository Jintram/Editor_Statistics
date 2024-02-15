
# See also document 
# /Users/m.wehrens/Documents/Labwork/Digital_journal_ELAB_BIOINF/2024-02_Thomas-Editors-ClinVar.docx
# for notes.

# See also directory:
# /Users/m.wehrens/Data/__resources/ClinVar

library(seqinr) # install.packages('seqinr')

#path_to_genome = "/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/GRCh38.p14/"
path_to_genome = "/Users/m.wehrens/Data_notbacked/references/"

################################################################################
# Some simple functions to convert DNA to (reverse) complement

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

################################################################################
# Load the applicable ClinVar table

ClinVar_Table_CM =
    read.table('/Users/m.wehrens/Data/__resources/ClinVar/clinvar_result_cardiomyopathy_pathogenic.txt', sep='\t', header=1)
View(ClinVar_Table_CM)

ClinVar_Table_CM_SNV = ClinVar_Table_CM[ClinVar_Table_CM$Variant.type=='single nucleotide variant',]

################################################################################
# I (manually) created a little overview of the chromosome files, just to be able
# to load them later easily.

mychromos = 
    read.table('/Users/m.wehrens/Documents/git_repos/Side_Projects/Editor_Statistics_Thomas/chromosome_files_overview.tsv', sep = '\t', header=1)

mychromo_files = mychromos$RefSeq
names(mychromo_files) = mychromos$Chromosome
mychromo_files

################################################################################
# Function to create a table with information

get_statistics_for_editor = function(the_editor, props, the_ClinVar_Table) {
    
    # Update table with which locations can be edited
    # First define the from and tos
    the_ClinVar_Table$MW_nt_WT   = sapply( the_ClinVar_Table$Canonical.SPDI, function(S) {strsplit( S , split=':')[[1]][3] })
    the_ClinVar_Table$MW_nt_var  = sapply( the_ClinVar_Table$Canonical.SPDI, function(S) {strsplit( S , split=':')[[1]][4] })
    
    # Now check where the editor can be applied
    the_ClinVar_Table$candidate_forward=
        (props[[the_editor]]$from == the_ClinVar_Table$MW_nt_var & the_ClinVar_Table$MW_nt_WT == props[[the_editor]]$to)
    the_ClinVar_Table$candidate_revcom=
        (props[[the_editor]]$from == revco[the_ClinVar_Table$MW_nt_var] & revco[the_ClinVar_Table$MW_nt_WT] == props[[the_editor]]$to)
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
        # table(the_ClinVar_Table$candidate_forward+the_ClinVar_Table$candidate_revcom)
    
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
    the_ClinVar_Table$MW_rowidx=1:nrow(the_ClinVar_Table)
    the_ClinVar_Table$PAM_region_fwd=NA
    the_ClinVar_Table$PAM_locations_fwd=NA
    the_ClinVar_Table$PAM_bystander_count_fwd=NA
    the_ClinVar_Table$PAM_region_rev=NA
    the_ClinVar_Table$PAM_locations_rev=NA
    the_ClinVar_Table$PAM_bystander_count_rev=NA
    the_ClinVar_Table$fwd_or_rev='none'
    the_ClinVar_Table$PAM_bystander_count_fwd_min=NA
    the_ClinVar_Table$PAM_bystander_count_rev_min=NA
    last_chromosome='none' # last_chromosome = '1'; current_chromosome='1'
    for (row_idx in 1:nrow(the_ClinVar_Table)) {
        # row_idx=1
        # row_idx=105
        # row_idx=37
        # row_idx=40
    
        # print(paste0(row_idx))
        
        if (the_ClinVar_Table$Canonical.SPDI[row_idx]=="") {seq_skipped=seq_skipped+1; next}
        
        current_chromosome = the_ClinVar_Table$GRCh38Chromosome[row_idx]
        
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
        
        if (the_ClinVar_Table$candidate_forward[row_idx]) {
            # print('testing forward')
            # FORWARD CASE
            
            current_SNV_loc = as.double(the_ClinVar_Table$GRCh38Location[row_idx])
            
            # save which strand this mutation could be edited
            # note that this assumes only one strand will ever work, which
            # might not be true for all editors..
            the_ClinVar_Table$fwd_or_rev[row_idx]='forward'
            
            # First obtain the region where potential PAM is located
            thePAMregion_coordinates = c(( current_SNV_loc + props[[the_editor]]$L_g - Y + 1 ), (current_SNV_loc + props[[the_editor]]$L_g + props[[the_editor]]$L_pam - props[[the_editor]]$X))
            thePAMregion = toupper(paste0( chromosome_data_current[[mychromo_files[current_chromosome]]][
                thePAMregion_coordinates[1]:thePAMregion_coordinates[2]]    , collapse=""))
            # Save that sequence
            the_ClinVar_Table$PAM_region_fwd[row_idx] = thePAMregion
    
            # now annotate whether it's editable given the PAM sequence
            PAM_locations = gregexpr(pattern = props[[the_editor]]$pamregexp, text = thePAMregion)[[1]][1]
            
            # Save hits; these are now relative to the PAM region
            the_ClinVar_Table$PAM_locations_fwd[row_idx] = 
                if (PAM_locations==-1) {NA} else {toString(PAM_locations)}
            
            # now count the number of bystander edits
            # This needs to be done for the region corresponding to the PAM hit
            if (PAM_locations != -1) {
                #for (current_PAM_pos in PAM_locations) {
                bystander_counts=
                    unlist(sapply(PAM_locations, function(current_PAM_pos) {
                        #current_PAM_pos=PAM_locations[1] # testing purposes
                        PAM_loc_start = thePAMregion_coordinates[1]+current_PAM_pos-1
                        guide_locations = (PAM_loc_start-1-props[[the_editor]]$L_g+1):(PAM_loc_start-1)
                        guide_nucleotides = chromosome_data_current[[mychromo_files[current_chromosome]]][guide_locations]
                        # paste0(guide_nucleotides, collapse="")
                        edit_target_sequence = guide_nucleotides[props[[the_editor]]$X:props[[the_editor]]$Y]
                        bystander_count = sum(sapply(edit_target_sequence, toupper)==props[[the_editor]]$from)
                            # Note that the ref sequence contains the wild type sequence.
                        return(bystander_count)
                    }))
                #}
                the_ClinVar_Table$PAM_bystander_count_fwd[row_idx] = toString(bystander_counts)
                the_ClinVar_Table$PAM_bystander_count_fwd_min[row_idx] = min(bystander_counts)
            }
        }
            
        
        if (the_ClinVar_Table$candidate_revcom[row_idx]) {
            
            # row_idx=36
            
            current_SNV_loc = as.double(the_ClinVar_Table$GRCh38Location[row_idx])
            
            # save which strand this mutation could be edited
            # note that this assumes only one strand will ever work, which
            # might not be true for all editors..
            the_ClinVar_Table$fwd_or_rev[row_idx]='reverse'
            
            # Get the region with potential PAM sequences
            thePAMregion_coordinates = c((current_SNV_loc - props[[the_editor]]$L_g - props[[the_editor]]$L_pam + props[[the_editor]]$X),
                                         ( current_SNV_loc - props[[the_editor]]$L_g + props[[the_editor]]$Y - 1 ))
    
            thePAMregion = get_revcom_sequence(
                    toupper(paste0( chromosome_data_current[[mychromo_files[current_chromosome]]][
                        (thePAMregion_coordinates[1]):(thePAMregion_coordinates[2])]    , collapse=""))
                )
            # Save that sequence
            the_ClinVar_Table$PAM_region_rev[row_idx] = thePAMregion
    
            
            # now annotate whether it's editable given the PAM sequence
            PAM_locations = gregexpr(pattern = props[[the_editor]]$pamregexp, text = thePAMregion)[[1]][1]
            
            # Save hits; these are now relative to the PAM region
            the_ClinVar_Table$PAM_locations_rev[row_idx] = 
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
                        guide_locations = (PAM_loc_start+1):(PAM_loc_start+1+props[[the_editor]]$L_g-1)
                        guide_nucleotides = get_revcom_charvector(
                            chromosome_data_current[[mychromo_files[current_chromosome]]][guide_locations])
                        # paste0(guide_nucleotides, collapse="")
                        edit_target_sequence = guide_nucleotides[props[[the_editor]]$X:props[[the_editor]]$Y]
                        bystander_count = sum(sapply(edit_target_sequence, toupper)==props[[the_editor]]$from)
                            # Note that the ref sequence contains the wild type sequence.
                        return(bystander_count)
                    }))
                #}
                the_ClinVar_Table$PAM_bystander_count_rev[row_idx] = toString(bystander_counts)
                the_ClinVar_Table$PAM_bystander_count_rev_min[row_idx] = min(bystander_counts)
            }
            
        } 
     
        if (row_idx %% 100 == 0) {
            #print('.')
            print(paste0(round(100*row_idx / nrow(the_ClinVar_Table),1), '% done'))
        }
           
        last_chromosome = current_chromosome
        
    }
    if (seq_skipped>0) { warning(paste0('Sequences had lacking SPDI, and where skipped! -- ',seq_skipped, ' sequences.')) }
 
    # Let's reorganize the table a little bit
    # Let's add common column for both forward or reverse edit
    # I didn't to this already above, because (a) this keeps 
    # the output a bit more organize, and (b) there might be
    # (prime) editors that can edit anything, so also
    # both rev and fwd; this would require an update of only 
    # the code below ..
    the_ClinVar_Table$PAM_region       = NA
    the_ClinVar_Table$PAM_locations    = NA
    the_ClinVar_Table$bystander_count  = NA
    the_ClinVar_Table$bystander_count_min = NA
    # fwd
    fwd_idxs = the_ClinVar_Table$fwd_or_rev=='forward'
    the_ClinVar_Table$PAM_region[fwd_idxs]       = the_ClinVar_Table$PAM_region_fwd[fwd_idxs]
    the_ClinVar_Table$PAM_locations[fwd_idxs]        = the_ClinVar_Table$PAM_locations_fwd[fwd_idxs]
    the_ClinVar_Table$bystander_count[fwd_idxs]  = the_ClinVar_Table$PAM_bystander_count_fwd[fwd_idxs]
    the_ClinVar_Table$bystander_count_min[fwd_idxs] = the_ClinVar_Table$PAM_bystander_count_fwd_min[fwd_idxs]
    # rev
    rev_idxs = the_ClinVar_Table$fwd_or_rev=='reverse'
    the_ClinVar_Table$PAM_region[rev_idxs]       = the_ClinVar_Table$PAM_region_rev[rev_idxs]
    the_ClinVar_Table$PAM_locations[rev_idxs]        = the_ClinVar_Table$PAM_locations_rev[rev_idxs]
    the_ClinVar_Table$bystander_count[rev_idxs]  = the_ClinVar_Table$PAM_bystander_count_rev[rev_idxs]
    the_ClinVar_Table$bystander_count_min[rev_idxs] = the_ClinVar_Table$PAM_bystander_count_rev_min[rev_idxs]
    
    # Also create a binary output
    the_ClinVar_Table$PAM_present = c('no','yes')[1+1*!is.na(the_ClinVar_Table$PAM_locations)]
    
    the_ClinVar_Table$editor = the_editor
    
    # now return the table of interest
    return(
        the_ClinVar_Table[, c('Name','Canonical.SPDI','editor','fwd_or_rev','PAM_region','PAM_present','PAM_locations','bystander_count','bystander_count_min','MW_rowidx')])
        # View(the_ClinVar_Table[, c('Name','Canonical.SPDI','fwd_or_rev','PAM_region','PAM_locations','bystander_count','MW_rowidx')])
}



################################################################################
# Editor properties

theprops = list()

## 'SAKKH_Abe9'
theprops$SAKKH_Abe9$from = 'A'
theprops$SAKKH_Abe9$to = 'G'
theprops$SAKKH_Abe9$pam = 'NNNRRT'
theprops$SAKKH_Abe9$pamregexp = '...[G|A][G|A]T'
theprops$SAKKH_Abe9$X = 7
theprops$SAKKH_Abe9$Y = 8
theprops$SAKKH_Abe9$L_g   = 22
theprops$SAKKH_Abe9$L_pam = nchar(props$SAKKH_Abe9$pam)

## 'SAKKH_Abe8'
theprops$SAKKH_Abe8$from = 'A'
theprops$SAKKH_Abe8$to = 'G'
theprops$SAKKH_Abe8$pam = 'NNNRRT'
theprops$SAKKH_Abe8$pamregexp = '...[G|A][G|A]T'
theprops$SAKKH_Abe8$X = 3
theprops$SAKKH_Abe8$Y = 14
theprops$SAKKH_Abe8$L_g   = 22
theprops$SAKKH_Abe8$L_pam = nchar(props$SAKKH_Abe8$pam)


## Now obtain statistics
the_Stats_Table=list()

the_Stats_Table$SAKKH_Abe8 =
    get_statistics_for_editor(the_editor = 'SAKKH_Abe8', props = theprops, the_ClinVar_Table = ClinVar_Table_CM_SNV)

the_Stats_Table$SAKKH_Abe9 =
    get_statistics_for_editor(the_editor = 'SAKKH_Abe9', props = theprops, the_ClinVar_Table = ClinVar_Table_CM_SNV)


################################################################################

final_df_editors = rbind(the_Stats_Table$SAKKH_Abe8, the_Stats_Table$SAKKH_Abe9)

final_df_editors$has_bystander = final_df_editors$bystander_count_min>0
final_df_editors$has_bystander_f = NA
final_df_editors$has_bystander_f[final_df_editors$has_bystander] = 'yes'
final_df_editors$has_bystander_f[!final_df_editors$has_bystander] = 'no'
final_df_editors$has_bystander_f = factor(final_df_editors$has_bystander_f, levels=c('yes','no',NA))

View(final_df_editors)

library(ggplot2)

ggplot(final_df_editors, aes(x=editor, fill=PAM_present, color=has_bystander))+
    geom_bar()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(final_df_editors[final_df_editors$PAM_present=='yes',], aes(x=editor, fill=has_bystander_f))+
    geom_bar()+
    scale_fill_manual(values = c('black','grey'))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    
    


