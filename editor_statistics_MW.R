
# See also document 
# /Users/m.wehrens/Documents/Project_files/Projects/SIDE_PROJECTS/2024-02_Thomas-Editors-ClinVar-notes.docx
# for notes.

# See also directory:
# /Users/m.wehrens/Data/__resources/ClinVar

######
# Downloading the ClinVar file:

#Downloaded ClinVar file from:
#https://www.ncbi.nlm.nih.gov/clinvar/
#https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
# 
# Note that we can also filter this database, e.g. looking for cardiomyopathy in the disease/phenotype field. https://www.ncbi.nlm.nih.gov/clinvar?term=cardiomyopathy%5BDisease%2FPhenotype%5D 
# This yields 60990 results, which can be further filtered for “pathogenic” only, yielding 2835 mutations.
# See the file 
# /Users/m.wehrens/Data/__resources/ClinVar/clinvar_result_cardiomyopathy_pathogenic.txt
# This is a slightly different but also convenient file format.
# Example how file looks:
# Name	Gene(s)	Protein change	Condition(s)	Accession	GRCh37Chromosome	GRCh37Location	GRCh38Chromosome	GRCh38Location	VariationID	AlleleID(s)	dbSNP ID	Canonical SPDI	Variant type	Molecular consequence	Germline classification	Germline date last evaluated	Germline review status	Somatic clinical impact	Somatic clinical impact date last evaluated	Somatic clinical impact review status	Oncogenicity classification	Oncogenicity date last evaluated	Oncogenicity review status	
# NM_022114.4(PRDM16):c.1573dup (p.Arg525fs)	PRDM16			VCV000060725	1	3328329 - 3328330	1	3411765 - 3411766	60725	75285	rs886041395	NC_000001.11:3411765:CCCCC:CCCCCC	Duplication	frameshift variant	Pathogenic	Jan 1, 2018	criteria provided, multiple submitters, no conflicts			

######
# Downloading the reference genome
#
#I downloaded GRCh38.p14 (latest version per 13.2.2024)
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
#https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
#
# Command I used to download the chromosome:
# curl --remote-name --remote-time https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
#
# THE REFERENCE GENOME WAS SPLIT INTO SINGLE FILES USING PYTHON --> CHECK MY WORK DISK FOR THESE FILES!

################################################################################

library(seqinr) # install.packages('seqinr')
library(stringr)

#path_to_genome = "/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/GRCh38.p14/"
path_to_genome = "/Users/m.wehrens/Data_notbacked/references/"

output_dir = '/Users/m.wehrens/Data/__other_analyses/ClinVar_Editors_Thomas/'

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

# I couldn't find how to return multiple partially overlapping matching
# in a sensible way, so this is a bit of a hack..

findallmatches_mw = function(thestring, thepattern) {
    hits=c()
    nchar_thepattern=nchar(thepattern)
    for (start_pos in 1:nchar(thestring)){
        if (1 == regexpr(pattern = thepattern, 
                         text = substr(   thestring, start = start_pos, stop = start_pos+nchar_thepattern-1   )))
            {hits=c(hits, start_pos)}
    }
    return(hits)
}

# findallmatches_mw('aabaaa','a')
# findallmatches_mw('aabaaa','aa')
# findallmatches_mw('aabaaa','a.a')


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

# Testing purposes
if (F) {
    the_ClinVar_Table = 
        ClinVar_Table_CM_SNV[ClinVar_Table_CM_SNV$Name %in% 
                                 c('NM_000256.3(MYBPC3):c.2827C>T (p.Arg943Ter)',
                                   'NM_001103.4(ACTN2):c.355G>A (p.Ala119Thr)',
                                   'NM_170707.4(LMNA):c.673C>T (p.Arg225Ter)'),]
    props=theprops
    the_editor='Sp_Cas9_Abe8'
}

get_statistics_for_editor = function(the_editor, props, the_ClinVar_Table) {
    
    # Update table with which locations can be edited
    # First define the from and tos
    # From the table, extract the wild type nucleotide identity
    the_ClinVar_Table$MW_nt_WT   = sapply( the_ClinVar_Table$Canonical.SPDI, function(S) {
        strsplit( S , split=':')[[1]][3] })
    # From the table, extract the variant nucleotide identity
    the_ClinVar_Table$MW_nt_var  = sapply( the_ClinVar_Table$Canonical.SPDI, function(S) {
        strsplit( S , split=':')[[1]][4] })
    
    # Now check where the editor can be applied
    # Note that independently of the strandedness of the gene, the editor can be applied to 
    # both strands
    the_ClinVar_Table$candidate_forward =
        (props[[the_editor]]$from == the_ClinVar_Table$MW_nt_var & # does the wild type nucleotide match the editor's converted nucleotide
             the_ClinVar_Table$MW_nt_WT == props[[the_editor]]$to) # does the variant wilde type nucleotide match the editor's "input" nucleotide
    the_ClinVar_Table$candidate_revcom=
        (props[[the_editor]]$from == revco[the_ClinVar_Table$MW_nt_var] & # as above, but reverse complement
             revco[the_ClinVar_Table$MW_nt_WT] == props[[the_editor]]$to) # as above, but revcom
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
    the_ClinVar_Table$guides_fwd = NA
    the_ClinVar_Table$annotated_sequence_list_fwd = NA
    the_ClinVar_Table$guides_rev = NA
    the_ClinVar_Table$annotated_sequence_list_rev = NA
    last_chromosome='none' # last_chromosome = '1'; current_chromosome='1'
    for (row_idx in 1:nrow(the_ClinVar_Table)) {
        # row_idx=1
        # row_idx=105
        # row_idx=37
        # row_idx=40
        # row_idx=2
    
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
            
            # This code is only here for testing/development to test specific cases
            # Convenient to identify the nucleotide of interest
            if (F) {
                # This just displays the surrounding sequence in the region of the mutation
                # In ape, it is
                # convenient to search for this piece:
                toupper(paste0( chromosome_data_current[[mychromo_files[current_chromosome]]][
                        ( current_SNV_loc - 10 ):(current_SNV_loc + 10)]    , collapse=""))
                # convenient to zoom in using only 6 nts:
                toupper(paste0( chromosome_data_current[[mychromo_files[current_chromosome]]][
                        ( current_SNV_loc - 3 ):(current_SNV_loc + 3)]    , collapse=""))    
            }
            
            # save which strand this mutation could be edited
            # note that this assumes only one strand will ever work, which
            # is probably true for all editors though..
            the_ClinVar_Table$fwd_or_rev[row_idx]='forward'
            
            # First obtain the region where potential PAM is located
            thePAMregion_coordinates = c(( current_SNV_loc + props[[the_editor]]$L_g - props[[the_editor]]$Y + 1 ), 
                                         (current_SNV_loc + props[[the_editor]]$L_g + props[[the_editor]]$L_pam - props[[the_editor]]$X))
            thePAMregion = toupper(paste0( chromosome_data_current[[mychromo_files[current_chromosome]]][
                thePAMregion_coordinates[1]:thePAMregion_coordinates[2]]    , collapse=""))
            # Save that sequence
            the_ClinVar_Table$PAM_region_fwd[row_idx] = thePAMregion
    
            # now annotate whether it's editable given the PAM sequence
            # bugfix! -- gregexpr doesn't return overlapping patterns!
            # PAM_locations = gregexpr(pattern = props[[the_editor]]$pamregexp, text = thePAMregion)[[1]][1]
            # # stri_locate_all_regex(str =thePAMregion, pattern =props[[the_editor]]$pamregexp, overlap=T) # doesn't work either
            # this does work, probably mutch slower ..
            PAM_locations = findallmatches_mw(thestring = thePAMregion, thepattern = props[[the_editor]]$pamregexp)
            
            # Save hits; these are now relative to the PAM region
            the_ClinVar_Table$PAM_locations_fwd[row_idx] = 
                if (is.null(PAM_locations)) {NA} else {toString(PAM_locations)}
            
            # now count the number of bystander edits
            # This needs to be done for the region corresponding to the PAM hit
            if (!is.null(PAM_locations)) {
                #bystander_counts=
                #    unlist(sapply(PAM_locations, function(current_PAM_pos) {
                bystander_count_list = c()
                guide_nucleotides_str_list = c()
                annotated_sequence_list = c()
                for (current_PAM_pos in PAM_locations) {
                        #current_PAM_pos=PAM_locations[1] # testing purposes
                        #current_PAM_pos=PAM_locations[2] # testing purposes
                        PAM_loc_start = thePAMregion_coordinates[1]+current_PAM_pos-1
                        guide_locations = (PAM_loc_start-1-props[[the_editor]]$L_g+1):(PAM_loc_start-1)
                        guide_nucleotides = chromosome_data_current[[mychromo_files[current_chromosome]]][guide_locations]
                        guide_nucleotides_str = paste0(guide_nucleotides, collapse="")
                        edit_target_sequence = guide_nucleotides[props[[the_editor]]$X:props[[the_editor]]$Y]
                        bystander_count = sum(sapply(edit_target_sequence, toupper)==props[[the_editor]]$from)
                            # Note that the ref sequence contains the wild type sequence.
                        
                        bystander_count_list   = c(bystander_count_list, bystander_count)
                        guide_nucleotides_str_list = c(guide_nucleotides_str_list, guide_nucleotides_str)
                        
                        # Create an annotated string for convenience
                        # where's the actual mutation that's targeted?
                        # merge strings
                        guide_left   = substr(guide_nucleotides_str, 1 , props[[the_editor]]$X-1 )
                        edit_target_sequence_str = paste0(edit_target_sequence, collapse="")
                        mut_in_tgt = nchar(edit_target_sequence_str) - current_PAM_pos + 1
                        guide_window = paste0(
                            substr(edit_target_sequence_str, 1, mut_in_tgt-1), 
                            toupper(substr(edit_target_sequence_str,  mut_in_tgt, mut_in_tgt)), # the targeted mutation depends on window shift, determined by PAM position
                            substr(edit_target_sequence_str,  mut_in_tgt+1, length(edit_target_sequence)),collapse="")
                        guide_right  = substr(guide_nucleotides_str, props[[the_editor]]$Y+1 , length(guide_nucleotides) )
                        PAM_ann      = substr(thePAMregion, current_PAM_pos, current_PAM_pos + props[[the_editor]]$L_pam - 1)
                        annotated_sequence = paste0('[',guide_left, '[', guide_window, ']', guide_right,']', PAM_ann)
                        annotated_sequence_list = c(annotated_sequence_list, annotated_sequence)
                }
                #}
                the_ClinVar_Table$PAM_bystander_count_fwd[row_idx] = toString(bystander_count_list)
                the_ClinVar_Table$PAM_bystander_count_fwd_min[row_idx] = min(bystander_count_list)
                the_ClinVar_Table$guides_fwd[row_idx] = toString(guide_nucleotides_str_list)
                the_ClinVar_Table$annotated_sequence_list_fwd[row_idx] = toString(annotated_sequence_list)
            }
        }
            
        
        if (the_ClinVar_Table$candidate_revcom[row_idx]) {
            # REVERSE CASE
            
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
            #PAM_locations = gregexpr(pattern = props[[the_editor]]$pamregexp, text = thePAMregion)[[1]][1]
            PAM_locations = findallmatches_mw(thestring = thePAMregion, thepattern = props[[the_editor]]$pamregexp)
            
            # Save hits; these are now relative to the PAM region
            the_ClinVar_Table$PAM_locations_rev[row_idx] = 
                if (is.null(PAM_locations)) {NA} else {toString(PAM_locations)}
    
            
            # now count the number of bystander edits
            # This needs to be done for the region corresponding to the PAM hit
            # THIS NEEDS TO BE UPDATED FOR REVERSE COORDINATES ... blergh
            if (!is.null(PAM_locations)) {
                
                bystander_count_list   = c()
                guide_nucleotides_str_list = c()
                annotated_sequence_list = c()
                for (current_PAM_pos in PAM_locations) {
                #bystander_counts=
                #    unlist(sapply(PAM_locations, function(current_PAM_pos) {
                        #current_PAM_pos=PAM_locations[1] # testing purposes
                        PAM_loc_start = thePAMregion_coordinates[2]-current_PAM_pos
                        guide_locations = (PAM_loc_start+1):(PAM_loc_start+1+props[[the_editor]]$L_g-1)
                        guide_nucleotides = tolower(get_revcom_charvector(
                            chromosome_data_current[[mychromo_files[current_chromosome]]][guide_locations]))
                        guide_nucleotides_str = paste0(guide_nucleotides, collapse="")
                        edit_target_sequence = guide_nucleotides[props[[the_editor]]$X:props[[the_editor]]$Y]
                        bystander_count = sum(sapply(edit_target_sequence, toupper)==props[[the_editor]]$from)
                            # Note that the ref sequence contains the wild type sequence.
                        #return(bystander_count)
                #    }))
                        bystander_count_list   = c(bystander_count_list, bystander_count)
                        guide_nucleotides_str_list = c(guide_nucleotides_str_list, guide_nucleotides_str)
                        
                        # THIS NEEDS TO BE UPDATED
                        # Create an annotated string for convenience
                        # where's the actual mutation that's targeted?
                        # merge strings
                        guide_left   = substr(guide_nucleotides_str, 1 , props[[the_editor]]$X-1 )
                        edit_target_sequence_str = paste0(edit_target_sequence, collapse="")
                        mut_in_tgt = nchar(edit_target_sequence_str) - current_PAM_pos + 1
                        guide_window = paste0(
                            substr(edit_target_sequence_str, 1, mut_in_tgt-1), 
                            toupper(substr(edit_target_sequence_str,  mut_in_tgt, mut_in_tgt)), # the targeted mutation depends on window shift, determined by PAM position
                            substr(edit_target_sequence_str,  mut_in_tgt+1, length(edit_target_sequence)),collapse="")
                        guide_right  = substr(guide_nucleotides_str, props[[the_editor]]$Y+1 , length(guide_nucleotides) )
                        PAM_ann      = substr(thePAMregion, current_PAM_pos, current_PAM_pos + props[[the_editor]]$L_pam - 1)
                        annotated_sequence = paste0('[',guide_left, '[', guide_window, ']', guide_right,']', PAM_ann)
                        annotated_sequence_list = c(annotated_sequence_list, annotated_sequence)
                }
                
                the_ClinVar_Table$PAM_bystander_count_rev[row_idx] = toString(bystander_count_list)
                the_ClinVar_Table$PAM_bystander_count_rev_min[row_idx] = min(bystander_count_list)
                the_ClinVar_Table$guides_rev[row_idx] = toString(guide_nucleotides_str_list)
                the_ClinVar_Table$annotated_sequence_list_rev[row_idx] = toString(annotated_sequence_list)
            }
            
        } 
     
        if (row_idx %% 100 == 0) {
            #print('.')
            print(paste0(round(100*row_idx / nrow(the_ClinVar_Table),1), '% done'))
        }
           
        last_chromosome = current_chromosome
        
    }
    if (seq_skipped>0) { warning(paste0('Sequences had lacking SPDI, and were skipped! -- ',seq_skipped, ' sequences.')) }
 
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
    the_ClinVar_Table$guides = NA
    the_ClinVar_Table$annotated_sequence_list = NA
    # fwd
    fwd_idxs = the_ClinVar_Table$fwd_or_rev=='forward'
    the_ClinVar_Table$PAM_region[fwd_idxs]       = the_ClinVar_Table$PAM_region_fwd[fwd_idxs]
    the_ClinVar_Table$PAM_locations[fwd_idxs]        = the_ClinVar_Table$PAM_locations_fwd[fwd_idxs]
    the_ClinVar_Table$bystander_count[fwd_idxs]  = the_ClinVar_Table$PAM_bystander_count_fwd[fwd_idxs]
    the_ClinVar_Table$bystander_count_min[fwd_idxs] = the_ClinVar_Table$PAM_bystander_count_fwd_min[fwd_idxs]
    the_ClinVar_Table$guides[fwd_idxs] = the_ClinVar_Table$guides_fwd[fwd_idxs]
    the_ClinVar_Table$annotated_sequence_list[fwd_idxs] = the_ClinVar_Table$annotated_sequence_list_fwd[fwd_idxs]
    # rev
    rev_idxs = the_ClinVar_Table$fwd_or_rev=='reverse'
    the_ClinVar_Table$PAM_region[rev_idxs]       = the_ClinVar_Table$PAM_region_rev[rev_idxs]
    the_ClinVar_Table$PAM_locations[rev_idxs]        = the_ClinVar_Table$PAM_locations_rev[rev_idxs]
    the_ClinVar_Table$bystander_count[rev_idxs]  = the_ClinVar_Table$PAM_bystander_count_rev[rev_idxs]
    the_ClinVar_Table$bystander_count_min[rev_idxs] = the_ClinVar_Table$PAM_bystander_count_rev_min[rev_idxs]
    the_ClinVar_Table$guides[rev_idxs] = the_ClinVar_Table$guides_rev[rev_idxs]
    the_ClinVar_Table$annotated_sequence_list[rev_idxs] = the_ClinVar_Table$annotated_sequence_list_rev[rev_idxs]
    
    
    # Also create a binary output
    the_ClinVar_Table$PAM_present = c('no','yes')[1+1*!is.na(the_ClinVar_Table$PAM_locations)]
    
    the_ClinVar_Table$editor = the_editor
    
    # now return the table of interest
    return(
        the_ClinVar_Table[, c('Name','Canonical.SPDI','editor','fwd_or_rev','PAM_region','PAM_present','PAM_locations','bystander_count','bystander_count_min','MW_rowidx','guides','annotated_sequence_list')])
        # View(the_ClinVar_Table[, c('Name','Canonical.SPDI','fwd_or_rev','PAM_region','PAM_locations','bystander_count','MW_rowidx')])
}



################################################################################
# Editor properties

theprops = list()

## 'SAKKH_Abe8'
theprops$SAKKH_Abe8$from = 'A'
theprops$SAKKH_Abe8$to = 'G'
theprops$SAKKH_Abe8$pam = 'NNNRRT'
theprops$SAKKH_Abe8$pamregexp = '...[G|A][G|A]T' 
theprops$SAKKH_Abe8$X = 3
theprops$SAKKH_Abe8$Y = 14
theprops$SAKKH_Abe8$L_g   = 22
theprops$SAKKH_Abe8$L_pam = nchar(theprops$SAKKH_Abe8$pam)
theprops$SAKKH_Abe8$dualorsingleAAV = 'single'

## 'eNmE-C-Abe8'
theprops$eNmE_C_Abe8$from = 'A'
theprops$eNmE_C_Abe8$to = 'G'
theprops$eNmE_C_Abe8$pam = 'NNNNCN'
theprops$eNmE_C_Abe8$pamregexp = '....C.'
theprops$eNmE_C_Abe8$X = 3
theprops$eNmE_C_Abe8$Y = 15
theprops$eNmE_C_Abe8$L_g   = 24
theprops$eNmE_C_Abe8$L_pam = nchar(theprops$eNmE_C_Abe8$pam)
theprops$eNmE_C_Abe8$dualorsingleAAV = 'single'

## 'SpRY-Abe8'
theprops$SpRY_Abe8$from = 'A'
theprops$SpRY_Abe8$to = 'G'
theprops$SpRY_Abe8$pam = 'NNN'
theprops$SpRY_Abe8$pamregexp = '...' # rather redundant but anyways
theprops$SpRY_Abe8$X = 3
theprops$SpRY_Abe8$Y = 11
theprops$SpRY_Abe8$L_g   = 20
theprops$SpRY_Abe8$L_pam = nchar(theprops$SpRY_Abe8$pam)
theprops$SpRY_Abe8$dualorsingleAAV = 'dual'

## 'Sp_Cas9_Abe8'
theprops$Sp_Cas9_Abe8$from = 'A'
theprops$Sp_Cas9_Abe8$to = 'G'
theprops$Sp_Cas9_Abe8$pam = 'NGG'
theprops$Sp_Cas9_Abe8$pamregexp = '.GG' # rather redundant but anyways
theprops$Sp_Cas9_Abe8$X = 3
theprops$Sp_Cas9_Abe8$Y = 11
theprops$Sp_Cas9_Abe8$L_g   = 20
theprops$Sp_Cas9_Abe8$L_pam = nchar(theprops$Sp_Cas9_Abe8$pam)
theprops$Sp_Cas9_Abe8$dualorsingleAAV = 'dual'

## New versions from Thomas

## 'SAKKH_Abe9'
theprops$SAKKH_Abe9$from = 'A'
theprops$SAKKH_Abe9$to = 'G'
theprops$SAKKH_Abe9$pam = 'NNNRRT'
theprops$SAKKH_Abe9$pamregexp = '...[G|A][G|A]T'
theprops$SAKKH_Abe9$X = 7
theprops$SAKKH_Abe9$Y = 8
theprops$SAKKH_Abe9$L_g   = 22
theprops$SAKKH_Abe9$L_pam = nchar(theprops$SAKKH_Abe9$pam)
theprops$SAKKH_Abe9$dualorsingleAAV = 'single'

# Not yet information about this one
# ## 'eNmE_C_Abe9'
# theprops$eNmE_C_Abe9$from = 'A'
# theprops$eNmE_C_Abe9$to = 'G'
# theprops$eNmE_C_Abe9$pam = 'NNNNCN'
# theprops$eNmE_C_Abe9$pamregexp = '....C.'
# theprops$eNmE_C_Abe9$X = 3
# theprops$eNmE_C_Abe9$Y = 15
# theprops$eNmE_C_Abe9$L_g   = 24
# theprops$eNmE_C_Abe9$L_pam = nchar(props$eNmE_C_Abe9$pam)
# theprops$eNmE_C_Abe9$dualorsingleAAV = 'single'

## 'SpRY_Abe9'
theprops$SpRY_Abe9$from = 'A'
theprops$SpRY_Abe9$to = 'G'
theprops$SpRY_Abe9$pam = 'NNN'
theprops$SpRY_Abe9$pamregexp = '...' # rather redundant but anyways
theprops$SpRY_Abe9$X = 5
theprops$SpRY_Abe9$Y = 6
theprops$SpRY_Abe9$L_g   = 20
theprops$SpRY_Abe9$L_pam = nchar(theprops$SpRY_Abe9$pam)
theprops$SpRY_Abe9$dualorsingleAAV = 'dual'

################################################################################
## Now obtain statistics


the_Stats_Table=list()

# Automatically loop over all editors
for (the_editor in names(theprops)) {
    
    # the_editor = "Sp_Cas9_Abe8"
    
    print('==============================')
    print(paste0('Working on editor ', the_editor))
    
    # And collect the statistics
    the_Stats_Table[[the_editor]] =
        get_statistics_for_editor(the_editor = the_editor, props = theprops, the_ClinVar_Table = ClinVar_Table_CM_SNV)
    
}

# format(Sys.time(), "%a %b %d %X %Y")

current_time=format(Sys.time(), "%Y-%b-%d_%X")
saveRDS(object = the_Stats_Table, file = paste0(output_dir, 'alleditors__the_Stats_Table__',current_time,'.Rds'))
    # the_Stats_Table = readRDS("/Users/m.wehrens/Data/__other_analyses/ClinVar_Editors_Thomas/alleditors__the_Stats_Table__2024-Mar-01_17-31-13.Rds")

################################################################################

# final_df_editors = rbind(the_Stats_Table$SAKKH_Abe8, the_Stats_Table$SAKKH_Abe9)

final_df_editors = Reduce( f = rbind, x = the_Stats_Table)

final_df_editors$has_bystander = final_df_editors$bystander_count_min>0
final_df_editors$has_bystander_f = NA
final_df_editors$has_bystander_f[final_df_editors$has_bystander] = 'yes'
final_df_editors$has_bystander_f[!final_df_editors$has_bystander] = 'no'
final_df_editors$has_bystander_f = factor(final_df_editors$has_bystander_f, levels=c('yes','no',NA))

# Just renaming for convenience when plotting
final_df_editors$bystanders = final_df_editors$has_bystander_f

View(final_df_editors)

library(ggplot2)

# ggplot(final_df_editors, aes(x=editor, fill=PAM_present, color=has_bystander))+
#     geom_bar()+
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p = ggplot(final_df_editors[final_df_editors$PAM_present=='yes',], aes(x=editor, fill=bystanders))+
        geom_bar()+
        scale_fill_manual(values = c('black','grey'))+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p
ggsave(plot=p, filename = paste0(output_dir, 'some_editor_stats.pdf'), width = 171/2, height = 171/2, units = 'mm', device = cairo_pdf)
    


