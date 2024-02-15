
# Script from
# Lavrov, A. V., Varenikov, G.G., and Skoblov, M.Y. (2020). Genome scale analysis of pathogenic variants targetable for single base editing. BMC Med. Genomics 13, 1–7. 10.1186/s12920-020-00735-8.



#function find only pathogenic or likely pathogenic mutations
#example preprocessing(read.xlsx(xlsx_file,sheet = 1), 'T>C', 'A>G', answer, 't', 'a')
#preprocessing[1] - data of mut_one
#preprocessing[2] - data of mut_two
#preprocessing[3] - numbers of all mutations in list, only mut_one and mut_two
#library("openxlsx")
#library("seqinr")




preprocessing <- function(mutations ,mut_one, mut_two, answer, nucleotide1, nucleotide2){
  data_need_mut <-  rbind(mutations[mutations['name'] == mut_one,],   #mut_one = mutation, which we want to edit
                          mutations[mutations['name'] == mut_two,])   #mut_two = complement mutation, which we want to edit
  
  answer = c(answer,nrow(data_need_mut))
  
  cic <-c()
  i <- 1
  while (i <= nrow(data_need_mut)){
    if ( grepl("athogenic", data_need_mut$clinSign[i]) == T   &&     grepl("Conflict", data_need_mut$clinSign[i]) == F  ){
      cic = c(cic, i)
    }
    i <- i + 1
  }
  data_need_mut = data_need_mut[cic, ]
  
  
  i <- 1
  while ( i <= nrow(data_need_mut)){
    data_need_mut$nucleotide[i] <- seq[[1]][data_need_mut$chromEnd[i]:data_need_mut$chromEnd[i]]
    i <- i + 1
  }
  
  
  data_need_mut["MedGen"] <- NA
  i <- 1
  while(i<=nrow(data_need_mut)) {
    place_medgen <- gregexpr(pattern ='MedGen:',data_need_mut$phenotype[i])
    list_medgen <- c(substr(data_need_mut$phenotype[i], place_medgen[[1]][1]+7, place_medgen[[1]][1]+14))
    g <- 2
    while(g<=length(place_medgen[[1]])) {
      list_medgen <- append(list_medgen, substr(data_need_mut$phenotype[i], place_medgen[[1]][g]+7, place_medgen[[1]][g]+14) )
      g <- g + 1
    }
    data_need_mut$MedGen[i] <- paste(list_medgen, collapse=" ")
    i <- i + 1
  }
  
  
  answer = c(answer, nrow(data_need_mut[data_need_mut['nucleotide'] == nucleotide1,]),nrow(data_need_mut[data_need_mut['nucleotide'] == nucleotide2,]) )
  return(list(data_need_mut[data_need_mut['nucleotide'] == nucleotide1,],data_need_mut[data_need_mut['nucleotide'] == nucleotide2,],answer))
  
}
  
  
  
#find only pathogenic mutations with PAM in the editing window
#'T>C'


  
find_PAM_C <- function(data, start, finish, PAM, type_PAM, len_pam, nucleotide, length_RNA ,answer){
  data['type_PAM'] <- type_PAM
  data['PAM_search'] <- NA
  data['PAM'] <- NA
  data['sequence'] <- NA
  data['frequency'] <- NA
  data['mut_location'] <- NA
  data['RNA'] <- NA
  data['PAM_seq'] <- NA
  
  
  i <- 1
  while ( i <= nrow(data)){
    start_wind_PAM =  data$chromEnd[i] - length_RNA + start
    finish_wind_PAM =  start_wind_PAM + len_pam + (finish - start - 1)
    data$PAMsearch[i] <- gsub(", ","",toString(seq[[1]][start_wind_PAM:finish_wind_PAM]))
    data$PAM[i] <- regexpr( PAM, data$PAMsearch[i])
    i <- i + 1
  }
  data <- data[data$PAM != -1,]
  answer = c(answer,nrow(data))
  paam = data
  
  data$PAM <- sapply(data$PAM, as.numeric)
  i <- 1
  while (i <= nrow(data)){
    zum <- 0
    #seque - if we find PAM, delete left part of seq to find more PAMs
    seque <- substr(data$PAMsearch[i], (data$PAM[i] + 1), nchar(data$PAMsearch[i]))
    if( regexpr(PAM, seque) != -1 ){
      while (nchar(seque) >= len_pam && regexpr(PAM, seque) != -1 ){
        zum <- zum + regexpr(PAM, seque)
        data$PAM[i] <- data$PAM[i] + regexpr(PAM, seque)
        seque <- substr(data$PAMsearch[i], (data$PAM[i] + 1), nchar(data$PAMsearch[i]))
        data[nrow(data) + 1,] = data[c(i), ]
      }
    }
    data$PAM[i] <- data$PAM[i] - zum
    i <-i+1
  }
  data = unique(data)
  answer = c(answer,nrow(data))
  
  i <- 1
  while ( i <= nrow(data)){
    loc_start = (data$chromEnd[i] - length_RNA + start) + (data$PAM[i] - 1)
    clothe = loc_start + length_RNA  - start
    open = clothe - (finish - start)
    data$sequence[i] <- gsub(", ","",toString(seq[[1]][open:clothe]))
    data$frequency[i] <- length(grep(nucleotide, data$sequence[i]))
    tolower(data$sequence[i])
    data$mut_location[i] <- length_RNA - (start - 1) - (data$PAM[i] - 1)
    i <- i + 1
  }
  if (nrow(data)>0){
    data = data[data['frequency'] == 0,]
  }
  else{
    data = data[c(0),]
  }
  answer = c(answer,length(unique(data$chromEnd)), nrow(data))
  
  i <- 1
  while (i <= nrow(data)){
    open <- (data$chromEnd[i] - length_RNA + start) + (data$PAM[i] - 1)
    clothe <- open + (length_RNA - 1)
    data$RNA[i] <- gsub(", ","",toString(seq[[1]][open:clothe]))
    data$PAM_seq[i] <- toupper(gsub(", ","",toString(seq[[1]][open:(open + (len_pam - 1) )])))
    i <- i + 1
  }
  data <- data[c(1,3,4,14,17,18,27,35,36,39,41,42,43)]
  return(list(data,answer,paam))
}


#find only pathogenic mutations with PAM in the editing window
#'A>G'


find_PAM_G <- function(data, start, finish, PAM, type_PAM, len_pam, nucleotide, length_RNA ,answer){
  data['type_PAM'] <- type_PAM
  data['PAM_search'] <- NA
  data['PAM'] <- NA
  data['sequence'] <- NA
  data['frequency'] <- NA
  data['mut_location'] <- NA
  data['RNA'] <- NA
  data['PAM_seq'] <- NA
  i <- 1
  while ( i <= nrow(data)){
    finish_wind_PAM = data$chromEnd[i] + length_RNA - start
    start_wind_PAM = finish_wind_PAM - len_pam - (finish - start - 1)
    data$PAMsearch[i] <- gsub(", ","",toString(seq[[1]][start_wind_PAM:finish_wind_PAM]))
    data$PAM[i] <- regexpr( PAM, data$PAMsearch[i])
    i = i + 1
  }
  data <- data[data$PAM != -1,]
  answer = c(answer,nrow(data))
  data$PAM <- sapply(data$PAM, as.numeric)
  paam = data
  
  i <- 1
  while (i <= nrow(data)){
    zum <- 0
    #seque - if we find PAM, delete left part of seq to find more PAMs
    seque <- substr(data$PAMsearch[i], (data$PAM[i] + 1), nchar(data$PAMsearch[i]))
    if( regexpr(PAM, seque) != -1 ){
      while (nchar(seque) >= len_pam && regexpr(PAM, seque) != -1 ){
        zum <- zum + regexpr(PAM, seque)
        data$PAM[i] <- data$PAM[i] + regexpr(PAM, seque)
        seque <- substr(data$PAMsearch[i], (data$PAM[i] + 1), nchar(data$PAMsearch[i]))
        data[nrow(data) + 1,] = data[c(i), ]
      }
    }
    data$PAM[i] <- data$PAM[i] - zum
    i <-i+1
  }
  data = unique(data)
  answer = c(answer,nrow(data))
  
  i <- 1
  while ( i <= nrow(data)){
    loc_start = (data$chromEnd[i] + length_RNA - start) - ( nchar(data$PAMsearch[i]) - data$PAM[i] - (len_pam - 1) ) 
    open = loc_start -  length_RNA + start
    clothe = open + (finish - start)
    data$sequence[i] <- gsub(", ","",toString(seq[[1]][open:clothe]))
    data$frequency[i] <- length(grep(nucleotide, data$sequence[i]))
    tolower(data$sequence[i])
    data$mut_location[i] <- (finish + 1) - data$PAM[i]
    i <- i + 1
  }
  if (nrow(data)){
    data = data[data['frequency'] == 0,]
  }
  else{
    data = data[c(0),]
  }
  answer = c(answer,length(unique(data$chromEnd)), nrow(data))
  
  i <- 1
  while (i <= nrow(data)){
    clothe <- (data$chromEnd[i] + length_RNA - start) - ( nchar(data$PAMsearch[i]) - data$PAM[i] - (len_pam - 1) ) 
    open <- clothe - length_RNA + 1
    data$RNA[i] <- gsub(", ","",toString(seq[[1]][open:clothe]))
    data$PAM_seq[i] <-  toupper(gsub(", ","",toString(seq[[1]][ (clothe - (len_pam - 1)) :clothe])))
    i <- i + 1
  }
  data <- data[c(1,3,4,14,17,18,27,35,36,39,41,42,43)]
  return(list(data,answer,paam))
}



variables = c(21) #c(1:22,'Y','X')
for (f in variables){
  fa_file = paste("chr", f,".fa", sep="")
  xlsx_file = paste("chr", f,".xlsx", sep="")
  mutations <- read.xlsx( xlsx_file, sheet = 1)
  seq = read.fasta(fa_file, seqtype="DNA", as.string = F)
  
  answer = c(nrow(mutations))
  results = data.frame(c('all_mutations', 'pathogenic', 'T', 'C' ,'T_with_PAM','allseq_T_with_PAM','T_with_PAM_and_freq0','all_T_with_PAM_and_freq0','C_with_PAM','allseq_C_with_PAM','C_with_PAM_and_freq0','all_C_with_PAM_and_freq0','number_of_mutations'))
  data_mut_G = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[1]] #only G mutations
  data_mut_C = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[2]] #only C mutations
  data_mut_A = preprocessing(mutations, 'G>A','C>T',answer,'g','c')[[1]] #only A mutations
  data_mut_T = preprocessing(mutations, 'G>A','C>T',answer,'g','c')[[2]] #only T mutations
  
  #APOBEC
  answer = c(nrow(mutations)) 
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first = find_PAM_G(data_mut_G,3,10,'[atgc]gg', 'APOBEC', 3,'c',23,answer)[[1]] 
  answer = find_PAM_G(data_mut_G,3,10,'[atgc]gg', 'APOBEC', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,3,10,'cc[atgc]', 'APOBEC', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,3,10,'cc[atgc]', 'APOBEC', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$APOBEC = answer                                                               #number found mutations at each step
  first_pam = find_PAM_G(data_mut_G,3,10,'[atgc]gg', 'APOBEC', 3,'c',23,answer)[[3]]   
  second_pam = find_PAM_C(data_mut_C,3,10,'cc[atgc]', 'APOBEC', 3,'g',23,answer)[[3]] 
  pam_data = rbind(first_pam,second_pam)
  
  #PmCDA1
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first = find_PAM_G(data_mut_G,1,5,'[atgc]gg', 'PmCDA1', 3,'c',23,answer)[[1]]
  answer = find_PAM_G(data_mut_G,1,5,'[atgc]gg', 'PmCDA1', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,1,5,'cc[atgc]', 'PmCDA1', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,1,5,'cc[atgc]', 'PmCDA1', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$PmCDA1 = answer
  first_pam = find_PAM_G(data_mut_G,1,5,'[atgc]gg', 'PmCDA1', 3,'c',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,1,5,'cc[atgc]', 'PmCDA1', 3,'g',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  #TadA
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'G>A','C>T',answer,'g','c')[[3]]
  first =  find_PAM_G(data_mut_A,4,9,'[atgc]gg','TadA',3,'a',23,answer)[[1]]
  answer = find_PAM_G(data_mut_A,4,9,'[atgc]gg','TadA',3,'a',23,answer)[[2]]
  second = find_PAM_C(data_mut_T,4,9,'cc[atgc]','TadA',3,'t',23,answer)[[1]]
  answer = find_PAM_C(data_mut_T,4,9,'cc[atgc]','TadA',3,'t',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$TadA = answer
  first_pam = find_PAM_G(data_mut_A,4,9,'[atgc]gg','TadA',3,'a',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_T,4,9,'cc[atgc]','TadA',3,'t',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #xCas9泡BE(NGV)
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'G>A','C>T',answer,'g','c')[[3]]
  first =  find_PAM_G(data_mut_A,4,8,'[atgc]g[acg]','xCas9泡BE(NGV)',3,'a',23,answer)[[1]]
  answer = find_PAM_G(data_mut_A,4,8,'[atgc]g[acg]','xCas9泡BE(NGV)',3,'a',23,answer)[[2]]
  second = find_PAM_C(data_mut_T,4,8,'[ctg]c[agc]','xCas9泡BE(NGV)',3,'t',23,answer)[[1]]
  answer = find_PAM_C(data_mut_T,4,8,'[ctg]c[agc]','xCas9泡BE(NGV)',3,'t',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$xCas9_ABE_NGV = answer
  first_pam = find_PAM_G(data_mut_A,4,8,'[atgc]g[acg]','xCas9泡BE(NGV)',3,'a',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_T,4,8,'[ctg]c[agc]','xCas9泡BE(NGV)',3,'t',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #xCas9泡BE(GAT)
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'G>A','C>T',answer,'g','c')[[3]]
  first =  find_PAM_G(data_mut_A,4,8,'gat','xCas9泡BE(GAT)',3,'a',23,answer)[[1]]
  answer = find_PAM_G(data_mut_A,4,8,'gat','xCas9泡BE(GAT)',3,'a',23,answer)[[2]]
  second = find_PAM_C(data_mut_T,4,8,'atc','xCas9泡BE(GAT)',3,'t',23,answer)[[1]]
  answer = find_PAM_C(data_mut_T,4,8,'atc','xCas9泡BE(GAT)',3,'t',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$xCas9_ABE_GAT = answer
  first_pam = find_PAM_G(data_mut_A,4,8,'gat','xCas9泡BE(GAT)',3,'a',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_T,4,8,'atc','xCas9泡BE(GAT)',3,'t',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #dCpf1-eBE
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first = find_PAM_G(data_mut_C,11,16,'[tgc]aaa', 'dCpf1-eBE', 4,'g',24,answer)[[1]]
  answer = find_PAM_G(data_mut_C,11,16,'[tgc]aaa', 'dCpf1-eBE', 4,'g',24,answer)[[2]]
  second = find_PAM_C(data_mut_G,11,16,'ttt[agc]', 'dCpf1-eBE', 4,'c',24,answer)[[1]]
  answer = find_PAM_C(data_mut_G,11,16,'ttt[agc]', 'dCpf1-eBE', 4,'c',24,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$dCpf1_eBE = answer
  first_pam = find_PAM_G(data_mut_C,11,16,'[tgc]aaa', 'dCpf1-eBE', 4,'g',24,answer)[[3]]
  second_pam = find_PAM_C(data_mut_G,11,16,'ttt[agc]', 'dCpf1-eBE', 4,'c',24,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #dCpf1-eBE-YE
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first = find_PAM_G(data_mut_C,13,15,'[tgc]aaa', 'dCpf1-eBE-YE', 4,'g',24,answer)[[1]]
  answer = find_PAM_G(data_mut_C,13,15,'[tgc]aaa', 'dCpf1-eBE-YE', 4,'g',24,answer)[[2]]
  second = find_PAM_C(data_mut_G,13,15,'ttt[agc]', 'dCpf1-eBE-YE', 4,'c',24,answer)[[1]]
  answer = find_PAM_C(data_mut_G,13,15,'ttt[agc]', 'dCpf1-eBE-YE', 4,'c',24,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$dCpf1_eBE_YE = answer
  first_pam = find_PAM_G(data_mut_C,13,15,'[tgc]aaa', 'dCpf1-eBE-YE', 4,'g',24,answer)[[3]]
  second_pam = find_PAM_C(data_mut_G,13,15,'ttt[agc]', 'dCpf1-eBE-YE', 4,'c',24,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  #A-BE3
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first = find_PAM_G(data_mut_G,4,9,'[atgc]gg', 'A-BE3', 3,'c',23,answer)[[1]]
  answer = find_PAM_G(data_mut_G,4,9,'[atgc]gg', 'A-BE3', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,4,9,'cc[atgc]', 'A-BE3', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,4,9,'cc[atgc]', 'A-BE3', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$A_BE3 = answer
  first_pam = find_PAM_G(data_mut_G,4,9,'[atgc]gg', 'A-BE3', 3,'c',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,4,9,'cc[atgc]', 'A-BE3', 3,'g',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  #Y-BE3
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first = find_PAM_G(data_mut_G,4,8,'[atgc]gg', 'Y-BE3', 3,'c',23,answer)[[1]]
  answer = find_PAM_G(data_mut_G,4,8,'[atgc]gg', 'Y-BE3', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,4,8,'cc[atgc]', 'Y-BE3', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,4,8,'cc[atgc]', 'Y-BE3', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$Y_BE3 = answer
  first_pam = find_PAM_G(data_mut_G,4,8,'[atgc]gg', 'Y-BE3', 3,'c',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,4,8,'cc[atgc]', 'Y-BE3', 3,'g',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #FE-BE3
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first = find_PAM_G(data_mut_G,5,7,'[atgc]gg', 'FE-BE3', 3,'c',23,answer)[[1]]
  answer = find_PAM_G(data_mut_G,5,7,'[atgc]gg', 'FE-BE3', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,5,7,'cc[atgc]', 'FE-BE3', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,5,7,'cc[atgc]', 'FE-BE3', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$FE_BE3 = answer
  first_pam = find_PAM_G(data_mut_G,5,7,'[atgc]gg', 'FE-BE3', 3,'c',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,5,7,'cc[atgc]', 'FE-BE3', 3,'g',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  
  #YEE-BE3
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first = find_PAM_G(data_mut_G,5,6,'[atgc]gg', 'YEE-BE3', 3,'c',23,answer)[[1]]
  answer = find_PAM_G(data_mut_G,5,6,'[atgc]gg', 'YEE-BE3', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,5,6,'cc[atgc]', 'YEE-BE3', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,5,6,'cc[atgc]', 'YEE-BE3', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$YEE_BE3 = answer
  first_pam = find_PAM_G(data_mut_G,5,6,'[atgc]gg', 'YEE-BE3', 3,'c',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,5,6,'cc[atgc]', 'YEE-BE3', 3,'g',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #xCas9烹E3(NGN)
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first = find_PAM_G(data_mut_G,5,6,'[atgc]g[atgc]', 'xCas9烹E3(NGN)', 3,'c',23,answer)[[1]]
  answer = find_PAM_G(data_mut_G,5,6,'[atgc]g[atgc]', 'xCas9烹E3(NGN)', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,5,6,'[atgc]c[atgc]', 'xCas9烹E3(NGN)', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,5,6,'[atgc]c[atgc]', 'xCas9烹E3(NGN)', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$xCas9_BE3_NGN = answer
  first_pam = find_PAM_G(data_mut_G,5,6,'[atgc]g[atgc]', 'xCas9烹E3(NGN)', 3,'c',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,5,6,'[atgc]c[atgc]', 'xCas9烹E3(NGN)', 3,'g',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #xCas9烹E3(GAW)
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first = find_PAM_G(data_mut_G,5,6,'ga[at]', 'xCas9烹E3(GAW)', 3,'c',23,answer)[[1]]
  answer = find_PAM_G(data_mut_G,5,6,'ga[at]', 'xCas9烹E3(GAW)', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,5,6,'[at]tc', 'xCas9烹E3(GAW)', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,5,6,'[at]tc', 'xCas9烹E3(GAW)', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$xCas9_BE3_GAW = answer
  first_pam = find_PAM_G(data_mut_G,5,6,'ga[at]', 'xCas9烹E3(GAW)', 3,'c',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,5,6,'[at]tc', 'xCas9烹E3(GAW)', 3,'g',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #BE-PLUS
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first = find_PAM_G(data_mut_G,5,16,'[atgc]gg', 'BE-PLUS', 3,'c',23,answer)[[1]]
  answer = find_PAM_G(data_mut_G,5,16,'[atgc]gg', 'BE-PLUS', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,5,16,'cc[atgc]', 'BE-PLUS', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,5,16,'cc[atgc]', 'BE-PLUS', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$BE_PLUS = answer
  first_pam = find_PAM_G(data_mut_G,5,16,'[atgc]gg', 'BE-PLUS', 3,'c',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,5,16,'cc[atgc]', 'BE-PLUS', 3,'g',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #APOBEC3A-Cas9
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first  = find_PAM_G(data_mut_G,5,9,'[atgc]gg', 'APOBEC3A-Cas9', 3,'c',23,answer)[[1]]
  answer = find_PAM_G(data_mut_G,5,9,'[atgc]gg', 'APOBEC3A-Cas9', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,5,9,'cc[atgc]', 'APOBEC3A-Cas9', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,5,9,'cc[atgc]', 'APOBEC3A-Cas9', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$APOBEC3A_Cas9 = answer
  first_pam = find_PAM_G(data_mut_G,5,9,'[atgc]gg', 'APOBEC3A-Cas9', 3,'c',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,5,9,'cc[atgc]', 'APOBEC3A-Cas9', 3,'g',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #EA3A-BE3(xCas9)-NGG
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first  = find_PAM_G(data_mut_G,4,8,'[atgc]gg', 'EA3A-BE3(xCas9)-NGG', 3,'c',23,answer)[[1]]
  answer = find_PAM_G(data_mut_G,4,8,'[atgc]gg', 'EA3A-BE3(xCas9)-NGG', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,4,8,'cc[atgc]', 'EA3A-BE3(xCas9)-NGG', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,4,8,'cc[atgc]', 'EA3A-BE3(xCas9)-NGG', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$EA3A_BE3_xCas9_NGG = answer
  first_pam = find_PAM_G(data_mut_G,4,8,'[atgc]gg', 'EA3A-BE3(xCas9)-NGG', 3,'c',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,4,8,'cc[atgc]', 'EA3A-BE3(xCas9)-NGG', 3,'g',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #EA3A-BE3(xCas9)-NGT
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first  = find_PAM_G(data_mut_G,4,8,'[atgc]gt', 'EA3A-BE3(xCas9)-NGT', 3,'c',23,answer)[[1]]
  answer = find_PAM_G(data_mut_G,4,8,'[atgc]gt', 'EA3A-BE3(xCas9)-NGT', 3,'c',23,answer)[[2]]
  second = find_PAM_C(data_mut_C,4,8,'ac[atgc]', 'EA3A-BE3(xCas9)-NGT', 3,'g',23,answer)[[1]]
  answer = find_PAM_C(data_mut_C,4,8,'ac[atgc]', 'EA3A-BE3(xCas9)-NGT', 3,'g',23,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$EA3A_BE3_xCas9_NGT = answer
  first_pam = find_PAM_G(data_mut_G,4,8,'[atgc]gt', 'EA3A-BE3(xCas9)-NGT', 3,'c',23,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,4,8,'ac[atgc]', 'EA3A-BE3(xCas9)-NGT', 3,'g',23,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #EA3A-BE3(VRQR)
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first  = find_PAM_G(data_mut_G,4,11,'[atgc]ga[atgc]', 'EA3A-BE3(VRQR)', 4,'c',24,answer)[[1]]
  answer = find_PAM_G(data_mut_G,4,11,'[atgc]ga[atgc]', 'EA3A-BE3(VRQR)', 4,'c',24,answer)[[2]]
  second = find_PAM_C(data_mut_C,4,11,'[atgc]tc[atgc]', 'EA3A-BE3(VRQR)', 4,'g',24,answer)[[1]]
  answer = find_PAM_C(data_mut_C,4,11,'[atgc]tc[atgc]', 'EA3A-BE3(VRQR)', 4,'g',24,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$EA3A_BE3_VRQR = answer
  first_pam = find_PAM_G(data_mut_G,4,11,'[atgc]ga[atgc]', 'EA3A-BE3(VRQR)', 4,'c',24,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,4,11,'[atgc]tc[atgc]', 'EA3A-BE3(VRQR)', 4,'g',24,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #YE1-VQR-Cas9
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first  = find_PAM_G(data_mut_G,5,6,'[atgc]ga[atgc]', 'YE1-VQR-Cas9', 4,'c',24,answer)[[1]]
  answer = find_PAM_G(data_mut_G,5,6,'[atgc]ga[atgc]', 'YE1-VQR-Cas9', 4,'c',24,answer)[[2]]
  second = find_PAM_C(data_mut_C,5,6,'[atgc]tc[atgc]', 'YE1-VQR-Cas9', 4,'g',24,answer)[[1]]
  answer = find_PAM_C(data_mut_C,5,6,'[atgc]tc[atgc]', 'YE1-VQR-Cas9', 4,'g',24,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$YE1_VQR_Cas9 = answer
  first_pam = find_PAM_G(data_mut_G,5,6,'[atgc]ga[atgc]', 'YE1-VQR-Cas9', 4,'c',24,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,5,6,'[atgc]tc[atgc]', 'YE1-VQR-Cas9', 4,'g',24,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #VQR-Cas9
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first  = find_PAM_G(data_mut_G,4,11,'[atgc]ga[atgc]', 'VQR-Cas9', 4,'c',24,answer)[[1]]
  answer = find_PAM_G(data_mut_G,4,11,'[atgc]ga[atgc]', 'VQR-Cas9', 4,'c',24,answer)[[2]]
  second = find_PAM_C(data_mut_C,4,11,'[atgc]tc[atgc]', 'VQR-Cas9', 4,'g',24,answer)[[1]]
  answer = find_PAM_C(data_mut_C,4,11,'[atgc]tc[atgc]', 'VQR-Cas9', 4,'g',24,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$VQR_Cas9 = answer
  first_pam = find_PAM_G(data_mut_G,4,11,'[atgc]ga[atgc]', 'VQR-Cas9', 4,'c',24,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,4,11,'[atgc]tc[atgc]', 'VQR-Cas9', 4,'g',24,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  #EQR-Cas9
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first  = find_PAM_G(data_mut_G,4,11,'[atgc]gag', 'EQR-Cas9', 4,'c',24,answer)[[1]]
  answer = find_PAM_G(data_mut_G,4,11,'[atgc]gag', 'EQR-Cas9', 4,'c',24,answer)[[2]]
  second = find_PAM_C(data_mut_C,4,11,'ctc[atgc]', 'EQR-Cas9', 4,'g',24,answer)[[1]]
  answer = find_PAM_C(data_mut_C,4,11,'ctc[atgc]', 'EQR-Cas9', 4,'g',24,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$EQR_Cas9 = answer
  first_pam = find_PAM_G(data_mut_G,4,11,'[atgc]gag', 'EQR-Cas9', 4,'c',24,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,4,11,'ctc[atgc]', 'EQR-Cas9', 4,'g',24,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  #VRER-Cas9
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first  = find_PAM_G(data_mut_G,2,11,'[atgc]gag', 'VRER-Cas9', 4,'c',24,answer)[[1]]
  answer = find_PAM_G(data_mut_G,2,11,'[atgc]gag', 'VRER-Cas9', 4,'c',24,answer)[[2]]
  second = find_PAM_C(data_mut_C,2,11,'ctc[atgc]', 'VRER-Cas9', 4,'g',24,answer)[[1]]
  answer = find_PAM_C(data_mut_C,2,11,'ctc[atgc]', 'VRER-Cas9', 4,'g',24,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$VRER_Cas9 = answer
  first_pam = find_PAM_G(data_mut_G,2,11,'[atgc]gag', 'VRER-Cas9', 4,'c',24,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,2,11,'ctc[atgc]', 'VRER-Cas9', 4,'g',24,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  
  #SaKKH-BE3
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first  = find_PAM_G(data_mut_G,4,12,'[atgc][atgc][atgc][ag][ag]t', 'SaKKH-BE3', 6,'c',26,answer)[[1]]
  answer = find_PAM_G(data_mut_G,4,12,'[atgc][atgc][atgc][ag][ag]t', 'SaKKH-BE3', 6,'c',26,answer)[[2]]
  second = find_PAM_C(data_mut_C,4,12,'a[tc][tc][atgc][atgc][atgc]', 'SaKKH-BE3', 6,'g',26,answer)[[1]]
  answer = find_PAM_C(data_mut_C,4,12,'a[tc][tc][atgc][atgc][atgc]', 'SaKKH-BE3', 6,'g',26,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$SaKKH_BE3 = answer
  first_pam = find_PAM_G(data_mut_G,4,12,'[atgc][atgc][atgc][ag][ag]t', 'SaKKH-BE3', 6,'c',26,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,4,12,'a[tc][tc][atgc][atgc][atgc]', 'SaKKH-BE3', 6,'g',26,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  #SaBE3
  answer = c(nrow(mutations))
  answer = preprocessing(mutations, 'T>C','A>G',answer,'t','a')[[3]]
  first  = find_PAM_G(data_mut_G,6,12,'[atgc][atgc]g[ag][ag]t', 'SaBE3', 6,'c',26,answer)[[1]]
  answer = find_PAM_G(data_mut_G,6,12,'[atgc][atgc]g[ag][ag]t', 'SaBE3', 6,'c',26,answer)[[2]]
  second = find_PAM_C(data_mut_C,6,12,'a[tc][tc]c[atgc][atgc]', 'SaBE3', 6,'g',26,answer)[[1]]
  answer = find_PAM_C(data_mut_C,6,12,'a[tc][tc]c[atgc][atgc]', 'SaBE3', 6,'g',26,answer)[[2]]
  full_table = unique(rbind(full_table,first,second))
  answer = c(answer,length(unique(unique(rbind(first,second))$chromEnd)))
  results$SaBE3 = answer
  first_pam = find_PAM_G(data_mut_G,6,12,'[atgc][atgc]g[ag][ag]t', 'SaBE3', 6,'c',26,answer)[[3]]
  second_pam = find_PAM_C(data_mut_C,6,12,'a[tc][tc]c[atgc][atgc]', 'SaBE3', 6,'g',26,answer)[[3]]
  pam_data = rbind(pam_data,first_pam,second_pam)
  
  
  library("xlsx")
  ress_freq = paste( '1statistic',f,".xlsx", sep="")
  ress_zero_full = paste(  f,"1.xlsx", sep="")
  write.xlsx(results, ress_freq, sheetName="1",  col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)
  write.xlsx(full_table, ress_zero_full, sheetName="1",  col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)
  re_full = paste(  f,"11.xlsx", sep="")
  write.xlsx(pam_data, re_full, sheetName="1",  col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)
  
  detach('package:xlsx')
  
}