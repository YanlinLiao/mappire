
#########Equation generator#########
#' create_bivalent_pairing
#' generate the format of bivalent pairing according to the ploidy level
#' @param ploidy crop ploidy level
#'
#' @return a vector of all possible pairing in two formats
#' In hexaploid, format1 is:
#'  "12|34|56" "12|35|46" "12|36|45" "13|26|45" "13|25|46" "13|24|56" "14|23|56" "16|23|45" "15|23|46"
#'  format2 is:
#'  "H1_H2|H3_H4|H5_H6" "H1_H2|H3_H5|H4_H6" "H1_H2|H3_H6|H4_H5" "H1_H3|H2_H6|H4_H5" "H1_H3|H2_H5|H4_H6" "H1_H3|H2_H4|H5_H6"
#'  "H1_H4|H2_H3|H5_H6" "H1_H6|H2_H3|H4_H5" "H1_H5|H2_H3|H4_H6"
#'
#' @export
#'
#' @examples
create_bivalent_pairing <- function(ploidy){
  allele <- seq(1,ploidy)

  library(combinat)
  gamete_matrix <- t(do.call(rbind,permn(allele)))

  ###format1
  pair_matrix_f1 <- do.call(cbind,lapply(seq(1,nrow(gamete_matrix),2),function(i){
    tmp <- gamete_matrix[i:(i+1),]
    as.matrix(sapply(1:ncol(tmp),function(j){
      paste0(sort(tmp[,j]),collapse = '')
    }))
  }))
  pair_matrix_f1 <- unique(t(apply(pair_matrix_f1,1,sort)))
  pairing_format1 <- apply(pair_matrix_f1, 1, paste, collapse='|')
  ###format2
  pair_matrix_f2 <- do.call(cbind,lapply(seq(1,nrow(gamete_matrix),2),function(i){
    tmp <- gamete_matrix[i:(i+1),]
    as.matrix(sapply(1:ncol(tmp),function(j){
      paste0(
        paste0('H',sort(tmp[,j]))
        ,collapse = '_')
    }))
  }))
  pair_matrix_f2 <- unique(t(apply(pair_matrix_f2,1,sort)))
  pairing_format2 <- apply(pair_matrix_f2, 1, paste, collapse='|')



  return(list('format1' = pairing_format1,
              'format2' = pairing_format2))
}


#' generate_recombi_nonrecombi
#' Generate the possible allele combinations together with their possibilities
#' @param ploidy crop ploidy level
#' @param hom_tmp a matrix of allele positions with their chromosomes, as the format of:
#'                     H1  H2
#'                [1,] "1" "2"
#'                [2,] "1" "2"
#' @return a list include two sublist with the recombination and non recombination
#' happened probability, as the format of:
#' $`1/2*s`
#' [1] "11" "22"
#'
#' $`1/2*r`
#' [1] "12" "21"
#'
#' @export
#'
#' @examples
generate_recombi_nonrecombi <- function(hom_tmp = tmp_P[,c(part1_hom1,part1_hom2)],
                                        ploidy){
  nonrecombi <- sapply(1:ncol(hom_tmp),function(j){paste0(hom_tmp[,j],collapse = '')})
  hom_tmp_reco <- hom_tmp; hom_tmp_reco[2,] <- hom_tmp_reco[2,c(2,1)]
  recombi <- sapply(1:ncol(hom_tmp_reco),function(j){paste0(hom_tmp_reco[,j],collapse = '')})

  res <- list('1/2*s' = nonrecombi,
              '1/2*r' = recombi)
  return(res)
}


#' generate all possible equation under polysomic inheritance
#'
#' @param P_genotype the genotype of marker P at one parent. It is a vector of characters, e.g. c("A","A","B","B")
#' @param Q_genotype the genotype of marker Q also at this parent. It is a vector of characters
#' @param generalized_equation generalized equation obtained from *disomic_equation_generator* according to the ploidy level
#'
#' @return a list include when parents with the given genotype, under each possible phasing, what is the
#' expected probability and number of recombination
#' @export
#'
#' @examples
generate_all_equation_polysomic <- function(P_genotype = strsplit('AABBBB','')[[1]],
                                  Q_genotype = strsplit('ABCDEF','')[[1]],
                                  generalized_equation){
  #load sufficient input
  ploidy <- generalized_equation$input$ploidy

  combs <- gtools::permutations(ploidy,ploidy)
  combs <- as.data.frame(unique(do.call(rbind,lapply(1:nrow(combs),function(i){
    Q_genotype[combs[i,]]
  }))))
  combs1 <- as.data.frame(apply(unique_hap(hap_mrkP = P_genotype,
                             hap_mrkQ = combs),2,as.character))
  if(ncol(combs1) == 1){
    combs1 <- t(combs1)
    rownames(combs1) <- NULL
  }
  equation <- lapply(1:nrow(combs1),function(z){
    equation <- simplify_equation(equation = generalized_equation$format1,
                                  comb_P = P_genotype,
                                  comb_Q = combs1[z,],
                                  ploidy = ploidy)
    equation$gsub_probability <- unlist(lapply(equation$probability,function(e){
      gsub(' ','',Deriv::Simplify(gsub('s','(1-r)',e)))
    }))
    equation
  })
  names(equation) <- sapply(1:nrow(combs1),function(z){nme <-
    paste0(combs1[z,],collapse = '')
  })
  return(equation)
}

#' generate all equation under simulated disomic inheritance when two marker genotype is given
#'
#' @param P_genotype the genotype of marker P at one parent. It is a vector of characters, e.g. c("A","A","B","B")
#' @param Q_genotype the genotype of marker Q also at this parent. It is a vector of characters
#' @param generalized_equation generalized equation obtained from *disomic_equation_generator* according to the ploidy level
#'
#' @return a list include when parents with the given genotype, under each possible phasing, what is the
#' expected probability and number of recombination
#' @export
#'
#' @examples
generate_all_equation_disomic <- function(P_genotype = strsplit('AABBBB','')[[1]],
                                          Q_genotype = strsplit('ABCDEF','')[[1]],
                                          generalized_equation){
  #load sufficient input
  ploidy <- generalized_equation$input$ploidy

  combs <- gtools::permutations(ploidy,ploidy)
  combs <- as.data.frame(unique(do.call(rbind,lapply(1:nrow(combs),function(i){
    Q_genotype[combs[i,]]
  }))))
  combs1 <- apply(unique_hap(hap_mrkP = P_genotype,
                             hap_mrkQ = combs),2,as.character)

  equation <- lapply(1:nrow(combs1),function(z){
    inheritance_equation <- lapply(names(generalized_equation$format1),function(p){
      each_inheritance <- generalized_equation$format1[[p]]
      equation <- simplify_equation(equation = each_inheritance,
                                    comb_P = P_genotype,
                                    comb_Q = combs1[z,],
                                    ploidy = ploidy)
      equation$gsub_probability <- unlist(lapply(equation$probability,function(e){
        gsub(' ','',Deriv::Simplify(gsub('s','(1-r)',e)))
      }))
      equation
    })

    names_list <- gsub('_','',gsub('H','',names(generalized_equation$format1)))
    names_list1 <- unlist(lapply(names_list,function(nme){
      suppressWarnings(paste0(plyr::mapvalues(strsplit(nme,'')[[1]],from = seq(1,ploidy),to = c(combs1[z,])),collapse = ''))
    }))
    names(inheritance_equation) <- names_list1
    inheritance_equation
  })
  names(equation) <- sapply(1:nrow(combs1),function(z){nme <-
    paste0(combs1[z,],collapse = '')
  })
  return(equation)
}



#' simplify_equation
#' simplify the equation according to specific marker type
#' the general matrix can be folded according to specific marker type
#' @param equation equation generated according to ploidy and pairing
#' @param comb_P a matrix including the genotype of one marker's one parent, e.g. matrix('A','A','B','C')
#' @param comb_Q a matrix including the genotype of the other marker's one parent, e.g. matrix('A','A','B','C')
#' @param ploidy crop ploidy level
#'
#' @return a dataframe including all
#' @export
#'
#' @examples
simplify_equation <- function(equation = tetraploid_equation,
                              comb_P = a_parents['P1',],
                              comb_Q = combs_b_P1[z,],
                              ploidy = 4){
  tmp_P <- as.matrix(rbind(comb_P,comb_Q))
  rownames(tmp_P) <- c('P','Q');colnames(tmp_P) <- paste0('H',seq(1,ploidy))

  LocP_gamete <- possible_gamete_parent(tmp = tmp_P,
                                        loc = 'P',ploidy = ploidy)
  LocQ_gamete <- possible_gamete_parent(tmp = tmp_P,
                                        loc = 'Q',ploidy = ploidy)

  recombination_matrix <- equation$recombination
  probability_matrix <- equation$probability
  rownames(recombination_matrix) <- rownames(probability_matrix) <- LocP_gamete
  colnames(recombination_matrix) <- colnames(probability_matrix) <- LocQ_gamete

  #simplify the recombination matrix
  class(recombination_matrix) <- 'table'; recombination_matrix <- as.data.frame(recombination_matrix)
  colnames(recombination_matrix) <- c('P','Q','type')
  recombination_matrix$P_Q <- paste0(recombination_matrix$P, '_',recombination_matrix$Q)
  #simplify the probability matrix
  class(probability_matrix) <- 'table'; probability_matrix <- as.data.frame(probability_matrix)
  colnames(probability_matrix) <- c('P','Q','probability')
  probability_matrix$P_Q <- paste0(probability_matrix$P, '_',probability_matrix$Q)

  combined_tmp <- do.call(rbind,lapply(unique(recombination_matrix$P_Q), function(x){
    recombination <- recombination_matrix[recombination_matrix$P_Q %in% x,]$type
    recombination <- ifelse(length(unique(recombination)) > 1,NA,recombination[1])

    probability <- as.character(probability_matrix[probability_matrix$P_Q %in% x,]$probability)
    probability <- gsub(' ','', Deriv::Simplify(paste0(probability,collapse = '+')))

    P <- strsplit(x,'_')[[1]][1];Q <- strsplit(x,'_')[[1]][2]
    data.frame(P,Q,recombination,probability)
  }))
  combined_tmp$probability <- as.character(combined_tmp$probability)

  return(combined_tmp)
}

#' find possible gamete of parent
#'
#' @param tmp a dataframe has two rows represent the genotype of each parent
#' @param loc a character, i.e. 'P' which is one of the rownames of *tmp*
#' @param ploidy crop ploidy level
#'
#' @return a vector including all possible gametes that can produced from the
#' specified parent
#' @export
#'
#' @examples
possible_gamete_parent <- function(tmp,
                                   loc = 'P',
                                   ploidy){
  gamete <- combn(tmp[loc,],ploidy/2)
  gamete <- sapply(1:ncol(gamete),function(j){
    paste0(sort(gamete[,j]),collapse = '')
  })
  return(gamete)
}


#' create_gamete_polysomic
#' generate the polysomic inheritance gametes
#'
#' @param parent one parent's multi-allelic marker, e.g. AAAC
#'
#' @return all possible gametes under polysomic inheritance
#' @export
#'
#' @examples
create_gamete_polysomic <- function(parent = P1,ploidy){ #parent:'AAAC'
  allele <- strsplit(parent,'')[[1]]
  gamete_matrix <- combn(allele,ploidy/2)
  gamete <- list(sapply(1:ncol(gamete_matrix),function(i){
    paste0(gamete_matrix[,i],collapse = '')
  }))
  names(gamete) <- parent
  frequency <- table(names(gamete))
  frequency <- frequency/sum(frequency)
  result_list <- gamete[names(frequency)]
  return(list('frequency' = frequency,
              'gamete' = result_list))
}


#' create_gamete_disomic
#' generate the disomic inheritance gametes under all situations
#'
#' @param parent one parent's multi-allelic marker, e.g. AAAC
#'
#' @return all possible gametes under disomic inheritance of all possibility situations
#' @export
#'
#' @examples
create_gamete_disomic <- function(parent = P1,ploidy){
  allele <- strsplit(parent,'')[[1]]
  pairing_nme <- create_bivalent_pairing(ploidy = ploidy)
  pairing_nme <- pairing_nme$format1

  #find the corresponding disomic pairing configuration
  pairing_new_nme <- unlist(lapply(pairing_nme,function(nme){
    paste0(plyr::mapvalues(strsplit(nme,'')[[1]],from = seq(1,ploidy),to = allele),collapse = '')
  }))
  names(pairing_new_nme) <- NULL

  result_list <- lapply(unique(pairing_new_nme),function(new_nme){
    gametes <- expand.grid(strsplit(strsplit(new_nme,'\\|')[[1]],''))
    gametes <- t(apply(gametes,1,sort))
    apply(gametes,1,paste,collapse = '')
  })
  names(result_list) <- unique(pairing_new_nme)

  frequency <- pairing_new_nme
  frequency <- table(frequency)
  frequency <- frequency/sum(frequency)

  return(list('frequency' = frequency,
              'gamete' = result_list))
}


######Infer offspring gametes######
#' found_probability
#' Find the possible gamete genotypes which contribute by each parents and their corresponding
#' probabilities.
#'
#' @param P1 parent1's genotype
#' @param P2 parent2's genotype
#' @param inheritance1 inheritance of P1, can be either 'Polysomic' or 'Disomic'
#' @param inheritance2 inheritance of P2, can be either 'Polysomic' or 'Disomic'
#' @param row genotype from one offspring, eg. "A" "A" "E" "F"
#' @param ind_nme the name of that offspring, eg. "F1_4"
#'
#' @return a dataframe contains the possible offspring genotype and their corresponding
#' probability, the colnames including:
#' 'individual P1_1 P1_2 P2_1 P2_2 probability   pairing pair_freq P1_gamete_freq P2_gamete_freq '
#'     ' F1_4    A    A    E    F           1 AAAB_EEEF         1            0.5            0.5 '
#' @export
#'
#' @examples
found_probability <- function(P1,
                              P2,
                              inheritance1,
                              inheritance2,
                              row = F1[i,],
                              ind_nme = F1_ind[i],
                              ploidy){
  if(inheritance1 == 'Polysomic'){
    P1_gamete <- create_gamete_polysomic(parent = P1,ploidy)
  }else{
    P1_gamete <- create_gamete_disomic(parent = P1,ploidy)
  }

  if(inheritance2 == 'Polysomic'){
    P2_gamete <- create_gamete_polysomic(parent = P2,ploidy)
  }else{
    P2_gamete <- create_gamete_disomic(parent = P2,ploidy)
  }
  result_temp <- list()
  for(parent1 in names(P1_gamete$gamete)){
    for(parent2 in names(P2_gamete$gamete)){
      P1_G <- P1_gamete$gamete[[parent1]]
      P2_G <- P2_gamete$gamete[[parent2]]

      P1_pass <- t(matrix(unlist(strsplit(P1_G,'')),ploidy/2))
      P2_pass <- t(matrix(unlist(strsplit(P2_G,'')),ploidy/2))

      result <- as.data.frame(do.call(rbind,lapply(1:nrow(P1_pass),function(i){
        do.call(rbind,lapply(1:nrow(P2_pass),function(j){
          possible <- c(as.character(P1_pass[i,]),as.character(P2_pass[j,]))
          if(all(sort(possible) == sort(row))) possible
        }))})))
      if(nrow(result) > 0){
        result <- plyr::count(result,vars = paste0('V',seq(1,ploidy)))
        result$individual <- ind_nme
        result$pairing <- paste0(parent1,'_',parent2)
        result$pair_freq <- P1_gamete$frequency[parent1] * P2_gamete$frequency[parent2]

        result$P1_freq <- sapply(1:nrow(result),function(r){
          table(P1_G)[paste0(sort(unlist(c(result[r,1:(ploidy/2)]))),collapse = '')]/length(P1_G)
        })
        result$P2_freq <- sapply(1:nrow(result),function(r){
          table(P2_G)[paste0(sort(unlist(c(result[r,(ploidy/2 + 1):(ploidy)]))),collapse = '')]/length(P2_G)
        })
        colnames(result) <- c(paste0('P1','_',seq(1,ploidy/2)),
                              paste0('P2','_',seq(1,ploidy/2)),
                              'probability','individual','pairing','pair_freq','P1_gamete_freq','P2_gamete_freq')
        result <- result[,c('individual',paste0('P1','_',seq(1,ploidy/2)),
                            paste0('P2','_',seq(1,ploidy/2)),'probability','pairing','pair_freq','P1_gamete_freq',
                            'P2_gamete_freq')]
        result$probability <- result$P1_gamete_freq * result$P2_gamete_freq
        result$probability <- result$probability/sum(result$probability)
        result <- apply(result,2,as.character)
        result_temp[[paste0(parent1,'_',parent2)]] <- result
      }
    }
  }
  result_tmp <- as.data.frame(do.call('rbind',result_temp),stringsAsFactors = FALSE)
  rownames(result_tmp) <- NULL
  return(result_tmp)
}



#' format_pairint_list
#'
#'
#' @param temp a dataframe contains the possible offspring genotype and their corresponding
#' probability.
#' @param F1_ind individual name of all F1 individuals
#'
#' @return a list for each paring, including the locus positions for marker P and Q,
#' and their probability of recombination, and the count from offspring genotypes.
#' @export
#'
#' @examples
format_pairint_list <- function(temp = poly_total,
                                F1_ind){
  final_list <- list()
  for(pair in unique(temp$pairing)){
    df <- temp[temp$pairing %in% pair,]
    not_appearing_ind <- setdiff(F1_ind,unique(df$individual))
    if(length(not_appearing_ind) > 0){
      df2 <- as.data.frame(matrix(ncol = ncol(df),nrow = length(not_appearing_ind)),stringsAsFactors = FALSE)
      colnames(df2) <- colnames(df)
      df2$individual <- not_appearing_ind
      df2$pairing <- pair
      df2$pair_freq <- unique(df$pair_freq)
      df <- rbind(df,df2)
    }
    final_list[[pair]] <- df
  }
  return(final_list)
}



######Markertype overview######
#' decide_MT
#' marker type identification
#'
#' @param parent parent genotype
#'
#' @return turn the parent genotype of multi-allelic marker to a number.
#' ABCD: 1_1_1_1
#' AABC: 1_1_2
#' AABB: 2_2
#' AAAB: 1_3
#' AAAA: 4
#' @export
#'
#' @examples
decide_MT <- function(parent = P1){
  allele <- strsplit(parent,'')[[1]]
  paste0(sort(table(allele)),collapse = "_")
}



#' find the unique haplotype configuration
#'
#' @param pair_comb a dataframe include two columns, first column include P1's haplotype (e.g. ABCD), second column include
#' P2's haplotype
#'
#' @return a same dataframe with the unique configuration (filter out the replicate configuration)
#' @export
#'
#' @examples
unique_hap <- function(hap_mrkP = a_parents['P1',],
                       hap_mrkQ = combs_b_P1){
  if(nrow(hap_mrkQ) > 1){
    hap_mrkQ <- as.data.frame(apply(hap_mrkQ,2,as.character),stringasfactors =FALSE)
  }
  hap_mrkQ$unique_hap <- sapply(1:nrow(hap_mrkQ),function(i){
    paste0(sort(paste0(hap_mrkP,hap_mrkQ[i,])),collapse = '_')
  })
  pair_comb_1 <- do.call(rbind,lapply(unique(hap_mrkQ$unique_hap),function(u){
    hap_mrkQ[hap_mrkQ$unique_hap %in% u,][1,] #always choose the first solution
  }))
  pair_comb_1$unique_hap <- NULL
  return(pair_comb_1)
}

#' identify the sharing between two parents
#'
#' @param mrk a string record one marker. e.g. 'AABC_AABE'
#'
#' @return a numeric about the percentage of sharing between two parents
#' @export
#'
#' @examples
sharing_identification <- function(mrk,ploidy){
  P1 <- strsplit(strsplit(mrk,'_')[[1]][1],'')[[1]]
  P2 <- strsplit(strsplit(mrk,'_')[[1]][2],'')[[1]]
  count_P1 <- table(P1)
  count_P2 <- table(P2)
  common <- intersect(P1,P2)

  sharing <- sum(unlist(sapply(common,function(c){
    min(count_P1[c],count_P2[c])
  })))/ploidy
  return(sharing)
}

######Polysomic x Disomic######
#' examine segregation inheritance
#'
#' @param possible_gamete expected possible gametes of one parent. It is generated either from
#' *create_gamete_disomic* or *create_gamete_polysomic*
#' @param parent 'P1' or 'P2'
#' @param probability_offspring probability offspring scores which are used for test.It is a dataframe contains:
#' + indivdual: individual name
#' + P1_1, P1_2: the gametes from P1 (not phased)
#' + P2_1, P2_2: the gametes from P2 (not phased)
#' + probability: the probability obtaining the gametes
#' + pairing: the pairing of two parents
#' + pair_freq: the frequency of this pairing from two parents
#' + P1_gamete_freq, P2_gamete_freq: Under this assumed pairing and inheritance, what is the expected
#' segregation ratio of obtaining this gamete.
#' @param seg_invalidrate the value used in binomial test
#' @param missing_nbr the number of missing individuals
#'
#' @return
#' @export
#'
#' @examples
check_segregation <- function(possible_gamete = P1_gamete_diso,
                              parent = 'P1',
                              probability_offspring = F1,
                              seg_invalidrate,
                              missing_nbr = missing_nbr){

  temp <- lapply(names(probability_offspring),function(i){
    # print(i)
    parent_pairing <-  strsplit(i,'_')[[1]][as.numeric(gsub('P','',parent))]
    seg_count <- possible_gamete$gamete[[parent_pairing]]
    seg_table <- table(seg_count)

    offspring_temp <- probability_offspring[[i]]
    valid_temp <- offspring_temp[complete.cases(offspring_temp),]
    valid_ind <- sum(as.numeric(valid_temp$probability))

    frequency <- seg_pro_NEW (F1 = offspring_temp,
                              parent = parent)

    if(valid_ind > 10){
      obs <- frequency[names(seg_table)]
      names(obs) <- names(seg_table)
      obs[is.na(obs)] <- 0

      invalidP <- pbinom(q= sum(frequency),   #nr of valid counts
                         size= length(unique(offspring_temp$individual)) - missing_nbr, #total individuals
                         prob=1 - seg_invalidrate)

      if (length(obs) == 1) {
        p_chi <- 1
      } else {
        if (sum(obs) == 0) {
          #all observations are in invalid classes
          p_chi <- 0.0
          invalidP <- 0.0
        } else {
          test <-  suppressWarnings(chisq.test(obs,p = seg_table/sum(seg_table)))
          p_chi <- test$p.value
        }
      }
      combined_P <- metap::sumlog(p =  c(invalidP,p_chi))

      list('chi_square'=round(p_chi,2),
           'binomial_p' = round(invalidP,2),
           'multiplied_p' = round(p_chi * invalidP,2),
           'combined_p' = round(combined_P$p,2),
           'observed_count' = obs,
           'expected_count' = seg_table,
           'invalid%' = round((length(unique(offspring_temp$individual)) - sum(obs))/length(unique(offspring_temp$individual)),2))
    }
  })
  names(temp) <- names(probability_offspring)
  return(temp)
}

#' convert the offspring information to required format
#'
#' @param F1 offspring data from poly. A dataframe with ploidy number's column
#' and individual numbers row. The genotype within each individual is not phased
#' @param parent 'P',the prefix of parent
#'
#' @return the counted frequency of offspring genotype
#' @export
#'
#' @examples
seg_pro_NEW <- function(F1 = offspring_temp,
                        parent = parent){
  parent_col <- paste0(parent,"_",seq(1,2))
  F1 <- F1[complete.cases(F1),] #remove missing individual
  F1[,parent_col] <- t(apply(F1[,parent_col],1,sort))
  F1 <- apply(F1,2,as.character)
  gametes <- unique(unlist(sapply(1:nrow(F1),function(i){
    paste0(sort(F1[i,parent_col]),collapse = '')
  })))
  frequency <- sapply(gametes,function(i){
    allele <- strsplit(i,'')[[1]]

    sum(ifelse(sapply(1:nrow(F1),function(j){
      sum(F1[j,parent_col] == sort(allele))
    }) == 2,1,0) * as.numeric(F1[,'probability']))

  })
  return(frequency)
} #for observarion


#######Linkage estimation#######
#' kronecker string
#' kronecker two matrix
#' @param matrix1 A vector or array.
#' @param matrix2 A vector or array.
#' @param symbol the symbol used to paste the string together, e.g. '*'
#'
#' @return An array A with dimensions dim(matrix1) * dim(matrix2).
#' @export
#'
#' @examples
kronecker_string <- function(matrix1 = g_matrix_P1,
                             matrix2 = g_matrix_P2,
                             symbol = '*'){
  class(matrix1) <- 'table';class(matrix2) <- 'table'
  matrix1 <- as.data.frame(matrix1);matrix2 <- as.data.frame(matrix2)

  matrix1$A_B <- paste0(matrix1$Var1,'_',matrix1$Var2); matrix2$A_B <- paste0(matrix2$Var1,'_',matrix2$Var2)
  rownames(matrix1) <- matrix1$A_B; rownames(matrix2) <- matrix2$A_B

  combined_matrix <- expand.grid(matrix1$A_B,matrix2$A_B)
  combined_matrix$Freq <- sapply(1:nrow(combined_matrix),function(i){
    paste0('(',as.character(matrix1[combined_matrix[i,]$Var1,'Freq']),')',
           symbol,
           '(',as.character(matrix2[combined_matrix[i,]$Var2,'Freq']),')')
  })

  parent1_name <- t(matrix(unlist(strsplit(as.character(combined_matrix$Var1),'_')),2))
  parent2_name <- t(matrix(unlist(strsplit(as.character(combined_matrix$Var2),'_')),2))
  combined_matrix$rowname <- paste0(parent1_name[,1],':',parent2_name[,1])
  combined_matrix$colname <- paste0(parent1_name[,2],':',parent2_name[,2])
  colnames(combined_matrix) <- c('P1','P2','probability','P','Q')
  combined_matrix$simplified_probability <- gsub(' ','',sapply(combined_matrix[,3],Deriv::Simplify))


  combined_table <- dataframe_to_contigency(tmp = combined_matrix,
                                            col1 = 'P',
                                            col2 = 'Q',
                                            count = 'simplified_probability')
  return(list('format1' = combined_matrix,
              'format2' = combined_table))
}

#' find the most likely phasing
#'
#' @param P_phasing phasing of marker p
#' @param Q_phasing phasing of marker Q
#' @param tmp_count a data
#' @param equation
#' @param eq_chosen
#' @param parent target parent. integer,'P1' or 'P2'
#'
#' @return
#' @export
#'
#' @examples
phasing_filter <- function(P_phasing =  P_parents,
                           Q_phasing = Q_parents,
                           tmp_count = tmp,
                           equation,
                           eq_chosen,
                           parent){

  #get the individual count of tmp
  if(eq_chosen == 1){
    P_phasing_parent <- P_phasing[parent,]
    Q_phasing_choices <- do.call(rbind,strsplit(names(equation),''))


    tmp_count[[parent]] <- paste0(tmp_count[[paste0(parent,'_P')]],'_',tmp_count[[paste0(parent,'_Q')]])
    tmp_count <- data.frame(cbind(unique(tmp_count[[parent]]),do.call(rbind,lapply(unique(tmp_count[[parent]]),function(m){
      sum(tmp_count[tmp_count[[parent]] == m,'count'])
    }))))
    colnames(tmp_count) <- c('phasing','count')
    tmp_count$count <- as.numeric(tmp_count$count)
    tmp_count <- tmp_count[order(tmp_count$count,decreasing = TRUE),]
  }else{
    P_phasing_parent <- Q_phasing[parent,]
    Q_phasing_choices <- do.call(rbind,strsplit(names(equation),''))


    tmp_count[[parent]] <- paste0(tmp_count[[paste0(parent,'_Q')]],'_',tmp_count[[paste0(parent,'_P')]])
    tmp_count <- data.frame(cbind(unique(tmp_count[[parent]]),do.call(rbind,lapply(unique(tmp_count[[parent]]),function(m){
      sum(tmp_count[tmp_count[[parent]] == m,'count'])
    }))))
    colnames(tmp_count) <- c('phasing','count')
    tmp_count$count <- as.numeric(tmp_count$count)
    tmp_count <- tmp_count[order(tmp_count$count,decreasing = TRUE),]
  }

  #look at the most abundant ones
  pairs_most <- strsplit(unlist(strsplit(tmp_count[1,1],'_')),'')
  first <- which(P_phasing_parent %in% pairs_most[[1]][1])
  second <- which(P_phasing_parent %in% pairs_most[[1]][2])
  pair_loc <- expand.grid(first,second)
  #filter for the repetitive
  pair_loc <- t(apply(pair_loc,1,sort))
  pair_loc <- pair_loc[!duplicated(pair_loc), ]
  pair_loc <- matrix(pair_loc,ncol = 2)
  #remove the ones with itself
  pair_loc <- pair_loc[which(unlist(
    lapply(1:nrow(pair_loc),function(t){length(unique(pair_loc[t,]))})
  ) > 1),,drop = FALSE]
  Q_phasing_choices <- Q_phasing_choices[which(unlist(lapply(1:nrow(Q_phasing_choices),function(j){
    check_all <- unlist(lapply(1:nrow(pair_loc),function(p){
      forward <- all(Q_phasing_choices[j,unlist(c(pair_loc[p,]))] == pairs_most[[2]])
      reverse <- all(Q_phasing_choices[j,unlist(c(pair_loc[p,]))] == pairs_most[[2]][c(2,1)])
      any(forward,reverse)
    }))
    any(check_all)
  }))),,drop = FALSE]

  chosen_phasing <- apply(Q_phasing_choices,1,paste0,collapse ='')
  return(chosen_phasing)
}

#' use maximumlikelihood method to estimate the rf
#'
#' @param equation maximumlikelihood function. Here is a string
#' @param plot if TRUE, how the function looks like will be plotted
#'
#' @return a list containing: 'rf'
#' @export
#'
#' @examples
maximum_likelihood_estimator <- function(equation = maximum_equation,
                                         plot){
  function_equation <- eval(parse(text = paste0("function(x)(",
                                                gsub('r','x',gsub(' ','',paste0(equation,collapse = '')))
                                                ,")")))

  if(plot){plot(function_equation)}

  maximumlikelihood1 <- optimize(function_equation,c(0, 0.1),maximum = TRUE, tol =1e-20)
  maximumlikelihood2 <- optimize(function_equation,c(0.1, 0.5),maximum = TRUE)
  if(maximumlikelihood1$objective == maximumlikelihood2$objective){
    rf <- maximumlikelihood2$maximum
  }else{
    rf <- c(maximumlikelihood1$maximum, maximumlikelihood2$maximum)[
      which.max(c(maximumlikelihood1$objective,maximumlikelihood2$objective))]
  }

  ###LOD score
  LOD <- log10(function_equation(rf)/function_equation(0.5))

  if(LOD < 0 | is.na(LOD)) LOD <- 0
  return(list('rf' = round(rf,3),
              'LOD' = round(LOD,3)))
}

######Evaluate results######
#' multiplot
#' Generating plots in one display window
#'
#' @param plotlist name of plots want to combine
#' @param cols number of plots want to display in columns
#' @param layout number of plots want to display in rows
#'
#' @return arranged plot displays
#' @export
#'
#' @examples
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' simulated_distance
#' Calculate the distance between real simulated marker positions and model estimated ones
#'
#' @param mapfile the marker position file for stimulated dataset
#' @param mrk1 the marker P in linkage estimation result
#' @param mrk2 the marker Q in linkage estimation result
#'
#' @return integer for the distance of locus 1 and locus 2
#' @export
#'
#' @examples
simulated_distance <- function(mapfile,
                               mrk1,
                               mrk2){
  loc1 <- mapfile[mapfile$marker %in% mrk1,'position']
  loc2 <- mapfile[mapfile$marker %in% mrk2,'position']
  return(abs(loc1 - loc2))
}

#' distance_to_rf
#' Convert chromosomal position to difference of recombination frequency
#'
#' @param dt integer for the distance of locus 1 and locus 2
#'
#' @return a number of converted recombination frequency
#' @export
#'
#' @examples
distance_to_rf <- function(dt){
  dist_Morgan <- dt/100
  r_true <- (1-exp(-2*dist_Morgan))/2
  return(r_true)
}


#' simulated_phasing
#' Extract real genotype phasing for each marker
#'
#' @param genfile the genotype phasing file for stimulated dataset
#' @param mrk1 the marker P name in linkage estimation result
#' @param mrk2 the marker Q name in linkage estimation result
#'
#' @return
#' @export
#'
#' @examples
simulated_phasing_gametes <- function(genfile,
                              mrk1,
                              mrk2){

  simulated_phase <- rbind(genfile[genfile$marker %in% mrk1,],
                           genfile[genfile$marker %in% mrk2,])
  return(sort(sapply(2:ncol(simulated_phase),function(j){
    paste0(simulated_phase[,j],collapse = '')
  })))
}

simulated_phasing <- function(genfile,
                              mrk1,
                              mrk2,
                              ploidy){
  simulated_phase <- rbind(genfile[genfile$marker %in% mrk1,],
                           genfile[genfile$marker %in% mrk2,])
  P1 <- simulated_phase[,paste0('P1_',seq(1,ploidy))]; P2 <- simulated_phase[,paste0('P2_',seq(1,ploidy))]

  P_phasing <- paste0(paste0(P1[1,order(unlist(c(P1[1,])))],collapse = ''),'_', paste0(P2[1,order(unlist(c(P2[1,])))],collapse = ''))
  Q_phasing <- paste0(paste0(P1[2,order(unlist(c(P1[1,])))],collapse = ''),'_', paste0(P2[2,order(unlist(c(P2[1,])))],collapse = ''))
  return(c(P_phasing,Q_phasing))
}



#' estimated_phasing
#' Extract estimated genotype phasing for each marker
#'
#' @param mrk1_phase string of marker P estimated phasing
#' @param mrk2_phase string of marker Q estimated phasing
#'
#' @return
#' @export
#'
#' @examples
estimated_phasing <- function(mrk1_phase,
                              mrk2_phase){
  mrk1_phase <- strsplit(as.character(gsub('\\|','',mrk1_phase)),'_')[[1]]
  mrk2_phase <- strsplit(as.character(gsub('\\|','',mrk2_phase)),'_')[[1]]

  phasing <- rbind(unlist(strsplit(mrk1_phase,'')),
                   unlist(strsplit(mrk2_phase,'')))
  return(sort(sapply(1:ncol(phasing),function(j){
    paste0(phasing[,j],collapse = '')
  })))
}

#' plot evaluated results
#'
#' @param result a dataframe include the linkage results and simulated r and phasing
#' @param filename saved filename
#'
#' @return
#' @export a figure named with specified filename. It includes four figures:
#' + r vs. LOD
#' + true r vs. estimated r (phasing correctness as groups)
#' + r vs. LOD (wrongly phased)
#' + r vs. LOD (correctly phased)
#'
#' @examples
plot_rf_comparison <- function(result,
                               filename){
  library(randomcoloR)
  library(ggplot2)
  result$simulated_phasing_correct <- as.character(result$simulated_phasing_correct)
  p1 <- ggplot(result, aes(x=r, y=LOD)) + geom_point()+
    xlab("r") +
    ylab("LOD") +
    labs(subtitle = paste0(nrow(result),' pairs'))+
    ggtitle("estimated: r vs LOD") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(plot.title = element_text(size= 25,face="bold"),plot.subtitle = element_text(size = 18))+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=20,face="bold"))+
    theme(legend.position = "none")

  p2 <-  ggplot(result, aes(x= simulated_rf, y= r, color= simulated_phasing_correct, group= simulated_phasing_correct)) +
    geom_point()+
    geom_abline(slope=1, intercept=0)+
    geom_point(size=1, aes(color=simulated_phasing_correct), fill="white") +
    xlab("simulated r") +
    ylab("estimated r") +
    xlim(0,0.5)+ylim(0,0.5)+
    ggtitle("r: estimated vs. simulated") +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(plot.title = element_text(size= 25,face="bold"))+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=20,face="bold"))+
    theme(legend.position = "bottom",legend.title= element_text(size = 15),legend.text = element_text(size = 15),
          legend.background = element_blank(),legend.box.background = element_rect(colour = "black")) +
    labs(color='correctly phased haplotypes')

  res_wrongly_phased <- result[result$simulated_phasing_correct != 8,]
  p3 <- ggplot(res_wrongly_phased, aes(x=r, y=LOD)) + geom_point()+
    xlab("r") +
    ylab("LOD") +
    xlim(0,0.5)+
    labs(subtitle = paste0(nrow(res_wrongly_phased),' pairs'))+
    ggtitle("Wrongly phased markers") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(plot.title = element_text(size= 25,face="bold"),plot.subtitle = element_text(size = 18))+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=20,face="bold"))+
    theme(legend.position = "none")


  res_correct_phased <- result[result$simulated_phasing_correct == 8,]
  p4 <- ggplot(res_correct_phased, aes(x=r, y=LOD)) + geom_point()+
    xlab("r") +
    ylab("LOD") +
    labs(subtitle = paste0(nrow(res_correct_phased),' pairs'))+
    ggtitle("Correctly phased markers") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(plot.title = element_text(size= 25,face="bold"),plot.subtitle = element_text(size = 18))+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=20,face="bold"))+
    theme(legend.position = "none")

  png(filename = paste0(filename,'.png'),width = 800,height = 600)
  multiplot(p1, p2,p3,p4, layout = matrix(seq(1,4),nrow=2, byrow=TRUE))
  dev.off()
}




######Matrix conversion######

#' dataframe_to_contigency
#' convert dataframe to contingency table
#'
#' @param tmp a dataframe
#' @param col1 the colname which will be converted in contingency table as row
#' @param col2 the colname which will be converted in contingency table as col
#' @param count the colname which will be used to fill in the contingency table
#'
#' @return contigency table after required conversion
#' @export
#'
#' @examples
dataframe_to_contigency <- function(tmp,
                                    col1 = 'A',
                                    col2 = 'B',
                                    count = 'count'){
  col1_uniq <- unique(as.character(tmp[,col1]))
  col2_uniq <- unique(as.character(tmp[,col2]))
  tmp_new <- matrix(nrow = length(col1_uniq),ncol = length(col2_uniq))
  rownames(tmp_new) <- col1_uniq
  colnames(tmp_new) <- col2_uniq
  for(i in 1:nrow(tmp)){
    tmp_new[tmp[i,col1],tmp[i,col2]] <- tmp[i,count]
  }
  return(tmp_new)
}


######Simulation######
#' Simulate multi_allelic markers by using PedigreeSim
#' the marker can be specified with different types. At each locus, it will include all marker types
#' @param chrnum numeric. Number of chromosomes
#' @param gap numeric. the gap length
#' @param length numeric. Total length of chromosome
#' @param ploidy numeric. crop ploidy level
#' @param filename string. filename which will be used to save files
#' @param popnum numeric. The number of F1 individuals
#' @param marker_type a vector of strings in the format of c("ABCD_EFGH","AABC_EEFG","AAAB_EEEF")
#' @param pairing numeric. 0 - polysomic; 1 - disomic
#' @param quadrivalents numeric. 0 - bivalent pairing; 1 - multivalent pairing
#'
#' @return
#' @export different types of files with prefix of your specified filename
#'
#' @examples
Multi_allelic_simulator_rep <- function(chrnum = 1,
                                        gap = 1,
                                        length = 99,
                                        ploidy=4,
                                        filename = "test",
                                        popnum = 200,
                                        marker_type = pair_1,
                                        pairing = 0,
                                        quadrivalents = 0){
  chromfile <- data.frame(chromosome = NULL, length = NULL, centromere = NULL, prefPairing = NULL, quadrivalents = NULL)
  TotalSNPset <- data.frame()
  TotalLinkageMap <- data.frame()
  Total_ped <- data.frame()
  for(i in 1:chrnum) {
    chrname <- as.character(as.roman(i))
    # markerpositions <- seq(0,(length(marker_type) - 1)*gap,gap)
    markerpositions <- seq(0,length,gap)
    posnames <- sprintf("%03d", markerpositions)
    ContigNme <- paste0("Cm", stringi::stri_rand_strings((length(marker_type) - 1)*gap/gap+1,5,'[0-9]'))
    contig_tmp <- expand.grid(ContigNme,posnames)
    markernames <- paste0(contig_tmp$Var1, "_", contig_tmp$Var2)


    SNPset_mat <- matrix(integer(), ncol = 2*ploidy)
    SNPset <- do.call(rbind,lapply(marker_type,function(MT){
      a <- unlist(strsplit(strsplit(MT,'_')[[1]],''))
      a_P1 <- a[1:ploidy];a_P2 <- a[(ploidy + 1):(ploidy*2)]
      c(a_P1[sample(length(a_P1))], a_P2[sample(length(a_P2))])#randomize the haplotype homologue
    }))
    SNPset_all <- replicate(length(posnames),SNPset)
    SNPset <- do.call(rbind,lapply(1:dim(SNPset_all)[3],function(i){SNPset_all[,,i]}))
    rownames(SNPset) <- markernames
    SNPset_mat <- rbind(SNPset_mat, SNPset)


    positions <- c(t(replicate(length(marker_type),markerpositions)))
    linkagemap <- data.frame(marker = rownames(SNPset_mat), chromosome = chrname, position = positions,
                             stringsAsFactors = FALSE)

    posorder <- order(linkagemap$position)

    sortedSNPset <- SNPset_mat[posorder,]
    sortedlinkagemap <- linkagemap[posorder,]

    sortedSNPset <- cbind(rownames(sortedSNPset), as.data.frame(sortedSNPset))

    colnames(sortedSNPset) <- c("marker", paste("P1", 1:ploidy, sep="_"), paste("P2", 1:ploidy, sep="_"))



    p <- format(round(pairing, 2), nsmall = 2)
    q <- format(round(quadrivalents, 2), nsmall = 2)

    ## Updates the file
    tempchromfile <- data.frame(chromosome = chrname, length = length, centromere = length/2, prefPairing = p, quadrivalents = q )

    chromfile <- rbind(chromfile, tempchromfile)
    TotalSNPset <- rbind(TotalSNPset,sortedSNPset)
    TotalLinkageMap <- rbind(TotalLinkageMap,sortedlinkagemap)
  }

  Total_ped <- rbind(data.frame(rbind(c("P1",NA,NA),c("P2",NA,NA))),
                     data.frame(t(sapply(1:popnum,function(i){c(paste0("F1_",i),"P1","P2")}))))
  colnames(Total_ped) <- c("Name","Parent1", "Parent2")

  write.table(Total_ped,paste(filename,c(".ped"),sep=""), quote = FALSE, sep = " ", row.names = FALSE)
  write.table(TotalSNPset, paste(filename,c(".gen"),sep=""), quote = FALSE, sep = " ", row.names = FALSE)
  write.table(TotalLinkageMap, paste(filename,c(".map"),sep=""), quote = FALSE, sep = " ", row.names = FALSE)
  write.table(chromfile, paste(filename,c(".chrom"),sep=""), quote = FALSE, sep = " ", row.names = FALSE)


  ## Creation of the .par file with usual defaults
  parfilename <- paste(filename,c(".par"),sep="")
  fileConn<-file(parfilename)
  writeLines(c(paste0("PLOIDY = ", ploidy), "MAPFUNCTION = HALDANE", "MISSING = NA", "HAPLOSTRUCT = myhaplostruct", paste("CHROMFILE = ",filename,".chrom",sep=""),
               "POPTYPE = F1", paste("POPSIZE = ", popnum, sep=""), paste("MAPFILE = ", filename, ".map",sep=""),  paste("PEDFILE = ", filename, ".ped",sep=""),
               paste("FOUNDERFILE = ", filename, ".gen", sep=""), paste("OUTPUT = ", filename, "_out", sep=""), paste0("NATURALPAIRING = ",pairing)
  ), fileConn)
  close(fileConn)

  message("Thank you. Please check your R working directory.")
}


#' Simulate polysomic x disomic dataset with multi-allelic markers
#'
#' @param filename string. filename which will be used to save files
#' @param marker_type a vector of strings in the format of c("ABCD_EFGH","AABC_EEFG","AAAB_EEEF")
#' @param pref_P1 numeric. 0 - polysomic; 1 - disomic (for parent 1)
#' @param quad_P1 numeric. 0 - bivalent pairing; 1 - multivalent pairing (for parent 1)
#' @param pref_P2 numeric. 0 - polysomic; 1 - disomic (for parent 2)
#' @param quad_P2  numeric. 0 - bivalent pairing; 1 - multivalent pairing (for parent 2)
#' @param ploidy numeric. crop ploidy level
#' @param popnum numeric. The number of F1 individuals
#' @param lg numeric. Number of chromosomes
#' @param gap numeric. the gap length
#' @param length numeric. Total length of chromosome
#'
#' @return different types of files with prefix of your specified filename
#' @export
#'
#' @examples
MultiAllelic_Inheritance_simulator <- function(filename = 'Pedtry',
                                               marker_type = c('ABCD_EFGH','AABB_EEFF'),
                                               pref_P1 = 0,
                                               quad_P1 = 1,
                                               pref_P2 = 0.5,
                                               quad_P2 = 0,
                                               ploidy = 4,
                                               popnum = 200,
                                               lg = 1,
                                               gap = 1,
                                               length = 99){
  filename_P1 <- paste0(filename,"_P1");filename_P2 <- paste0(filename,"_P2")
  Multi_allelic_simulator_rep(chrnum = lg,
                              gap = gap,
                              length = length,
                              ploidy=ploidy,
                              filename = filename_P1,
                              popnum = popnum,
                              marker_type = marker_type,
                              pairing = pref_P1,
                              quadrivalents = quad_P1)
  ###Genfile
  P1_gen <- read.table(paste0(filename_P1,'.gen'),header = TRUE,stringsAsFactors = FALSE,colClasses = c("character"))
  P2_gen <- P1_gen
  P2_gen[,2:(ploidy+1)] <- matrix('Z',nrow = nrow(P1_gen),ncol = ploidy)
  P1_gen[,(ploidy+2):(2*ploidy+1)] <- matrix('Z',nrow = nrow(P1_gen),ncol = ploidy)
  write.table(P1_gen,paste(filename_P1,c(".gen"),sep=""), quote = FALSE, sep = " ", row.names = FALSE)
  write.table(P2_gen,paste(filename_P2,c(".gen"),sep=""), quote = FALSE, sep = " ", row.names = FALSE)
  ###Chromfile
  P2_chrom <- read.table(paste0(filename_P1,'.chrom'),header = TRUE,stringsAsFactors = FALSE,colClasses = c("character"))
  P2_chrom[,4] <- as.character(pref_P2);P2_chrom[,5] <- as.character(quad_P2)
  write.table(P2_chrom,paste(filename_P2,c(".chrom"),sep=""), quote = FALSE, sep = " ", row.names = FALSE)
  ###par file
  fileConn<- file(paste0(filename_P2,'.par'))
  writeLines(c(paste("PLOIDY =",ploidy),
               "MAPFUNCTION = HALDANE",
               'HAPLOSTRUCT = myhaplostruct',
               "MISSING = NA",
               paste0("CHROMFILE = ",filename_P2,".chrom"),
               "POPTYPE = F1",
               paste0("POPSIZE = ", popnum),
               paste0("MAPFILE = ", filename_P1, ".map"),
               paste0("PEDFILE = ", filename_P1, ".ped"),
               paste0("FOUNDERFILE = ", filename_P2,".gen"),
               paste0("OUTPUT = ", filename_P2,"_out"),
               paste0("NATURALPAIRING = ",0)), fileConn)
  close(fileConn)

  run.PedigreeSim(paste0(filename_P1,'.par'))
  run.PedigreeSim(paste0(filename_P2,'.par'))


  ######Combine allele dose file
  P1dose <-  readDatfile(paste0(filename_P1,'_out_allAlleledose.dat'))
  P2dose <-  readDatfile(paste0(filename_P2,'_out_allAlleledose.dat'))

  F1dose <- rbind(P1dose[!P1dose$allele %in% 'Z',],
                  P2dose[!P2dose$allele %in% 'Z',])
  F1dose <- F1dose[order(F1dose$marker),]

  write.table(F1dose,
              paste0(filename,"_out_alleledose.dat"), sep = "\t", row.names=FALSE)

  ######combine founder allele
  P1founder <- as.matrix(read.table(paste0(filename, "_P1_out_founderalleles.dat"),
                                    header = TRUE, stringsAsFactors = FALSE,
                                    sep = "\t",row.names=1))
  P2founder <- as.matrix(read.table(paste0(filename, "_P2_out_founderalleles.dat"),
                                    header = TRUE, stringsAsFactors = FALSE,
                                    sep = "\t",row.names=1))
  P1founder.X <- P1founder[,-seq(2*ploidy)]
  P2founder.X <- P2founder[,-seq(2*ploidy)]
  ## Generate new colnames for P2founder.1
  popNames <- colnames(F1dose)[5:ncol(F1dose)]
  P1founder.1 <- P1founder.X[,sort(c(matrix(1:ncol(P1founder.X),ncol=ploidy,byrow = TRUE)[,1:(ploidy/2)]))]
  P2founder.1 <- P2founder.X[,sort(c(matrix(1:ncol(P2founder.X),ncol=ploidy,byrow = TRUE)[,(ploidy/2 + 1):(ploidy)]))]
  ## P1founder.1 numbering is ok - starts from 0 as expected. P2founder.1 numbering is not right and must be adjusted.
  d <- ploidy - ploidy
  P2founder.1 <- P2founder.1 + d
  colnames(P2founder.1) <- sort(c(sapply((ploidy/2 + 1):((ploidy + ploidy)/2),
                                         function(n) paste(popNames,n,sep="_"))))
  founder.out <- cbind(P1founder.1,P2founder.1)
  founder.out <- founder.out[,order(colnames(founder.out))]
  ###Add parental columns back###
  P2parcols <- do.call(rbind,lapply(1:nrow(P2founder), function(cl) ploidy:(ploidy + ploidy - 1)))
  colnames(P2parcols) <- paste0("P2_",1:ploidy)
  founder.out <- cbind("marker" = rownames(P1founder),P1founder[,1:ploidy],P2parcols,founder.out)
  write.table(founder.out, paste0(filename,"_out_founderalleles.dat"), sep = "\t", row.names=FALSE)

  ######combine gen file
  F1gen.out <- cbind(P1_gen[,1:(ploidy+1)],P2_gen[,(ploidy+2):(2*ploidy+1)])
  write.table(F1gen.out, paste0(filename,".gen"), quote = FALSE, sep = " ", row.names = FALSE)

  ######combine hsa and hsb file
  P1hsa <- read.table(paste0(filename, "_P1_out.hsa"),header=FALSE,stringsAsFactors = FALSE)
  P1hsb <- read.table(paste0(filename, "_P1_out.hsb"),header=FALSE,stringsAsFactors = FALSE)
  P2hsa <- read.table(paste0(filename, "_P2_out.hsa"),header=FALSE,stringsAsFactors = FALSE)
  P2hsb <- read.table(paste0(filename, "_P2_out.hsb"),header=FALSE,stringsAsFactors = FALSE)
  #if they are not in the same length, make them the same
  asjust_length <- function(tmp1, tmp2){
    combine_list <- list(tmp1,tmp2)
    colLen <- unlist(lapply(combine_list,ncol))
    if(diff(colLen) != 0){
      chosen <- which.min(colLen);NTchosen <- which.max(colLen)
      NA_matrix <- t(matrix(rep(NA,nrow(combine_list[[chosen]]) * abs(diff(colLen))),abs(diff(colLen))))
      chosen_tmp <- cbind(combine_list[[chosen]],NA_matrix)
      colnames(chosen_tmp) <- colnames(combine_list[[NTchosen]])
      combine_list[[chosen]] <- chosen_tmp
    }
    return(list(combine_list[[1]],combine_list[[2]]))
  } #add NA to the one with less columns
  hsa_list <- asjust_length(tmp1 = P1hsa,
                            tmp2 = P2hsa)
  P1hsa <- hsa_list[[1]]; P2hsa <- hsa_list[[2]]
  hsb_list <- asjust_length(tmp1 = P1hsb,
                            tmp2 = P2hsb)
  P1hsb <- hsb_list[[1]]; P2hsb <- hsb_list[[2]]

  ## Get rid of P2 from the P1 files, P1 from the P2 files, and combine...
  P1hsa.hd <- P1hsa[1:(ploidy*lg),]; P2hsa.hd <- P2hsa[(ploidy*lg+1):(2*ploidy*lg),]
  P1hsb.hd <- P1hsb[1:(ploidy*lg),]; P2hsb.hd <- P2hsb[(ploidy*lg+1):(2*ploidy*lg),]

  P1hsa <- P1hsa[-(1:(2*ploidy*lg)),]; P2hsa <- P2hsa[-(1:(2*ploidy*lg)),]
  P1hsb <- P1hsb[-(1:(2*ploidy*lg)),]; P2hsb <- P2hsb[-(1:(2*ploidy*lg)),]

  ## Add difference in ploidy (d) to P2hsa
  P2hsa.hd[,ploidy] <- P2hsa.hd[,ploidy] + d
  P2hsa[,(ploidy + 1):ncol(P2hsa)] <- P2hsa[,(ploidy + 1):ncol(P2hsa)] + d

  P1set <- which(rep(c(rep(T,ploidy/2),rep(F,ploidy/2)),popnum*lg))
  P2set <- which(rep(c(rep(F,ploidy/2),rep(T,ploidy/2)),popnum*lg))

  P12hsa <- rbind(P1hsa[P1set,],P2hsa[P2set,])
  P12hsb <- rbind(P1hsb[P1set,],P2hsb[P2set,])
  reord <- order(P12hsa[,1],P12hsa[,2])

  P12hsa <- rbind(P1hsa.hd,P2hsa.hd,P12hsa[reord,])
  P12hsb <- cbind("",rbind(P1hsb.hd,P2hsb.hd,P12hsb[reord,])) #Space in first column is in original .hsb files

  write.table(P12hsa, paste0(filename,"_out.hsa"), sep = "\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
  write.table(P12hsb, paste0(filename,"_out.hsb"), sep = "\t", quote = FALSE, row.names=FALSE, col.names = FALSE)

  ######Write the .ped file (can use either of the parental .ped files)
  pedfl <- read.table(paste0(filename,"_P1.ped"), header = TRUE, stringsAsFactors = FALSE)
  write.table(pedfl, paste0(filename,"_out.ped"), quote = FALSE, row.names = FALSE, sep = "\t")

  ######Generate the genotypes file
  P1genotypes <- as.matrix(read.table(paste0(filename, "_P1_out_genotypes.dat"),
                                      header = TRUE, stringsAsFactors = FALSE,
                                      sep = "\t",row.names=1,colClasses = c("character")))
  P2genotypes <- as.matrix(read.table(paste0(filename, "_P2_out_genotypes.dat"),
                                      header = TRUE, stringsAsFactors = FALSE,
                                      sep = "\t",row.names=1,colClasses = c("character")))

  P1genotypes.parental <- P1genotypes[,1:ploidy]; P2genotypes.parental <- P2genotypes[,(ploidy+1):(2*ploidy)]
  P1genotypes <- P1genotypes[,-seq(2*ploidy)]; P2genotypes <- P2genotypes[,-seq(2*ploidy)]

  if(any(rownames(P1genotypes) != rownames(P2genotypes))) stop("Error in combining genotype.dat files!")

  P1set <- which(rep(c(rep(T,ploidy/2),rep(F,ploidy/2)),popnum))
  P2set <- which(rep(c(rep(F,ploidy/2),rep(T,ploidy/2)),popnum))

  colnames(P2genotypes)[P2set] <- paste0(substr(colnames(P2genotypes[,P2set]),1,
                                                nchar(colnames(P2genotypes[,P2set]))-1),(ploidy/2+1):(ploidy/2+ploidy/2))

  P12genotypes <- cbind(P1genotypes[,P1set],P2genotypes[,P2set])
  P12genotypes <- cbind(P1genotypes.parental,P2genotypes.parental,P12genotypes[,order(colnames(P12genotypes))])

  write.table(cbind("marker" = rownames(P12genotypes),P12genotypes),
              paste0(filename,"_out_genotypes.dat"), sep = "\t", row.names=FALSE)


  ## Tidy-up - delete the intermediate files:
  file.remove(paste0(filename,"_P1.ped"))
  file.remove(paste0(filename,"_P1.gen"));file.remove(paste0(filename,"_P2.gen"))
  file.remove(paste0(filename,"_P1.par"));file.remove(paste0(filename,"_P2.par"))
  file.remove(paste0(filename,"_P1_out.hsa"));file.remove(paste0(filename,"_P2_out.hsa"))
  file.remove(paste0(filename,"_P1_out.hsb"));file.remove(paste0(filename,"_P2_out.hsb"))
  # file.remove(paste0(filename,"_P1_out.ped"));file.remove(paste0(filename,"_P2_out.ped"))
  file.remove(paste0(filename,"_P1_out_founderalleles.dat"));file.remove(paste0(filename,"_P2_out_founderalleles.dat"))
  file.remove(paste0(filename,"_P1_out_alleledose.dat"));file.remove(paste0(filename,"_P2_out_alleledose.dat"))
  file.remove(paste0(filename,"_P1_out_allAlleledose.dat"));file.remove(paste0(filename,"_P2_out_allAlleledose.dat"))
  file.remove(paste0(filename,"_P1_out_genotypes.dat"));file.remove(paste0(filename,"_P2_out_genotypes.dat"))

  P1_map <- read.table(paste0(filename_P1,'.map'),header = TRUE,stringsAsFactors = FALSE,colClasses = c("character"))
  write.table(P1_map,paste(filename,c(".map"),sep=""), quote = FALSE, sep = " ", row.names = FALSE)
  file.remove(paste0(filename,"_P1.map"))
  print("Thank you. Please check your R working directory.")
}

run.PedigreeSim <- function(parfile, path.to.PedigreeSim="PedigreeSim.jar") {
  ##To execute (e.g.), type: run.PedigreeSim("myparfile.par")
  ps <- system2(command = "java",
                args = c("-jar",
                         path.to.PedigreeSim,
                         parfile),
                stdout = TRUE,
                stderr = TRUE)
} #run.PedigreeSim()


##Inheritance check
filter_inheritance <- function(res,
                               plot = TRUE,
                               inheritance_P1 = 'Poly',
                               inheritance_P2 = 'Poly',
                               threshold = 0.1){

  res_P1 <- res[res$parent %in% 'P1',]
  res_P2 <- res[res$parent %in% 'P2',]

  ketp1 <- as.character(res_P1[res_P1[[inheritance_P1]] > threshold,]$mrk)
  ketp2 <- as.character(res_P1[res_P1[[inheritance_P2]] > threshold,]$mrk)

  if(plot){
    plot(res$Poly, ylab = 'inheritance scores',main = 'Check segregation',pch = 16)
    abline(h = threshold,col = 'red')
  }

  kept <- intersect(ketp1, ketp2)
  writeLines(paste0('When threshold is ', threshold,','))
  writeLines(paste0(length(kept),' markers were kept'))
  writeLines(paste0(round((1 - length(kept)/nrow(res_P1)) * 100,2),'% (',
                    length(setdiff(res$mrk, kept)),' markers) were filtered, and they are:'))
  res <- res[res$mrk %in% kept,]
  return(res)
}

######Clustering######
#'
#'
#' @param linkage_res a list include input and output. In the input, it has: ploidy, Mrk_P, and Mrk_Q.
#' In the output, it includes the linkage estimation result. The linkage estimation result is a list.
#' Each list's name is the marker pair. In each list, it is a dataframe include: marker_P, marker_Q, rf, LOD, phasing_P, phasing_Q
#'
#' @return a dataframe include: marker_P, marker_Q, rf, LOD for all markers
#' @export
#'
#' @examples
linkagelist_to_linkage <- function(linkage_res){
  linkage <- do.call(rbind,lapply(names(linkage_res$output),function(pairs){
    marker_a <- strsplit(pairs,' x ')[[1]][1]
    marker_b <- strsplit(pairs,' x ')[[1]][2]
    r <- unique(linkage_res$output[[pairs]]$r)
    LOD <- unique(linkage_res$output[[pairs]]$LOD)
    if(length(LOD) > 1){
      LOD <- max(LOD)
    }
    data.frame(marker_a,marker_b,r,LOD)
  }))
  return(linkage)
}



######Phasing######
#' code_converter
#' according to possible MT in your dataset and make a list of all possible phasing code
#' @param MT a dataframe record each markers' marker type of two parents (MT_P1, MT_P2),the
#' numer of shared haplotypes between two parents, and the unique haplotype
#' @param ploidy numeric. crop ploidy level
#'
#' @return a list of dataframe. In each dataframe, it is two markers all possible pairing and the
#' corresponding code
#' @export
#'
#' @examples
code_converter <- function(MT,
                           ploidy = 4){
  all <- strsplit(MT$parents,'_')
  P1_possi <- unlist(lapply(all,function(t){t[1]}))
  P2_possi <- unlist(lapply(all,function(t){t[2]}))
  all_possi <- unique(rbind(expand.grid(P1_possi,P1_possi),expand.grid(P2_possi,P2_possi)))
  temp_list <- list()
  for(i in 1:nrow(all_possi)){
    letters_P <- sort(strsplit(as.character(all_possi[i,1]),'')[[1]]);  letters_Q <- sort(strsplit(as.character(all_possi[i,2]),'')[[1]])
    combs <- gtools::permutations(ploidy,ploidy)
    combs_P <- as.data.frame(unique(do.call(rbind,lapply(1:nrow(combs),function(i){
      letters_P[combs[i,]]
    }))))
    combs_Q <- as.data.frame(unique(do.call(rbind,lapply(1:nrow(combs),function(i){
      letters_Q[combs[i,]]
    }))))

    pairwise <- expand.grid(seq(1,nrow(combs_P)),seq(1,nrow(combs_Q)))

    unique_P <- unique(letters_P)
    unique_Q <- unique(letters_Q)
    all_loc <- expand.grid(unique_P,unique_Q)
    colnames(all_loc) <- c('P','Q')


    temp <- do.call(rbind,lapply(1:nrow(pairwise),function(i){
      P <- combs_P[pairwise[i,1],]
      Q <- combs_Q[pairwise[i,2],]


      score <- sapply(1:nrow(all_loc),function(j){
        P_loc <- which(P == as.character(all_loc[j,'P']))
        Q_loc <- which(Q == as.character(all_loc[j,'Q']))
        length(intersect(P_loc,Q_loc))
      })

      data.frame('P' = paste0(P,collapse = ''),
                 'Q' = paste0(Q,collapse = ''),
                 'score' = paste0(score,collapse = ''))
    }))

    temp_table <- dataframe_to_contigency(temp,'P','Q','score')
    temp_list[[paste0(paste0(letters_P,collapse = ''),'x',
                      paste0(letters_Q,collapse = ''))]] <- temp_table
  }
  return(temp_list)
}


#' find the order to change fix to observed
#'
#' @param ploidy crop ploidy level
#' @param fix dataframe/vector/numeric. It need to include the same element as observed
#' @param observed dataframe/vector/numeric. It need to include the same element as fix
#'
#' @return a dataframe of all possible ordering to make fix as observed
#' @export
#'
#' @examples
FindOrder <- function(ploidy,
                      fix =  fixed[1,],
                      observed = phasing[1,]){
  combs <- gtools::permutations(ploidy,ploidy)
  combs_letters_P1 <- do.call(rbind,lapply(1:nrow(combs),function(i){
    e <- fix[combs[i,]]
    colnames(e) <- NULL;rownames(e) <- NULL
    e <- as.matrix(e)
    if(nrow(e) > 1){
      t(e)
    }else{
      e
    }
  }))
  identified_order <- which(apply(combs_letters_P1,1,paste0,collapse = '') == paste0(observed,collapse = ''))
  return(combs[identified_order,,drop = FALSE])
}


#' find the possible phasing
#'
#' @param ele dataframe. it is a part of all pooled linkage information. It includes rows:
#' marker_P, marker_Q, rf, LOD, phasing_P, phasing_Q, probability, pair, phase_pair, score
#' - marker_P, marker_Q: marker name
#' - rf: estimated recombination frequency
#' - LOD: estimated LOD score
#' - phasing_P, phasing_Q: phasing of marker_P, marker_Q
#' - probability: probability of having this phasing, estimated from linkage estimation function
#' - pair: put marker_P and marker_Q together to name this marker pair
#' - score: phasing code, which can refer to the unique phasing code table
#' It is generated from 'pool_linkage'
#' @param m the markername of the marker that we want to look at
#' @param p_list a list of vectors In each list, it is each markers' all possible phasing that is calculated
#' from its marker type. It is generated from 'prepare_phasing_list'
#' @param temp_list a list of dataframe. In each dataframe, it is two markers all possible pairing and the
#' corresponding code.This is generated from 'code_converted'
#'
#' @return a string which is the possible phasing
#' @export
#'
#' @examples
find_possibility <- function(ele,
                             m = other_marker,
                             p_list,
                             temp_list){
  if(nrow(ele) == 1){
    ele1 <- ele[1,]
    target <- which(ele1[,c(1,2)] == m);rest <- setdiff(c(1,2),target)
    target <- c('P','Q')[target];rest <- c('P','Q')[rest]
    rest_marker <- as.character(ele1[,paste0('marker_',rest)])

    phase_target <- p_list[[m]]
    temp_table <- temp_list[[ele1$phase_pair]]
    if(target == 'P'){
      possibility <- names(which(temp_table[phase_target,] == ele1$score))
      if(nrow(temp_table) == 1){
        possibility <- colnames(temp_table[phase_target,,drop = FALSE] == ele1$score)
      }
    }else{
      possibility <- rownames(temp_table)[which(temp_table[,phase_target,drop = FALSE] == ele1$score)]
      if(nrow(temp_table) == 1){
        possibility <- rownames(temp_table[,phase_target,drop = FALSE] == ele1$score)
      }
    }

  }else{
    possibility <- c()
    for(i in 1:nrow(ele)){
      ele1 <- ele[i,]
      target <- which(ele1[,c(1,2)] == m);rest <- setdiff(c(1,2),target)
      target <- c('P','Q')[target];rest <- c('P','Q')[rest]
      rest_marker <- as.character(ele1[,paste0('marker_',rest)])

      phase_target <- p_list[[m]]
      temp_table <- temp_list[[ele1$phase_pair]]
      if(target == 'P'){
        possi <- names(which(temp_table[phase_target,] == ele1$score))
      }else{
        possi <- names(which(temp_table[,phase_target] == ele1$score))
      }
      possibility <- c(possibility,possi)
    }
  }
  return(unique(possibility))
}


#' find the most significant marker with the 1st and 2nd marker
#'
#' @param marker_1 the first marker's name
#' @param marker_2 the second marker's name
#' @param pairs_tmp a dataframe which include two columns. each row is an occured marker pairs in the linkage
#' estimation. Two columns are the name of markers in this pair
#' @param temp dataframe. it is a part of all pooled linkage information. It includes rows:
#' marker_P, marker_Q, rf, LOD, phasing_P, phasing_Q, probability, pair, phase_pair, score
#' - marker_P, marker_Q: marker name
#' - rf: estimated recombination frequency
#' - LOD: estimated LOD score
#' - phasing_P, phasing_Q: phasing of marker_P, marker_Q
#' - probability: probability of having this phasing, estimated from linkage estimation function
#' - pair: put marker_P and marker_Q together to name this marker pair
#' - score: phasing code, which can refer to the unique phasing code table
#' It is generated from 'pool_linkage'
#' @param markers_pool vector of characters which include the markernames which need to be phased
#'
#' @return the markername of the marker which has the most significant linkage with 1st and 2nd marker
#' @export
#'
#' @examples
significant_marker <- function(marker_1 = start_marker,
                               marker_2 = other_marker,
                               pairs_tmp,
                               temp,
                               markers_pool){
  pairs_1 <-  names(which(apply(pairs_tmp ==marker_1,1,sum) > 0))
  temp1 <- temp[temp$pair %in% pairs_1,]
  temp1$markers <- gsub(' x ','',gsub(marker_1,'',temp1$pair))
  temp1 <- do.call(rbind,lapply(unique(temp1$markers),function(m){
    LOD <- max(temp1[temp1$markers %in% m,]$LOD)
    data.frame('marker' = m,LOD)
  }))
  rownames(temp1) <- temp1$marker


  pairs_2 <-  names(which(apply(pairs_tmp ==marker_2,1,sum) > 0))
  temp2 <- temp[temp$pair %in% pairs_2,]
  temp2$markers <- gsub(' x ','',gsub(marker_2,'',temp2$pair))
  temp2 <- do.call(rbind,lapply(unique(temp2$markers),function(m){
    LOD <- max(temp2[temp2$markers %in% m,]$LOD)
    data.frame('marker' = m,LOD)
  }))
  rownames(temp2) <- temp2$marker

  allmrks <- unique(c(rownames(temp1),rownames(temp2)))
  allmrks <- intersect(allmrks,markers_pool)
  temp1 <- temp1[allmrks,];temp2 <- temp2[allmrks,]
  scores <- temp1$LOD * temp2$LOD

  return(allmrks[which.max(scores)])
}


#' generate the all possible phasing of all markers
#'
#' @param ploidy numeric. crop ploidy level
#' @param poly formatted list of each haploblock (temp, parent)
#'
#' @return a list of vectors In each list, it is each markers' all possible phasing that is calculated
#' from its marker type
#' @export
#'
#' @examples
prepare_phasing_list <- function(ploidy,
                                 poly){
  P1_list <- P2_list <-list()
  for(marker in names(poly)){
    P1_letters <- sort(as.character(poly[[marker]]$parent_info[1,]))
    P2_letters <- sort(as.character(poly[[marker]]$parent_info[2,]))

    combs <- gtools::permutations(ploidy,ploidy)
    combs_P1 <- as.data.frame(unique(do.call(rbind,lapply(1:nrow(combs),function(i){
      P1_letters[combs[i,]]
    }))))
    combs_P2 <- as.data.frame(unique(do.call(rbind,lapply(1:nrow(combs),function(i){
      P2_letters[combs[i,]]
    }))))

    combs_P1 <- apply(unique(combs_P1),1,paste0,collapse = '')
    combs_P2 <- apply(unique(combs_P2),1,paste0,collapse = '')
    P1_list[[marker]] <- as.character(combs_P1)
    P2_list[[marker]] <- as.character(combs_P2)
  }
  return(list('P1' = P1_list,
              'P2' = P2_list))
}


#' according to linkage result to find all possible pairs
#'
#' @param linkage_res a list include input and output. In the input, it has: ploidy, Mrk_P, and Mrk_Q.
#' In the output, it includes the linkage estimation result. The linkage estimation result is a list.
#' Each list's name is the marker pair. In each list, it is a dataframe include: marker_P, marker_Q, rf, LOD, phasing_P, phasing_Q
#'
#' @return a list include each markername, and each markers' all possible linked pairs
#' @export
#'
#' @examples
prepare_pairwise_order <- function(linkage_res = Linkage_all$output){
  linkage <- do.call(rbind,lapply(names(linkage_res),function(pairs){
    marker_a <- strsplit(pairs,' x ')[[1]][1]
    marker_b <- strsplit(pairs,' x ')[[1]][2]
    r <- unique(linkage_res[[pairs]]$r)
    LOD <- unique(linkage_res[[pairs]]$LOD)
    if(length(LOD) > 1){
      LOD <- max(LOD)
    }
    data.frame(marker_a,marker_b,r,LOD)
  }))
  markers <- names(sort(table(unlist(c(linkage[linkage$LOD > 3,c(1,2)]))),decreasing = TRUE))
  markers <- c(markers,setdiff(names(linkage_res),markers))
  linkage_used <- linkage_res
  linkage_reformat <- do.call(rbind,strsplit(names(linkage_used),' x '))
  all_used <- which(linkage_reformat[,1] %in% markers & linkage_reformat[,2] %in% markers)
  #make the order of the list
  all_possible_pairs <- names(linkage_used[all_used])
  sort_order <- list()
  i <- 1
  pairs_tmp <- do.call(rbind,strsplit(all_possible_pairs, ' x '))
  rownames(pairs_tmp) <- all_possible_pairs
  while(length(all_possible_pairs) != 0){
    # print(i)
    chosen_pairs <- names(which(apply(pairs_tmp == markers[[i]],1,sum) > 0))
    sort_order[[markers[[i]]]] <- intersect(chosen_pairs,all_possible_pairs)
    all_possible_pairs <- setdiff(all_possible_pairs,chosen_pairs)
    i <- i+1
  }
  return(sort_order)
}


#' pool all linkage results together
#'
#' @param poly formatted list of each haploblock (temp, parent)
#' @param sort_order a list include each markername, and each markers' all possible linked pairs.
#' It is generated from 'prepare_pairwise_order'
#' @param temp_list a list of dataframe. In each dataframe, it is two markers all possible pairing and the
#' @param linkage_res a list include input and output. In the input, it has: ploidy, Mrk_P, and Mrk_Q.
#' In the output, it includes the linkage estimation result. The linkage estimation result is a list.
#' Each list's name is the marker pair. In each list, it is a dataframe include: marker_P, marker_Q, rf, LOD, phasing_P, phasing_Q
#' @return dataframe. it is a part of all pooled linkage information. It includes rows:
#' marker_P, marker_Q, rf, LOD, phasing_P, phasing_Q, probability, pair, phase_pair, score
#' - marker_P, marker_Q: marker name
#' - rf: estimated recombination frequency
#' - LOD: estimated LOD score
#' - phasing_P, phasing_Q: phasing of marker_P, marker_Q
#' - probability: probability of having this phasing, estimated from linkage estimation function
#' - pair: put marker_P and marker_Q together to name this marker pair
#' - score: phasing code, which can refer to the unique phasing code table
#' @export
#'
#' @examples
pool_linkage <- function(poly,
                         sort_order,
                         temp_list,
                         linkage_res){
  P1_all <- P2_all <- data.frame()
  for(mrk in names(poly)){
    group_P1 <- c();group_P2 <- c()
    multiple_solution_p1 <- multiple_solution_p2 <- c()

    for(pair in sort_order[[mrk]]){
      chosen <- linkage_res[[pair]]
      rest_marker <- setdiff(unique(c(as.character(chosen$marker_P),as.character(chosen$marker_Q))),mrk)
      P_phasing <- strsplit(as.character(chosen$phasing_P),'_');Q_phasing <- strsplit(as.character(chosen$phasing_Q),'_')
      # decide_P1 <- unlist(lapply(1:length(P_phasing),function(i){
      #   P_phasing[[i]][1] == Q_phasing[[i]][1]
      # }))
      # decide_P2 <- unlist(lapply(1:length(Q_phasing),function(i){
      #   P_phasing[[i]][2] == Q_phasing[[i]][2]
      # }))
      # if(any(decide_P1 == TRUE)){
      #   group_P1 <- c(group_P1,rest_marker)
      # }
      # if(any(decide_P2 == TRUE)){
      #   group_P2 <- c(group_P2,rest_marker)
      # }

      # if(length(unique(unlist(lapply(Q_phasing,function(l){l[1]}))))>1){
      #   multiple_solution_p1 <- c(multiple_solution_p1,rest_marker)
      # }
      # if(length(unique(unlist(lapply(Q_phasing,function(l){l[2]}))))>1){
      #   multiple_solution_p2 <- c(multiple_solution_p2,rest_marker)
      # }

      P_phasing1 <- do.call(rbind,P_phasing);Q_phasing1 <- do.call(rbind,Q_phasing)

      P1_phasing <- cbind(chosen[,1:4],P_phasing1[,1],Q_phasing1[,1]);P2_phasing <- cbind(chosen[,1:4],P_phasing1[,2],Q_phasing1[,2])
      colnames(P1_phasing)[5:6] <- colnames(P2_phasing)[5:6] <- c('phasing_P','phasing_Q')

      P1_phasing$probability <- 1/nrow(P1_phasing); P2_phasing$probability <- 1/nrow(P2_phasing)



      # P1_phasing$chosen <- ifelse(decide_P1,1,0); P2_phasing$chosen <- ifelse(decide_P2,1,0)
      P1_phasing$pair <- paste0(P1_phasing$marker_P,' x ',P1_phasing$marker_Q)
      P2_phasing$pair <- paste0(P2_phasing$marker_P,' x ',P2_phasing$marker_Q)

      colnames(P1_phasing)[3] <- 'rf';colnames(P2_phasing)[3] <- 'rf'
      P1_all <- rbind(P1_all,P1_phasing)
      P2_all <- rbind(P2_all,P2_phasing)
    }
  }
  P1_all <- cbind(P1_all,do.call(rbind,lapply(1:nrow(P1_all),function(i){
    P_phase <- paste0(sort(strsplit(P1_all[i,'phasing_P'],'')[[1]]),collapse = '')
    Q_phase <- paste0(sort(strsplit(P1_all[i,'phasing_Q'],'')[[1]]),collapse = '')
    phase_pair <- paste0(P_phase,'x',Q_phase)
    temp_table <- temp_list[[paste0(P_phase,'x',Q_phase)]]
    score <- temp_table[P1_all[i,'phasing_P'],P1_all[i,'phasing_Q']]
    data.frame(phase_pair,score)
  })))
  P2_all <- cbind(P2_all,do.call(rbind,lapply(1:nrow(P2_all),function(i){
    P_phase <- paste0(sort(strsplit(P2_all[i,'phasing_P'],'')[[1]]),collapse = '')
    Q_phase <- paste0(sort(strsplit(P2_all[i,'phasing_Q'],'')[[1]]),collapse = '')
    phase_pair <- paste0(P_phase,'x',Q_phase)
    temp_table <- temp_list[[paste0(P_phase,'x',Q_phase)]]
    score <- temp_table[P2_all[i,'phasing_P'],P2_all[i,'phasing_Q']]
    data.frame(phase_pair,score)
  })))
  return(list('P1' = P1_all,
              'P2' = P2_all))
}



######IBD estimation######
plot_IBD <- function(IBD = proba_res3,
                     ploidy = ploidy,
                     ind = individual){
  max_position_chromosome <- max(IBD$position)
  min_position_chromsome <- min(IBD$position)


  LG <- colnames(IBD)[1:(2*ploidy)]
  #define the chromosome
  max_position_chromosome_each <- max_position_chromosome
  min_position_chromsome_each <- min_position_chromsome
  min <- 0
  max <- round(max_position_chromosome + round(max_position_chromosome/20))
  axis_number <- seq(min,max,by=ceiling((max-min)/10))

  plot(NULL, ylim = c(0,round(max_position_chromosome+round(max_position_chromosome/10))),
       xlim = c(0,ploidy*2 + 1), yaxt = "n",
       bty = "n",
       xaxt = "n",
       xlab = '', ylab ="", main = NULL)
  axis(side=2,at=axis_number,cex.axis = 2)

  for(j in seq_along(LG)){
    chrm.width <- 0.08
    chrm.offset <- 0.1 + j - 1
    symbols(x= chrm.offset, y=min_position_chromsome_each, circles = chrm.width, bg="white", add=TRUE, inches = FALSE)
    symbols(x= chrm.offset, y=max_position_chromosome_each, circles = chrm.width, bg="white", add=TRUE, inches = FALSE)
    rect((chrm.offset - chrm.width), (min_position_chromsome_each - chrm.width) , (chrm.offset + chrm.width), (max_position_chromosome_each + chrm.width) , col = "white")

    colors <- obtain_color_IBD(score = IBD[,j])

    #draw chromosome out
    for(i in 1:length(IBD[,j])){
      haplotype_size <- 1 #represent 1 cM
      rect((chrm.offset - chrm.width), (IBD[i,'position']- haplotype_size) ,
           (chrm.offset + chrm.width), (IBD[i,'position'] + haplotype_size) , col = colors[i])
    }
  }
  axis(1, at=  0.1 + seq_along(LG) - 1, labels= LG, cex.lab = 1.5, cex.axis = 1.5)
  title(xlab = "LG", ylab = "Position(cM)", main = ind,cex.main = 2,cex.lab = 1.5)
  par(new = FALSE)
}

obtain_color_IBD <- function(score = IBD$P1_1){

  if("colorspace" %in% rownames(installed.packages())){
    library(colorspace)
  }else{
    install.packages('colorspace')
    library(colorspace)
  }
  #set the color matrix
  maxColorValue <- length(seq(0,1,0.01))
  palette <- colorspace::sequential_hcl(maxColorValue)
  color <- unlist(sapply(score,function(i){
    seq_choice <- round(seq(0,1,0.01),2)
    #-1: disomic - blue
    #1: polysomic - red
    if(round(i,2) %in% seq_choice){
      palette[which(as.numeric(i) == seq_choice)]
    }else{
      print(i)
    }
  }))
  color
}



