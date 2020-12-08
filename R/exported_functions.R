######Import data######
#' readDatfile
#' Importing the data with R function read.table.
#' @param file the name of the file which the data are to be read from.
#' @param header a logical value indicating whether the file contains the names of the variables as its first line.
#' @param sep the field separator character.
#' @param check.names logical. If TRUE then the names of the variables in the data frame are checked to
#' ensure that they are syntactically valid variable names.
#'
#' @return a matrix contain the dosage information for all individuals,
#' according to their markers and alleles.
#' @export
#'
#' @examples
readDatfile <- function (file, header = TRUE, sep = "\t", check.names = FALSE,
                         ...) {
  return(read.table(file = file, header = header, sep = sep,
                    check.names = check.names, ...))
}


#' format_results
#' format result from PedigreeSim format to Multi-allelic mode inheritance format
#'
#' @param dosages Pedigreesim simulated multi-allelic marker format. Markers as row, individuals as column.
#' @param ploidy crop ploidy level
#' @param Parent1 the column name of parent 1
#' @param Parent2 the column name of parent 2
#' @param ncores number of cores used for running, suggested to use -2 cores of your computers
#'
#' @return formatted list of each haploblock (temp, parent)
#' @export
#'
#' @examples
format_results <- function(dosages = poly_dosages,ploidy = 4,Parent1,Parent2,ncores){
  mrk_name <- as.character(unique(dosages$marker))

  #run parallel function set up
  library(foreach)
  library(doParallel)
  win <- Sys.info()["sysname"] == "Windows"
  if (win) {
    cl <- parallel::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)
  } else {
    doParallel::registerDoParallel(cores = ncores)
  }

  A <- foreach::foreach(i = seq_along(mrk_name), .multicombine=FALSE,
                        .export = NULL) %dopar% {
                          mrk <- mrk_name[i]
                          temp <- dosages[dosages$marker %in% mrk,]
                          #identify the NA
                          NA_col <- names(which(colSums(is.na(temp)) > 0))
                          non_NA <- names(which(colSums(is.na(temp)) == 0))

                          #non NA part
                          temp_nonNA <- temp[,non_NA]
                          class <- as.character(temp_nonNA$allele)
                          temp_convert <- as.data.frame(t(sapply(3:ncol(temp_nonNA),function(j){
                            unlist(sapply(1:nrow(temp_nonNA),function(i){
                              if(!is.na(temp_nonNA[i,j])){
                                replicate(temp_nonNA[i,j],class[i])
                              }
                            }))
                          })))
                          rownames(temp_convert) <- colnames(temp)[3:ncol(temp_nonNA)]
                          #missing part
                          temp_missing <- matrix(nrow = length(NA_col),ncol = ploidy)
                          rownames(temp_missing) <- NA_col
                          colnames(temp_missing) <- colnames(temp_convert)
                          #fill possibly not all missing
                          if(length(NA_col) > 0){
                            if(colSums(!is.na(temp[,NA_col,drop = FALSE])) > 0){
                              writeLines(paste0('Please check ',mrk))
                            }
                          }

                          temp_res <- rbind(temp_convert,temp_missing)
                          rownames(temp_res) <- c(non_NA[3:ncol(temp_nonNA)],NA_col)
                          temp_res <- temp_res[colnames(dosages)[3:ncol(dosages)],]


                          format_list <- list()
                          format_list[["parent_info"]] <- temp_res[c(Parent1,Parent2),]
                          format_list[["temp"]] <- temp_res
                          format_list
                        }
  names(A) <- mrk_name
  return(A)
}

######Extract MT & Sharing######
#' Have an overview on marker type and the sharing of haplotypes between two parents
#'
#' @param poly
#' @param ploidy crop ploidy level
#' @param plot it can be either True/False. If true, the overview graph will be generated
#'
#' @return a dataframe record each markers' marker type of two parents (MT_P1, MT_P2),the
#' numer of shared haplotypes between two parents, and the unique haplotype
#'
#' @export
#'
#' @examples
MarkerType_Overview <- function(poly,ploidy,
                                plot = TRUE){
  res <- do.call(rbind,lapply(names(poly),function(mrk){
    parent <- poly[[mrk]]$parent_info
    parent <- apply(parent,2,as.character)
    P1_MT <- decide_MT(paste0(parent[1,],collapse = ''));P2_MT <- decide_MT(paste0(parent[2,],collapse = ''))

    P1 <- c(parent[1,]);P2 <- c(parent[2,])

    count_P1 <- table(P1)
    count_P2 <- table(P2)
    common <- intersect(P1,P2)

    sharing <- sum(unlist(sapply(common,function(c){
      min(count_P1[c],count_P2[c])
    })))/ploidy

    unique_hap <- sum(table(c(P1,P2)) == 1)
    data.frame(mrk,P1_MT,P2_MT,sharing,unique_hap)
  }))
  res$sharing <- res$sharing * ploidy
  res$MT <- paste0(apply(res[,2:3],1,paste,collapse = ' x '),'_Shar',
                   res$sharing)
  parents <- sapply(poly,function(m){
    #identify all the possible parents genotype
    paste0(apply(apply(m$parent_info,1,sort),2,paste,collapse = ''),collapse = '_')
  })
  parents <- data.frame('marker' = names(poly),parents)
  parents <- cbind(parents,res[,2:6])
  rownames(parents) <- NULL
  if(plot){
    library(ggplot2)
    P1_plot <- as.data.frame(table(res$P1_MT))
    pic1 <- ggplot(P1_plot, aes(x=Var1 , y=Freq, fill=Var1)) +
      geom_bar(stat="identity")+
      xlab("Marker type") +
      ylab("Nbr of markers") +
      labs(subtitle = paste0(sum(P1_plot$Freq),' markers'))+
      ggtitle("P1") +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(plot.title = element_text(size= 25,face="bold"),plot.subtitle = element_text(size = 18))+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"))+
      theme(legend.position = "none")

    P2_plot <- as.data.frame(table(res$P2_MT))
    pic2 <- ggplot(P2_plot, aes(x=Var1 , y=Freq, fill=Var1)) +
      geom_bar(stat="identity")+
      xlab("Marker type") +
      ylab("Nbr of markers") +
      labs(subtitle = paste0(sum(P1_plot$Freq),' markers'))+
      ggtitle("P2") +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(plot.title = element_text(size= 25,face="bold"),plot.subtitle = element_text(size = 18))+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"))+
      theme(legend.position = "none")

    Sharing <- as.data.frame(table(res$sharing))
    pic3 <- ggplot(Sharing, aes(x=Var1 , y=Freq, fill=Var1)) +
      geom_bar(stat="identity")+
      xlab("Sharing haplotypes") +
      ylab("Nbr of markers") +
      labs(subtitle = paste0(sum(P1_plot$Freq),' markers'))+
      ggtitle("P1 & P2") +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(plot.title = element_text(size= 25,face="bold"),plot.subtitle = element_text(size = 18))+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"))+
      theme(legend.position = "none")

    Summary <- as.data.frame(table(res$MT))
    colnames(Summary) <- c('MarkerType', 'Frequency')
    pic4 <- ggplot(Summary, aes(x=MarkerType , y=Frequency, fill=MarkerType)) +
      geom_bar(stat="identity")+
      xlab("Marker type") +
      ylab("Nbr of markers") +
      labs(subtitle = paste0(sum(P1_plot$Freq),' markers'))+
      ggtitle("P1 & P2") +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(plot.title = element_text(size= 25,face="bold"),plot.subtitle = element_text(size = 18))+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"))+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())

    multiplot(pic1, pic2,pic3,pic4, layout = matrix(seq(1,4),nrow=2, byrow=TRUE))
    print(Summary)
  }
  return(parents)
}


#########Equation generator#########
#' Generate all possible equations that will be used in your datasets
#' Here it made use of the generalized equation and your own dataset to create all possible equations that probably used later.
#' @param P1_inheritance 'polysomic' or 'disomic'. It will depend on this to generate different equations
#' @param P2_inheritance 'polysomic' or 'disomic'. It will depend on this to generate different equations
#' @param MT a dataframe include markername, parent genotype, P1 and P2 marker type, the number of
#' shared haplotypes, the number of unique haplotypes, combined marker type
#' @param ploidy crop ploidy level
#'
#' @return a list include:
#'  - equations: equations for P1 and P2 separately
#'  - input: P1 and P2's inheritance, ploidy level, and marker type (MT)
#' @export
#'
#' @examples
all_equations_generator <- function(P1_inheritance = 'Polysomic',
                                    P2_inheritance = 'Disomic',
                                    MT = MT,
                                    ploidy){
  #take the unique ones
  all_pairs <- as.character(unique(MT$parents))

  #find the possible combinations of two parents in P1 and P2 separately
  mrk_P1 <- unlist(lapply(strsplit(all_pairs,'_'),function(each){each[1]}))
  P1_mrkcomb <- rbind(t(combn(mrk_P1,2)),replicate(2,mrk_P1))
  mrk_P2 <- unlist(lapply(strsplit(all_pairs,'_'),function(each){each[2]}))
  P2_mrkcomb <- rbind(t(combn(mrk_P2,2)),replicate(2,mrk_P2))

  #obtain the generalized function
  if(any(c(P1_inheritance,P2_inheritance) %in% 'Polysomic')){
    polysomic_equation <- bivalent_equation_generator(ploidy = ploidy)
  }
  if(any(c(P1_inheritance,P2_inheritance) %in% 'Disomic')){
    disomic_equation <- disomic_equation_generator(ploidy = ploidy)
  }

  #all equations P1
  if(P1_inheritance == 'Polysomic'){
    all_equation_P1 <- lapply(1:nrow(P1_mrkcomb),function(i){
      generate_all_equation_polysomic(P_genotype = strsplit(P1_mrkcomb[i,1],'')[[1]],
                                      Q_genotype = strsplit(P1_mrkcomb[i,2],'')[[1]],
                                      polysomic_equation)
    })
    names(all_equation_P1) <- sapply(1:nrow(P1_mrkcomb),function(i){paste0(P1_mrkcomb[i,],collapse = '_')})
  }else{#disomic inheritance
    all_equation_P1 <- lapply(1:nrow(P1_mrkcomb),function(i){
      generate_all_equation_disomic(P_genotype = strsplit(P1_mrkcomb[i,1],'')[[1]],
                                    Q_genotype = strsplit(P1_mrkcomb[i,2],'')[[1]],
                                    disomic_equation)
    })
    names(all_equation_P1) <- sapply(1:nrow(P1_mrkcomb),function(i){paste0(P1_mrkcomb[i,],collapse = '_')})
  }
  #all equations P2
  if(P2_inheritance == 'Polysomic'){
    all_equation_P2 <- lapply(1:nrow(P2_mrkcomb),function(i){
      generate_all_equation_polysomic(P_genotype = strsplit(P2_mrkcomb[i,1],'')[[1]],
                                      Q_genotype = strsplit(P2_mrkcomb[i,2],'')[[1]],
                                      polysomic_equation)
    })
    names(all_equation_P2) <- sapply(1:nrow(P2_mrkcomb),function(i){paste0(P2_mrkcomb[i,],collapse = '_')})
  }else{#disomic inheritance
    all_equation_P2 <- lapply(1:nrow(P2_mrkcomb),function(i){
      generate_all_equation_disomic(P_genotype = strsplit(P2_mrkcomb[i,1],'')[[1]],
                                    Q_genotype = strsplit(P2_mrkcomb[i,2],'')[[1]],
                                    disomic_equation)
    })
    names(all_equation_P2) <- sapply(1:nrow(P2_mrkcomb),function(i){paste0(P2_mrkcomb[i,],collapse = '_')})
  }
  #prepare to write the output
  equations <- list('P1' = all_equation_P1,
                    'P2' = all_equation_P2)

  input <- list('P1_inheritance'= P1_inheritance,
                'P2_inheritance' = P2_inheritance,
                'ploidy' = ploidy,
                'MT' = MT)

  return(list('equations' = equations,
              'input'= input))
}
#'bivalent_equation_generator
#' Generate equation according to ploidy level under assumed bivalent pairing
#' @param ploidy crop ploidy level
#'
#' @return a list include two formats of results (probability, recombination) and include
#' format 1 is contingency table. format 2 is melt table
#' @export
#'
#' @examples
bivalent_equation_generator <- function(ploidy){
  allele <- seq(1,ploidy)
  gametes <- apply(combn(allele,ploidy/2),2,paste,collapse='')

  #creates an empty table
  probability_table <- matrix(0,nrow = length(gametes),ncol = length(gametes))
  rownames(probability_table) <- colnames(probability_table) <- gametes

  #basic pairwise markers fixed phasing
  tmp <- t(replicate(2,as.character(seq(1,ploidy))))
  colnames(tmp) <- paste0('H',seq(1,ploidy))

  #found bivalent configuration
  pairing_method <- create_bivalent_pairing(ploidy)
  pairing_method <- pairing_method$format2

  #basic pairwise markers fixed phasing
  tmp <- t(replicate(2,as.character(seq(1,ploidy))))
  colnames(tmp) <- paste0('H',seq(1,ploidy))

  #fill the R,S in:
  pairing_gamete <- list();res_all <- data.frame()
  for(p in pairing_method){
    ###generalized all ploidy model
    parts <- strsplit(p,'\\|')[[1]]
    parts_recombi <- list()
    for(pi in 1:length(parts)){
      hom1 <- strsplit(parts[pi],'_')[[1]][1];hom2 <- strsplit(parts[p],'_')[[1]][2]
      parts_recombi[[paste0('part',pi)]] <- generate_recombi_nonrecombi(hom_tmp = tmp[,strsplit(parts[pi],'_')[[1]]],                                                                 ploidy = ploidy)
    }
    x <- replicate(length(parts_recombi),names(parts_recombi$part1))
    gametes_reformation <- apply(expand.grid(split(x, rep(1:ncol(x), each = nrow(x)))),2,as.character)
    colnames(gametes_reformation) <- names(parts_recombi)


    unique_nme <-  unique(unlist(lapply(1:nrow(gametes_reformation), function(i){
      paste0(sort(gametes_reformation[i,]),collapse = '_')
    })))
    res <- vector(mode = "list", length = length(unique_nme))
    names(res) <- unique_nme

    res_df <- data.frame()
    for(i in 1:nrow(gametes_reformation)){
      each <- gametes_reformation[i,]
      possible_choice <- lapply(names(each),function(l){
        parts_recombi[[l]][[as.character(each[l])]]
      })
      allcombs <- apply(expand.grid(possible_choice),2,as.character)

      df <- data.frame();possible_gametes <- c()
      for(n in 1:nrow(allcombs)){
        cob <- sapply(apply(allcombs[n,,drop = FALSE],2,strsplit,''),"[[", 1)
        cob1 <- apply(cob,1,paste,collapse='')

        #generate recombination
        count <- sum(each %in% '1/2*r') #check how many recombination is there

        #generate probability
        locs <- apply(apply(cob,1,sort),2,paste,collapse='')
        locP <- locs[1]
        locQ <- locs[2]
        probability <- paste0('1/',length(pairing_method),'*',
                              paste0(each,collapse = '*'))
        probability <- gsub(' ','',Deriv::Simplify(gsub('s','(1-r)',probability)))
        probability_table[locs[1],locs[2]] <- paste0(probability_table[locs[1],locs[2]],
                                                     '+',probability)

        df <- rbind(df,data.frame(locP,locQ,
                                  probability,count))
        possible_gametes <- c(possible_gametes,paste0(locP,'_',locQ))
      }
      res[[paste0(sort(each),collapse = '_')]] <- c(res[[paste0(sort(each),collapse = '_')]],possible_gametes)
      res_df <- rbind(res_df,df)
    }
     res_df$pairing <- p
    res_all <- rbind(res_all,res_df)
    pairing_gamete[[p]] <- res
  }
  #simplify it
  probability_table <- gsub(' ','',apply(probability_table, c(1,2), Deriv::Simplify))

  #generate count_table
  library(dplyr);library(reshape2)
  count_matrix <- res_all %>%
    group_by(locP,locQ) %>%
    summarise_at(vars('count'), mean)
  count_matrix1 <- dataframe_to_contigency(as.data.frame(count_matrix),
                                           'locP','locQ','count')
  count_matrix1 <- count_matrix1[gametes,gametes]

  format1 <- list('recombination' = count_matrix1,
                  'probability' = probability_table)

  # format2 <-
  format2 <- cbind(melt(count_matrix1),melt(probability_table)[,3])
  colnames(format2) <- c('markerP','markerQ','recombination','probability')
  input <- list('ploidy'=ploidy)

  return(list('format1' = format1,
              'format2' = format2,
              'input' = input))
}

#' generate equation according to ploidy level under assumed disomic inheritance
#'
#' @param ploidy crop ploidy level
#'
#' @return  a list include two formats of results (probability, recombintion) and include
#' format 1 is contigency table. format 2 is melt table. Each format is a list include all po
#' possible pairing's result.
#' @export
#'
#' @examples
disomic_equation_generator <- function(ploidy){
  allele <- seq(1,ploidy)
  gametes <- apply(combn(allele,ploidy/2),2,paste,collapse='')

  #basic pairwise markers fixed phasing
  tmp <- t(replicate(2,as.character(seq(1,ploidy))))
  colnames(tmp) <- paste0('H',seq(1,ploidy))

  #found bivalent configuration
  pairing_method <- create_bivalent_pairing(ploidy)
  pairing_method <- pairing_method$format2

  #basic pairwise markers fixed phasing
  tmp <- t(replicate(2,as.character(seq(1,ploidy))))
  colnames(tmp) <- paste0('H',seq(1,ploidy))

  #fill the R,S in:
  pairing_gamete <- list();res_all <- data.frame()
  for(p in pairing_method){
    ###generalized all ploidy model
    parts <- strsplit(p,'\\|')[[1]]
    parts_recombi <- list()
    for(pi in 1:length(parts)){
      hom1 <- strsplit(parts[pi],'_')[[1]][1];hom2 <- strsplit(parts[p],'_')[[1]][2]
      parts_recombi[[paste0('part',pi)]] <- generate_recombi_nonrecombi(hom_tmp = tmp[,strsplit(parts[pi],'_')[[1]]],                                                                 ploidy = ploidy)
    }
    x <- replicate(length(parts_recombi),names(parts_recombi$part1))
    gametes_reformation <- apply(expand.grid(split(x, rep(1:ncol(x), each = nrow(x)))),2,as.character)
    colnames(gametes_reformation) <- names(parts_recombi)


    unique_nme <-  unique(unlist(lapply(1:nrow(gametes_reformation), function(i){
      paste0(sort(gametes_reformation[i,]),collapse = '_')
    })))
    res <- vector(mode = "list", length = length(unique_nme))
    names(res) <- unique_nme

    res_df <- data.frame()
    for(i in 1:nrow(gametes_reformation)){
      each <- gametes_reformation[i,]
      possible_choice <- lapply(names(each),function(l){
        parts_recombi[[l]][[as.character(each[l])]]
      })
      allcombs <- apply(expand.grid(possible_choice),2,as.character)

      df <- data.frame();possible_gametes <- c()
      for(n in 1:nrow(allcombs)){
        cob <- sapply(apply(allcombs[n,,drop = FALSE],2,strsplit,''),"[[", 1)
        cob1 <- apply(cob,1,paste,collapse='')

        #generate recombination
        # count <- decide_recombintion(cob1[1],cob1[2]) #method 1: from the haplotype
        count <- sum(each %in% '1/2*r') #method 2: check how many recombination is there

        #generate probability
        locs <- apply(apply(cob,1,sort),2,paste,collapse='')
        locP <- locs[1]
        locQ <- locs[2]
        probability <- paste0(paste0(each,collapse = '*'))
        probability <- gsub(' ','',Deriv::Simplify(gsub('s','(1-r)',probability)))

        df <- rbind(df,data.frame(locP,locQ,
                                  probability,count))
        possible_gametes <- c(possible_gametes,paste0(locP,'_',locQ))
      }
      res[[paste0(sort(each),collapse = '_')]] <- c(res[[paste0(sort(each),collapse = '_')]],possible_gametes)
      res_df <- rbind(res_df,df)
    }
    res_df$pairing <- p
    res_all <- rbind(res_all,res_df)
    pairing_gamete[[p]] <- res
  }

  library(dplyr);library(reshape2)
  #creates an empty table
  probability_table <- matrix(0,nrow = length(gametes),ncol = length(gametes))
  rownames(probability_table) <- colnames(probability_table) <- gametes

  format1 <- lapply(pairing_method,function(p){
    default_proba <- probability_table
    tmp <- dataframe_to_contigency(res_all[res_all$pairing %in% p,],'locP','locQ','probability')
    default_proba[rownames(tmp),colnames(tmp)] <- tmp

    default_count <- probability_table
    tmp <- dataframe_to_contigency(res_all[res_all$pairing %in% p,],'locP','locQ','count')
    default_count[rownames(tmp),colnames(tmp)] <- tmp

    list('probability' = default_proba,
         'recombination' = default_count)
  })
  names(format1) <- pairing_method

  format2 <- lapply(pairing_method,function(p){
    each_inheritance <- format1[[p]]
    format2 <- cbind(melt(each_inheritance$recombination),melt(each_inheritance$probability)[,3])
    colnames(format2) <- c('markerP','markerQ','recombination','probability')
    format2
  })
  names(format2) <- pairing_method

  input <- list('ploidy'=ploidy)

  return(list('format1' = format1,
              'format2' = format2,
              'input' = input))
}

######Infer offspring gametes probability######
#' offspring_probability_estimator
#' generate the probability for offspring's genotype under specified inheritance
#'
#' @param F1 dataframe of offspring's multi-allelic marker genotype
#' @param P1 parent1's genotype
#' @param P2 parent2's genotype
#' @param ploidy crop ploidy level
#' @param P1_inheritance inheritance of P1, can be either 'Polysomic' or 'Disomic'
#' @param P2_inheritance inheritance of P2, can be either 'Polysomic' or 'Disomic'
#'
#' @return a dataframe include each offspring's genotype probability under different mode of inheritance
#' @export
#'
#' @examples
offspring_probability_estimator <- function(F1 = offspring_score,
                                            P1,
                                            P2,
                                            P1_inheritance,
                                            P2_inheritance,
                                            ploidy = 4){


  F1_ind <- rownames(F1)
  F1 <- apply(F1,2,as.character)
  total_list <- format_pairint_list(temp = do.call(rbind,lapply(1:nrow(F1),function(i){
    found_probability (P1,
                       P2,
                       inheritance1 = P1_inheritance,
                       inheritance2 = P2_inheritance,
                       row = F1[i,],
                       ploidy = ploidy,
                       ind_nme = F1_ind[i])
  })),
  F1_ind = F1_ind)
  return(total_list)
}

#' Infer offspring inheritance
#'
#' @param poly formatted list include all individuals and parents genotype. It is a list of all marker's names
#' Within each marker, it includes:
#' + parental info: a dataframe include two parents' genotype
#' + temp: a dataframe include alll individuals' genotype
#' All the information here is not phased.
#' @param ploidy crop ploidy level
#' @param P1_inheritance either 'Polysomic' or 'Disomic'
#' @param P2_inheritance either 'Polysomic' or 'Disomic'
#'
#' @return It includes two lists: input and output,
#' The 'input' include: P1_inheritance, P2_inheritance, ploidy. The output include a list of all markers. Within each marker,
#' it include the missing_nbr(estimated number of missing individuals), and all possible pairing configuration, what is the probability
#' of each gametes. It is a dataframe contains:
#' + indivdual: individual name
#' + P1_1, P1_2: the gametes from P1 (not phased)
#' + P2_1, P2_2: the gametes from P2 (not phased)
#' + probability: the probability obtaining the gametes
#' + pairing: the pairing of two parents
#' + pair_freq: the frequency of this pairing from two parents
#' + P1_gamete_freq, P2_gamete_freq: Under this assumed pairing and inheritance, what is the expected
#' segregation ratio of obtaining this gamete.
#' @export
#'
#' @examples
Infer_offspring_inheritance <- function(poly,
                                        ploidy,
                                        P1_inheritance,
                                        P2_inheritance){
  offspring_score_list <- lapply(names(poly),function(mrk){
    parents <- apply(poly[[mrk]]$parent_info,2,as.character)
    F1 <- poly[[mrk]]$temp; F1 <- F1[-c(1,2),]  #remove parents

    pair <- paste(do.call(paste0, as.data.frame(parents)),collapse = '_')

    F1_infer <- offspring_probability_estimator(F1 = F1,
                                                P1 = strsplit(pair,'_')[[1]][1],
                                                P2 = strsplit(pair,'_')[[1]][2],
                                                P1_inheritance = P1_inheritance,
                                                P2_inheritance = P2_inheritance,
                                                ploidy = ploidy)
    missing_nbr <- sum(!complete.cases(F1))
    list('F1_inference' = F1_infer,
         missing_nbr = missing_nbr)
  })
  names(offspring_score_list) <- names(poly)
  input <- list('P1_inheritance' = P1_inheritance,
                'P2_inheritance' = P2_inheritance,
                'ploidy' = ploidy)
  return(list('output' = offspring_score_list,
              'input' = input))
}


######Polysomic x Disomic######
#' examine which pairing configuration the parent follow disomic inheritance goes
#' This function aims to help the parents who follows disomic inheritance to decide whith pairing
#' it follows exactly. This choice is made by using the F1 segregation ratio.
#' @param MT a dataframe generated from *MarkerType_Overview*. It needs to include the columns:
#' marker (marker name), parents (parents genotype), P1_MT (P1 marker type), P2_MT (P2_marker type),
#' sharing (the number of shared haplotypes between two parents), unique_hap, MT (combined marker type)
#' @param seg_invalidrate the value used in binomial test
#' @param offspring_score_list a list generated from *Infer_offspring_inheritance*. It needs to include two lists: input and output,
#' The 'input' include: P1_inheritance, P2_inheritance, ploidy. The output include a list of all markers. Within each marker,
#' it include the missing_nbr(estimated number of missing individuals), and all possible pairing configuration, what is the probability
#' of each gametes. It is a dataframe contains:
#' + indivdual: individual name
#' + P1_1, P1_2: the gametes from P1 (not phased)
#' + P2_1, P2_2: the gametes from P2 (not phased)
#' + probability: the probability obtaining the gametes
#' + pairing: the pairing of two parents
#' + pair_freq: the frequency of this pairing from two parents
#' + P1_gamete_freq, P2_gamete_freq: Under this assumed pairing and inheritance, what is the expected
#' segregation ratio of obtaining this gamete.
#' @return a list include
#' + output: a data frame include the marker name and the pairing each marker follows of the target parent
#' + target_parent: which parent follows disomic inheritance
#' + inheritance scores: to determine the pairing configuration, mode of inheritance study was performed,
#' what is the inheritance scores and parameters in it. This is a combined list include each marker, under each
#' pairing configuration, their: chi-square value, binomial_p, multiplied_p, combined_p, observed_count,
#' expected count, and invalid percentage.
#' @export
#'
#' @examples
pairing_examinator <- function(MT = MT,
                               seg_invalidrate = 0.05,
                               offspring_score_list = offspring_score_list){

  P1_inheritance <- offspring_score_list$input$P1_inheritance
  P2_inheritance <- offspring_score_list$input$P2_inheritance
  ploidy <- offspring_score_list$input$ploidy
  target_parent <- c('P1','P2')[which(c(P1_inheritance,P2_inheritance) %in% 'Disomic')]
  if(length(target_parent) > 1){
    stop('Both parents follow disomic inheritance and there is no need to perform this step')
  }else if(length(target_parent) == 0){
    stop('Both parents follow polysomic inheritance and there is no need to perform this step')
  }else{
    writeLines(paste0(target_parent,' follows disomic inheritance. Here it decide which
                    pairing it follows for each marker'))

    #estimate inheritance score
    inheritance <- lapply(names(offspring_score_list$output),function(mrk){
      F1 <- offspring_score_list$output[[mrk]]$F1_inference
      missing_nbr <- offspring_score_list$output[[mrk]]$missing_nbr

      #check segregation and obtain the P value
      pair <- MT[MT$marker %in% mrk,]$parents
      target_parent_genotype <- strsplit(pair,'_')[[1]][as.numeric(gsub('P','',target_parent))]
      other_parent_genotype <- strsplit(pair,'_')[[1]][setdiff(c(1,2),
                                                               as.numeric(gsub('P','',target_parent)))]

      target_parent_gametes <- create_gamete_disomic(parent = target_parent_genotype,
                                                     ploidy = ploidy)
      other_parent_gametes <- create_gamete_polysomic(parent = other_parent_genotype,
                                                      ploidy = ploidy)


      inheritance_scores <- check_segregation(possible_gamete = target_parent_gametes,
                                              parent = target_parent,
                                              probability_offspring = F1,
                                              seg_invalidrate,
                                              missing_nbr = missing_nbr)
      inheritance_scores
    })
    names(inheritance) <- names(offspring_score_list$output)

    #find the solution
    pairing_chosen <- do.call(rbind,lapply(names(inheritance),function(mrk){
      combined_scores <- sapply(inheritance[[mrk]],function(m){
        m$combined_p
      })
      c(mrk,names(which.max(combined_scores)))
    }))
    colnames(pairing_chosen) <- c('marker','pairing')

    return(list('output' = pairing_chosen,
                'target_parent' = target_parent,
                'inheritance_scores' = inheritance))
  }

}

#' based on tested pairing to choose to offspring score
#'
#' @param offspring_score_list It includes two lists: input and output,
#' The 'input' include: P1_inheritance, P2_inheritance, ploidy. The output include a list of all markers. Within each marker,
#' it include the missing_nbr(estimated number of missing individuals), and all possible pairing configuration, what is the probability
#' of each gametes. It is a dataframe contains:
#' + indivdual: individual name
#' + P1_1, P1_2: the gametes from P1 (not phased)
#' + P2_1, P2_2: the gametes from P2 (not phased)
#' + probability: the probability obtaining the gametes
#' + pairing: the pairing of two parents
#' + pair_freq: the frequency of this pairing from two parents
#' + P1_gamete_freq, P2_gamete_freq: Under this assumed pairing and inheritance, what is the expected
#' segregation ratio of obtaining this gamete.
#' @param pairing a list include
#' + output: a data frame include the marker name and the pairing each marker follows of the target parent
#' + target_parent: which parent follows disomic inheritance
#' + inheritance scores: to determine the pairing configuration, mode of inheritance study was performed,
#' what is the inheritance scores and parameters in it. This is a combined list include each marker, under each
#' pairing configuration, their: chi-square value, binomial_p, multiplied_p, combined_p, observed_count,
#' expected count, and invalid percentage.
#'
#' @return It includes two lists: input and output,
#' The 'input' include: P1_inheritance, P2_inheritance, ploidy. The output include a list of all markers. Within each marker,
#' it include the missing_nbr(estimated number of missing individuals), and all possible pairing configuration, what is the probability
#' of each gametes. It is a dataframe contains:
#' + indivdual: individual name
#' + P1_1, P1_2: the gametes from P1 (not phased)
#' + P2_1, P2_2: the gametes from P2 (not phased)
#' + probability: the probability obtaining the gametes
#' + pairing: the pairing of two parents
#' + pair_freq: the frequency of this pairing from two parents
#' + P1_gamete_freq, P2_gamete_freq: Under this assumed pairing and inheritance, what is the expected
#' segregation ratio of obtaining this gamete.
#' Please note here only the chosen pairing dataframe is kept
#' @export
#'
#' @examples
choose_offspringscore <- function(offspring_score_list,
                                  pairing){
  suggested_pairing <- as.data.frame(pairing$output)

  offspring_score_list1 <- lapply(names(offspring_score_list$output),function(mrk){
    tmp <- offspring_score_list$output[[mrk]]$F1_inference
    chosen_pairing <- suggested_pairing[suggested_pairing$marker %in% mrk,]$pairing

    output1 <- list( tmp[[chosen_pairing]]); names(output1) <- chosen_pairing
    list('F1_inference' = output1,
         'missing_nbr' = offspring_score_list$output[[mrk]]$missing_nbr)
  })

  names(offspring_score_list1) <- names(offspring_score_list$output)

  input <- list('P1_inheritance' = offspring_score_list$input$P1_inheritance,
                'P2_inheritance' = offspring_score_list$input$P2_inheritance,
                'ploidy' = offspring_score_list$input$ploidy)

  return(list('input' = input,
              'output' = offspring_score_list1))
}






######Linkage estimation#######
#' perform linkage anlaysis in chosen marker type
#'
#' @param Mrk_P one marker, it need to be in the format of "AABB_EEFF"
#' @param Mrk_Q another marker, it need to be in the format of "AABB_EEFF"
#' @param offspring_score_list It includes two lists: input and output,
#' The 'input' include: P1_inheritance, P2_inheritance, ploidy. The output include a list of all markers. Within each marker,
#' it include the missing_nbr(estimated number of missing individuals), and all possible pairing configuration, what is the probability
#' of each gametes. It is a dataframe contains:
#' + indivdual: individual name
#' + P1_1, P1_2: the gametes from P1 (not phased)
#' + P2_1, P2_2: the gametes from P2 (not phased)
#' + probability: the probability obtaining the gametes
#' + pairing: the pairing of two parents
#' + pair_freq: the frequency of this pairing from two parents
#' + P1_gamete_freq, P2_gamete_freq: Under this assumed pairing and inheritance, what is the expected
#' segregation ratio of obtaining this gamete.
#' @param poly formatted list include all individuals and parents genotype. It is a list of all marker's names
#' Within each marker, it includes:
#' + parental info: a dataframe include two parents' genotype
#' + temp: a dataframe include alll individuals' genotype
#' All the information here is not phased.
#' @param cores number of cores used to run the analysis. It is suggested to use less than your computer capacity. To
#' check how many cores you computer has, you can use *detectCores()*
#' @param all_equations a list include:
#'  - equations: equations for P1 and P2 separately
#'  - input: P1 and P2's inheritance, ploidy level, and marker type (MT)
#' @param pairing a list include
#' + output: a data frame include the marker name and the pairing each marker follows of the target parent
#' + target_parent: which parent follows disomic inheritance
#' + inheritance scores: to determine the pairing configuration, mode of inheritance study was performed,
#' what is the inheritance scores and parameters in it. This is a combined list include each marker, under each
#' pairing configuration, their: chi-square value, binomial_p, multiplied_p, combined_p, observed_count,
#' expected count, and invalid percentage.
#'
#' @return a list include input and output. In the input, it has: ploidy, Mrk_P, and Mrk_Q.
#' In the output, it includes the linkage estimation result
#' @export
#'
#' @examples
# Haplotype_linkage <- function(Mrk_P = combs[c,][1],
#                               Mrk_Q = combs[c,][2],
#                               offspring_score_list = offspring_score_list1,
#                               poly,
#                               cores,
#                               all_equations,
#                               pairing){
#   time1 <- Sys.time()
#   ###Step 1:make the tempset for chosen marker type###
#   ploidy <- offspring_score_list$input$ploidy
#   share_haplotypes_P <- sharing_identification(Mrk_P)* ploidy
#   share_haplotypes_Q <- sharing_identification(Mrk_Q)* ploidy
#   #determine the share haplotypes is useful because it will help to choose which method
#   #used in linkage estimation
#   share_nbr <- ifelse(any(c(share_haplotypes_P,share_haplotypes_Q) > 0),TRUE,FALSE) #True - corrected p MML
#
#   M_P <- paste0(paste0(sapply(strsplit(Mrk_P,'_')[[1]],decide_MT),collapse = ' x '),
#                 '_Shar',share_haplotypes_P)
#   M_Q <- paste0(paste0(sapply(strsplit(Mrk_Q,'_')[[1]],decide_MT),collapse = ' x '),
#                 '_Shar',share_haplotypes_Q)
#
#   if(M_P == M_Q){
#     marker_P <- marker_Q <- MT[MT$MT %in% M_P,]$marker
#     tempset <- t(combn(marker_P,2))
#   }else{
#     marker_P <- MT[MT$MT %in% M_P,]$marker
#     marker_Q <- MT[MT$MT %in% M_Q,]$marker
#     tempset <- apply(expand.grid(marker_P,marker_Q),2,as.character)
#     tempset <- tempset[tempset[,1] != tempset[,2],]  #remove the ones pair with itself
#   }
#   colnames(tempset) <- c('marker_P','marker_Q')
#
#
#
#   # combined_equation <-  kronecker_string(matrix1 = dataframe_to_contigency(Q_P1_equation$ABCC,'P','Q','probability'),
#   #                                        matrix2 = dataframe_to_contigency(Q_P2_equation$ABEF,'P','Q','probability'),
#   #                                        symbol = '*')
#
#   #run parallel function set up
#   library(foreach)
#   library(doParallel)
#   win <- Sys.info()["sysname"] == "Windows"
#   if (win) {
#     cl <- parallel::makeCluster(cores)
#     doSNOW::registerDoSNOW(cl)
#   } else {
#     doParallel::registerDoParallel(cores = cores)
#   }
#
#   ######Step 3: perform linkage estimation
#   individuals <- rownames(poly[[1]]$temp)[-c(1,2)]
#
#   Linkage_res <- foreach::foreach(i = 1:200, .combine = rbind, .inorder = F,
#                                   .export = c('recombination_estimator',
#                                               'maximum_likelihood_estimator')) %dopar% {
#
#                                                 marker_P <- as.character(tempset[i,'marker_P']); marker_Q <- as.character(tempset[i,'marker_Q'])
#
#
#                                                 ######Step 2: define genotype
#                                                 P_parents <- apply(poly[[marker_P]]$parent_info,2,as.character)
#                                                 Q_parents <- apply(poly[[marker_Q]]$parent_info,2,as.character)
#                                                 rownames(P_parents) <- rownames(Q_parents) <- c('P1','P2')
#                                                 phasing_P <- paste(do.call(paste0, as.data.frame(P_parents)),collapse = '_')
#
#                                                 ######Step 3: find the corresponding equation
#                                                 all_equation_P1 <- all_equations$equations$P1; all_equation_P2 <- all_equations$equations$P2
#                                                 Q_P1_equation <- all_equation_P1[[paste0(paste0(P_parents[1,],collapse = ''),'_',paste0(Q_parents[1,],collapse = ''))]]
#                                                 Q_P2_equation <- all_equation_P2[[paste0(paste0(P_parents[2,],collapse = ''),'_',paste0(Q_parents[2,],collapse = ''))]]
#
#
#                                                 F1_infer_P <- offspring_score_list$output[[marker_P]]$F1_inference[[1]]
#                                                 F1_infer_Q <- offspring_score_list$output[[marker_Q]]$F1_inference[[1]]
#
#                                                 ##for each offspring, we calculate the gamete formation table
#                                                 offspring_probability_table <- do.call(rbind,lapply(individuals, function(ind){
#                                                   F1_infer_P <- F1_infer_P[complete.cases(F1_infer_P),]
#                                                   F1_infer_Q <- F1_infer_Q[complete.cases(F1_infer_Q),]
#                                                   if((ind %in% F1_infer_P$individual) & ind %in% (F1_infer_Q$individual)){
#                                                     mrk_P <- F1_infer_P[F1_infer_P$individual %in% ind,c(paste0('P1_',seq(1,ploidy/2)),
#                                                                                                          paste0('P2_',seq(1,ploidy/2)),
#                                                                                                          'probability'),drop = FALSE]
#                                                     mrk_Q <- F1_infer_Q[F1_infer_Q$individual %in% ind,c(paste0('P1_',seq(1,ploidy/2)),
#                                                                                                          paste0('P2_',seq(1,ploidy/2)),
#                                                                                                          'probability'),drop = FALSE]
#                                                     P1_P <- apply(apply(mrk_P[,paste0('P1_',seq(1,ploidy/2)),drop = FALSE],1,sort),2,paste,collapse = '')
#                                                     P2_P <- apply(apply(mrk_P[,paste0('P2_',seq(1,ploidy/2)),drop = FALSE],1,sort),2,paste,collapse = '')
#                                                     P_tmp <- data.frame(cbind(P1_P,P2_P,'probability' = mrk_P$probability),stringsAsFactors = FALSE)
#                                                     P_tmp$probability <- as.numeric(P_tmp$probability)
#
#                                                     P1_Q <- apply(apply(mrk_Q[,paste0('P1_',seq(1,ploidy/2)),drop = FALSE],1,sort),2,paste,collapse = '')
#                                                     P2_Q <- apply(apply(mrk_Q[,paste0('P2_',seq(1,ploidy/2)),drop = FALSE],1,sort),2,paste,collapse = '')
#                                                     Q_tmp <- data.frame(cbind(P1_Q,P2_Q,'probability' = mrk_Q$probability),stringsAsFactors = FALSE)
#                                                     Q_tmp$probability <- as.numeric(Q_tmp$probability)
#
#                                                     tmp <- do.call(rbind,lapply(1:nrow(P_tmp),function(i){
#                                                       do.call(rbind,lapply(1:nrow(Q_tmp),function(j){cbind(P_tmp[i,c(1,2)],Q_tmp[j,c(1,2)],P_tmp[i,3] * Q_tmp[j,3])}))}))
#                                                     colnames(tmp)[5] <- 'count'
#                                                     tmp$individual <- ind
#                                                     tmp
#                                                   }
#                                                 }))
#
#
#                                                 # #If we want to calculate rf from both parents together
#                                                 # #create gamete count table of offspring
#                                                 # P_gamete <- gamete_combinations_per_marker(P_parents)
#                                                 # Q_gamete <- gamete_combinations_per_marker(Q_parents)
#                                                 # P1_P2_count <- matrix(nrow = length(P_gamete),ncol = length(Q_gamete),0)
#                                                 # rownames(P1_P2_count) <- P_gamete; colnames(P1_P2_count) <- Q_gamete
#                                                 # for(o in 1:nrow(offspring_probability_table)){
#                                                 #   offspring_pro <- offspring_probability_table[o,]
#                                                 #   P1_P2_count[paste0(offspring_pro$P1_P,':',offspring_pro$P2_P),
#                                                 #               paste0(offspring_pro$P1_Q,':',offspring_pro$P2_Q)] <-
#                                                 #     P1_P2_count[paste0(offspring_pro$P1_P,':',offspring_pro$P2_P),
#                                                 #                 paste0(offspring_pro$P1_Q,':',offspring_pro$P2_Q)] + offspring_pro$probability
#                                                 # }
#
#                                                 #count the R
#                                                 R_P1 <- recombination_estimator(equation = Q_P1_equation,
#                                                                                 tmp = offspring_probability_table,
#                                                                                 share_nbr = share_nbr,
#                                                                                 pairing = pairing,
#                                                                                 mrk = marker_Q,
#                                                                                 parent = 'P1')
#                                                 R_P2 <- recombination_estimator(equation = Q_P2_equation,
#                                                                                 tmp = offspring_probability_table,
#                                                                                 share_nbr = share_nbr,
#                                                                                 pairing = pairing,
#                                                                                 mrk = marker_Q,
#                                                                                 parent = 'P2')
#
#
#
#                                                 P1_chosen <- R_P1[R_P1$rf %in% min(R_P1$rf),]
#                                                 P2_chosen <- R_P2[R_P2$rf %in% min(R_P2$rf),]
#
#                                                 if(nrow(P1_chosen) > 1){
#                                                   P1_chosen <- P1_chosen[1,]
#                                                 }
#                                                 if(nrow(P2_chosen) > 1){
#                                                   P2_chosen <- P2_chosen[2,]
#                                                 }
#
#                                                 Total <- length(individuals)*ploidy
#                                                 R_count <- (P1_chosen$rf+P2_chosen$rf) * Total/2
#                                                 NR_count <- Total - R_count
#                                                 r <- R_count/Total
#                                                 LOD <- round(log10((2^Total) * ((1 -r)^NR_count) * ((r)^R_count)),2)
#
#                                                 if(!is.null(pairing)){
#                                                   suggested_pairing <- as.data.frame(pairing$output)
#                                                   target_parent <- pairing$target_parent
#                                                   phasing_P <- as.character(suggested_pairing[suggested_pairing$marker %in% marker_P,]$pairing)
#                                                   # chosen_pairing <- strsplit(chosen_pairing,'_')[[1]][as.numeric(gsub('P','',target_parent))]
#                                                 }else{
#                                                   phasing_P <- paste0(apply(P_parents,1,paste,collapse = ''),collapse = '_')
#                                                 }
#
#                                                 phasing_Q <-  paste0(P1_chosen$nme,'_',P2_chosen$nme)
#
#                                                 data.frame(marker_P,marker_Q,r,LOD,phasing_P,phasing_Q)
#                                               }
#
#   input <- list('ploidy' = ploidy,
#                 'marker_P' = Mrk_P,
#                 'marker_Q' = Mrk_Q)
#
#   time2 <- Sys.time()
#   writeLines(paste0('It took totally ',time2 - time1,' to obtain the results'))
#
#   return(list('input' = input,
#               'output' = Linkage_res))
# }
Haplotype_linkage <- function(Mrk_P = combs[c,][1],
                              Mrk_Q = combs[c,][2],
                              offspring_score_list = offspring_score_list1,
                              poly,
                              MT,
                              cores,
                              all_equations,
                              pairing){
  time1 <- Sys.time()
  ###Step 1:make the tempset for chosen marker type###
  ploidy <- offspring_score_list$input$ploidy
  share_haplotypes_P <- sharing_identification(Mrk_P)* ploidy
  share_haplotypes_Q <- sharing_identification(Mrk_Q)* ploidy
  #determine the share haplotypes is useful because it will help to choose which method
  #used in linkage estimation
  share_nbr <- ifelse(any(c(share_haplotypes_P,share_haplotypes_Q) > 0),TRUE,FALSE) #True - corrected p MML

  M_P <- paste0(paste0(sapply(strsplit(Mrk_P,'_')[[1]],decide_MT),collapse = ' x '),
                '_Shar',share_haplotypes_P)
  M_Q <- paste0(paste0(sapply(strsplit(Mrk_Q,'_')[[1]],decide_MT),collapse = ' x '),
                '_Shar',share_haplotypes_Q)

  if(M_P == M_Q){
    marker_P <- marker_Q <- MT[MT$MT %in% M_P,]$marker
    tempset <- t(combn(marker_P,2))
  }else{
    marker_P <- MT[MT$MT %in% M_P,]$marker
    marker_Q <- MT[MT$MT %in% M_Q,]$marker
    tempset <- apply(expand.grid(marker_P,marker_Q),2,as.character)
    tempset <- tempset[tempset[,1] != tempset[,2],]  #remove the ones pair with itself
  }
  colnames(tempset) <- c('marker_P','marker_Q')



  # combined_equation <-  kronecker_string(matrix1 = dataframe_to_contigency(Q_P1_equation$ABCC,'P','Q','probability'),
  #                                        matrix2 = dataframe_to_contigency(Q_P2_equation$ABEF,'P','Q','probability'),
  #                                        symbol = '*')

  #run parallel function set up
  library(foreach)
  library(doParallel)
  win <- Sys.info()["sysname"] == "Windows"
  if (win) {
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
  } else {
    doParallel::registerDoParallel(cores = cores)
  }

  ######Step 3: perform linkage estimation
  individuals <- rownames(poly[[1]]$temp)[-c(1,2)]

  Linkage_res <- foreach::foreach(i = 1:nrow(tempset), .combine = rbind, .inorder = F,
                                  .export = c('recombination_estimator',
                                              'maximum_likelihood_estimator',
                                              'kronecker_string',
                                              'dataframe_to_contigency')) %dopar% {

                                                marker_P <- as.character(tempset[i,'marker_P']); marker_Q <- as.character(tempset[i,'marker_Q'])


                                                ######Step 2: define genotype
                                                P_parents <- apply(poly[[marker_P]]$parent_info,2,as.character)
                                                Q_parents <- apply(poly[[marker_Q]]$parent_info,2,as.character)
                                                rownames(P_parents) <- rownames(Q_parents) <- c('P1','P2')
                                                phasing_P <- paste(do.call(paste0, as.data.frame(P_parents)),collapse = '_')

                                                ######Step 3: find the corresponding equation
                                                all_equation_P1 <- all_equations$equations$P1; all_equation_P2 <- all_equations$equations$P2
                                                phasing_P1 <- paste0(paste0(P_parents[1,],collapse = ''),'_',paste0(Q_parents[1,],collapse = ''))
                                                phasing_P1_r <- paste0(paste0(Q_parents[1,],collapse = ''),'_',paste0(P_parents[1,],collapse = ''))
                                                phasing_P2 <- paste0(paste0(P_parents[2,],collapse = ''),'_',paste0(Q_parents[2,],collapse = ''))
                                                phasing_P2_r <- paste0(paste0(Q_parents[2,],collapse = ''),'_',paste0(P_parents[2,],collapse = ''))

                                                P1_eq_chosen <- which(c(phasing_P1,phasing_P1_r) %in% names(all_equation_P1))[1]
                                                P2_eq_chosen <- which(c(phasing_P2,phasing_P2_r) %in% names(all_equation_P2))[1]

                                                Q_P1_equation <- all_equation_P1[[c(phasing_P1,phasing_P1_r)[P1_eq_chosen]]]
                                                Q_P2_equation <- all_equation_P2[[c(phasing_P2,phasing_P2_r)[P2_eq_chosen]]]


                                                F1_infer_P <- offspring_score_list$output[[marker_P]]$F1_inference[[1]]
                                                F1_infer_Q <- offspring_score_list$output[[marker_Q]]$F1_inference[[1]]

                                                ##for each offspring, we calculate the gamete formation table
                                                offspring_probability_table <- do.call(rbind,lapply(individuals, function(ind){
                                                  F1_infer_P <- F1_infer_P[complete.cases(F1_infer_P),]
                                                  F1_infer_Q <- F1_infer_Q[complete.cases(F1_infer_Q),]
                                                  if((ind %in% F1_infer_P$individual) & ind %in% (F1_infer_Q$individual)){
                                                    mrk_P <- F1_infer_P[F1_infer_P$individual %in% ind,c(paste0('P1_',seq(1,ploidy/2)),
                                                                                                         paste0('P2_',seq(1,ploidy/2)),
                                                                                                         'probability'),drop = FALSE]
                                                    mrk_Q <- F1_infer_Q[F1_infer_Q$individual %in% ind,c(paste0('P1_',seq(1,ploidy/2)),
                                                                                                         paste0('P2_',seq(1,ploidy/2)),
                                                                                                         'probability'),drop = FALSE]
                                                    P1_P <- apply(apply(mrk_P[,paste0('P1_',seq(1,ploidy/2)),drop = FALSE],1,sort),2,paste,collapse = '')
                                                    P2_P <- apply(apply(mrk_P[,paste0('P2_',seq(1,ploidy/2)),drop = FALSE],1,sort),2,paste,collapse = '')
                                                    P_tmp <- data.frame(cbind(P1_P,P2_P,'probability' = mrk_P$probability),stringsAsFactors = FALSE)
                                                    P_tmp$probability <- as.numeric(P_tmp$probability)

                                                    P1_Q <- apply(apply(mrk_Q[,paste0('P1_',seq(1,ploidy/2)),drop = FALSE],1,sort),2,paste,collapse = '')
                                                    P2_Q <- apply(apply(mrk_Q[,paste0('P2_',seq(1,ploidy/2)),drop = FALSE],1,sort),2,paste,collapse = '')
                                                    Q_tmp <- data.frame(cbind(P1_Q,P2_Q,'probability' = mrk_Q$probability),stringsAsFactors = FALSE)
                                                    Q_tmp$probability <- as.numeric(Q_tmp$probability)

                                                    tmp <- do.call(rbind,lapply(1:nrow(P_tmp),function(i){
                                                      do.call(rbind,lapply(1:nrow(Q_tmp),function(j){cbind(P_tmp[i,c(1,2)],Q_tmp[j,c(1,2)],P_tmp[i,3] * Q_tmp[j,3])}))}))
                                                    colnames(tmp)[5] <- 'count'
                                                    tmp$individual <- ind
                                                    tmp
                                                  }
                                                }))

                                                #if need to combine two parents together, what is the possible phasing
                                                P_common <- sum(P_parents[1,] %in% P_parents[,2]); Q_common <- sum(Q_parents[1,] %in% Q_parents[,2])

                                                if(any(c(P_common, Q_common) %in% ploidy)){# if it is all sharing, then follow this method
                                                  phasing_Q_possibilities <- apply(expand.grid(names(Q_P1_equation),names(Q_P2_equation)),2,as.character)


                                                  res <- do.call(rbind,lapply(1:nrow(phasing_Q_possibilities),function(pha){
                                                    chosen_phasing <- paste0(c(phasing_Q_possibilities[pha,]),collapse = 'x')
                                                    combined_equation <-  kronecker_string(matrix1 = dataframe_to_contigency(Q_P1_equation[[phasing_Q_possibilities[pha,1]]],'P','Q','probability'),
                                                                                           matrix2 = dataframe_to_contigency(Q_P2_equation[[phasing_Q_possibilities[pha,2]]],'P','Q','probability'),
                                                                                           symbol = '*')
                                                    equation <- combined_equation$format1
                                                    tmp <- offspring_probability_table

                                                    rownames(equation) <- paste0(equation$P,'|',equation$Q)

                                                    all_individuals <- sapply(individuals,function(ind){
                                                      df <- tmp[tmp$individual %in% ind,]
                                                      df$probability <- equation[paste0(df$P1_P,':',df$P2_P,'|',
                                                                                        df$P1_Q,':',df$P2_Q),'simplified_probability']

                                                      paste0('(',gsub(' ','',Deriv::Simplify(paste0(paste0(df$count, '*', df$probability),collapse = '+'))),
                                                             ')')
                                                    })

                                                    table_all_individuals <- table(all_individuals)
                                                    multi_value <- round(length(individuals)/length(table_all_individuals)) * 2
                                                    all_individuals1 <- sapply(1:length(table_all_individuals),function(ins){
                                                      paste0('(',names(table_all_individuals)[ins],'^',table_all_individuals[ins],'*(10^',multi_value,'))')
                                                    })
                                                    all_individuals2 <-  paste0(all_individuals1,collapse = '*')

                                                    probability_maximumlikelihood_approach <- maximum_likelihood_estimator(equation = all_individuals2,plot = TRUE)
                                                    rf <- round(probability_maximumlikelihood_approach$rf,3)


                                                    data.frame(chosen_phasing,rf)
                                                  }))
                                                  res_chosen <- res[res$rf %in% min(res$rf),]

                                                  if(nrow(res_chosen) > 1){
                                                    res_chosen <- res_chosen[1,]
                                                  }

                                                  Total <- length(individuals)*ploidy
                                                  R_count <- res_chosen$rf * Total
                                                  NR_count <- Total - R_count
                                                  r <- R_count/Total
                                                  LOD <- round(log10((2^Total) * ((1 -r)^NR_count) * ((r)^R_count)),2)
                                                  phasing_Q <- res_chosen$chosen_phasing
                                                  phasing_Q <- gsub('x','_',phasing_Q)

                                                }else{
                                                  ###Separating parents method
                                                  R_P1 <- recombination_estimator(equation = Q_P1_equation,
                                                                                  tmp = offspring_probability_table,
                                                                                  share_nbr = share_nbr,
                                                                                  pairing = pairing,
                                                                                  mrk_P = marker_P,
                                                                                  mrk_Q = marker_Q,
                                                                                  parent = 'P1',
                                                                                  eq_chosen = P1_eq_chosen)
                                                  R_P2 <- recombination_estimator(equation = Q_P2_equation,
                                                                                  tmp = offspring_probability_table,
                                                                                  share_nbr = share_nbr,
                                                                                  pairing = pairing,
                                                                                  mrk_P = marker_P,
                                                                                  mrk_Q = marker_Q,
                                                                                  parent = 'P2',
                                                                                  eq_chosen = P2_eq_chosen)
                                                  P1_chosen <- R_P1[R_P1$rf %in% min(R_P1$rf),]
                                                  P2_chosen <- R_P2[R_P2$rf %in% min(R_P2$rf),]

                                                  if(nrow(P1_chosen) > 1){
                                                    P1_chosen <- P1_chosen[1,]
                                                  }
                                                  if(nrow(P2_chosen) > 1){
                                                    P2_chosen <- P2_chosen[2,]
                                                  }

                                                  Total <- length(individuals)*ploidy
                                                  R_count <- (P1_chosen$rf+P2_chosen$rf) * Total/2
                                                  NR_count <- Total - R_count
                                                  r <- R_count/Total
                                                  LOD <- round(log10((2^Total) * ((1 -r)^NR_count) * ((r)^R_count)),2)

                                                  if(!is.null(pairing)){
                                                    suggested_pairing <- as.data.frame(pairing$output)
                                                    target_parent <- pairing$target_parent
                                                    phasing_P <- as.character(suggested_pairing[suggested_pairing$marker %in% marker_P,]$pairing)
                                                  }else{
                                                    phasing_P <- paste0(apply(P_parents,1,paste,collapse = ''),collapse = '_')
                                                  }

                                                  phasing_Q <-  paste0(P1_chosen$phasing,'_',P2_chosen$phasing)


                                                }


                                                data.frame(marker_P,marker_Q,r,LOD,phasing_P,phasing_Q)
                                              }

  input <- list('ploidy' = ploidy,
                'marker_P' = Mrk_P,
                'marker_Q' = Mrk_Q)

  time2 <- Sys.time()
  writeLines(paste0('It took totally ',time2 - time1,' to obtain the results'))

  return(list('input' = input,
              'output' = Linkage_res))
}


#' estimate recombination frequency using MML per marker pair
#'
#' @param equation the corresponding equation of the chosen marker
#' @param tmp corresponding offspring probability table
#' @param share_nbr True/False. When it is True, it share more than one haplotype
#' @param pairing a list include
#' + output: a data frame include the marker name and the pairing each marker follows of the target parent
#' + target_parent: which parent follows disomic inheritance
#' + inheritance scores: to determine the pairing configuration, mode of inheritance study was performed,
#' what is the inheritance scores and parameters in it. This is a combined list include each marker, under each
#' pairing configuration, their: chi-square value, binomial_p, multiplied_p, combined_p, observed_count,
#' expected count, and invalid percentage.
#' @param mrk the marker name of marker Q at this comparison
#' @param parent the target parent, 'P1' or 'P2'
#'
#' @return a dataframe include all phasing situations' estimated rf
#' @export
#'
#' @examples
# recombination_estimator <- function(equation = Q_P2_equation,
#                                     tmp = offspring_probability_table,
#                                     share_nbr,
#                                     pairing,
#                                     mrk = marker_Q,
#                                     parent = 'P2'){
#   recombintion_rf <- do.call(rbind,lapply(names(equation),function(nme){
#     tmp1 <- equation[[nme]]
#     if(!is.data.frame(tmp1)){#if it is a list then we need to choose the most likely pairing configuration
#       suggested_pairing <- as.data.frame(pairing$output)
#       target_parent <- pairing$target_parent
#
#       chosen_pairing <- as.character(suggested_pairing[suggested_pairing$marker %in% mrk,]$pairing)
#       chosen_pairing <- strsplit(chosen_pairing,'_')[[1]][as.numeric(gsub('P','',target_parent))]
#       chosen_pairing <- paste0(sort(strsplit(chosen_pairing,'\\|')[[1]]),collapse = '|')
#
#       part1 <- apply(apply(do.call(rbind,strsplit(unlist(lapply(strsplit(names(tmp1),'\\|'), `[[`, 1)),'')),1,sort),2,paste,collapse = '')
#       part2 <- apply(apply(do.call(rbind,strsplit(unlist(lapply(strsplit(names(tmp1),'\\|'), `[[`, 2)),'')),1,sort),2,paste,collapse = '')
#
#       sorted_names <- apply(apply(data.frame(part1,part2),1,sort),2,paste,collapse = '|')
#
#       tmp1 <- tmp1[[which(sorted_names %in% chosen_pairing)[1]]]
#     }
#
#     if(share_nbr){
#       ###Solution 2: maximum likelihood function:use probability from individual to correct per individual###
#       rownames(tmp1) <- paste0(tmp1$P,'_',tmp1$Q)
#       individuals <- unique(tmp$individual)
#       all_individuals <- sapply(individuals,function(ind){
#         df <- tmp[tmp$individual %in% ind,c(paste0(parent,c('_P','_Q')),'count')]
#         df$probability <- tmp1[paste0(df[[paste0(parent,'_P')]],'_',df[[paste0(parent,'_Q')]]),'gsub_probability']
#         paste0('(',gsub(' ','',Deriv::Simplify(paste0(paste0(df$count, '*', df$probability),collapse = '+'))),
#                ')')
#       })
#       table_all_individuals <- table(all_individuals)
#
#       all_individuals <- sapply(1:length(table_all_individuals),function(ins){
#         paste0('(',names(table_all_individuals)[ins],'^',table_all_individuals[ins],')')
#       })
#
#       all_individuals1 <- paste0(all_individuals,collapse = '*')
#       probability_maximumlikelihood_approach <- maximum_likelihood_estimator(all_individuals1,plot = FALSE)
#       rf <- round(probability_maximumlikelihood_approach$rf,2)
#     }else{
#       ###Solution 1: maximum likelihood function###
#       tmp1$offspring <- sapply(1:nrow(tmp1),function(g){
#         sum(offspring_probability_table[offspring_probability_table[[paste0(parent,'_P')]] %in% as.character(tmp1[g,'P']) &
#                                           offspring_probability_table[[paste0(parent,'_Q')]] %in% as.character(tmp1[g,'Q']),'count'])
#       })
#       colnames(tmp1) <- c('P','Q','recombination','probability','gsub_probability','count')
#       non_zero <- tmp1[tmp1$count > 0,]
#
#       library(dplyr)
#       non_zero1 <- non_zero %>%
#         group_by(gsub_probability) %>%
#         summarise_at(vars('count'), sum)
#       non_zero1 <- as.data.frame(non_zero1)
#
#       maximum_equation <- paste0(sapply(1:nrow(non_zero1),function(h){
#         paste0('((',non_zero1[h,'gsub_probability'],')^',non_zero1[h,'count'],')')
#       }),collapse = '*')
#       maximumlikelihood_approach <- maximum_likelihood_estimator(maximum_equation,
#                                                                  plot = FALSE)
#       rf <- round(maximumlikelihood_approach$rf,2)
#     }
#     data.frame(nme,rf)
#   }))
#
#   tmp1 <- equation[[1]]
#   if(!is.data.frame(tmp1)){
#     recombintion_rf$nme <- sapply(names(equation),function(nme){
#       tmp1 <- equation[[nme]]
#       if(is.list(tmp1)){#if it is a list then we need to choose the most likely pairing configuration
#         suggested_pairing <- as.data.frame(pairing$output)
#         target_parent <- pairing$target_parent
#
#         chosen_pairing <- as.character(suggested_pairing[suggested_pairing$marker %in% mrk,]$pairing)
#         chosen_pairing <- strsplit(chosen_pairing,'_')[[1]][as.numeric(gsub('P','',target_parent))]
#         chosen_pairing <- paste0(sort(strsplit(chosen_pairing,'\\|')[[1]]),collapse = '|')
#
#         part1 <- apply(apply(do.call(rbind,strsplit(unlist(lapply(strsplit(names(tmp1),'\\|'), `[[`, 1)),'')),1,sort),2,paste,collapse = '')
#         part2 <- apply(apply(do.call(rbind,strsplit(unlist(lapply(strsplit(names(tmp1),'\\|'), `[[`, 2)),'')),1,sort),2,paste,collapse = '')
#
#         sorted_names <- apply(apply(data.frame(part1,part2),1,sort),2,paste,collapse = '|')
#
#
#         names(tmp1)[[which(sorted_names %in% chosen_pairing)[1]]]
#       }
#     })
#   }
#   return(recombintion_rf)
# }
recombination_estimator <- function(equation = Q_P2_equation,
                                    tmp = offspring_probability_table,
                                    share_nbr,
                                    pairing,
                                    mrk_P = marker_P,
                                    mrk_Q = marker_Q,
                                    parent = 'P2',
                                    eq_chosen){
  rf_MML <- function(equation,
                     tmp,
                     parent,
                     share_nbr = share_nbr,
                     eq_chosen){
    # if(share_nbr){
    ###Solution 2: maximum likelihood function:use probability from individual to correct per individual###
    rownames(equation) <- paste0(equation$P,'_',equation$Q)
    individuals <- unique(tmp$individual)
    all_individuals <- sapply(individuals,function(ind){
      df <- tmp[tmp$individual %in% ind,c(paste0(parent,c('_P','_Q')),'count')]

      if(eq_chosen == 2){
        df$probability <- equation[paste0(df[[paste0(parent,'_Q')]],'_',df[[paste0(parent,'_P')]]),'gsub_probability']
      }else{
        df$probability <- equation[paste0(df[[paste0(parent,'_P')]],'_',df[[paste0(parent,'_Q')]]),'gsub_probability']
      }

      paste0('(',gsub(' ','',Deriv::Simplify(paste0(paste0(df$count, '*', df$probability),collapse = '+'))),
             ')')
    })
    table_all_individuals <- table(all_individuals)
    all_individuals <- sapply(1:length(table_all_individuals),function(ins){
      paste0('(',names(table_all_individuals)[ins],'^',table_all_individuals[ins],')')
    })

    all_individuals1 <- paste0(all_individuals,collapse = '*')
    probability_maximumlikelihood_approach <- maximum_likelihood_estimator(equation = all_individuals1,plot = FALSE)
    rf_MMLp <- round(probability_maximumlikelihood_approach$rf,2)

    # }else{
    ###Solution 1: maximum likelihood function###
    equation$offspring <- sapply(1:nrow(equation),function(g){
      if(eq_chosen == 2){
        sum(tmp[tmp[[paste0(parent,'_Q')]] %in% as.character(equation[g,'P']) &
                  tmp[[paste0(parent,'_P')]] %in% as.character(equation[g,'Q']),'count'])
      }else{
        sum(tmp[tmp[[paste0(parent,'_P')]] %in% as.character(equation[g,'P']) &
                  tmp[[paste0(parent,'_Q')]] %in% as.character(equation[g,'Q']),'count'])
      }


    })
    colnames(equation) <- c('P','Q','recombination','probability','gsub_probability','count')
    non_zero <- equation[equation$count > 0,]

    library(dplyr)
    non_zero1 <- non_zero %>%
      group_by(gsub_probability) %>%
      summarise_at(vars('count'), sum)
    non_zero1 <- as.data.frame(non_zero1)

    maximum_equation <- paste0(sapply(1:nrow(non_zero1),function(h){
      paste0('((',non_zero1[h,'gsub_probability'],')^',non_zero1[h,'count'],')')
    }),collapse = '*')
    maximumlikelihood_approach <- maximum_likelihood_estimator(maximum_equation,
                                                               plot = FALSE)
    rf_MML <- round(maximumlikelihood_approach$rf,2)
    # }



    rf <- min(c(rf_MMLp,rf_MML))[1]
    return(rf)
  }

  #if there is disomic pairing, need to choose the inheritance decided pairing
  if(!is.null(pairing)){
    suggested_pairing <- as.data.frame(pairing$output)
    target_parent <- pairing$target_parent

    if(any(parent == target_parent)){#if there is disomic inheritance
      chosen_Q <- as.character(suggested_pairing[suggested_pairing$marker %in% mrk_Q,]$pairing)
      chosen_P <- as.character(suggested_pairing[suggested_pairing$marker %in% mrk_P,]$pairing)
      Q_c <- strsplit(chosen_Q,'_')[[1]][as.numeric(gsub('P','',parent))]
      P_c <- strsplit(chosen_P,'_')[[1]][as.numeric(gsub('P','',parent))]

      Q_c <- paste0(sort(strsplit(Q_c,'\\|')[[1]]),collapse = '|')
      P_c <- paste0(sort(strsplit(P_c,'\\|')[[1]]),collapse = '|')

      if(eq_chosen == 1){
        chosen_equation <- equation[[P_c]][[Q_c]]
      }else{
        chosen_equation <- equation[[Q_c]][[P_c]]
      }


      recombintion_rf <- rf_MML(equation = chosen_equation,
                                tmp = tmp,
                                parent = parent,
                                share_nbr = share_nbr,
                                eq_chosen = eq_chosen)
      estimated_rf <- data.frame('rf' = recombintion_rf,
                                 'phasing' = P_c)

    }else{
      estimated_rf <- data.frame(sapply(names(equation),function(eq){
        chosen_equation <- equation[[eq]]
        recombintion_rf <- rf_MML(equation = chosen_equation,
                                  tmp = tmp,
                                  parent = parent,
                                  share_nbr = share_nbr,
                                  eq_chosen = eq_chosen)
      }))
      colnames(estimated_rf) <- 'rf'
      estimated_rf$phasing <- rownames(estimated_rf)
      rownames(estimated_rf) <- NULL
    }
  }else{
    estimated_rf <- data.frame(sapply(names(equation),function(eq){
      chosen_equation <- equation[[eq]]
      recombintion_rf <- rf_MML(equation = chosen_equation,
                                tmp = tmp,
                                parent = parent,
                                share_nbr = share_nbr,
                                eq_chosen = eq_chosen)
    }))
    colnames(estimated_rf) <- 'rf'
    estimated_rf$phasing <- rownames(estimated_rf)
    rownames(estimated_rf) <- NULL
  }


  return(estimated_rf)
}


######Evaluate results######
#' evaluate estimated linkage result
#'
#' @param result linkage estimation results. A dataframe include:
#' marker_P, marker_Q, r, LOD, phasing_P, phasing_Q
#' @param genfile simulated .gen file used for check
#' @param mapfile simulated .map file used for check
#' @param plot TRUE/FALSE. If True, a plot will be saved in your folder 'Plot/'
#' @param filename filename that the figure will be saved
#'
#' @return a dataframe include the linkage results and simulated r and phasing
#' @export a figure
#'
#' @examples
results_evaluation <- function(result = Linkage_res$output,
                               genfile = genfile,
                               mapfile = mapfile,
                               plot =TRUE,
                               filename){
  result$simulated_rf <- unlist(lapply(1:nrow(result),function(i){
    distance_to_rf(simulated_distance(mapfile = mapfile,
                                      # mrk1 = 'Cm01452_000',
                                      # mrk2 = 'Cm16186_050'))
                                      mrk1 =  result[i,'marker_P'],
                                      mrk2 = result[i,'marker_Q']))

  }))
  result$simulated_phasing_correct <- unlist(lapply(1:nrow(result), function(i){
    sim_phase <- simulated_phasing(genfile = genfile,
                                   mrk1 = result[i,'marker_P'],
                                   mrk2 = result[i,'marker_Q'])
    esti_phase <- estimated_phasing(mrk1_phase = result[i,]$phasing_P,
                                    mrk2_phase = result[i,]$phasing_Q)
    sum(sim_phase == esti_phase)
  }))
  if(plot){
    dir.create('Plot_1/')
    plot_rf_comparison(result = result,
                       filename = paste0('Plot_1/',filename))
  }
  writeLines(paste0('Figures are generated and stored as: ',
                    paste0('Plot_1/',filename)))
  return(result)
}














