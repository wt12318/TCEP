
#' Sliding peptide sequence
#'
#' @param pep_seq, Character vector, peptide sequence
#' @param req_len, numberic, the length of resulting sub-sequence
#'
#' @return Character vector of sub-sequence
#' @export
#'
#' @examples split_pep(c("GHAHKVPRRLLKAA", "SLYNTVATLY"),8)
split_pep <- function(pep_seq, req_len){
  req_len <- req_len - 1
  pep_res <- sapply(pep_seq,
                    function(x){
                      t <- sapply(c(1:nchar(x)),
                                  function(y){
                                    if (nchar(x)-req_len < y){
                                      NA
                                      }else{
                                        substr(x,y,y+req_len)
                                        }
                                    })
                      t <- t[!is.na(t)]
                    },simplify = FALSE)
  pep_res <- data.frame(
    seq_num = rep((1:length(pep_seq)),times=lengths(pep_res)),
    peptide = pep_res %>% unlist() %>% unname()
  )
  pep_res$length <- nchar(pep_res$peptide)
  pep_res <- pep_res %>%
    group_by(seq_num) %>% mutate(start=row_number()) %>%
    ungroup() %>%
    mutate(end = start + req_len)
  return(pep_res)
}

#' Predict proteasomal processing using NetChop
#'
#' @param pep, character vector of peptide sequence
#' @param temp_dir, temp dir
#' @param netchop_path, installed IEDB Netchop path
#'
#' @return Dataframe containing 3 column: "amino_acid", "prediction_score", "sequence_id"
#' @export
#'
#' @examples netchop_processing(c("GHAHKVPRRLLKAA","SLYNTVATLY"),"~/tmp/","~/software/netchop/")
netchop_processing <- function(pep,temp_dir,netchop_path){
  temp_list <- as.list(pep)
  names(temp_list) <- paste0("seq",c(1:length(temp_list)))

  t <- mapply(function(x,y){
    seqinr::write.fasta(x,names = y,
                        file.out = paste0(temp_dir,"/a.",y))
  },temp_list,names(temp_list))
  command_run <- sapply(names(temp_list),function(x){
    paste0(netchop_path, "predict.py --method netchop ",
           paste0(temp_dir,"/a.",x)," -n > ",paste0(temp_dir,"/b.",x))
  })
  ##run
  mess <- lapply(command_run,function(x){system(x)})
  tmp <- lapply(names(temp_list),function(x){
    read.table(paste0(temp_dir,"/b.",x),skip = 1)
  }) %>%
    bind_rows() %>%
    select(-V1)
  colnames(tmp) <- c("amino_acid", "prediction_score", "sequence_id")
  return(tmp)
}

#' seq2neo help function to get ic50 and tap
#'
#' @param pep_len peptide length
#' @param allele hla alleles
#' @param peptide peptide sequences
#' @param tmp_dir tmp dir
#' @param netchop_path the path of neochop
#' @param client_path the path of IEDB tools
#'
#' @return dataframe containg TAP and IC50
#' @export
#'
#' @examples seq2neo_help(pep_len=c(8,9),allele=c("HLA-A*01:01","HLA-A*01:02"),
#'                        tmp_dir="~/tmp/", netchop_path="~/software/netchop/",
#'                        peptide=c("SLYNTVATLY","GHAHKVPR"),
#'                        client_path = "~/software/mhc_i/src/")
seq2neo_help <- function(pep_len, allele, peptide, tmp_dir, netchop_path, client_path){
  tap_pre <- TCAP:::general_mhcbinding(hla_type = "I", length = pep_len,
                                        allele = gsub("[*]","",allele),pre_method = "NetCTLpan",
                                        method_type="Processing",
                                        peptide = peptide,
                                        tmp_dir=tmp_dir,netchop_path = netchop_path)
  tap_pre <- tap_pre %>% select(peptide,allele,tap_pre)
  ic50 <- TCAP:::general_mhcbinding(hla_type = "I", length = pep_len,
                                        allele = allele, pre_method = "netmhcpan_ba",
                                        method_type="Binding",
                                        peptide = peptide,
                                        tmp_dir= tmp_dir,client_path = client_path)
  ic50 <- ic50 %>% select(peptide,allele,ic50)
  ic50 <- left_join(
    ic50 %>% mutate(index=paste(peptide,gsub("[*]","",allele),sep = "_")),
    tap_pre %>% mutate(index=paste(peptide,allele,sep = "_")) %>% select(index,tap_pre)
  )
  ic50 <- ic50 %>% select(peptide,allele,ic50,tap_pre)
  ic50$allele <- gsub("[*]","",ic50$allele)
  return(ic50)
}

#' @title Predicting HLA peptide binding
#'
#' @description  This is the wrapped function for IEDB commond tools.
#' @param peptide A character vector of input protein sequence.
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}.
#' @param length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-14
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#' @param tmp_dir Character, the temp dir
#' @param hla_type HLA type, I or II
#' @param method_type, which type prediction method used, could be "Binding", "Processing" or "Immuno"
#' @param mhcflurry_type, is calculating wt peptide for mhcflurry
#' @param netchop_path, the path of Netchop
#' @param Immuno_IEDB_path, the path of IEDB Class I Immunogenicity tool
#' @param Immuno_Deepimmuno_path, the deepimmuno-cnn.py script path
#' @param Deepimmuno_env, the conda envrionment of Deepimmuno
#' @param MixMHCpred_path, the intalled Mixmhcpred path
#' @param PRIME_path, PRIME tool path
#' @param seq2neo_env, the conda env of seq2neo
#' @param seq2neo_path, the of `immuno_Prediction` dir of Seq2Neo
#' @export
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).

mhcbinding_client <- function(client_path,
                              peptide,
                              allele=NULL,
                              length,
                              pre_method,tmp_dir,hla_type,
                              method_type,
                              mhcflurry_type="mt",
                              mhcflurry_env="mhcflurry-env",
                              mhcnuggets_env="mhcnuggets",
                              netchop_path,Immuno_IEDB_path,
                              Immuno_Deepimmuno_path,Deepimmuno_env,
                              MixMHCpred_path,PRIME_path,seq2neo_env,seq2neo_path){
  input_len <- nchar(peptide)
  max_len <- max(input_len)
  min_len <- min(input_len)
  if(any(as.numeric(length)> input_len)){
    warning("Some of input peptides are shorter than the predicted core length specified by the user \n",immediate. = T)
  }

  if (! dir.exists(tmp_dir)){
    dir.create(tmp_dir)
  }

  temp_dir <- tmp_dir
  file.create(paste0(temp_dir,"/a"))
  length <- as.numeric(length)

  res <- vector("list",length(length)*length(allele))
  k <- 1
  for (i in seq_along(allele)){
    for (j in seq_along(length)){
      if (hla_type == "I"){
        available_length <- available_len(pre_method,allele[i])
      }else{
        available_length <- c(11:30)
      }
      if(length[j] > max_len | !(length[j] %in% available_length)){
        next
      }

      ## 如果有一个序列小于指定的长度就会报错，因此去掉小于指定长度的序列
      ship_lines <- 0
      filter_pep <- peptide[which(input_len >= length[j])]
      if (pre_method %in% c("mhcflurry","mhcnuggets","NetCTLpan","IEDB","DeepImmuno","PRIME2.0","Seq2Neo-CNN")){
        if (pre_method == "mhcflurry"){
          if (mhcflurry_type == "wt"){
            pep_res <- data.frame(peptide = filter_pep)
            pep_res$allele <- allele[i]
          }else{
            pep_res <- split_pep(filter_pep,req_len = length[j])
            pep_res$allele <- allele[i]##add seq_num and length
          }
          write.csv(pep_res,file = paste0(temp_dir,"/a"),row.names = FALSE)
        }

        if (pre_method %in% c("mhcnuggets","IEDB", "DeepImmuno","PRIME2.0")){
          ##mhcnuggets,IEDB,DeepImmuno
          pep_res <- split_pep(filter_pep,req_len = length[j])
          pep_res_pep <- data.frame(peptide = pep_res$peptide)
          if (pre_method %in% c("DeepImmuno")){
            pep_res_pep$allele <- allele[i]
            write.table(pep_res_pep,file = paste0(temp_dir,"/a"),sep = ",",col.names = F,row.names = F,quote = F)
          }else{
            write.table(pep_res_pep,file = paste0(temp_dir,"/a"),sep = "\t",col.names = F,row.names = F,quote = F)
          }

        }

        if (pre_method == "NetCTLpan"){
          temp_list <- as.list(filter_pep)
          names(temp_list) <- paste0("seq",c(1:length(temp_list)))
          pep_res <- split_pep(filter_pep,req_len = length[j])
          t <- mapply(function(x,y){
            seqinr::write.fasta(x,names = y,
                                file.out = paste0(temp_dir,"/a.",y))
          },temp_list,names(temp_list))
        }

        if (pre_method == "Seq2Neo-CNN"){
          pep_res <- split_pep(filter_pep,req_len = length[j])
          tap_ic50 <- TCAP::seq2neo_help(pep_len=length[j],allele=allele[i],
                                               tmp_dir=temp_dir, netchop_path=netchop_path,
                                               peptide=filter_pep,
                                               client_path = client_path)
          write.table(tap_ic50,file = paste0(temp_dir,"/a"),sep = ",",
                      row.names = FALSE, col.names = TRUE,quote = FALSE)
        }

      }else{
        temp_list <- as.list(filter_pep)
        names(temp_list) <- paste0("seq",c(1:length(temp_list)))
        seqinr::write.fasta(temp_list,names = names(temp_list),file.out = paste0(temp_dir,"/a"))
      }


      file.create(paste0(temp_dir,"/b"))
      if (method_type == "Binding"){
        if (hla_type == "II"){
          if (pre_method == "mhcnuggets"){
            system(paste0("touch ",temp_dir,"/pre.py"))
            writeLines(paste0("from mhcnuggets.src.predict import predict\npredict(",
                              paste0("class_='",hla_type,"',"),
                              paste0("peptides_path='",normalizePath(temp_dir),"/a'",","),
                              paste0("mhc='",allele[i],"')")),
                       paste0(temp_dir,"/pre.py"))
            command_run <- paste0(paste0("conda run -n ",mhcnuggets_env," python "),
                                  paste0(temp_dir,"/pre.py")," > ",paste0(temp_dir,"/b"))
          }else{
            command_run <- paste0(client_path,'mhc_II_binding.py ',pre_method," ",allele[i]," ",paste0(temp_dir,"/a "),length[j]," ",
                                  ' > ',paste0(temp_dir,"/b"))
          }
        }else{
          if (pre_method %in% c("mhcflurry","mhcnuggets")){
            ##add mhcnuggets method
            if (pre_method == "mhcflurry"){
              command_run <- paste0(paste0("conda run -n ",mhcflurry_env," mhcflurry-predict "),
                                    paste0(temp_dir,"/a"), " > ",paste0(temp_dir,"/b"))
            }
            if (pre_method == "mhcnuggets"){
              system(paste0("touch ",temp_dir,"/pre.py"))
              writeLines(paste0("from mhcnuggets.src.predict import predict\npredict(",
                                paste0("class_='",hla_type,"',"),
                                paste0("peptides_path='",normalizePath(temp_dir),"/a'",","),
                                paste0("mhc='",allele[i],"')")),
                         paste0(temp_dir,"/pre.py"))
              command_run <- paste0(paste0("conda run -n ",mhcnuggets_env," python "),
                                    paste0(temp_dir,"/pre.py")," > ",paste0(temp_dir,"/b"))
            }
          }else{
            command_run <- paste0(client_path,'predict_binding.py ',pre_method," ",allele[i]," ",length[j]," ",paste0(temp_dir,"/a"),
                                  ' > ',paste0(temp_dir,"/b"))
          }
        }
      }

      if (method_type == "Processing"){
        ##netctlpan only for HLA-I
        command_run <- sapply(names(temp_list),function(x){
          paste0(netchop_path,
                 paste0("predict.py --method netctlpan -a ",allele[i]," -l ",length[j]," -n "),
                 paste0(temp_dir,"/a.",x)," > ",paste0(temp_dir,"/b.",x))
        })
        }


      if (method_type == "Immuno"){
        if (pre_method == "IEDB"){
          ##add iedb method
          command_run <- paste0(Immuno_IEDB_path,
                                "predict_immunogenicity.py --allele=", allele[i]," ",
                                paste0(temp_dir,"/a"), " > ",paste0(temp_dir,"/b"))
        }
        if (pre_method == "PRIME2.0"){
          ##add PRIME2.0
          command_run <- paste0(PRIME_path,
                                "/PRIME -i ", paste0(temp_dir,"/a"),
                                " -o ",paste0(temp_dir,"/b")," -a ",allele[i],
                                " -mix ",MixMHCpred_path)
        }
        if (pre_method == "DeepImmuno"){
          ##add DeepImmuno
          command_run <- paste0("cd ",Immuno_Deepimmuno_path,"; conda run -n ",Deepimmuno_env," python ",
                                "deepimmuno-cnn.py --mode 'multiple' --intdir ",
                                paste0(temp_dir,"/a")," --outdir ",temp_dir)
        }
        if (pre_method == "Seq2Neo-CNN"){
          seq2neo_path <- normalizePath(seq2neo_path)
          system(paste0("touch ",temp_dir,"/pre.py"))
          writeLines(paste0("import sys\nsys.path.append('",seq2neo_path,"')\n",
                            "from _cnn import *\n",
                            "import pandas as pd\n",
                            "a = pd.read_csv('",paste0(normalizePath(temp_dir),"/a"),"')\n",
                            "file_process(a,'",normalizePath(temp_dir),"')"),
                     paste0(temp_dir,"/pre.py"))
          command_run <- paste0(paste0("conda run -n ",seq2neo_env," python "),
                                paste0(temp_dir,"/pre.py"))
        }
      }

      cat("Predicting using ",pre_method,"\n")

      if (pre_method == "NetCTLpan"){
        mess <- lapply(command_run,function(x){system(x)})
      }else{
        mess <- system(command_run)
      }

      if (pre_method == "mhcflurry"){
        tmp <- read.table(paste0(temp_dir,"/b"),header = T,skip = 4,sep = ",")
      }else if(pre_method == "mhcnuggets"){
        tmp <- data.table::fread(paste0(temp_dir,"/b"),skip = "peptide,ic50",data.table = FALSE)
        tmp <- left_join(tmp,pep_res)
        tmp$allele <- allele[i]
      }else if (pre_method == "NetCTLpan"){
        tmp <- lapply(names(temp_list),function(x){
          read.table(paste0(temp_dir,"/b.",x),skip = 1)
        }) %>%
          bind_rows() %>%
          select(-V1)
        colnames(tmp) <- c("peptide","mhc_pre","tap_pre","cleavage_pre","combined_score_%rank")
        tmp <- left_join(tmp,pep_res)
        tmp$allele <- allele[i]
      }else if (pre_method == "IEDB") {
        tmp <- data.table::fread(paste0(temp_dir,"/b"),skip = "peptide,length,score",data.table = FALSE)
        tmp <- left_join(tmp,pep_res %>% select(-length))
        tmp$allele <- allele[i]
      }else if (pre_method == "DeepImmuno"){
        tmp <- read.table(paste0(temp_dir,"/deepimmuno-cnn-result.txt"),header = T)
        tmp <- left_join(tmp,pep_res)
        tmp$allele <- allele[i]
      }else if (pre_method == "PRIME2.0") {
        tmp <- data.table::fread(paste0(temp_dir,"/b"),data.table = FALSE)
        tmp <- tmp[,1:4]
        colnames(tmp)[1:4] <- c("peptide","PRIME_rank","PRIME_score","MixMHCpred_rank")
        tmp <- left_join(tmp,pep_res)
        tmp$allele <- allele[i]
      }else if (pre_method == "Seq2Neo-CNN") {
        tmp <- read.csv(paste0(temp_dir,"/cnn_results.csv"))
        tmp <- tmp %>% select(-HLA,-pseudosequence)
        colnames(tmp)[1] <- "peptide"
        tmp <- left_join(tmp,pep_res)
        tmp$allele <- allele[i]
      }else {
        tmp <- read.table(paste0(temp_dir,"/b"),header = T)
      }
      if (pre_method == "consensus"){
        suppressWarnings(tmp <- tmp %>% mutate(across(7:13,as.numeric)))
      }
      res[[k]] <- tmp
      k <- k + 1
    }
  }


  if(file.exists(paste0(temp_dir,"/a"))){
    file.remove(paste0(temp_dir,"/a"))
  }

  if(file.exists(paste0(temp_dir,"/b"))){
    file.remove(paste0(temp_dir,"/b"))
  }

  # if(temp_dir != tempdir()){
  #   unlink(temp_dir,recursive = T)
  # }
  res <- dplyr::bind_rows(res)
  return(res)
}




#' @title Predicting MHC-I binding based on IEDB local tools.
#'
#' @description  This is the wrapped function for IEDB commond tools.
#' @param peptide A character vector of input protein sequence.
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}. For `Netchop`, allele is not needed.
#' @param length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-14
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#' @param tmp_dir Character, the temp dir
#' @param method_type, which type prediction method used, could be "Binding", "Processing" or "Immuno"
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples
#' test <- mhcIbinding_client(peptide = "SLYNTVATLY",
#'                            allele = "HLA-A*01:01",length = "8",
#'                            pre_method = "netmhcpan_el",client_path="~/software/mhc_i/src/",tmp_dir=tempdir())
mhcIbinding_client <- function(client_path,
                               peptide = c("GHAHKVPRRLLKAA","SLYNTVATLY"),
                               allele = c("HLA-A*01:01","HLA-A*03:01"),
                               length ,
                               pre_method = c("ann","comblib_sidney2008","consensus",
                                              "netmhccons","netmhcpan_ba","netmhcpan_el",
                                              "netmhcstabpan","pickpocket","IEDB_recommended",
                                              "smm","smmpmbec","mhcflurry","mhcnuggets",
                                              "IEDB","PRIME2.0","DeepImmuno",
                                              "Seq2Neo-CNN","Netchop","NetCTLpan"),
                               tmp_dir=tempdir(),
                               method_type = c("Binding","Processing","Immuno"),
                               mhcflurry_type="mt",
                               mhcflurry_env="mhcflurry-env",mhcnuggets_env="mhcnuggets",
                               netchop_path,
                               Immuno_IEDB_path,Immuno_Deepimmuno_path,Deepimmuno_env,
                               MixMHCpred_path,PRIME_path,
                               seq2neo_env,seq2neo_path){
  length <- match.arg(as.character(length),
                      choices = available_len(pre_method,allele),
                      several.ok=T)
  pre_method <- match.arg(pre_method)
  method_type <- match.arg(method_type)
  res <- TCAP::mhcbinding_client(client_path=client_path,peptide=peptide,
                                        allele=allele,length=length,pre_method=pre_method,tmp_dir=tmp_dir,
                                        hla_type = "I",
                                        method_type=method_type,mhcflurry_type=mhcflurry_type,
                                        mhcflurry_env=mhcflurry_env,mhcnuggets_env=mhcnuggets_env,
                                        netchop_path=netchop_path,Immuno_IEDB_path=Immuno_IEDB_path,
                                        Immuno_Deepimmuno_path=Immuno_Deepimmuno_path,
                                        Deepimmuno_env=Deepimmuno_env,
                                        MixMHCpred_path=MixMHCpred_path,PRIME_path=PRIME_path,
                                       seq2neo_env=seq2neo_env,seq2neo_path=seq2neo_path)
  return(res)
}


#' @title Predicting MHC-II binding based on IEDB local tools.
#'
#' @description  This is the wrapped function for IEDB commond tools.
#' @param peptide A character vector of input protein sequence.
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}.
#' @param length A numeric or character vector, indicating the length for which to make predictions. For MHC-II, the length can be 11-30
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#' @param tmp_dir Character, the temp dir
#' @param method_type, which type prediction method used, for HLA-II , it can only be "Binding".
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples
#' test <- mhcIIbinding_client(peptide = "GHAHKVPRRLLKAAR",
#'                            allele = "DRB1*15:01",length = "14",
#'                            pre_method = "netmhciipan_ba",client_path="~/software/mhc_ii/",tmp_dir=tempdir())
mhcIIbinding_client <- function(client_path,
                                peptide = c("GHAHKVPRRLLKAAR","LKAADASADADGSGSGSGSG"),
                                allele = c("DRB1*15:01","DPB1*04:01"),
                                length = c(14,15),
                                pre_method = c("comblib", "consensus3",
                                               "IEDB_recommended", "netmhciipan_el",
                                               "netmhciipan_ba","nn_align",
                                               "smm_align","sturniolo","mhcnuggets"),
                                tmp_dir=tempdir(),
                                method_type = c("Binding"),
                                mhcnuggets_env="mhcnuggets"){
  length <- match.arg(as.character(length),
                      choices = as.character(seq(11,30)),
                      several.ok=T)
  pre_method <- match.arg(pre_method)
  method_type <- match.arg(method_type)
  res <- TCAP::mhcbinding_client(client_path=client_path,peptide=peptide,
                                        allele=allele,length=length,pre_method=pre_method,tmp_dir=tmp_dir,
                                        hla_type  = "II",method_type=method_type,
                                        mhcnuggets_env=mhcnuggets_env)
  return(res)
}
