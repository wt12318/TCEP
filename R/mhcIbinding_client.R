#' @title Predicting MHC peptide binding
#'
#' @description  This is the wrapped function for IEDB commond tools.
#' @param peptide A character vector of input protein sequence.
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}.
#' @param length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-14
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#' @param tmp_dir Character, the temp dir
#' @param mhc_type MHC type, I or II
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).

mhcbinding_client <- function(client_path,
                              peptide,
                              allele,
                              length,
                              pre_method,tmp_dir,mhc_type){
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
      if(length[j] > max_len){
        next
      }

      ## 如果有一个序列小于指定的长度就会报错，因此去掉小于指定长度的序列
      if (length[j] < min_len){
        temp_list <- as.list(peptide)
        names(temp_list) <- paste0("seq",c(1:length(peptide)))
        seqinr::write.fasta(temp_list,names = names(temp_list),file.out = paste0(temp_dir,"/a"))
      }else{
        temp_list <- as.list(peptide[which(input_len >= length[j])])
        names(temp_list) <- paste0("seq",c(1:length(temp_list)))
        seqinr::write.fasta(temp_list,names = names(temp_list),file.out = paste0(temp_dir,"/a"))
      }


      file.create(paste0(temp_dir,"/b"))
      if (mhc_type == "II"){
        command_run <- paste0(client_path,'mhc_II_binding.py ',pre_method," ",allele[i]," ",paste0(temp_dir,"/a "),length[j]," ",
                              ' > ',paste0(temp_dir,"/b"))
      }else{
        command_run <- paste0(client_path,'predict_binding.py ',pre_method," ",allele[i]," ",length[j]," ",paste0(temp_dir,"/a"),
                              ' > ',paste0(temp_dir,"/b"))
      }
      message("Predicting using local IEDB tools ... \n")
      mess <- system(command_run)
      tmp <- read.table(paste0(temp_dir,"/b"),header = T)
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

  if(temp_dir != tempdir()){
    unlink(temp_dir,recursive = T)
  }
  res <- dplyr::bind_rows(res)
  return(res)
}




#' @title Predicting MHC-I binding based on IEDB local tools.
#'
#' @description  This is the wrapped function for IEDB commond tools.
#' @param peptide A character vector of input protein sequence.
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}.
#' @param length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-14
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#' @param tmp_dir Character, the temp dir
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
                               length = c(8,9),
                               pre_method = c("ann","comblib_sidney2008","consensus",
                                              "netmhccons","netmhcpan_ba","netmhcpan_el",
                                              "netmhcstabpan","pickpocket","IEDB_recommended",
                                              "smm","smmpmbec"),tmp_dir=tempdir()){
  length <- match.arg(as.character(length),
                      choices = c("8","9", "10", "11", "12", "13", "14"),
                      several.ok=T)
  pre_method <- match.arg(pre_method)
  res <- MHCbinding:::mhcbinding_client(client_path=client_path,peptide=peptide,
                                        allele=allele,length=length,pre_method=pre_method,tmp_dir=tmp_dir,
                                        mhc_type = "I")
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
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples
#' test <- mhcIbinding_client(peptide = "GHAHKVPRRLLKAAR",
#'                            allele = "DRB1*15:01",length = "14",
#'                            pre_method = "recommended",client_path="~/software/mhc_i/src/",tmp_dir=tempdir())
mhcIIbinding_client <- function(client_path,
                                peptide = c("GHAHKVPRRLLKAAR","LKAADASADADGSGSGSGSG"),
                                allele = c("DRB1*15:01","DPB1*04:01"),
                                length = c(14,15),
                                pre_method = c("IEDB_recommended", "consensus3", "NetMHCIIpan", "nn_align", "smm_align",
                                               "comblib", "sturniolo"),tmp_dir=tempdir()){
  length <- match.arg(as.character(length),
                      choices = as.character(seq(11,30)),
                      several.ok=T)
  pre_method <- match.arg(pre_method)
  res <- MHCbinding:::mhcbinding_client(client_path=client_path,peptide=peptide,
                                        allele=allele,length=length,pre_method=pre_method,tmp_dir=tmp_dir,
                                        mhc_type = "II")
  return(res)
}
