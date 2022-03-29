#' @title Predicting MHC-I binding based on IEDB local tools.
#'
#' @description  This is the wrapped function for IEDB commond tools.
#' @param peptide A character vector of input protein sequence.
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}.
#' @param length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-14
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#'
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples
#' test <- mhcIbinding_client(peptide = "SLYNTVATLY",
#'                            allele = "HLA-A*01:01",length = "8",
#'                            pre_method = "netmhcpan_el",client_path="~/software/mhc_i/src/")
mhcIbinding_client <- function(client_path,
                               peptide = c("GHAHKVPRRLLKAAR","LKAADASADADGSGSGSGSG"),
                               allele = c("HLA-A*01:01","HLA-A*03:01"),
                               length = c(8,9),
                               pre_method = c("ann","comblib_sidney2008","consensus",
                                              "netmhccons","netmhcpan_ba","netmhcpan_el",
                                              "netmhcstabpan","pickpocket","recommended",
                                              "smm","smmpmbec")){

  input_len <- nchar(peptide)
  max_len <- max(input_len)
  if(any(as.numeric(length)> input_len)){
    warning("Some of input peptides are shorter than the predicted core length specified by the user \n",immediate. = T)
  }

  temp_dir <- tempdir()
  file.create(paste0(temp_dir,"/a"))
  temp_list <- as.list(peptide)
  names(temp_list) <- paste0("seq",c(1:length(peptide)))
  seqinr::write.fasta(temp_list,names = names(temp_list),file.out = paste0(temp_dir,"/a"))

  length <- match.arg(as.character(length),
                      choices = c("8","9", "10", "11", "12", "13", "14"),
                      several.ok=T)
  length <- as.numeric(length)
  pre_method <- match.arg(pre_method)

  res <- vector("list",length(length)*length(allele))
  k <- 1
  for (i in seq_along(allele)){
    for (j in seq_along(length)){
      if(length[j] > max_len){
        next
      }
      file.create(paste0(temp_dir,"/b"))
      command_run <- paste0(client_path,'predict_binding.py ',pre_method," ",allele[i]," ",length[j]," ",paste0(temp_dir,"/a"),
                            ' > ',paste0(temp_dir,"/b"))
      message("Predicting using local IEDB tools ... \n")
      mess <- system(command_run)
      tmp <- read.table(paste0(temp_dir,"/b"),header = T)
      res[[k]] <- tmp
      k <- k + 1
    }
  }
  file.remove(paste0(temp_dir,"/a"),paste0(temp_dir,"/b"))
  res <- dplyr::bind_rows(res)
  return(res)
}
