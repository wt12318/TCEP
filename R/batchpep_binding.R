#' Predict binding affinity of peptides with MHC in batch mode.
#'
#' @param pep_file Character, the path of peptide file needed to be processed, each row refers to a peptide sequence.
#' @param mhc_type MHC class, can be MHC-I or MHC-II.
#' @param pep_length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-15, for MHC-II, the length can be 11-30 or asis (take the length of input sequence as the peptide length)
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#'
#'
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples res <- batchpep_binding(pep_file=system.file("extdata", "random.pep", package = "MHCbinding"),
#'                                   mhc_type="MHC-I",pep_length=c(9,10),allele="HLA-A*01:01",pre_method="ann")
batchpep_binding <- function(pep_file,mhc_type,pep_length,allele,pre_method){
  pep <- read.table(pep_file)
  pep <- pep %>% mutate(seq_num=row_number())
  pep1 <- paste(pep$V1, sep= ",")
  res <- MHCbinding:::general_mhcbinding(mhc_type=mhc_type,length=pep_length,
                                         allele=allele,pre_method=pre_method,peptide=pep1)
  pep <- left_join(
    pep,res
  )
  colnames(pep)[1] <- "query_pep"
  return(pep)
}

