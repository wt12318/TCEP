#' @title Show available predicted peptide length for sepcific method and allele for MHC-I
#'
#' @param pre_method Character, which method used to do prodiction
#' @param pre_allele Character, which MHC-I allele need to be predicted
#'
#' @return A numberic vector
#' @export
#'
#' @examples
#' available_len("comblib_sidney2008","HLA-A*30:01")
available_len <- function(pre_method,pre_allele){
  dt <- MHCbinding::available_lens
  dt <- dt %>%
    filter(MHC %in% pre_allele & method %in% pre_method)
  if (pre_method == "DeepImmuno"){
    return(c(9,10))
  }else{
   return(unique(dt$PeptideLength))
  }
}
