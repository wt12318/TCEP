#' General MCHbinding
#'
#' @param mhc_type MHC class
#' @param ... Other para passed to mhcIbing or mhcIIbind
#'
#' @return Dataframe
#'
#' @examples
#' general_mhcbinding(mhc_type = "MHC-I",peptide = "SLYNTVATLY",
#'                   allele = "HLA-A*01:01",length = "8",
#'                    pre_method = "netmhcpan_el")

general_mhcbinding <- function(mhc_type=c("MHC-I","MHC-II"),...){
  mhc_type <- match.arg(mhc_type)
  if (mhc_type == "MHC-I"){
    res <- mhcIbinding(...)
    return(res)
  }else{
    res <- mhcIIbinding(...)
    return(res)
  }
}
