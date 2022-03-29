#' General MCHbinding
#'
#' @param mhc_type MHC class
#' @param get_method The way to predict, can be api or client
#' @param client_path The path of local IEDB tools, used when setting get_method as client
#' @param ... Other para passed to mhcIbing or mhcIIbind
#'
#' @return Dataframe
#'
#' @examples
#' general_mhcbinding(mhc_type = "MHC-I",peptide = "SLYNTVATLY",
#'                   allele = "HLA-A*01:01",length = "8",
#'                    pre_method = "netmhcpan_el")

general_mhcbinding <- function(get_method=c("api","client"),mhc_type=c("MHC-I","MHC-II"),client_path,...){
  mhc_type <- match.arg(mhc_type)
  get_method <- match.arg(get_method)
  if (mhc_type == "MHC-I"){
    if (get_method == "api"){
      res <- mhcIbinding(...)
    }else {
      res <- mhcIbinding_client(client_path=client_path,...)
    }
    return(res)
  }else{
    if (get_method == "api"){
      res <- mhcIIbinding(...)
    }else {
      res <- mhcIIbinding_client(client_path=client_path,...)
    }
    return(res)
  }
}
