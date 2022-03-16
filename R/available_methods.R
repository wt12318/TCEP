#' @title Show available prediction methods for MHC-I or MHC-II
#'
#' @param get_methods Character, which method to used for predicting MHC-peptide binding, can be api, refering to IEDB API and client, refering to local client
#' @param pre_type Character, which MHC type need to be predicted, can be MHC-I or MHC-II
#'
#' @return A dataframe (tibble) contains two columns of method and corresponding version.
#' @export
#'
#' @examples
#' available_methods("api","MHC-I")
available_methods <- function(get_methods=c("api","client"),
                              pre_type=c("MHC-I","MHC-II")){
  get_methods <- match.arg(get_methods)
  pre_type <- match.arg(pre_type)
  if (get_methods == "api"){
    if(pre_type == "MHC-I"){
      print(dplyr::as_tibble(MHCbinding::mhcIbinding_api_methods))
      return(MHCbinding::mhcIbinding_api_methods[,1])
    }else{
      print(dplyr::as_tibble(MHCbinding::mhcIIbinding_api_methods))
      return(MHCbinding::mhcIIbinding_api_methods[,1])
    }
  }else{
    NULL ##To do
  }
}
