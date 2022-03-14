#' Show avaliable prediction methods for MHC-I or MHC-II
#'
#' @param get_methods
#' @param pre_type
#'
#' @return
#' @export
#'
#' @examples
available_methods <- function(get_methods=c("api","client"),
                              pre_type=c("MHC-I","MHC-II")){
  get_methods <- match.arg(get_methods)
  pre_type <- match.arg(pre_type)
  if (get_methods == "api"){
    if(pre_type == "MHC-I"){
      print(dplyr::as_tibble(MHCbinding::mhcIbinding_api_methods))
    }else{
      print(dplyr::as_tibble(MHCbinding::mhcIIbinding_api_methods))
    }
  }else{
    NULL ##To do
  }
}
