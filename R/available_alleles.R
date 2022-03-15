#' @title Show available alleles for MHC-I or MHC-II and corresponding methods.
#'
#' @param pre_method Character, which method to used for predicting MHC-peptide binding. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#' @param pre_type Character, which MHC type need to be predicted, can be MHC-I or MHC-II
#'
#' @return A dataframe (tibble) contains one column of available alleles.
#' @export
#'
#' @examples
#' available_alleles("MHC-I","ann")
available_alleles <- function(pre_type = c("MHC-I","MHC-II"),
                              pre_method){
  pre_type <- match.arg(pre_type)
  if (pre_type == "MHC-I"){
    pre_method <- match.arg(pre_method,
                            choices = c("ann","comblib_sidney2008","consensus",
                                        "netmhccons","netmhcpan_ba","netmhcpan_el",
                                        "netmhcstabpan","pickpocket","recommended",
                                        "smm","smmpmbec"))
    #print(mhcIallele[mhcIallele$method==pre_method,"alleles"])
    return(mhcIallele[mhcIallele$method==pre_method,"alleles"])
  }else{
    pre_method <- match.arg(pre_method,
                            choices = c("recommended","consensus","netmhciipan",
                                        "nn_align","smm_align","comblib","tepitope"))
    #print(mhcIIallele[mhcIIallele$method==pre_method,"alleles"])
    return(mhcIIallele[mhcIIallele$method==pre_method,"alleles"])
  }
}
