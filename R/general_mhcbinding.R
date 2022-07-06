#' General MCHbinding
#'
#' @param hla_type HLA class
#' @param client_path The path of local IEDB tools, used when setting get_method as client
#' @param ... Other para passed to mhcIbing or mhcIIbind
#'
#' @return Dataframe
#'
#' @examples
#' general_mhcbinding(hla_type = "I",peptide = "SLYNTVATLY",
#'                    allele = "HLA-A*01:01",length = "8",
#'                    pre_method = "netmhcpan_el",
#'                    client_path="~/software/mhc_i/src/",tmp_dir=tmp())

general_mhcbinding <- function(hla_type=c("I","II"),
                               client_path,peptide,allele,length,
                               pre_method,tmp_dir,method_type,mhcflurry_type="mt"){
  hla_type <- match.arg(hla_type)
  if (hla_type == "I"){
    res <- mhcIbinding_client(client_path=client_path,tmp_dir=tmp_dir,
                              peptide=peptide,allele =allele,
                              length =length,pre_method =pre_method,
                              method_type = method_type,mhcflurry_type=mhcflurry_type)
  }else{
    res <- mhcIIbinding_client(client_path=client_path,tmp_dir=tmp_dir,
                               peptide=peptide,allele =allele,
                               length =length,pre_method =pre_method,
                               method_type = method_type)
  }
  return(res)
}
