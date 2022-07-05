#' @title Show available alleles for MHC-I or MHC-II and corresponding methods.
#'
#' @param pre_type Character, what kind of method to use, could be "Binging", "Processing", or "Immuno"
#' @param pre_method Character, which method to used for predicting MHC-peptide binding. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#' @param HLA_type, Character, HLA type need to be predicted, can be I or II (Processing and Immuno can only be I).
#' @return A dataframe (tibble) contains one column of available alleles.
#' @export
#'
#' @examples
#' available_alleles("MHC-I","ann")
available_alleles <- function(pre_type = c("Binding","Processing","Immuno"),
                              HLA_type = c("I","II"),
                              pre_method){
  pre_type <- match.arg(pre_type)
  if (pre_type == "Binding"){
    HLA_type <- match.arg(HLA_type)
    if (HLA_type == "I"){
      pre_method <- match.arg(pre_method,
                              choices = c("ann","comblib_sidney2008","consensus",
                                          "netmhccons","netmhcpan_ba","netmhcpan_el",
                                          "netmhcstabpan","pickpocket","IEDB_recommended",
                                          "smm","smmpmbec","mhcnuggets","mhcflurry"))
      #print(mhcIallele[mhcIallele$method==pre_method,"alleles"])
      alleles <- mhcIallele$alleles[mhcIallele$method==pre_method]
    }else{
      pre_method <- match.arg(pre_method,
                              choices = c("comblib", "consensus3", "IEDB_recommended", "netmhciipan_el",
                                          "netmhciipan_ba","nn_align", "smm_align","sturniolo","mhcnuggets"))
      #print(mhcIIallele[mhcIIallele$method==pre_method,"alleles"])
      alleles <- mhcIIallele_client$alleles[mhcIIallele_client$method==pre_method]
    }
  }else if (pre_type == "Processing"){
    pre_method <- match.arg(pre_method,
                            choices = c("Netchop","NetCTLpan"))
    if (HLA_type == "II"){
      cat(crayon::red("NetCTLpan only for HLA-I","\n"))
    }
    alleles <- NetCTLpan_alleles$Allele
  }else {
    pre_method <- match.arg(pre_method,
                            choices = c("IEDB","PRIME2.0","DeepImmuno","Seq2Neo-CNN"))
    if (HLA_type == "II"){
      cat(crayon::red("Immuno only for HLA-I","\n"))
    }
    alleles <- immuno_alleles$alleles[immuno_alleles$methods==pre_method]
  }
  return(alleles)
}
