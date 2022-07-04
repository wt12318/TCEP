#' @title Show available prediction methods for MHC-I or MHC-II
#'
#' @param get_methods Character, which method to used for predicting MHC-peptide binding, can be api, refering to IEDB API and client, refering to local client
#' @param pre_type Character, which MHC type need to be predicted, can be MHC-I or MHC-II
#'
#' @return A dataframe (tibble) contains two columns of method and corresponding version.
#' @export
#'
#' @examples
#' available_methods("MHC-I","Binding")
available_methods <- function(pre_type=c("Binding","Processing","Immuno")){

  pre_type <- match.arg(pre_type)
  if (pre_type == "Binding"){
    dt_mhci <- MHCbinding::mhcIbinding_api_methods
    dt_mhci$method[9] <- "IEDB_recommended"
    cat(crayon::green("HLA-I antige binding methods: ","\n"))
    print(dplyr::as_tibble(dt_mhci))
    methods_mhci <- c(dt_mhci[,1],"mhcflurry","mhcnuggets")

    dt_mhcii <- data.frame(
      method = c("comblib", "consensus3", "IEDB_recommended", "netmhciipan_el",
                 "netmhciipan_ba","nn_align", "smm_align","sturniolo"),
      version = c("1.0","2.22","2.22","4.0","4.0","2.3","1.1","1.0")
    )
    cat(crayon::green("HLA-II antige binding methods: ","\n"))
    print(dplyr::as_tibble(dt_mhcii))
    methods_pre <- c(methods_mhci,dt_mhcii$method)
  }else if (pre_type == "Processing"){
    methods_pre <- c("Netchop","NetCTLpan")
    cat(crayon::green("Antigen processing methods: ","\n"))
    print(methods_pre)
  }else {
    methods_pre <- c("IEDB","NetTepi","PRIME","iPred","DeepImmuno","Seq2Neo-CNN")
    cat(crayon::green("Antigen immunogenicity methods: ","\n"))
    print(methods_pre)
  }
  return(invisible(methods_pre))
}
