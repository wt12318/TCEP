#' @title Show available prediction methods for HLA-antigen binding, antigen processing or Immunogenicity.
#'
#' @param pre_type Character, which kind of predicted method, could be Binding, Processing, or Immuno
#'
#' @return A dataframe (tibble) contains two columns of method and corresponding version.
#' @export
#'
#' @examples
#' available_methods("Binding")
available_methods <- function(pre_type=c("Binding","Processing","Immuno")){

  pre_type <- match.arg(pre_type)
  if (pre_type == "Binding"){
    dt_mhci <- TCAP::mhcIbinding_api_methods
    dt_mhci$method[9] <- "IEDB_recommended"
    cat(crayon::green("HLA-I antigen binding methods: ","\n"))
    print(dplyr::as_tibble(bind_rows(
      dt_mhci,
      data.frame(method = c("mhcflurry","mhcnuggets"),
                 version = c("2.0.0","2.3"))
    )))
    methods_mhci <- c(dt_mhci[,1],"mhcflurry","mhcnuggets")

    dt_mhcii <- data.frame(
      method = c("comblib", "consensus3", "IEDB_recommended", "netmhciipan_el",
                 "netmhciipan_ba","nn_align", "smm_align","sturniolo","mhcnuggets"),
      version = c("1.0","2.22","2.22","4.0","4.0","2.3","1.1","1.0","2.3")
    )
    cat(crayon::green("HLA-II antigen binding methods: ","\n"))
    print(dplyr::as_tibble(dt_mhcii))
    methods_pre <- list(mhc_i=methods_mhci,
                        mhc_ii=dt_mhcii$method)
  }else if (pre_type == "Processing"){
    methods_pre <- c("Netchop","NetCTLpan")
    cat(crayon::green("Antigen processing methods: ","\n"))
    print(methods_pre)
  }else {
    methods_pre <- c("IEDB","PRIME2.0","DeepImmuno","Seq2Neo-CNN","Immuno-GNN")
    cat(crayon::green("Antigen immunogenicity methods: ","\n"))
    print(methods_pre)
  }
  return(invisible(methods_pre))
}
