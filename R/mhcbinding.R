#' Predicting MHC binding
#'
#' @param peptide
#' @param allele
#' @param length
#' @param pre_method
#' @param get_method
#'
#' @return
#' @export
#'
#' @examples
#' test <- mhcIbinding(peptide = "SLYNTVATLY",
#'                     allele = "HLA-A*01:01",length = "8",
#'                     pre_method = "netmhcpan_el",get_method = "api")
mhcIbinding <- function(peptide = "SLYNTVATLYCVHQRIDV",
                        allele = "HLA-A*01:01",
                        length = c(	8, 9, 10, 11, 12, 13, 14, 15),
                       pre_method = c("ann","comblib_sidney2008","consensus",
                                      "netmhccons","netmhcpan_ba","netmhcpan_el",
                                      "netmhcstabpan","pickpocket","recommended",
                                      "smm","smmpmbec"),
                       get_method = c("api","client")){
  length <- match.arg(length)
  if (length(length) != length(allele)){
    stop("The length for which to make predictions must the paired with alleles")
  }
  length <- paste(as.character(length),collapse = ",")
  allele <- paste(allele,collapse = ",")
  pre_method <- match.arg(pre_method)
  get_method <- match.arg(get_method)

  temp_file <- tempfile()
  file.create(temp_file)
  command_run <- paste0('curl --data "method=',pre_method,'&sequence_text=',
                        peptide,'&allele=',allele,'&length=',length,'" ',
                        'http://tools-cluster-interface.iedb.org/tools_api/mhci/',
                        ' > ',temp_file)
  message("Retrieving data from server ...")
  mess <- try(system(command_run))
  if (mess == 0){
    message("Succeed !")
  }else{
    stop("Failed retrieving", immediate. = TRUE)
  }
  res <- read.table(temp_file,header = T)
  file.remove(temp_file)
  return(res)
}

