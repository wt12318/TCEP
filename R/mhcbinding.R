#' Predicting MHC-I binding
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
mhcIbinding <- function(peptide = c("GHAHKVPRRLLKAAR","LKAADASADADGSGSGSGSG"),
                        allele = c("HLA-A*01:01","HLA-A*03:01"),
                        length = c(8,9),
                       pre_method = c("ann","comblib_sidney2008","consensus",
                                      "netmhccons","netmhcpan_ba","netmhcpan_el",
                                      "netmhcstabpan","pickpocket","recommended",
                                      "smm","smmpmbec")){
  if (length(peptide) != 1){
    ##To submit multiple sequences at a time,
    ##escape the special characters in a fasta-formatted sequence with URI codes
    pep_number <- map_chr(c(1:length(peptide)),~gsub("number",.x,"%3Epeptidenumber%0A"))
    pep_con <- c(paste0(peptide[1:length(peptide)-1],"%0A"),peptide[length(peptide)])
    peptide <- paste0(pep_number,pep_con,collapse = "")
  }
  length <- match.arg(as.character(length),
                      choices = c("8","9", "10", "11", "12", "13", "14", "15"),
                      several.ok=T)
  if (length(length) != length(allele)){
    stop("The number of length of peptide for which to make predictions must the paired with alleles")
  }
  length <- paste(as.character(length),collapse = ",")
  allele <- paste(allele,collapse = ",")
  pre_method <- match.arg(pre_method)

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

    for (i in 1:3){
      warning(paste0("Failed retrieving, retrying ",i," times"))
      mess1 <- try(system(command_run))
      if (mess1 == 0){
        message("Succeed !")
        break
      }
    }
    if (j == 10){
      stop("Failed retrieving, stop", immediate. = TRUE)
    }
  }
  res <- read.table(temp_file,header = T)
  file.remove(temp_file)
  return(res)
}

#' Predicting MHC-II binding
#'
#' @param peptide
#' @param allele
#' @param length
#' @param pre_method
#'
#' @return
#' @export
#'
#' @examples
mhcIIbinding <- function(peptide = c("GHAHKVPRRLLKAAR","LKAADASADADGSGSGSGSG"),
                        allele = c("HLA-DRB1*01:01","HLA-A*03:01"),
                        length = c(8,9),
                        pre_method = c("ann","comblib_sidney2008","consensus",
                                       "netmhccons","netmhcpan_ba","netmhcpan_el",
                                       "netmhcstabpan","pickpocket","recommended",
                                       "smm","smmpmbec")){
  if (length(peptide) != 1){
    ##To submit multiple sequences at a time,
    ##escape the special characters in a fasta-formatted sequence with URI codes
    pep_number <- map_chr(c(1:length(peptide)),~gsub("number",.x,"%3Epeptidenumber%0A"))
    pep_con <- c(paste0(peptide[1:length(peptide)-1],"%0A"),peptide[length(peptide)])
    peptide <- paste0(pep_number,pep_con,collapse = "")
  }
  length <- match.arg(as.character(length),
                      choices = c("8","9", "10", "11", "12", "13", "14", "15"),
                      several.ok=T)
  if (length(length) != length(allele)){
    stop("The number of length of peptide for which to make predictions must the paired with alleles")
  }
  length <- paste(as.character(length),collapse = ",")
  allele <- paste(allele,collapse = ",")
  pre_method <- match.arg(pre_method)

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

    for (i in 1:3){
      warning(paste0("Failed retrieving, retrying ",i," times"))
      mess1 <- try(system(command_run))
      if (mess1 == 0){
        message("Succeed !")
        break
      }
    }
    if (j == 10){
      stop("Failed retrieving, stop", immediate. = TRUE)
    }
  }
  res <- read.table(temp_file,header = T)
  file.remove(temp_file)
  return(res)
}
