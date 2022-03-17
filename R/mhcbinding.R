#' @title Predicting MHC-I binding based on IEDB API
#'
#' @description  This is the wrapped function for IEDB API, the full document can refer to http://tools.iedb.org/main/tools-api/
#' @param peptide A character vector of input protein sequence.
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}.
#' @param length A numeric or character vector, indicating the length for which to make predictions, must be paired with HLA allele. Foe MHC-I, the length can be 8-15
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#'
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples
#' test <- mhcIbinding(peptide = "SLYNTVATLY",
#'                     allele = "HLA-A*01:01",length = "8",
#'                     pre_method = "netmhcpan_el")
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

  #For alleles, a comma-separated list of the alleles for which to make predictions.
  #This list gets paired with the length list, so there must be a corresponding length for each allele.

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
  message("Retrieving data from server ... \n")
  mess <- try(system(command_run))
  if (mess == 0){
    message("Succeed ! \n")
  }else{

    for (i in 1:10){
      warning(paste0("Failed retrieving, retrying ",i," times \n"),immediate. =TRUE)
      mess1 <- try(system(command_run))
      if (mess1 == 0){
        message("Succeed ! \n")
        break
      }
    }
    if (i == 10){
      stop("Failed retrieving, stop \n", immediate. = TRUE)
    }
  }
  res <- read.table(temp_file,header = T)
  file.remove(temp_file)
  return(res)
}

#' @title Predicting MHC-II binding based on IEDB API
#'
#' @description  This is the wrapped function for IEDB API, the full document can refer to http://tools.iedb.org/main/tools-api/
#' @param peptide A character vector of input protein sequence.
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}.
#' @param length A numeric or character vector, indicating the length for which to make predictions. For MHC-II, the length can be 11-30 or asis (take the length of input sequence as the peptide length)
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#'
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples
#' dt <- mhcIIbinding(pre_method = "nn_align")
mhcIIbinding <- function(peptide = c("GHAHKVPRRLLKAAR"),
                        allele = c("DRB1*15:01"),
                        length = 15,
                        pre_method = c("recommended","consensus","netmhciipan",
                                       "nn_align","smm_align","comblib","tepitope")){
  if (length(peptide) != 1){
    ##To submit multiple sequences at a time,
    ##escape the special characters in a fasta-formatted sequence with URI codes
    pep_number <- map_chr(c(1:length(peptide)),~gsub("number",.x,"%3Epeptidenumber%0A"))
    pep_con <- c(paste0(peptide[1:length(peptide)-1],"%0A"),peptide[length(peptide)])
    peptide <- paste0(pep_number,pep_con,collapse = "")
  }
  length <- match.arg(as.character(length),
                      choices = c("11","12","13","14","15","16","17","18","19",
                                  "20","21","22","23","24","25","26","27","28","29","30","asis"),
                      several.ok=T)
  ##MHC-II not need paired alleles and lengths

  length <- paste(as.character(length),collapse = ",")
  allele <- paste(allele,collapse = ",")
  pre_method <- match.arg(pre_method)

  temp_file <- tempfile()
  file.create(temp_file)
  command_run <- paste0('curl --data "method=',pre_method,'&sequence_text=',
                        peptide,'&allele=',allele,'&length=',length,'" ',
                        'http://tools-cluster-interface.iedb.org/tools_api/mhcii/',
                        ' > ',temp_file)
  message("Retrieving data from server ... \n")
  mess <- try(system(command_run))
  if (mess == 0){
    message("Succeed ! \n")
  }else{

    for (i in 1:10){
      warning(paste0("Failed retrieving, retrying ",i," times \n"),immediate. =TRUE)
      mess1 <- try(system(command_run))
      if (mess1 == 0){
        message("Succeed ! \n")
        break
      }
    }
    if (i == 10){
      stop("Failed retrieving, stop \n", immediate. = TRUE)
    }
  }
  res <- read.table(temp_file,header = T,row.names = NULL)
  res$length <- nchar(res$peptide)
  file.remove(temp_file)
  return(res)
}
