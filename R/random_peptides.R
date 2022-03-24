#' Generate random peptide sequences
#'
#' @description This function samples amino acid to form random peptide sequences from 20 common amino acids.
#' @param len, the length of peptide
#' @param n, the number needed peptides
#'
#' @return A dataframe with one column containing the generated random peptides.
#' @export
#'
#' @examples random_peptides(len=9,n=10)
random_peptides <- function(len,n){
  all_code <- "GALMFWKQESPVICYHRNDT"
  all_code <- strsplit(all_code,split = "") %>% unlist()

  dt <- data.frame(index=1:n,pep=NA)
  dt$pep <- sapply(c(1:nrow(dt)),function(x){sample(all_code,len) %>% paste0(.,collapse = "")})
  dt <- dt %>% select(-index)
}
