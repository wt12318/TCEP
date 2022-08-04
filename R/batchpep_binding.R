#' Predict binding affinity of peptides with MHC in batch mode.
#' @param pep_file Character, the path of peptide file needed to be processed, each row refers to a peptide sequence.
#' @param hla_type HLA class, can be HLA-I or HLA-II.
#' @param pep_length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-15, for MHC-II, the length can be 11-30 or asis (take the length of input sequence as the peptide length)
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#' @param client_path The path of local IEDB tools, used when setting get_method as client
#' @param tmp_dir Character, the temp dir
#' @param mhcflurry_env, the installed conda environment of mhcflurry, default is "mhcflurry-env"
#' @param mhcnuggets_env, the installed conda environment of mhcnuggets, default is "mhcnuggets"
#' @param netchop_path, the installed netchop path
#' @param Immuno_IEDB_path, the installed IEDB immunogenicity tool
#' @param Immuno_Deepimmuno_path, the deepimmuno-cnn.py script path
#' @param Deepimmuno_env, the conda envrionment of Deepimmuno
#' @param MixMHCpred_path, the intalled Mixmhcpred path
#' @param PRIME_path, PRIME tool path
#' @param seq2neo_env, the conda env of seq2neo
#' @param seq2neo_path, the of `immuno_Prediction` dir of Seq2Neo
#'
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples
#' test <- batchpep_binding(pep_file=system.file("extdata", "random.pep", package = "TCEP"),
#'                          hla_type = "I",pep_length = c(9),allele = c("HLA-A*01:01"),
#'                          pre_method = "netmhcpan_el",tmp_dir=tempdir(),
#'                          num_thread=1,method_type = "Binding",client_path="~/software/mhc_i/src/")
batchpep_binding <- function(pep_file,hla_type,pep_length,allele,
                             pre_method,method_type,client_path,tmp_dir,
                             mhcflurry_env="mhcflurry-env",mhcnuggets_env="mhcnuggets",
                             netchop_path,Immuno_IEDB_path,Immuno_Deepimmuno_path,Deepimmuno_env,
                             MixMHCpred_path,PRIME_path,seq2neo_env,seq2neo_path){
  pep <- read.table(pep_file)
  pep_seq <- rep(pep$V1,length(pep_length))

  pep_dt <- data.frame(pep_seq=pep_seq,pre_len=rep(pep_length,each=nrow(pep)))
  pep_dt <- pep_dt[which(nchar(pep_dt$pep_seq) >= as.numeric(pep_dt$pre_len)),]

  pep_dt <- pep_dt %>%
    dplyr::group_by(pre_len) %>%
    dplyr::mutate(seq_num=row_number()) %>%
    dplyr::mutate(index=paste(pre_len,seq_num,sep = ":")) %>%
    as.data.frame()

  pep_length <- unique(pep_dt$pre_len)
  pre_res <- vector("list",length = length(pep_length))
  names(pre_res) <- pep_length
  for (i in seq_along(pre_res)){
    pep <- pep_dt[pep_dt$pre_len == names(pre_res)[i],"pep_seq"]
    pre_res[[i]] <- TCEP:::general_mhcbinding(hla_type = hla_type,
                                                    length = pep_length[i],
                                                    allele = allele,pre_method = pre_method,
                                                    method_type=method_type,
                                                    peptide = pep,client_path = client_path,
                                                    tmp_dir=tmp_dir,mhcflurry_type="mt",
                                                    mhcflurry_env=mhcflurry_env,
                                                    mhcnuggets_env=mhcnuggets_env,netchop_path=netchop_path,
                                                    Immuno_IEDB_path=Immuno_IEDB_path,
                                                    Immuno_Deepimmuno_path=Immuno_Deepimmuno_path,
                                                    Deepimmuno_env=Deepimmuno_env,
                                                    MixMHCpred_path=MixMHCpred_path,PRIME_path=PRIME_path,
                                                    seq2neo_env=seq2neo_env,seq2neo_path=seq2neo_path)
  }

  pre_res <- dplyr::bind_rows(pre_res)
  pre_res <- pre_res %>% dplyr::mutate(index=paste(length,seq_num,sep = ":"))
  pep <- left_join(
    pre_res,pep_dt %>% dplyr::select(pep_seq,index)
  ) %>% dplyr::select(-index)
  colnames(pep)[which(colnames(pep) == "pep_seq")] <- "query_pep"
  return(pep)
}

