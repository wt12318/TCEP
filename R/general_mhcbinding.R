#' General MCHbinding
#'
#' @param hla_type HLA class
#' @param client_path The path of local IEDB tools, used when setting get_method as client
#' @param ... Other para passed to mhcIbing or mhcIIbind
#'
#' @return Dataframe
#'
#' @examples
#'test <- TCEP:::general_mhcbinding(hla_type = "I", length = c(10),
#'                                       allele = c("HLA-A*01:01"),pre_method = "Seq2Neo-CNN",
#'                                       method_type="Immuno",
#'                                       peptide = c("NVDTHPGSGK","QTSEKALLRR"),
#'                                       tmp_dir="~/tmp/",netchop_path = "~/software/netchop/",
#'                                       client_path = "~/software/mhc_i/src/",
#'                                       seq2neo_env = "DeepImmuno",
#'                                       seq2neo_path = "~/software/Seq2Neo/seq2neo/function/immuno_Prediction/")

general_mhcbinding <- function(hla_type=c("I","II"),
                               client_path,peptide,allele,length,
                               pre_method,tmp_dir,method_type,mhcflurry_type="mt",
                               mhcflurry_env="mhcflurry-env",
                               mhcnuggets_env="mhcnuggets",netchop_path,
                               Immuno_IEDB_path,Immuno_Deepimmuno_path,Deepimmuno_env,
                               MixMHCpred_path,PRIME_path,
                               seq2neo_env,seq2neo_path){
  hla_type <- match.arg(hla_type)
  if (hla_type == "I"){
    res <- mhcIbinding_client(client_path=client_path,tmp_dir=tmp_dir,
                              peptide=peptide,allele =allele,
                              length =length,pre_method =pre_method,
                              method_type = method_type,mhcflurry_type=mhcflurry_type,
                              mhcflurry_env=mhcflurry_env,mhcnuggets_env=mhcnuggets_env,
                              netchop_path=netchop_path,Immuno_IEDB_path=Immuno_IEDB_path,
                              Immuno_Deepimmuno_path=Immuno_Deepimmuno_path,Deepimmuno_env=Deepimmuno_env,
                              MixMHCpred_path = MixMHCpred_path, PRIME_path = PRIME_path,
                              seq2neo_env=seq2neo_env,seq2neo_path=seq2neo_path)
  }else{
    res <- mhcIIbinding_client(client_path=client_path,tmp_dir=tmp_dir,
                               peptide=peptide,allele =allele,
                               length =length,pre_method =pre_method,
                               method_type = method_type,
                               mhcnuggets_env=mhcnuggets_env)
  }
  return(res)
}
