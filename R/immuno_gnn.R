#' Predict immunogenicity score of peptide using Immuno-GNN model
#'
#' @param input_type, multiple mode for inputing file containing peptide and hla columns, or single mode for
#'                    one pair of peptide-hla
#' @param input_file, csv file path for multiple mode of input_type
#' @param pep, peptides for single mode of input_type
#' @param hla, single hla allele for single mode of input_type
#' @param immuno_gnn_env, the conda environment of Immuno-GNN
#' @param immuno_gnn_path, the dir path of Immuno_gnn.py, aaindex1_pca.csv , pMHCDataset and last_model.pth
#' @param num_thread, the thread should be used for prediction
#' @param temp_dir, the path of temp dir
#' @return A dataframe containing pep, hla and predicted immunogenicity
#' @export
#'
#' @examples immuno_gnn(pep = c("GHAHKVPRRLLKAA","SLYNTVATLY"),
#'                      hla = "HLA-A01:01", immuno_gnn_env="dl_gpu",
#'                      immuno_gnn_path="./data-raw/",num_thread=20,temp_dir=tempdir())
immuno_gnn <- function(input_type="single", input_file, pep, hla,
                       immuno_gnn_env, immuno_gnn_path, num_thread, temp_dir){
  input_type <- match.arg(arg = input_type, choices = c("single","multiple"))
  pseudu <- TCEP::pseudu
  if (input_type == "multiple"){
    pep_file <- read.csv(input_file)
  }else{
    pep_file <- data.frame(pep=pep,HLA_Allele=hla)
  }
  pep_file <- pep_file %>%
    rowwise() %>%
    mutate(hla_seq=unique(pseudu$V2[pseudu$V1==HLA_Allele])) %>%
    ungroup() %>%
    mutate(type=0)
  write.csv(pep_file,file = paste0(temp_dir,"/pep_file.csv"),quote = FALSE,row.names = FALSE)
  system(paste0("conda run -n ",immuno_gnn_env," python ",
                normalizePath(immuno_gnn_path),"/Immuno_gnn.py ",
                "-a ",normalizePath(immuno_gnn_path),"/aaindex1_pca.csv ",
                "-m ",normalizePath(immuno_gnn_path),"/last_model.pth ",
                "-t ",num_thread," -i ",paste0(temp_dir,"/pep_file.csv"),
                " -o ",temp_dir))
  res <- read.csv(paste0(temp_dir,"/pred_res.csv"))
  res <- res %>% select(-X,-label,-pred,-hla_seq,-type)
  colnames(res)[3] <- "immunogenicity"
  return(res)
}

