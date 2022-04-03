#' Convert TXT file to mutated peptides
#'
#' @description This function extract neo-peptides produced from the mutation recorded in the TXT file using the scripts provided by annovar.The txt file is actually
#' a annovar input file, and the first five space- or tab- delimited columns represent chromosome, start position, end position, the reference nucleotides and the observed nucleotides
#' @param annovar_path Character, the install path of annovar.
#' @param txt_path Character, the path of TXT file needed to be processed.
#' @param genome_version Character, which genome build version of the TXT file, can be hg19 or hg38
#' @param tmp_dir Character, the temp dir
#' @details The command run is like: annotate_variation.pl -out out_dir -build hg19 anaovar_input.txt humandb/ -hgvs -dbtype ensGene
#'          coding_change.pl common_driver.exonic_variant_function humandb/hg19_ensGene.txt humandb/hg19_ensGeneMrna.fa -includesnp -onlyAltering --outfile coding.txt

#' @return A dataframe containg wildtype and the mutated protein sequence.
#' @export
#'
#' @examples
#' pep <- txt2pep(annovar_path = annovar_path,txt_path = system.file("extdata", "test_avinput.txt", package = "MHCbinding"),
#'                genome_version = "hg19",tmp_dir=tempdir())
txt2pep <- function(annovar_path,txt_path,
                    genome_version=c("hg19","hg38"),tmp_dir){
  genome_version <- match.arg(genome_version)
  mut <- data.table::fread(txt_path,data.table = F)
  mut <- mut[,1:5]
  colnames(mut) <- c("chr","start","end","ref","alt")

  temp_dir <- tmp_dir
  comm_annotate <- paste0(annovar_path,"/annotate_variation.pl -out ",temp_dir,"/test -build ",genome_version," ",txt_path,
                          " ",annovar_path,"/humandb/ -hgvs -dbtype ensGene")
  system(comm_annotate)
  comm_coding <- paste0(annovar_path,"/coding_change.pl ",temp_dir,"/test.exonic_variant_function ",annovar_path,"/humandb/",
                        genome_version,"_ensGene.txt ",annovar_path,"/humandb/",genome_version,"_ensGeneMrna.fa ",
                        "-includesnp -onlyAltering --alltranscript --outfile ",temp_dir,"/coding.txt")
  system(comm_coding)

  dt <- seqinr::read.fasta(file = paste0(temp_dir,"/coding.txt"),seqtype = "AA",as.string = TRUE,whole.header = T)

  variants <- names(dt)[seq(2,length(dt),by=2)]
  pos_alter <- stringr::str_extract_all(variants,"p[.].+ p") %>% gsub(" p","",.)
  lines <- stringr::str_extract_all(variants,"line.+ ENST") %>% gsub(" ENST","",.) %>%
    gsub("line","",.) %>% as.numeric()
  ENST <- stringr::str_extract_all(variants,"line.+ ENST[0-9]+") %>% gsub("line.+ ","",.)##extract transcripts

  mut_need <- mut[lines,]
  pep_seq_mt <- dt[seq(2,length(dt),by=2)] %>% as.character()
  pep_seq_wt <- dt[seq(1,length(dt),by=2)] %>% as.character()
  mut_need$seq_wt <- pep_seq_wt
  mut_need$seq_mt <- pep_seq_mt
  mut_need$pos_alter <- pos_alter
  mut_need$transcript <- ENST
  files_exist <- c(list.files(temp_dir,pattern = "test",full.names = T),
                   list.files(temp_dir,pattern = "coding.txt",full.names = T))##remove all temp files
  file.remove(files_exist)
  return(mut_need)
}

#' Convert TXT file to mutated peptides and extract proper range of amino acids
#'
#' @description This function use \code{\link{txt2pep}} to convert TXT file to mutated peptides, and use \code{\link{extractSeq}} to extract proper range of amino acids.
#' @param annovar_path Character, the install path of annovar.
#' @param txt_path Character, the path of TXT file needed to be converted.
#' @param genome_version Character, which genome build version of the TXT file, can be hg19 or hg38
#' @param len Numeric, the epitope length needed to predicted.
#' @param tmp_dir Character, the temp dir
#' @return A dataframe containg wildtype and the mutated protein sequence, plus exacted sequence from mutated peptides.
#' @export
#'
#' @examples
#' txt2seq(annovar_path = annovar_path,txt_path = system.file("extdata", "test_avinput.txt", package = "MHCbinding"),
#'         genome_version = "hg19",len = 9,tmp_dir=tempdir())
txt2seq <- function(annovar_path,txt_path,
                    genome_version=c("hg19","hg38"),len,tmp_dir){

  pep <- txt2pep(annovar_path = annovar_path,txt_path = txt_path,
                 genome_version = genome_version,tmp_dir=tmp_dir)

  pep <- pep %>%
    mutate(pos = stringr::str_extract(pos_alter,"[0-9]+") %>% as.numeric()) %>%
    rowwise() %>%
    mutate(indel=ifelse(ref != "-" & alt != "-" & nchar(ref) == 1 & nchar(alt) == 1,"FALSE","TRUE")) %>%
    as.data.frame()
  ext_seqs_mt <- mapply(extractSeq,seq=pep$seq_mt,pos=as.numeric(pep$pos),len=as.numeric(len),indel=pep$indel) %>% unname()
  #ext_seqs_wt <- mapply(extractSeq,seq=pep$seq_wt,pos=pep$pos,len=len,indel=pep$indel) %>% unname()
  ##remove the last star
  a <-  substr(ext_seqs_mt,nchar(ext_seqs_mt),nchar(ext_seqs_mt))=="*"
  #b <-  substr(ext_seqs_wt,nchar(ext_seqs_wt),nchar(ext_seqs_wt))=="*"
  ext_seqs_mt[a] <- gsub("\\*","",ext_seqs_mt[a])
  #ext_seqs_wt[b] <- gsub("\\*","",ext_seqs_wt[b])

  #pep$ext_seqs_wt <- ext_seqs_wt
  pep$ext_seqs_mt <- ext_seqs_mt
  return(pep)
}

#' Predict peptide MHC binding based on TXT file.
#'
#' @description This function use \code{\link{txt2seq}} to exacted mutated "new peptide" from TXT file, and predicted the binding affinity of these peptide with specific MHC allele using \code{\link{general_mhcbinding}}
#' @param get_method The way to predict, can be api or client.
#' @param annovar_path Character, the install path of annovar.
#' @param txt_path Character, the path of TXT file needed to be converted.
#' @param genome_version Character, which genome build version of the VCF file, can be hg19 or hg38
#' @param mhc_type MHC class, can be MHC-I or MHC-II.
#' @param pep_length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-15, for MHC-II, the length can be 11-30 or asis (take the length of input sequence as the peptide length)
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#' @param client_path The path of local IEDB tools, used when setting get_method as client
#' @param tmp_dir Character, the temp dir
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples
#' test <- txt2binding(get_method="api",annovar_path = "~/software/annovar/",txt_path = system.file("extdata", "test_avinput.txt", package = "MHCbinding"),
#'                     genome_version = "hg19",mhc_type = "MHC-I",pep_length = c(9,10),
#'                     allele = c("HLA-A*01:01", "HLA-A*03:01"),pre_method = "ann",tmp_dir=tempdir())
txt2binding <- function(get_method=c("api","client"),annovar_path,txt_path,
                        genome_version=c("hg19","hg38"),mhc_type,pep_length,allele,
                        pre_method,client_path,tmp_dir){

  if (! dir.exists(tmp_dir)){
    dir.create(tmp_dir)
  }
  get_method <- match.arg(get_method)
  res <- vector("list",length = length(pep_length))
  names(res) <- pep_length
  for (i in seq_along(res)){
    res[[i]] <- txt2seq(annovar_path = annovar_path,txt_path = txt_path,genome_version = genome_version,
                        len = pep_length[i],tmp_dir=tmp_dir)
  }
  res <- dplyr::bind_rows(res,.id = "predicted_length")
  res <- res[which(nchar(res$ext_seqs_mt) >= as.numeric(res$predicted_length)),]

  pep_length <- unique(res$predicted_length)
  pre_res <- vector("list",length = length(pep_length))
  names(pre_res) <- pep_length
  for (i in seq_along(pre_res)){
    pep <- res[res$predicted_length == names(pre_res)[i],"ext_seqs_mt"]
    pre_res[[i]] <- MHCbinding:::general_mhcbinding(get_method = get_method,mhc_type = mhc_type, length = pep_length[i],
                                                    allele = allele,pre_method = pre_method, peptide = pep,client_path = client_path,
                                                    tmp_dir=tmp_dir)
  }

  pre_res <- dplyr::bind_rows(pre_res)
  res <- res %>%
    group_by(predicted_length) %>%
    mutate(seq_num=row_number()) %>%
    mutate(index=paste(predicted_length,seq_num,sep = ":"))
  pre_res <- pre_res %>% mutate(index=paste(length,seq_num,sep = ":"))
  pre_res <- left_join(
    pre_res %>% rename(pep_start=start,pep_end=end),
    res %>% ungroup() %>% select(chr,start,end,ref,alt,index,ext_seqs_mt,pos_alter,transcript)
  ) %>% select(-index,-seq_num) %>% select(chr,start,end,ref,alt,ext_seqs_mt,pos_alter,transcript,everything())

  return(pre_res)
}

