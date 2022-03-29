#' Convert MAF file to mutated peptides
#'
#' @description This function extract neo-peptides produced from the mutation recorded in the MAF file using the scripts provided by annovar.
#' @param annovar_path Character, the install path of annovar.
#' @param maf_path Character, the path of MAF file needed to be converted.
#' @param need_allsamples Logical. Whether need all samples when multi-sample MAF file is supplied.
#' @param need_samples Character, the needed sample when multi-sample MAF file is supplied.
#'
#' @details The command run is like: annotate_variation.pl -out out_dir -build hg19 anaovar_input.txt humandb/ -hgvs -dbtype ensGene
#'          coding_change.pl common_driver.exonic_variant_function humandb/hg19_ensGene.txt humandb/hg19_ensGeneMrna.fa -includesnp -onlyAltering --outfile coding.txt

#' @return A dataframe containg wildtype and the mutated protein sequence.
#' @export
#'
#' @examples
#' pep <- maf2pep(annovar_path = annovar_path,maf_path = system.file("extdata", "test.maf", package = "MHCbinding"),
#'                need_allsamples = TRUE)
maf2pep <- function(annovar_path,maf_path,
                    need_allsamples=FALSE,need_samples){
  need_allsamples <- match.arg(as.character(need_allsamples),choices = c("TRUE","FALSE"))
  maf <- data.table::fread(maf_path,data.table = F,skip = "Hugo_Symbol")
  genome_version <- if (unique(maf$NCBI_Build) == "GRCh37"){"hg19"} else {"hg38"}
  if (need_allsamples){
    annovar_input <- maf %>%
      select(Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2,Tumor_Sample_Barcode)
  }else{
    annovar_input <- maf %>%
      select(Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2,Tumor_Sample_Barcode) %>%
      filter(Tumor_Sample_Barcode %in% need_samples)
  }

  temp_file <- tempfile()
  temp_dir <- tempdir()
  file.create(temp_file)

  write.table(annovar_input,file = temp_file,col.names = F,row.names = F,sep = "\t",quote = F)
  comm_annotate <- paste0(annovar_path,"/annotate_variation.pl -out ",temp_dir,"/test -build ",genome_version," ",temp_file,
                          " ",annovar_path,"/humandb/ -hgvs -dbtype ensGene")
  system(comm_annotate)
  comm_coding <- paste0(annovar_path,"/coding_change.pl ",temp_dir,"/test.exonic_variant_function ",annovar_path,"/humandb/",
                        genome_version,"_ensGene.txt ",annovar_path,"/humandb/",genome_version,"_ensGeneMrna.fa ",
                        "-includesnp -onlyAltering --outfile ",temp_dir,"/coding.txt")
  system(comm_coding)

  dt <- seqinr::read.fasta(file = paste0(temp_dir,"/coding.txt"),seqtype = "AA",as.string = TRUE,whole.header = T)

  variants <- names(dt)[seq(2,length(dt),by=2)]
  pos_alter <- stringr::str_extract_all(variants,"p[.].+ p") %>% gsub(" p","",.)
  lines <- stringr::str_extract_all(variants,"line.+ ENST") %>% gsub(" ENST","",.) %>%
    gsub("line","",.) %>% as.numeric()

  mut_need <- annovar_input[lines,]
  pep_seq_mt <- dt[seq(2,length(dt),by=2)] %>% as.character()
  pep_seq_wt <- dt[seq(1,length(dt),by=2)] %>% as.character()
  mut_need$seq_wt <- pep_seq_wt
  mut_need$seq_mt <- pep_seq_mt
  mut_need$pos_alter <- pos_alter
  files_exist <- c(list.files(temp_dir,pattern = gsub(paste0(temp_dir,"/"),"",temp_file),full.names = T),
                   list.files(temp_dir,pattern = "test",full.names = T),
                   list.files(temp_dir,pattern = "coding.txt",full.names = T))##remove all temp files
  file.remove(files_exist)
  return(mut_need)
}


#' Convert MAF file to mutated peptides and extract proper range of amino acids
#'
#' @description This function use \code{\link{maf2pep}} to convert MAF file to mutated peptides, and use \code{\link{extractSeq}} to extract proper range of amino acids.
#' @param annovar_path Character, the install path of annovar.
#' @param maf_path Character, the path of MAF file needed to be converted.
#' @param need_allsamples Logical. Whether need all samples when multi-sample MAF file is supplied.
#' @param need_samples Character, the needed sample when multi-sample MAF file is supplied.
#' @param len Numeric, the epitope length needed to predicted.
#'
#' @return A dataframe containg wildtype and the mutated protein sequence, plus exacted sequence from mutated peptides.
#' @export
#'
#' @examples
#' maf2seq(annovar_path = annovar_path,maf_path = system.file("extdata", "test.maf", package = "MHCbinding"),
#'         need_allsamples = TRUE,len = 9)
maf2seq <- function(annovar_path,maf_path,need_allsamples=TRUE,need_samples,len){

  pep <- maf2pep(annovar_path = annovar_path,maf_path = maf_path,
                 need_allsamples = need_allsamples,
                 need_samples = need_samples)

  colnames(pep)[1:6] <- c("chr","start","end","ref","alt","tumor_sample_barcode")[1:6]
  pep <- pep %>%
    mutate(pos = stringr::str_extract(pos_alter,"[0-9]+") %>% as.numeric()) %>%
    rowwise() %>%
    mutate(indel=ifelse(ref != "-" & alt != "-" & nchar(ref) == 1 & nchar(alt) == 1,"FALSE","TRUE")) %>%
    as.data.frame()
  ext_seqs_mt <- mapply(extractSeq,seq=pep$seq_mt,pos=as.numeric(pep$pos),len=as.numeric(len),indel=pep$indel) %>% unname()
  ##remove the last star
  a <-  substr(ext_seqs_mt,nchar(ext_seqs_mt),nchar(ext_seqs_mt))=="*"
  ext_seqs_mt[a] <- gsub("\\*","",ext_seqs_mt[a])
  pep$ext_seqs_mt <- ext_seqs_mt
  return(pep)
}

#' Predict peptide MHC binding based on MAF file.
#'
#' @description This function use \code{\link{maf2seq}} to exacted mutated "new peptide" from MAF file, and predicted the binding affinity of these peptide with specific MHC allele using \code{\link{general_mhcbinding}}
#' @param get_method The way to predict, can be api or client.
#' @param annovar_path Character, the install path of annovar.
#' @param maf_path Character, the path of MAF file needed to be converted.
#' @param need_allsamples Logical. Whether need all samples when multi-sample MAF file is supplied.
#' @param need_samples Character, the needed sample when multi-sample MAF file is supplied.
#' @param mhc_type MHC class, can be MHC-I or MHC-II.
#' @param pep_length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-15, for MHC-II, the length can be 11-30 or asis (take the length of input sequence as the peptide length)
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#' @param client_path The path of local IEDB tools, used when setting get_method as client
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples
#' test <- maf2binding(get_method="api",annovar_path = "~/software/annovar/",maf_path = system.file("extdata", "test.maf", package = "MHCbinding"),
#'                     need_allsamples = TRUE,mhc_type = "MHC-I",pep_length = c(9,10),
#'                     allele = c("HLA-A*01:01", "HLA-A*03:01"),pre_method = "ann")

maf2binding <- function(get_method=c("api","client"),annovar_path,maf_path,need_allsamples=FALSE,need_samples,
                        mhc_type,pep_length,allele,pre_method,client_path){
  get_method <- match.arg(get_method)
  res <- vector("list",length = length(pep_length))
  names(res) <- pep_length
  for (i in seq_along(res)){
    res[[i]] <- maf2seq(annovar_path = annovar_path,maf_path = maf_path,
                        need_allsamples = need_allsamples,need_samples = need_samples,len = pep_length[i])
  }
  res <- dplyr::bind_rows(res,.id = "predicted_length")
  res <- res[which(nchar(res$ext_seqs_mt) >= as.numeric(res$predicted_length)),]

  pre_res <- vector("list",length = length(pep_length))
  names(pre_res) <- pep_length
  for (i in seq_along(pre_res)){
    pep <- res[res$predicted_length == names(pre_res)[i],"ext_seqs_mt"]
    pre_res[[i]] <- MHCbinding:::general_mhcbinding(get_method = get_method,mhc_type = mhc_type, length = pep_length[i],
                                                    allele = allele,pre_method = pre_method, peptide = pep,client_path = client_path)
  }

  pre_res <- dplyr::bind_rows(pre_res)
  res <- res %>%
    group_by(predicted_length) %>%
    mutate(seq_num=row_number()) %>%
    mutate(index=paste(predicted_length,seq_num,sep = ":"))
  pre_res <- pre_res %>% mutate(index=paste(length,seq_num,sep = ":"))
  pre_res <- left_join(
    pre_res %>% rename(pep_start=start,pep_end=end),
    res %>% ungroup() %>% select(chr,start,end,ref,alt,index,ext_seqs_mt)
  )%>% select(-index)

  return(pre_res)
}
