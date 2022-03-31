
#' Convert vcf file format into input file format of annovar.
#'
#' @description This function called convert2annovar.pl script provided by \href{https://annovar.openbioinformatics.org/en/latest/}{Annovar} to convert VCF to input file of annovar.
#' @param annovar_path Character, the install path of annovar.
#' @param vcf_path Character, the path of VCF file needed to be converted.
#' @param out_file Character, the path of converted file.
#' @param need_allsamples Logical. Whether need all samples when multi-sample VCF file is supplied.
#' @param need_samples Character, the needed sample when multi-sample VCF file is supplied.
#' @details The command run is like :convert2annovar.pl -format vcf4 example/ex2.vcf -outfile ex2.avinput -allsample -withfreq -include
#' @return NULL
#' @export
#'
#' @examples
#' temp_file <- tempfile()
#' file.create(temp_file)
#' vcf2annova(annovar_path = annovar_path,vcf_path = system.file("extdata", "test_grch38.vcf", package = "MHCbinding"),
#'            out_file = temp_file,need_allsamples = FALSE,
#'            need_samples = "TUMOR")
#'
vcf2annova <- function(annovar_path,vcf_path,out_file,
                      need_allsamples=TRUE,need_samples){

  need_allsamples <- match.arg(as.character(need_allsamples),choices = c("TRUE","FALSE"))

  if (need_allsamples){
    commond <- paste0(annovar_path,"/convert2annovar.pl -format vcf4 ",vcf_path," -outfile ",out_file," -allsample -withfreq")
    system(command = commond)
  }else{
    commond <- paste0(annovar_path,"/convert2annovar.pl -format vcf4 ",vcf_path," -outfile ",out_file," -allsample")
    system(command = commond)
    need_samples_file <- paste0("cat ",out_file,".",need_samples,".avinput >> ",out_file)
    system(need_samples_file)
  }
}


#' Convert VCF file to mutated peptides
#'
#' @description This function extract neo-peptides produced from the mutation recorded in the VCF file using the scripts provided by annovar.
#' @param annovar_path Character, the install path of annovar.
#' @param vcf_path Character, the path of VCF file needed to be converted.
#' @param genome_version Character, which genome build version of the VCF file, can be hg19 or hg38
#' @param need_allsamples Logical. Whether need all samples when multi-sample VCF file is supplied.
#' @param need_samples Character, the needed sample when multi-sample VCF file is supplied.
#'
#' @details The command run is like: annotate_variation.pl -out out_dir -build hg19 anaovar_input.txt humandb/ -hgvs -dbtype ensGene
#'          coding_change.pl common_driver.exonic_variant_function humandb/hg19_ensGene.txt humandb/hg19_ensGeneMrna.fa -includesnp -onlyAltering --outfile coding.txt

#' @return A dataframe containg wildtype and the mutated protein sequence.
#' @export
#'
#' @examples
#' pep <- vcf2pep(annovar_path = annovar_path,vcf_path = system.file("extdata", "test_grch38.vcf", package = "MHCbinding"),
#'                genome_version = "hg38",need_allsamples = FALSE,
#'                need_samples = "TUMOR")
vcf2pep <- function(annovar_path,vcf_path,
                    genome_version=c("hg19","hg38"),need_allsamples=TRUE,need_samples){
  need_allsamples <- match.arg(as.character(need_allsamples),choices = c("TRUE","FALSE"))
  genome_version <- match.arg(genome_version)
  temp_file <- tempfile()
  temp_dir <- tempdir()
  file.create(temp_file)
  vcf2annova(annovar_path = annovar_path,vcf_path = vcf_path,out_file = temp_file,need_allsamples = need_allsamples,
             need_samples = need_samples)

  mut <- read.table(temp_file)
  mut <- mut[,1:5]
  colnames(mut) <- c("chr","start","end","ref","alt")
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

  mut_need <- mut[lines,]
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

#' Extract sepcific sub-squence of peptide
#'
#' @description Given a peptide, position of mutated amino acid and needed epitope length, exacted the proper range of amino acids that used for silde window prediction of epitope.
#' @param seq Character, the peptide sequence.
#' @param pos Numeric, the position of mutated amino acid.
#' @param len Numeric, the epitope length needed to predicted.
#' @param indel Logical, whether the type of mutation is indel.
#'
#' @return Character, the exacted sub-peptide sequence
#' @export
#'
#' @examples
#' extractSeq("GHAHKVPRRLLKAAR",pos=2,len=8,indel=FALSE)
#'
extractSeq <- function(seq,pos,len,indel=FALSE){
  indel <- match.arg(as.character(indel),choices = c("TRUE","FALSE"))
  if (pos < len){
    ext_seq <- substr(seq,1,(pos+len-1))
  }else if ((nchar(seq)-pos)<(len-1)){
    ext_seq <- substr(seq,(pos-len+1),nchar(seq))
  }else{
    ext_seq <- substr(seq,(pos-len+1),(pos+len-1))
  }

  if (indel){
    if (pos < len){
      ext_seq <- substr(seq,1,nchar(seq))
    }else{
      ext_seq <- substr(seq,(pos-len+1),nchar(seq))
    }
  }

  return(ext_seq)
}

#' Convert VCF file to mutated peptides and extract proper range of amino acids
#'
#' @description This function use \code{\link{vcf2pep}} to convert VCF file to mutated peptides, and use \code{\link{extractSeq}} to extract proper range of amino acids.
#' @param annovar_path Character, the install path of annovar.
#' @param vcf_path Character, the path of VCF file needed to be converted.
#' @param genome_version Character, which genome build version of the VCF file, can be hg19 or hg38
#' @param need_allsamples Logical. Whether need all samples when multi-sample VCF file is supplied.
#' @param need_samples Character, the needed sample when multi-sample VCF file is supplied.
#' @param len Numeric, the epitope length needed to predicted.
#'
#' @return A dataframe containg wildtype and the mutated protein sequence, plus exacted sequence from mutated peptides.
#' @export
#'
#' @examples
#' vcf2seq(annovar_path = annovar_path,vcf_path = system.file("extdata", "test_grch38.vcf", package = "MHCbinding"),
#'         genome_version = "hg38",
#'         need_allsamples = FALSE,need_samples = "TUMOR",len = 9)
vcf2seq <- function(annovar_path,vcf_path,
                    genome_version=c("hg19","hg38"),need_allsamples=TRUE,need_samples,len){

  pep <- vcf2pep(annovar_path = annovar_path,vcf_path = vcf_path,
                 genome_version = genome_version,need_allsamples = need_allsamples,
                 need_samples = need_samples)

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

#' Predict peptide MHC binding based on VCF file.
#'
#' @description This function use \code{\link{vcf2seq}} to exacted mutated "new peptide" from VCF file, and predicted the binding affinity of these peptide with specific MHC allele using \code{\link{general_mhcbinding}}
#' @param get_method The way to predict, can be api or client.
#' @param annovar_path Character, the install path of annovar.
#' @param vcf_path Character, the path of VCF file needed to be converted.
#' @param genome_version Character, which genome build version of the VCF file, can be hg19 or hg38
#' @param need_allsamples Logical. Whether need all samples when multi-sample VCF file is supplied.
#' @param need_samples Character, the needed sample when multi-sample VCF file is supplied.
#' @param mhc_type MHC class, can be MHC-I or MHC-II.
#' @param pep_length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-15, for MHC-II, the length can be 11-30 or asis (take the length of input sequence as the peptide length)
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}
#' @param pre_method Character, indicating the prediction method. Available methods for MHC-I or MHC-II can be obtained by \code{\link{available_methods}}
#' @param client_path The path of local IEDB tools, used when setting get_method as client
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples
#' test <- vcf2binding(get_method="api",annovar_path = "~/software/annovar/",vcf_path = system.file("extdata", "test_grch38.vcf", package = "MHCbinding"),
#'                     genome_version = "hg38",need_allsamples = FALSE,need_samples = "TUMOR",mhc_type = "MHC-I",pep_length = c(9,10),
#'                     allele = c("HLA-A*01:01", "HLA-A*03:01"),pre_method = "ann")
vcf2binding <- function(get_method=c("api","client"),annovar_path,vcf_path,
                        genome_version=c("hg19","hg38"),need_allsamples=FALSE,need_samples,mhc_type,pep_length,allele,
                        pre_method,client_path){
  get_method <- match.arg(get_method)
  res <- vector("list",length = length(pep_length))
  names(res) <- pep_length
  for (i in seq_along(res)){
    res[[i]] <- vcf2seq(annovar_path = annovar_path,vcf_path = vcf_path,genome_version = genome_version,
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
    res %>% ungroup() %>%  select(chr,start,end,ref,alt,index,ext_seqs_mt)
  ) %>% select(-index)
  return(pre_res)
}

