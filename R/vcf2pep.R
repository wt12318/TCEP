

##convert2annovar.pl -format vcf4 example/ex2.vcf -outfile ex2.avinput -allsample -withfreq -include
vcf2annova <- function(annovar_path,vcf_path,out_file,
                      need_allsamples=TRUE,need_samples){

  need_allsamples <- match.arg(as.character(need_allsamples),choices = c("TRUE","FALSE"))

  if (need_allsamples){
    commond <- paste0(annovar_path,"/convert2annovar.pl -format vcf4 ",vcf_path," -outfile ",out_file," -allsample -withfreq")
    system(command = commond)
  }else{
    commond <- paste0(annovar_path,"/convert2annovar.pl -format vcf4 ",vcf_path," -outfile ",out_file," -allsample")
    system(command = commond)
    need_samples_file <- paste0("cat ",out_file,".",need_samples,".avinput > ",out_file)
    system(need_samples_file)
  }
}



#annotate_variation.pl -out out_dir -build hg19 anaovar_input.txt humandb/ -hgvs -dbtype ensGene
#coding_change.pl common_driver.exonic_variant_function humandb/hg19_ensGene.txt humandb/hg19_ensGeneMrna.fa -includesnp -onlyAltering --outfile coding.txt
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
  ext_seqs <- mapply(extractSeq,seq=pep$seq_mt,pos=pep$pos,len=len,indel=pep$indel) %>% unname()

  ##remove the last star
  a <-  substr(ext_seqs,nchar(ext_seqs),nchar(ext_seqs))=="*"
  ext_seqs[a] <- gsub("\\*","",ext_seqs[a])

  pep$ext_seq <- ext_seqs
  return(pep)
}

vcf2binding <- function(annovar_path,vcf_path,
                        genome_version=c("hg19","hg38"),need_allsamples=FALSE,need_samples,mhc_type,pep_length,allele,
                        pre_method){
  res <- vector("list",length = length(pep_length))
  names(res) <- pep_length
  for (i in seq_along(res)){
    res[[i]] <- vcf2seq(annovar_path = annovar_path,vcf_path = vcf_path,genome_version = genome_version,
                        need_allsamples = need_allsamples,need_samples = need_samples,len = pep_length[i])
  }
  res <- dplyr::bind_rows(res,.id = "predicted_length")

  pre_res <- vector("list",length = length(pep_length))
  names(pre_res) <- pep_length
  for (i in seq_along(pre_res)){
    pep <- res[res$predicted_length == names(pre_res)[i],"ext_seq"]
    pre_res[[i]] <- MHCbinding:::general_mhcbinding(mhc_type = mhc_type, length = pep_length[i],
                                                    allele = allele,pre_method = pre_method, peptide = pep)
  }

  tt <- dplyr::bind_rows(pre_res)
  res <- res %>% mutate(index = paste(predicted_length,map(table(res$predicted_length) %>% unname,~seq(1,.x)) %>% unlist(),sep = ":"))
  tt <- tt %>% mutate(index=paste(length,seq_num,sep = ":"))
  tt <- left_join(
    tt %>% rename(pep_start=start,pep_end=end),
    res %>% select(chr,start,end,ref,alt,index)
  ) %>% select(-index)

}

annovar_path = "~/software/annovar/"
vcf_path = "inst/extdata/test_grch38.vcf"
genome_version = "hg38"
need_allsamples=FALSE
need_samples = "TUMOR"
mhc_type = "MHC-I"
pep_length = c(9,10)
allele = c("HLA-A*01:01", "HLA-A*03:01")
pre_method = "ann"

