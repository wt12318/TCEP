#' Convert TXT file to mutated peptides
#'
#' @description This function extract neo-peptides produced from the mutation recorded in the TXT file using the scripts provided by annovar.The txt file is actually
#' a annovar input file, and the first five space- or tab- delimited columns represent chromosome, start position, end position, the reference nucleotides and the observed nucleotides
#' @param annovar_path Character, the install path of annovar.
#' @param txt_path Character, the path of TXT file needed to be processed.
#' @param genome_version Character, which genome build version of the TXT file, can be hg19 or hg38
#' @param tmp_dir Character, the temp dir
#' @param num_thread specify the number of threads to be used in annotation
#' @details The command run is like: annotate_variation.pl -out out_dir -build hg19 anaovar_input.txt humandb/ -hgvs -dbtype ensGene
#'          coding_change.pl common_driver.exonic_variant_function humandb/hg19_ensGene.txt humandb/hg19_ensGeneMrna.fa -includesnp -onlyAltering --outfile coding.txt

#' @return A dataframe containg wildtype and the mutated protein sequence.
#' @export
#'
#' @examples
#' pep <- txt2pep(annovar_path = annovar_path,txt_path = system.file("extdata", "test_avinput.txt", package = "TCEP"),
#'                genome_version = "hg19",tmp_dir=tempdir(),num_thread=1)
txt2pep <- function(annovar_path,txt_path,
                    genome_version=c("hg19","hg38","mm10","mm9"),tmp_dir,num_thread){
  genome_version <- match.arg(genome_version)
  mut <- data.table::fread(txt_path,data.table = F)
  mut <- mut[,1:5]
  colnames(mut) <- c("chr","start","end","ref","alt")

  temp_dir <- tmp_dir
  db <- ifelse(grepl("mm",genome_version),"/mm10db/","/humandb/")
  comm_annotate <- paste0(annovar_path,"/table_annovar.pl ",txt_path," ",annovar_path,db,
                          " -build ",genome_version,
                          " --outfile ",temp_dir,"/myanno",
                          ' -protocol refGene -operation g  --codingarg "-includesnp -onlyAltering"',
                          " --thread ",num_thread)
  system(comm_annotate)

  dt <- try(seqinr::read.fasta(file = paste0(temp_dir,"/myanno.refGene.fa"),seqtype = "AA",as.string = TRUE,
                               whole.header = T))
  if ('try-error' %in% class(dt)){
    return(NULL)
  }
  variants <- names(dt)[seq(2,length(dt),by=2)]
  variants <- variants[grepl("protein-altering",variants)]
  if (length(variants) == 0){
    return(NULL)
  }else{
    pos_alter <- stringr::str_extract_all(variants,"p[.].+ protein-altering") %>% gsub(" protein-altering","",.)
    cdna <- stringr::str_extract_all(variants,"c[.].+ p[.]") %>% gsub(" p[.]","",.)
    lines <- stringr::str_extract_all(variants,"line.+ NM_") %>% gsub(" NM_","",.) %>%
      gsub("line","",.) %>% as.numeric()
    ENST <- stringr::str_extract_all(variants,"line.+ NM_[0-9]+") %>% gsub("line.+ ","",.)##extract transcripts

    mut_need <- mut[lines,]
    pep_seq_mt <- dt[sapply(variants,function(x){which(names(dt)==x)}) %>% unname()] %>% as.character()
    pep_seq_wt <- dt[(sapply(variants,function(x){which(names(dt)==x)}) %>% unname())-1] %>% as.character()
    mut_need$seq_wt <- pep_seq_wt
    mut_need$seq_mt <- pep_seq_mt
    mut_need$pos_alter <- pos_alter
    mut_need$cdna <- cdna
    mut_need$transcript <- ENST
    files_exist <- c(list.files(temp_dir,pattern = "myanno",full.names = T))##remove all temp files
    file.remove(files_exist)
    return(mut_need)
  }
}

#' Convert TXT file to mutated peptides and extract proper range of amino acids
#'
#' @description This function use \code{\link{txt2pep}} to convert TXT file to mutated peptides, and use \code{\link{extractSeq}} to extract proper range of amino acids.
#' @param annovar_path Character, the install path of annovar.
#' @param txt_path Character, the path of TXT file needed to be converted.
#' @param genome_version Character, which genome build version of the TXT file, can be hg19 or hg38
#' @param len Numeric, the epitope length needed to predicted.
#' @param tmp_dir Character, the temp dir
#' @param num_thread specify the number of threads to be used in annotation
#' @return A dataframe containg wildtype and the mutated protein sequence, plus exacted sequence from mutated peptides.
#' @export
#'
#' @examples
#' txt2seq(annovar_path = annovar_path,txt_path = system.file("extdata", "test_avinput.txt", package = "TCEP"),
#'         genome_version = "hg19",len = 9,tmp_dir=tempdir(),num_thread=1)
txt2seq <- function(annovar_path,txt_path,
                    genome_version=c("hg19","hg38"),len,tmp_dir,num_thread){

  pep <- txt2pep(annovar_path = annovar_path,txt_path = txt_path,
                 genome_version = genome_version,tmp_dir=tmp_dir,num_thread=num_thread)

  if (is.null(pep)){
    return(NULL)
  }else{
    pep <- pep %>%
      dplyr::mutate(pos = stringr::str_extract(pos_alter,"[0-9]+") %>% as.numeric()) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(indel=ifelse(ref != "-" & alt != "-" & nchar(ref) == 1 & nchar(alt) == 1,"FALSE","TRUE")) %>%
      as.data.frame()
    ext_seqs_mt <- mapply(TCEP::extractSeq,seq=pep$seq_mt,pos=as.numeric(pep$pos),
                          len=as.numeric(len),indel=pep$indel) %>% unname()
    ext_seqs_wt <- mapply(TCEP::extractSeq,seq=pep$seq_wt,pos=as.numeric(pep$pos),
                          len=as.numeric(len),indel=pep$indel) %>% unname()
    ##remove the last star
    a <-  substr(ext_seqs_mt,nchar(ext_seqs_mt),nchar(ext_seqs_mt))=="*"
    b <-  substr(ext_seqs_wt,nchar(ext_seqs_wt),nchar(ext_seqs_wt))=="*"
    ext_seqs_mt[a] <- gsub("\\*","",ext_seqs_mt[a])
    ext_seqs_wt[b] <- gsub("\\*","",ext_seqs_wt[b])

    pep$ext_seqs_wt <- ext_seqs_wt
    pep$ext_seqs_mt <- ext_seqs_mt
    return(pep)
  }
}

#' Predict peptide MHC binding based on TXT file.
#'
#' @description This function use \code{\link{txt2seq}} to exacted mutated "new peptide" from TXT file, and predicted the binding affinity of these peptide with specific MHC allele using \code{\link{general_mhcbinding}}
#' @param annovar_path Character, the install path of annovar.
#' @param txt_path Character, the path of TXT file needed to be converted.
#' @param genome_version Character, which genome build version of the VCF file, can be hg19 or hg38
#' @param hla_type HLA class, can be I or II.
#' @param pep_length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-15, for MHC-II, the length can be 11-30 or asis (take the length of input sequence as the peptide length)
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}
#' @param pre_method Character, indicating the prediction method. Available methods for HLA-I or HLA-II can be obtained by \code{\link{available_methods}}
#' @param method_type, which type prediction method used, could be "Binding", "Processing" or "Immuno"
#' @param client_path The path of local IEDB tools
#' @param tmp_dir Character, the temp dir
#' @param num_thread specify the number of threads to be used in annotation
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
#' @return A dataframe contains the predicted IC50 and precentile rank (if available).
#' @export
#'
#' @examples
#' test <- txt2binding(annovar_path = "~/software/annovar/annovar/",
#'                     txt_path = system.file("extdata", "test_avinput.txt", package = "TCEP"),
#'                     genome_version = "hg19",hla_type = "I",pep_length = c(9),
#'                     allele = c("HLA-A*01:01"),pre_method = "netmhcpan_el",tmp_dir=tempdir(),
#'                     num_thread=1,method_type = "Binding",client_path="~/software/mhc_i/src/")
txt2binding <- function(annovar_path,txt_path,
                        genome_version=c("hg19","hg38","mm10","mm9"),hla_type,pep_length,allele,
                        pre_method,method_type,client_path,tmp_dir,num_thread,
                        mhcflurry_env="mhcflurry-env",mhcnuggets_env="mhcnuggets",
                        netchop_path,Immuno_IEDB_path,Immuno_Deepimmuno_path,Deepimmuno_env,
                        MixMHCpred_path,PRIME_path,seq2neo_env,seq2neo_path){

  if (! dir.exists(tmp_dir)){
    dir.create(tmp_dir,recursive = TRUE)
  }
  res <- vector("list",length = length(pep_length))
  names(res) <- pep_length
  for (i in seq_along(res)){
    res[[i]] <- txt2seq(annovar_path = annovar_path,txt_path = txt_path,
                        genome_version = genome_version,
                        len = pep_length[i],tmp_dir=tmp_dir,num_thread=num_thread)
  }
  res <- dplyr::bind_rows(res,.id = "predicted_length")
  if(nrow(res)==0){
    message("There is no peptides!")
    return(NULL)
  }else{
    res <- res[which(nchar(res$ext_seqs_mt) >= as.numeric(res$predicted_length)),]
    res <- res %>%
      dplyr::select(-seq_mt,-seq_wt) %>%
      dplyr::group_by(predicted_length,chr,start,end,ref,alt,indel,ext_seqs_mt,ext_seqs_wt) %>%
      summarise(cdna=paste(cdna,collapse = ","),
                transcript=paste(transcript,collapse = ","),
                pos_alter=paste(pos_alter,collapse = ",")) %>%
      ungroup()

    pep_length <- unique(res$predicted_length)
    pre_res_mt <- vector("list",length = length(pep_length))
    names(pre_res_mt) <- pep_length
    for (i in seq_along(pre_res_mt)){
      pep_mt <- res[res$predicted_length == names(pre_res_mt)[i],"ext_seqs_mt"]
      pre_res_mt[[i]] <- TCEP:::general_mhcbinding(hla_type = hla_type, length = names(pre_res_mt)[i],
                                                         allele = allele,pre_method = pre_method,method_type=method_type,
                                                         peptide = pep_mt$ext_seqs_mt,client_path = client_path,
                                                         tmp_dir=tmp_dir,mhcflurry_type="mt",
                                                         mhcflurry_env=mhcflurry_env,
                                                         mhcnuggets_env=mhcnuggets_env,netchop_path=netchop_path,
                                                         Immuno_IEDB_path=Immuno_IEDB_path,
                                                         Immuno_Deepimmuno_path=Immuno_Deepimmuno_path,
                                                         Deepimmuno_env=Deepimmuno_env,
                                                         MixMHCpred_path=MixMHCpred_path,PRIME_path=PRIME_path,
                                                         seq2neo_env=seq2neo_env,seq2neo_path=seq2neo_path)
    }

    pre_res_mt <- dplyr::bind_rows(pre_res_mt)

    res <- res %>%
      dplyr::group_by(predicted_length) %>%
      dplyr::mutate(seq_num=row_number()) %>%
      dplyr::mutate(index=paste(predicted_length,seq_num,sep = ":"))
    pre_res_mt <- pre_res_mt %>% dplyr::mutate(index=paste(length,seq_num,sep = ":"))
    pre_res_mt <- left_join(
      pre_res_mt %>% dplyr::rename(pep_start=start,pep_end=end),
      res %>% dplyr::ungroup()
      %>% dplyr::select(chr,start,end,ref,alt,index,ext_seqs_mt,ext_seqs_wt,pos_alter,cdna,transcript)
    ) %>%
      dplyr::select(-index,-seq_num) %>%
      dplyr::select(chr,start,end,ref,alt,
                    ext_seqs_mt,ext_seqs_wt,pos_alter,cdna,transcript,everything())
    pre_res_mt$peptide_wt <- substr(pre_res_mt$ext_seqs_wt,pre_res_mt$pep_start,pre_res_mt$pep_end)

    pre_res_wt <- vector("list",length = length(unique(pre_res_mt$length)))
    names(pre_res_wt) <- unique(pre_res_mt$length)
    for (i in seq_along(pre_res_wt)){
      wt_pep_dt <- pre_res_mt %>%
        filter(length == names(pre_res_wt)[i])
      pre_res_wt[[i]] <- TCEP:::general_mhcbinding(hla_type = hla_type, length = names(pre_res_wt)[i],
                                                         allele = unique(wt_pep_dt$allele),
                                                         pre_method = pre_method,
                                                         peptide = unique(wt_pep_dt$peptide_wt),
                                                         client_path = client_path,
                                                         method_type=method_type,
                                                         tmp_dir=tmp_dir,mhcflurry_type="wt",
                                                         mhcflurry_env=mhcflurry_env,
                                                         mhcnuggets_env=mhcnuggets_env,
                                                         netchop_path=netchop_path,Immuno_IEDB_path=Immuno_IEDB_path,
                                                         Immuno_Deepimmuno_path=Immuno_Deepimmuno_path,
                                                         Deepimmuno_env=Deepimmuno_env,
                                                         MixMHCpred_path=MixMHCpred_path,PRIME_path=PRIME_path,
                                                         seq2neo_env=seq2neo_env,seq2neo_path=seq2neo_path)
    }
    pre_res_wt <- bind_rows(pre_res_wt)

    if (hla_type == "I"){
      pre_res_wt <- pre_res_wt %>%
        dplyr::select(allele,peptide,res_cols[pre_method][[1]]) %>%
        dplyr::distinct_all(.keep_all = T) %>%
        dplyr::mutate(index=paste(allele,peptide,sep = ":")) %>%
        dplyr::select(index,res_cols[pre_method][[1]]) %>%
        dplyr::rename_with(function(x){paste0("wt_",x)},res_cols[pre_method][[1]])
    }else{
      pre_res_wt <- pre_res_wt %>%
        dplyr::select(allele,peptide,res_cols_ii[pre_method][[1]]) %>%
        dplyr::distinct_all(.keep_all = T) %>%
        dplyr::mutate(index=paste(allele,peptide,sep = ":")) %>%
        dplyr::select(index,res_cols_ii[pre_method][[1]]) %>%
        dplyr::rename_with(function(x){paste0("wt_",x)},res_cols_ii[pre_method][[1]])
    }

    pre_res <- left_join(
      pre_res_mt %>% dplyr::mutate(index=paste(allele,peptide_wt,sep = ":")),
      pre_res_wt
    ) %>% dplyr::select(-index)
    return(pre_res)
  }
}

