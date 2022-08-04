#' Convert MAF file to mutated peptides
#'
#' @description This function extract neo-peptides produced from the mutation recorded in the MAF file using the scripts provided by annovar.
#' @param annovar_path Character, the install path of annovar.
#' @param maf_path Character, the path of MAF file needed to be converted.
#' @param need_allsamples Logical. Whether need all samples when multi-sample MAF file is supplied.
#' @param need_samples Character, the needed sample when multi-sample MAF file is supplied.
#' @param num_thread specify the number of threads to be used in annotation
#' @details The command run is like: annotate_variation.pl -out out_dir -build hg19 anaovar_input.txt humandb/ -hgvs -dbtype ensGene
#'          coding_change.pl common_driver.exonic_variant_function humandb/hg19_ensGene.txt humandb/hg19_ensGeneMrna.fa -includesnp -onlyAltering --outfile coding.txt

#' @return A dataframe containg wildtype and the mutated protein sequence.
#' @export
#'
#' @examples
#' pep <- maf2pep(annovar_path = annovar_path,maf_path = system.file("extdata", "test.maf", package = "TCEP"),
#'                need_allsamples = TRUE,num_thread=1)
maf2pep <- function(annovar_path,maf_path,
                    need_allsamples=FALSE,need_samples,num_thread){
  need_allsamples <- match.arg(as.character(need_allsamples),choices = c("TRUE","FALSE"))
  maf <- data.table::fread(maf_path,data.table = F,skip = "Hugo_Symbol")
  genome_version <- if (unique(maf$NCBI_Build) == "GRCh37"){"hg19"} else {"hg38"}
  if (need_allsamples){
    annovar_input <- maf %>%
      dplyr::select(Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2,Tumor_Sample_Barcode)
  }else{
    annovar_input <- maf %>%
      dplyr::select(Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2,Tumor_Sample_Barcode) %>%
      dplyr::filter(Tumor_Sample_Barcode %in% need_samples)
  }

  temp_file <- tempfile()
  temp_dir <- tempdir()
  file.create(temp_file)

  write.table(annovar_input,file = temp_file,col.names = F,row.names = F,sep = "\t",quote = F)
  comm_annotate <- paste0(annovar_path,"/table_annovar.pl ",temp_file," ",annovar_path,"/humandb/",
                          " -build ",genome_version,
                          " --outfile ",temp_dir,"/myanno",
                          ' -protocol refGene -operation g  --codingarg "-includesnp -onlyAltering"',
                          " --thread ",num_thread)
  system(comm_annotate)

  dt <- seqinr::read.fasta(file = paste0(temp_dir,"/myanno.refGene.fa"),
                           seqtype = "AA",as.string = TRUE,whole.header = T)

  variants <- names(dt)[seq(2,length(dt),by=2)]
  variants <- variants[grepl("protein-altering",variants)]
  pos_alter <- stringr::str_extract_all(variants,"p[.].+ protein-altering") %>% gsub(" protein-altering","",.)
  cdna <- stringr::str_extract_all(variants,"c[.].+ p[.]") %>% gsub(" p[.]","",.)
  lines <- stringr::str_extract_all(variants,"line.+ NM_") %>% gsub(" NM_","",.) %>%
    gsub("line","",.) %>% as.numeric()
  ENST <- stringr::str_extract_all(variants,"line.+ NM_[0-9]+") %>% gsub("line.+ ","",.)##extract transcripts

  mut_need <- annovar_input[lines,]
  pep_seq_mt <- dt[sapply(variants,function(x){which(names(dt)==x)}) %>% unname()] %>% as.character()
  pep_seq_wt <- dt[(sapply(variants,function(x){which(names(dt)==x)}) %>% unname())-1] %>% as.character()
  mut_need$seq_wt <- pep_seq_wt
  mut_need$seq_mt <- pep_seq_mt
  mut_need$pos_alter <- pos_alter
  mut_need$cdna <- cdna
  mut_need$transcript <- ENST
  files_exist <- c(list.files(temp_dir,pattern = gsub(paste0(temp_dir,"/"),"",temp_file),full.names = T),
                   list.files(temp_dir,pattern = "myanno",full.names = T))##remove all temp files
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
#' @param num_thread specify the number of threads to be used in annotation
#' @return A dataframe containg wildtype and the mutated protein sequence, plus exacted sequence from mutated peptides.
#' @export
#'
#' @examples
#' maf2seq(annovar_path = annovar_path,maf_path = system.file("extdata", "test.maf", package = "TCEP"),
#'         need_allsamples = TRUE,len = 9,num_thread=1)
maf2seq <- function(annovar_path,maf_path,need_allsamples=TRUE,need_samples,len,num_thread){

  pep <- maf2pep(annovar_path = annovar_path,maf_path = maf_path,
                 need_allsamples = need_allsamples,
                 need_samples = need_samples,num_thread=num_thread)

  colnames(pep)[1:6] <- c("chr","start","end","ref","alt","tumor_sample_barcode")[1:6]
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

#' Predict peptide MHC binding based on MAF file.
#'
#' @description This function use \code{\link{maf2seq}} to exacted mutated "new peptide" from MAF file, and predicted the binding affinity of these peptide with specific MHC allele using \code{\link{general_mhcbinding}}
#' @param annovar_path Character, the install path of annovar.
#' @param maf_path Character, the path of MAF file needed to be converted.
#' @param need_allsamples Logical. Whether need all samples when multi-sample MAF file is supplied.
#' @param need_samples Character, the needed sample when multi-sample MAF file is supplied.
#' @param hla_type HLA class, can be HLA-I or HLA-II.
#' @param pep_length A numeric or character vector, indicating the length for which to make predictions. For MHC-I, the length can be 8-15, for MHC-II, the length can be 11-30 or asis (take the length of input sequence as the peptide length)
#' @param allele A character vector of HLA alleles, available alleles for specific method can be obtained by \code{\link{available_alleles}}
#' @param pre_method Character, indicating the prediction method. Available methods for HLA-I or HLA-II can be obtained by \code{\link{available_methods}}
#' @param method_type, which type prediction method used, could be "Binding", "Processing" or "Immuno"
#' @param client_path The path of local IEDB tools, used when setting get_method as client
#' @param num_thread specify the number of threads to be used in annotation
#' @param tmp_dir temp dir
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
#'test <- maf2binding(annovar_path = "~/software/annovar/annovar/",need_allsamples=TRUE,
#'                    maf_path = system.file("extdata", "test.maf", package = "TCEP"),
#'                    hla_type = "I",pep_length = c(9,10),
#'                    allele = c("HLA-A*0203"),pre_method = "DeepImmuno",tmp_dir=tempdir(),
#'                    num_thread=1,method_type = "Immuno",
#'                    Deepimmuno_env = "DeepImmuno",
#'                    Immuno_Deepimmuno_path = ""~/software/DeepImmuno/"")
maf2binding <- function(annovar_path,maf_path,
                        need_allsamples=FALSE,need_samples,
                        hla_type,pep_length,allele,
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
    res[[i]] <- maf2seq(annovar_path = annovar_path,maf_path = maf_path,
                        need_allsamples = need_allsamples,need_samples = need_samples,len = pep_length[i],
                        num_thread=num_thread)
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
