

##convert2annovar.pl -format vcf4 example/ex2.vcf -outfile ex2.avinput -allsample -withfreq -include
vcf2input <- function(annovar_path,vcf_path,out_file,
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
#coding_change.pl common_driver.exonic_variant_function humandb/hg19_ensGene.txt humandb/hg19_ensGeneMrna.fa -includesnp --alltranscript -onlyAltering --outfile coding.txt
vcf2pep <- function(annovar_path,vcf_path,pep_len,
                    genome_version=c("hg19","hg38")){
  genome_version <- match.arg(genome_version)
  temp_file <- tempfile()
  temp_dir <- tempdir()
  file.create(temp_file)
  vcf2input(annovar_path = annovar_path,vcf_path = vcf_path,out_file = temp_file,need_allsamples = FALSE,need_samples = "TUMOR")

  mut <- read.table(temp_file)
  colnames(mut) <- c("chr","start","end","ref","alt","zygosity_status", "genotype_quality","read_depth")
  comm_annotate <- paste0(annovar_path,"/annotate_variation.pl -out ",temp_dir,"/test -build ",genome_version," ",temp_file,
                          " ",annovar_path,"/humandb/ -hgvs -dbtype ensGene")
  system(comm_annotate)
  comm_coding <- paste0(annovar_path,"/coding_change.pl ",temp_dir,"/test.exonic_variant_function ",annovar_path,"/humandb/",
                        genome_version,"_ensGene.txt ",annovar_path,"/humandb/",genome_version,"_ensGeneMrna.fa ",
                        "-includesnp --alltranscript -onlyAltering --outfile ",temp_dir,"/coding.txt")
  system(comm_coding)

  dt <- seqinr::read.fasta(file = paste0(temp_dir,"/coding.txt"),seqtype = "AA",as.string = TRUE,whole.header = T)
  file.remove(temp_file)

  variants <- names(dt)[seq(2,length(dt),by=2)]
  pos_num <- stringr::str_extract_all(variants,"p[.].+ p") %>% gsub(" p","",.)
  pos <- stringr::str_extract_all(pos_num,"[0-9]+") %>% unlist()
  aa <- pos_num %>% gsub("[0-9]","",.) %>% gsub("p[.]","",.)
  pep_seq_mt <- dt[seq(2,length(dt),by=2)] %>% as.character()
  pep_seq_wt <- dt[seq(1,length(dt),by=2)] %>% as.character()

  anno <- data.frame(len=rep(pep_len,times=length(pos)),
                     pos = rep(pos,each=length(pep_len)),
                     aa = rep(aa,each=length(pep_len)))
  anno <- anno %>%
    mutate(pos = as.numeric(pos)) %>%
    mutate(left=pos-(len-1),right=pos+(len-1)) %>%
    rowwise() %>%
    mutate(seq_mt=substr(pep_seq_mt,left,right),
           seq_wt=substr(pep_seq_wt,left,right))

}
