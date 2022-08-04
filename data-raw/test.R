test <- txt2binding(annovar_path = "~/software/annovar/annovar/",
                    txt_path = system.file("extdata", "test_avinput.txt", package = "TCEP"),
                    genome_version = "hg19",hla_type = "I",pep_length = c(9),
                    allele = c("HLA-A01:01"),pre_method = "NetCTLpan",tmp_dir=tempdir(),
                    num_thread=1,method_type = "Processing",netchop_path = "~/software/netchop/",
                    PRIME_path = "~/software/PRIME-2.0/",MixMHCpred_path = "~/software/MixMHCpred/MixMHCpred")

system(paste0("touch ",tmp_dir,"/pre.py"))
writeLines("from mhcnuggets.src.predict import predict\npredict(class_='I', peptides_path='/home/wt/tmp/a',mhc='HLA-A02:01')",
           paste0(tmp_dir,"/pre.py"))
system(paste0("conda run -n mhcnuggets python ",paste0(tmp_dir,"/pre.py")," > ",paste0(tmp_dir,"/b")))

tt <- TCEP:::general_mhcbinding(hla_type = "I", length = c(9,8),
                                allele = c("HLA-A*01:01","HLA-A*01:02"),pre_method = "netmhcpan_ba",
                                method_type="Binding",
                                peptide = c("SLYNTVATLY","GHAHKVPR"),
                                tmp_dir="~/tmp/",client_path = "~/software/mhc_i/src/")
tt1 <- TCEP:::general_mhcbinding(hla_type = "I", length = c(9,8),
                                      allele = c("HLA-A01:01","HLA-A01:02"),pre_method = "NetCTLpan",
                                      method_type="Processing",
                                      peptide = c("SLYNTVATLY","GHAHKVPR"),
                                      tmp_dir="~/tmp/",netchop_path = "~/software/netchop/")
tt <- seq2neo_help(pep_len=c(10),allele=c("HLA-A*01:01"),
                   tmp_dir="~/tmp/", netchop_path="~/software/netchop/",
                   peptide=c("NVDTHPGSGK"),
                   client_path = "~/software/mhc_i/src/")
tt1 <- seq2neo_help(pep_len=c(10),allele=c("HLA-A*01:01"),
                   tmp_dir="~/tmp/", netchop_path="~/software/netchop/",
                   peptide=c("NVDAHPGSGK"),
                   client_path = "~/software/mhc_i/src/")
test <- txt2binding(annovar_path = "~/software/annovar/annovar/",
                    txt_path = system.file("extdata", "test_avinput.txt", package = "TCEP"),
                    genome_version = "hg19",hla_type = "I",pep_length = c(9,10),
                    allele = c("HLA-A*01:01","HLA-A*01:03"),pre_method = "Seq2Neo-CNN",tmp_dir=tempdir(),
                    num_thread=1,method_type = "Immuno",
                    netchop_path = "~/software/netchop/",client_path = "~/software/mhc_i/src/",
                    seq2neo_env = "DeepImmuno",
                    seq2neo_path = "~/software/Seq2Neo/seq2neo/function/immuno_Prediction/"
                    )
test <- maf2binding(annovar_path = "~/software/annovar/annovar/",need_allsamples=TRUE,
                    maf_path = system.file("extdata", "test.maf", package = "TCEP"),
                    hla_type = "I",pep_length = c(9,10),
                    allele = c("HLA-A*01:01","HLA-A*01:03"),pre_method = "Seq2Neo-CNN",tmp_dir=tempdir(),
                    num_thread=1,method_type = "Immuno",
                    netchop_path = "~/software/netchop/",client_path = "~/software/mhc_i/src/",
                    seq2neo_env = "DeepImmuno",
                    seq2neo_path = "~/software/Seq2Neo/seq2neo/function/immuno_Prediction/")
test <- vcf2binding(annovar_path = "~/software/annovar/annovar/",
                   vcf_path = system.file("extdata", "test_grch38.vcf", package = "TCEP"),
                   hla_type = "I",pep_length = c(9,10),genome_version = "hg38",
                   need_allsamples = FALSE,need_samples = "TUMOR",
                    allele = c("HLA-A*01:01","HLA-A*01:03"),pre_method = "Seq2Neo-CNN",tmp_dir=tempdir(),
                    num_thread=1,method_type = "Immuno",
                    netchop_path = "~/software/netchop/",client_path = "~/software/mhc_i/src/",
                    seq2neo_env = "DeepImmuno",
                    seq2neo_path = "~/software/Seq2Neo/seq2neo/function/immuno_Prediction/")

test <- batchpep_binding(pep_file=system.file("extdata", "random.pep", package = "TCEP"),
                    hla_type = "I",pep_length = c(9,10),
                    allele = c("HLA-A*01:01","HLA-A*01:03"),pre_method = "Seq2Neo-CNN",tmp_dir=tempdir(),
                    num_thread=1,method_type = "Immuno",
                    netchop_path = "~/software/netchop/",client_path = "~/software/mhc_i/src/",
                    seq2neo_env = "DeepImmuno",
                    seq2neo_path = "~/software/Seq2Neo/seq2neo/function/immuno_Prediction/")

tt1 <- TCEP:::general_mhcbinding(hla_type = "I", length = c(10),
                                       allele = c("HLA-A*01:01"),pre_method = "Seq2Neo-CNN",
                                       method_type="Immuno",
                                       peptide = c("NVDTHPGSGK","QTSEKALLRR"),
                                       tmp_dir="~/tmp/",netchop_path = "~/software/netchop/",
                                       client_path = "~/software/mhc_i/src/",
                                       seq2neo_env = "DeepImmuno",
                                       seq2neo_path = "~/software/Seq2Neo/seq2neo/function/immuno_Prediction/")
