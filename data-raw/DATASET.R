## code to prepare `DATASET` dataset goes here
## MHC-I methods
dt <- read.table("clipboard",sep = "\t")
colnames(dt) <- c("method","version")
mhcIbinding_api_methods <- dt
usethis::use_data(mhcIbinding_api_methods,overwrite = T)
## MHC-II methods
dt <- read.table("clipboard",sep = "\t")
colnames(dt) <- c("method","version")
mhcIIbinding_api_methods <- dt
usethis::use_data(mhcIIbinding_api_methods,overwrite = T)

### MHC-I allele
for (i in mhcIbinding_api_methods$method){
  i <- gsub(" \\(netmhcpan\\)","",i)
  if (i == "recommended"){i <- "IEDB_recommended"}
  file.create(paste0("~/tmp/mhciallele_",i))
  comm <- paste0("~/software/mhc_i/src/predict_binding.py ",i,paste0(" mhc > ~/tmp/mhciallele_",i))
  mess <- try(system(comm))
}
files <- list.files("~/tmp",full.names = T,pattern = "mhciallele_")
res <- purrr::map(files,
                  function(x){
                    dt <- data.table::fread(x,blank.lines.skip = T,data.table = F) %>% as_tibble()
                    dt <- dt[1:(nrow(dt)-2),]
                    dt$method <- gsub("/home/data/t040201/tmp/","",x)
                    return(dt)})
res <- bind_rows(res)
colnames(res)[1] <- "alleles"
res$method <- gsub("mhciallele_","",res$method)
mhcIallele <- res
mhcIallele <- mhcIallele %>%
  filter(alleles == "human") %>%
  distinct(MHC,method,.keep_all = T) %>%
  select(-PeptideLength)
mhcIallele <- mhcIallele %>%
  select(-alleles) %>%
  rename(alleles=MHC)
usethis::use_data(mhcIallele,overwrite = T)

mhcIallele <- mhcIallele
load("~/TCEP/data-raw/mhcnuggest_alleles.Rdata")
mhcIallele <- bind_rows(
  mhcIallele,
  data.frame(alleles = mhcnuggest$V1,
             method = c(rep("mhcnuggest",nrow(mhcnuggest))))
)

mhcflurry <- read.csv("data-raw/class1_pseudosequences.csv")
mhcflurry <- mhcflurry %>%
  filter(grepl("HLA",allele))
mhcIallele <- bind_rows(
  mhcIallele,
  data.frame(alleles = mhcflurry$allele,
             method = c(rep("mhcflurry",nrow(mhcflurry))))
)
usethis::use_data(mhcIallele,overwrite = T)

mhcIallele <- mhcIallele
mhcIallele$method[which(mhcIallele$method == "mhcnuggest")] <- "mhcnuggets"

mhcIallele <- bind_rows(
  mhcIallele,
  data.frame(alleles = pseudu$V1 %>% unique(),
             method = "Immuno-GNN")
)
usethis::use_data(mhcIallele,overwrite = T)
###MHC-II
for (i in mhcIIbinding_api_methods$method){
  i <- gsub(" \\(netmhcii\\)","",i)
  file.create(paste0("~/tmp/mhciiallele_",i))
  comm <- paste0('curl --data "method=',i,'" http://tools-cluster-interface.iedb.org/tools_api/mhcii/ > ',paste0("~/tmp/mhciiallele_",i))
  mess <- try(system(comm))
  if (mess == 0){
    message("Succeed !")
  }else{

    for (j in 1:10){
      warning(paste0("Failed retrieving, retrying ",j," times"))
      mess1 <- try(system(comm))
      if (mess1 == 0){
        message("Succeed !")
        break
      }
    }
    if (j == 10){
      stop("Failed retrieving, stop", immediate. = TRUE)
    }
  }
}
files <- list.files("~/tmp",full.names = T,pattern = "mhciiallele")
res <- purrr::map(files,
                  function(x){
                    dt <- data.table::fread(x,blank.lines.skip = T,data.table = F) %>% as_tibble()
                    dt <- dt[1:(nrow(dt)-2),]
                    dt$method <- gsub("/home/wt/tmp/","",x)
                    return(dt)})
res <- bind_rows(res)
colnames(res)[1] <- "alleles"
res$method <- gsub("mhciiallele_","",res$method)
mhcIIallele <- res
usethis::use_data(mhcIIallele)

###alleles for MHCII client
netmhcpan_ba <- mhcIIallele %>% filter(method=="netmhciipan") %>% mutate(method="netmhciipan_ba")
netmhcpan_el <- mhcIIallele %>% filter(method=="netmhciipan") %>% mutate(method="netmhciipan_el")
dt <- data.table::fread("~/software/mhc_ii/all_alleles",nrows = 14,skip = 3,header = F)
consensus3 <- data.frame(alleles=c(dt[,2:7] %>% unlist() %>% unname(),dt$V8[1]),method="consensus3")

dt <- data.table::fread("~/software/mhc_ii/all_alleles",nrows = 11,skip = 17,header = F)
comblib <- data.frame(alleles=c(dt$V2,dt$V3[1:5]),method="comblib")

dt <- data.table::fread("~/software/mhc_ii/all_alleles",nrows = 11,skip = 31,header = F)
smm_align <- data.frame(alleles=c(dt[,2:3] %>% unlist() %>% unname(),dt$V4[1:7]),method="smm_align")

dt <- data.table::fread("~/software/mhc_ii/all_alleles",nrows = 29,skip = 45,header = F)
nn_align <- data.frame(alleles=c(dt[,2:3] %>% unlist() %>% unname(),dt$V4[1:3]),method="nn_align")

dt <- data.table::fread("~/software/mhc_ii/all_alleles",nrows = 15,skip = 77,header = F)
sturniolo <- data.frame(alleles=c(dt[,2:4] %>% unlist() %>% unname(),dt$V5[1:6]),method="sturniolo")

mhcIIallele_client <- bind_rows(netmhcpan_ba,netmhcpan_el,consensus3,comblib,smm_align,nn_align,sturniolo)

IEDB_recommended <- mhcIIallele %>% filter(method=="recommended") %>% mutate(method="IEDB_recommended")
mhcIIallele_client <- bind_rows(mhcIIallele_client,IEDB_recommended)
usethis::use_data(mhcIIallele_client,overwrite = T)

mhcIIallele_client <- mhcIIallele_client %>% select(-`method == "IEDB_recommended"`)
mhcIIallele_client <- bind_rows(
  mhcIIallele_client,
  data.frame(alleles = mhcnuggest2$V1,
             method = c(rep("mhcnuggest",nrow(mhcnuggest2))))
)
usethis::use_data(mhcIIallele_client,overwrite = T)
mhcIIallele_client <- mhcIIallele_client
mhcIIallele_client$method[which(mhcIIallele_client$method == "mhcnuggest")] <- "mhcnuggets"
###pep file
dt <- random_peptides(len = 12,n=10)
write.table(dt,file = "inst/extdata/random.pep",row.names = F,col.names = F,quote = F)

###the allele for method has support length
methods <- available_methods("client","MHC-I")
res <- vector("list",11)
for (i in methods){
  command_run <- paste0('~/software/mhc_i/src/predict_binding.py ',
                        i," mhc",' > ',paste0("data-raw/",i))
  system(command_run)
  dt <- data.table::fread(paste0("data-raw/",i),data.table = F,skip = 1)
  dt$method <- i
  res[[i]] <- dt
}
res <- bind_rows(res)
res <- res %>% filter(Species=="human")
available_lens <- res
usethis::use_data(available_lens)
available_lens <- available_lens
available_lens <- bind_rows(
  available_lens,
  data.frame(Species="human",MHC=rep(available_alleles("Processing","I","NetCTLpan"),each=4),
             PeptideLength=rep(c(8:11),length(available_alleles("Processing","I","NetCTLpan"))),
             method="NetCTLpan")
)
available_lens <- bind_rows(
  available_lens,
  data.frame(Species="human",MHC=rep(available_alleles("Binding","I","mhcflurry"),each=11),
             PeptideLength=rep(c(5:15),length(available_alleles("Binding","I","mhcflurry"))),
             method="mhcflurry")
)
available_lens <- bind_rows(
  available_lens,
  data.frame(Species="human",MHC=rep(available_alleles("Binding","I","mhcnuggets"),each=11),
             PeptideLength=rep(c(5:15),length(available_alleles("Binding","I","mhcnuggets"))),
             method="mhcnuggets")
)
available_lens <- bind_rows(
  available_lens,
  data.frame(Species="human",MHC=rep(available_alleles("Immuno","I","IEDB"),each=1),
             PeptideLength=rep(c(9),length(available_alleles("Immuno","I","IEDB"))),
             method="IEDB")
)
available_lens <- bind_rows(
  available_lens,
  data.frame(Species="human",MHC=rep(available_alleles("Immuno","I","DeepImmuno"),each=2),
             PeptideLength=rep(c(9:10),length(available_alleles("Immuno","I","DeepImmuno"))),
             method="DeepImmuno")
)
available_lens <- bind_rows(
  available_lens,
  data.frame(Species="human",MHC=rep(available_alleles("Immuno","I","Seq2Neo-CNN"),each=4),
             PeptideLength=rep(c(8:11),length(available_alleles("Immuno","I","Seq2Neo-CNN"))),
             method="Seq2Neo-CNN")
)
available_lens <- bind_rows(
  available_lens,
  data.frame(Species="human",MHC=rep(available_alleles("Immuno","I","PRIME2.0"),each=7),
             PeptideLength=rep(c(8:14),length(available_alleles("Immuno","I","PRIME2.0"))),
             method="PRIME2.0")
)
usethis::use_data(available_lens,overwrite = TRUE)

available_lens <- available_lens
available_lens <- bind_rows(
  available_lens,
  data.frame(Species="human",MHC=rep(available_alleles("Immuno","I","Immuno-GNN"),each=15),
             PeptideLength=rep(c(6:20),length(available_alleles("Immuno","I","Immuno-GNN"))),
             method="Immuno-GNN")
)
usethis::use_data(available_lens,overwrite = TRUE)
###各种方法的列名
#MHC-I
smmpmbec <- c("ic50","rank")
smm <- c("ic50","rank")
IEDB_recommended <- c("score","rank")
pickpocket <- c("ic50","rank")
netmhcstabpan <- c("ic50","rank")
netmhcpan_el <- c("score","rank")
netmhcpan_ba <- c("ic50","rank")
netmhccons <- c("ic50","rank")
comblib_sidney2008 <- c("score","rank")
ann <- c("ic50","rank")
consensus <- colnames(consensus_method_res)[16:22]
res_cols <- list(smmpmbec,smm,IEDB_recommended,pickpocket,
                 netmhcstabpan,netmhcpan_el,netmhcpan_ba,netmhccons,
                 comblib_sidney2008,ann,consensus)
names(res_cols) <- c("smmpmbec","smm","IEDB_recommended","pickpocket",
                     "netmhcstabpan","netmhcpan_el","netmhcpan_ba","netmhccons",
                     "comblib_sidney2008","ann","consensus")
usethis::use_data(res_cols)
load("~/TCEP/data/res_cols.rda")
res_cols["mhcflurry"][[1]] <- colnames(tmp)[3:7]
res_cols["mhcnuggets"][[1]] <- "ic50"
usethis::use_data(res_cols, overwrite = TRUE)

res_cols <- res_cols
res_cols["NetCTLpan"][[1]] <- colnames(tmp)[2:5]
usethis::use_data(res_cols, overwrite = TRUE)

res_cols <- res_cols
res_cols["IEDB"][[1]] <- colnames(tmp)[3]
usethis::use_data(res_cols, overwrite = TRUE)

res_cols <- res_cols
res_cols["DeepImmuno"][[1]] <- colnames(tt)[3]
usethis::use_data(res_cols, overwrite = TRUE)

res_cols <- res_cols
res_cols["PRIME2.0"][[1]] <- c("PRIME_rank","PRIME_score","MixMHCpred_rank")
usethis::use_data(res_cols, overwrite = TRUE)

res_cols <- res_cols
res_cols["Seq2Neo-CNN"][[1]] <- colnames(tmp)[2:4]
usethis::use_data(res_cols, overwrite = TRUE)
##MHC-II
method <- available_methods("client","MHC-II")
pre_method <- method[8]
allele <- available_alleles("MHC-II",pre_method)[1,] %>% unlist()
test <- txt2binding(get_method="client",annovar_path = "~/software/annovar/",
                    txt_path = system.file("extdata", "test_avinput.txt", package = "TCEP"),
                    genome_version = "hg19",mhc_type = "MHC-II",pep_length = c(14),
                    allele = allele,pre_method = pre_method,
                    tmp_dir=tempdir(),num_thread=1,client_path = "~/software/mhc_ii/")
comblib <- c("ic50","percentile_rank","adjusted_rank")
consensus3 <- c("consensus_percentile_rank","adjusted_consensus_percentile_rank","comblib_core",
                "comblib_score","comblib_rank","adjusted_comblib_rank","smm_align_core",
                "smm_align_ic50","smm_align_rank",
                "adjusted_smm_align_rank","nn_align_core","nn_align_ic50",
                "nn_align_rank","adjusted_nn_align_rank","sturniolo_core",
                "sturniolo_score","sturniolo_rank","adjusted_sturniolo_rank" )
IEDB_recommended <- IEDB_rem_mhcii
netmhciipan_el <- c("score","percentile_rank")
netmhciipan_ba <- c("ic50","percentile_rank")
nn_align <- c("ic50","percentile_rank","adjusted_rank")
smm_align <- c("ic50","percentile_rank","adjusted_rank")
sturniolo <- c("score","percentile_rank","adjusted_rank")
res_cols_ii <- list(comblib,consensus3,IEDB_recommended,netmhciipan_el,
                    netmhciipan_ba,nn_align,smm_align,sturniolo)
names(res_cols_ii) <- c("comblib","consensus3","IEDB_recommended","netmhciipan_el",
                     "netmhciipan_ba","nn_align","smm_align","sturniolo")
usethis::use_data(res_cols_ii,overwrite = TRUE)
res_cols_ii <- res_cols_ii
res_cols_ii["mhcnuggets"][[1]] <- "ic50"
###netChop alleles
netChop_alleles <- readRDS("~/TCEP/data-raw/netChop_alleles.rds")
colnames(netChop_alleles) <- "Allele"
NetCTLpan_alleles <- netChop_alleles
usethis::use_data(NetCTLpan_alleles)

###immuno alleles
iedb_imm_alleles <- readRDS("~/TCEP/data-raw/iedb_imm_alleles.rds")
prime_alleles <- read.table("data-raw/PRIME_alleles.txt")
deepimmuno <- read.table("data-raw/hla2paratopeTable_aligned.txt",header = T)
seq2neo_alleles <- mhcIallele[mhcIallele$method=="netmhcpan_el","alleles"] %>%
  unlist() %>% unname()

immuno_alleles <- data.frame(methods = c(rep("IEDB",length(iedb_imm_alleles)),
                                         rep("PRIME",nrow(prime_alleles)),
                                         rep("DeepImmuno",nrow(deepimmuno)),
                                         rep("Seq2Neo-CNN",length(seq2neo_alleles))),
                             alleles = c(iedb_imm_alleles,prime_alleles$V1,
                                         deepimmuno$HLA,
                                         seq2neo_alleles))
usethis::use_data(immuno_alleles,overwrite = TRUE)
immuno_alleles <- immuno_alleles
immuno_alleles$methods[which(immuno_alleles$methods=="PRIME")] <- "PRIME2.0"

immuno_alleles <- immuno_alleles
immuno_alleles <- bind_rows(
  immuno_alleles,
  data.frame(
    methods = "Immuno-GNN",
    alleles = pseudu$V1 %>% unique())
)
usethis::use_data(immuno_alleles,overwrite = T)
##Repitope
devtools::install_github("masato-ogishi/Repitope",lib="/home/data/R/R-4.1.0/library")
library(Repitope,lib.loc = "/home/data/R/R-4.1.0/library")
##下载数据 https://data.mendeley.com/datasets/sydw5xnxpt/1
fragLibDT <- fst::read_fst("~/tmp/Repitope/FragmentLibrary_TCRSet_Public_RepitopeV3.fst", as.data.table=T)
featureDF_MHCI <- fst::read_fst("~/tmp/Repitope/FeatureDF_MHCI_Weighted.10000_RepitopeV3.fst", as.data.table=T)
featureDF_MHCII <- fst::read_fst("~/tmp/Repitope/FeatureDF_MHCII_Weighted.10000_RepitopeV3.fst", as.data.table=T)

res_MHCI <- EpitopePrioritization(
  featureDF=featureDF_MHCI[Peptide%in%MHCI_Human$Peptide,],
  metadataDF=MHCI_Human[,.(Peptide,Immunogenicity)],
  peptideSet=c("GHAHKVPRRLLKAA", "SLYNTVATLY"),
  peptideLengthSet=8:15,
  fragLib=fragLibDT,
  aaIndexIDSet="all",
  fragLenSet=3:8,
  fragDepth=10000,
  fragLibType="Weighted",
  featureSet=MHCI_Human_MinimumFeatureSet,
  seedSet=1:5,
  coreN=parallel::detectCores(logical=F),
  outDir="./Output"  ## Intermediate and final output files will be stored under this directory
)

##
pseudu  <- data.table::fread("~/software/neo_pre/netMHCpan-4.1/data/MHC_pseudo.dat",header = F,fill=TRUE,data.table = F)
pseudu <- pseudu %>% filter(grepl("HLA",V1)) %>% filter(grepl(":",V1))
usethis::use_data(pseudu)

