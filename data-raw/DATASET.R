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

###各种方法的列名
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
