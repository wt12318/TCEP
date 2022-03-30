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

###pep file
dt <- random_peptides(len = 12,n=10)
write.table(dt,file = "inst/extdata/random.pep",row.names = F,col.names = F,quote = F)



