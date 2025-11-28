library(httr)
library(jsonlite)
library(dplyr)
library(data.table)
library(edgeR)
library(limma)
library(tidyr)
library(stringr)
library(R.utils)
library(GEOquery)

source("~/Documents/Publications/webGEO/script/helper_functions.R")

setwd("~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/")

getOption("timeout")
options(timeout = max(1000, getOption("timeout")))

###########################################################
###             NCBI-generated RNAseq Count             ###
###########################################################

projects = read.csv("data/geo_studies_with_ncbi_generated_rnaseq_counts.csv", header = T, stringsAsFactors = F)
accessions <- projects$Accession

# geo_accession, title, source_name_ch1, description, relation ..., characteristics, ...
# sample_id, title, source_name, characteristics, description, srx_accession, biosample_accession

organize_geo_sample_metadata <- function(metadata) {
  
  characteristics_cols <- names(metadata)[grepl(":ch1", names(metadata))]
  characteristics <- metadata %>% select(geo_accession, all_of(characteristics_cols)) %>%
    rename_with(~ gsub(":ch1$", "", .x))
  
  relation_cols <- names(metadata)[grepl("^relation", names(metadata))]
  
  biosample_accession <- NA
  sra_accession <- NA
  
  if (length(relation_cols)>0) {
    for (col in relation_cols) {
      if (all(grepl("^BioSample", metadata[[col]]))) {
        biosample_accession <- str_extract(metadata[[col]], "SAMN[0-9]+")
      } else if (all(grepl("^SRA", metadata[[col]]))) {
        sra_accession <- str_extract(metadata[[col]], "SRX[0-9]+")
      }
    }
  }
  
  relation <- metadata %>% select(geo_accession) %>%
    mutate(biosample_accession = biosample_accession,
           sra_accession = sra_accession)
  
  description_cols <- names(metadata)[grepl("^description", names(metadata))]
  if (length(description_cols) > 0) {
    description <- metadata %>% select(geo_accession, all_of(description_cols)) %>%
      mutate(description_combined = apply(select(., all_of(description_cols)), 1, function(row) {
        paste(na.omit(row), collapse = "\n")
      })) %>%
      select(-all_of(description_cols))
  } else {
    description <- metadata %>% select(geo_accession) %>%
      mutate(description_combined = NA)
  }

  metadata <- metadata %>% mutate(sample_id = geo_accession) %>% 
    select(sample_id, title, source_name_ch1, geo_accession) %>%
    rename(source_name = source_name_ch1)
  
  ### description column is not included in the final metadata table
  # metadata_list <- list(metadata, characteristics, description, relation)
  # metadata <- Reduce(function(x, y) merge(x, y, by = "geo_accession"), metadata_list) %>%
  #   relocate(geo_accession, .before = biosample_accession) %>%
  #   rename(description = description_combined)
  
  metadata_list <- list(metadata, characteristics, relation)
  metadata <- Reduce(function(x, y) merge(x, y, by = "geo_accession"), metadata_list) %>%
    relocate(geo_accession, .before = biosample_accession)
  
  return(metadata)
}

annot = fread("data/Human.GRCh38.p13.annot.tsv.gz", sep = "\t", header = T, stringsAsFactors = F)
annot <- annot %>% select("GeneID", "Symbol")
annot

ncbi_data_ingestion <- function(gse_acc, annot) {

  print("Step1: get sample metadata")
  
  metadata <- parseSOFTMetaFun(accession=gse_acc, dest.dir = "data/geo/metadata/")
  saveRDS(metadata, file = sprintf("data/geo/metadata/%s_Sample_Metadata.raw.rds",gse_acc))
  
  metadata <- organize_geo_sample_metadata(metadata)
  saveRDS(metadata, file = sprintf("data/geo/metadata/%s_Sample_Metadata.rds",gse_acc))
  
  Sys.sleep(1)
  
  print("Step2: get raw count data")
  
  url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=%s&format=file&file=%s_raw_counts_GRCh38.p13_NCBI.tsv.gz", 
                 gse_acc, gse_acc)
  
  count_data <- fread(url, sep = "\t", header = TRUE)
  
  count_data <- count_data %>%
    left_join(annot, by = "GeneID") %>%
    relocate(GeneID, Symbol) %>%
    as.data.frame()
  
  write.csv(count_data, file=sprintf("data/geo/count/%s_Count_Data.ncbi.csv",gse_acc),
            quote = F, row.names = F)
  
  saveRDS(count_data, file = sprintf("data/geo/count/%s_Count_Data.ncbi.rds",gse_acc))
  
  Sys.sleep(1)
  
  print("Step3: data normalization")
  
  dge <- DGEList(counts = count_data %>% select(-GeneID, -Symbol), genes = annot)
  dge <- calcNormFactors(dge, method = 'TMM')
  
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 0.1)
  logcpm <- data.frame(annot, logcpm)
  
  # write.csv(logcpm, file = sprintf("data/geo/normalized/%s_Normalized_Data.log2cpm.ncbi.csv",gse_acc), 
  #           row.names = FALSE, quote=F)
  
  saveRDS(logcpm, file = sprintf("data/geo/normalized/%s_Normalized_Data.log2cpm.ncbi.rds",gse_acc))
}


#gse_acc <- "GSE262419"
#data <- ncbi_data_ingestion(gse_acc = gse_acc, annot = annot)


completed <- gsub("_Normalized_Data.log2cpm.ncbi.rds", "", list.files("data/geo/normalized/"))
completed
remaining <- accessions[which(!accessions %in% completed)]
remaining

for (gse_acc in remaining[-(1:9)]) {
  message(gse_acc)
  
  data <- ncbi_data_ingestion(gse_acc = gse_acc, annot = annot)
}

# GSE243919
# GSE178270
# GSE186121
# GSE171360
# GSE152978
# GSE130644
# GSE79864
# GSE115572
# GSE77466
# GSE85839




###########################################################
###                recount3 RNAseq Count                ###
###########################################################

library(recount3)

recount3_cache_rm()

# projects = read.csv("data/recount3_projects_with_gse_accession.csv", header = T, stringsAsFactors = F)

annotation_options("human")
annotation_options("mouse")

## Obtain all available projects
projects <- rbind(
  recount3::available_projects("human"),
  recount3::available_projects("mouse")
)

# ## Locate the project raw files at the gene level using the default annotation
# projects$gene <- apply(projects, 1, function(x)
#   locate_url(
#     project = x["project"],
#     project_home = x["project_home"],
#     type = "gene",
#     organism = x["organism"],
#     annotation = ifelse(x["organism"]=="human", "gencode_v29", "gencode_v23")
#   ))
# 
# projects

## Locate the project raw files at the exon level using the default annotation
# projects$exon <- apply(projects, 1, function(x)
#   locate_url(
#     project = x["project"],
#     project_home = x["project_home"],
#     type = "exon",
#     organism = x["organism"],
#     annotation = ifelse(x["organism"]=="human", "gencode_v29", "gencode_v23")
#   ))

# ## Locate the project raw exon-exon junction files
# projects <-
#   cbind(projects, do.call(rbind, apply(projects, 1, function(x) {
#     x <-
#       locate_url(
#         project = x["project"],
#         project_home = x["project_home"],
#         type = "jxn",
#         organism = x["organism"]
#       )
#     res <- data.frame(t(x))
#     colnames(res) <-
#       paste0("jxn_", gsub("^.*\\.", "", gsub("\\.gz", "", colnames(res))))
#     return(res)
#   })))
# rownames(projects) <- NULL

i <- 1

proj_info <- projects[i,]
proj_info

View(projects)

annotation <- ifelse(proj_info$organism=="human", "gencode_v29", "gencode_v23")

rse_gene <- create_rse(proj_info, 
                       annotation = annotation, 
                       type = "gene")

counts <- compute_read_counts(rse_gene)
counts

metadata <- data.frame(rse_gene@colData, check.names = F)
metadata

all(colnames(counts) == rownames(metadata))


completed <- gsub("_Normalized_Data.log2cpm.recount3.rds", "", list.files("data/recount3/normalized/"))
#completed
remaining <- which(!projects$project %in% completed)
remaining

for (i in remaining) {
  print(i)
  
  proj_info <- projects[i,]
  project <- proj_info$project
  
  annotation <- ifelse(proj_info$organism=="human", "gencode_v29", "gencode_v23")
  
  rse_gene <- create_rse(proj_info, 
                         annotation = annotation, 
                         type = "gene")
  
  count_data <- compute_read_counts(rse_gene)
  
  idx <- which(colSums(count_data)>0)
  print(length(idx))
  print(ncol(count_data))
  count_data <- count_data[,idx]
  
  count_data <- data.frame(EnsemblID=gsub("\\.\\d+", "", rownames(count_data)), count_data, stringsAsFactors = F, check.names = F)
  
  metadata <- data.frame(rse_gene@colData, check.names = F)
  
  if (! all(colnames(count_data)[-1] == rownames(metadata))) {
    print("sample names not consistent")
  #  break
  }
  
  saveRDS(metadata, file = sprintf("data/recount3/metadata/%s_Sample_Metadata.recount3.rds",project))
  
  write.csv(count_data, file=sprintf("data/recount3/count/%s_Count_Data.recount3.csv", project),
            quote = F, row.names = F)
  
  saveRDS(count_data, file = sprintf("data/recount3/count/%s_Count_Data.recount3.rds", project))
  
  dge <- DGEList(counts = count_data %>% select(-EnsemblID))
  dge <- calcNormFactors(dge, method = 'TMM')
  
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 0.1)
  logcpm <- data.frame(EnsemblID=count_data$EnsemblID, logcpm, stringsAsFactors = F)
  
  # write.csv(logcpm, file = sprintf("data/geo/normalized/%s_Normalized_Data.log2cpm.ncbi.csv",gse_acc), 
  #           row.names = FALSE, quote=F)
  
  saveRDS(logcpm, file = sprintf("data/recount3/normalized/%s_Normalized_Data.log2cpm.recount3.rds", project))
  
}

# 218 695, 1138, 3359, 3581, 3657 3755 4171 4628 5194 6098 6312 7097 7987 8743 9011 10245 0 or 1 sample
# 3068, 5105, 10856 12103 16944 too much samples
# 8797 8802 8820 10288 10314 10321 10342 10350 10351 10871 10898 10947 13482 13483 14262 14458 14475 14494 14524 14533 15477 15533 16473 16511 17945 18652 18677 Not all 'samples' are present in the 'counts_file'


###########################################################
###                      Microarray                     ###
###########################################################

### Affymetrix

# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
# R Source Package: O

# Affymetrix Human Exon 1.0 ST Array: pd.huex10st.hs.gencodeg
# Affymetrix Human Gene 2.0 ST Array: pd.huex20st.hs.gencodeg
# Affymetrix Human Transcriptome Array 2.0: pd.hta20.hs.gencodeg
# Affymetrix Human Genome U133A Array: pd.hgu133a.hs.gencodeg
# Affymetrix Human Genome U133 Plus 2.0 Array: pd.hgu133plus2.hs.gencodeg
# Affymetrix Human Genome U133A 2.0 Array: pd.hgu133a2.hs.gencodeg

# library(pd.hg.u133.plus.2) # from Bioconductor (NOT IN USE)

install.packages("http://mbni.org/customcdf/25.0.0/gencodeg.download/pd.hgu133plus2.hs.gencodeg_25.0.0.tar.gz",
                 repos = NULL, type = "source")

pkgs <- list.files("~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/rpackages", full.names = T)

install.packages(pkgs, repos = NULL, type = "source")

# for (pkg in pkgs) {
#   install.packages(pkg, repos = NULL, type = "source")
# }

# pd_pkgs <- rownames(installed.packages())
# pd_pkgs <- pd_pkgs[grepl("^pd", pd_pkgs)]

# Load each package using library()
# for(pkg in pd_pkgs) {
#   library(pkg, character.only = TRUE)
# }


library(pd.huex10st.hs.gencodeg)
library(pd.hgu133a.hs.gencodeg)
library(pd.hgu133plus2.hs.gencodeg)
library(pd.hgu133a2.hs.gencodeg)

library(GEOquery)

platforms <- read.csv("~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/gse_platforms_summary_table.csv",
                     header = T, stringsAsFactors = F)
platforms

View(platforms)
sort(table(platforms$Contact), decreasing = T)

idx <- which(grepl("Homo sapiens|Mus musculus", platforms$Taxonomy, ignore.case = T) & grepl("affy", platforms$Title, ignore.case = T))
platforms.affy <- platforms[idx,]

platforms.affy

extract_brackets <- function(text_vec) {
  sapply(text_vec, function(x) {
    m <- regmatches(x, regexpr("\\[.*?\\]", x))
    if (length(m) == 0 || m == "") {
      return("")
    } else {
      # Remove the brackets
      return(gsub("\\[|\\]", "", m))
    }
  })
}

chips <- extract_brackets(platforms.affy$Title)


stats <- table(chips)
sort(stats, decreasing = T)[31:60]
View(stats)

platforms.affy$Chip[chips == "HG-U133_Plus_2"] <- "HG-U133 Plus 2.0"
platforms.affy$Chip[chips == "HuEx-1_0-st"] <- "Human Exon 1.0 ST"
platforms.affy$Chip[chips == "HuGene-1_0-st"] <- "Human Gene 1.0 ST"
platforms.affy$Chip[chips == "MoGene-1_0-st"] <- "Mouse Gene 1.0 ST"
platforms.affy$Chip[chips == "HuGene-2_0-st"] <- "Human Gene 2.0 ST"
platforms.affy$Chip[chips == "MoEx-1_0-st"] <- "Mouse Exon 1.0 ST"
platforms.affy$Chip[chips == "HTA-2_0"] <- "HTA 2.0"
platforms.affy$Chip[chips == "Mouse430_2"] <- "MG-430 2.0"
platforms.affy$Chip[chips == "MoGene-2_0-st"] <- "Mouse Gene 2.0 ST"
platforms.affy$Chip[chips == "HuGene-2_1-st"] <- "Human Gene 2.1 ST"
platforms.affy$Chip[chips == "HT_MG-430_PM"] <- "HT MG-430 PM"
platforms.affy$Chip[chips == "HuGene-1_1-st"] <- "Human Gene 1.1 ST"
platforms.affy$Chip[chips == "MoGene-2_1-st"] <- "Mouse Gene 2.1 ST"
platforms.affy$Chip[chips == "PrimeView"] <- "Human PrimeView"
platforms.affy$Chip[chips == "HG-U219"] <- "HG-U219"
platforms.affy$Chip[chips == "MoGene-1_1-st"] <- "Mouse Gene 1.1 ST"
platforms.affy$Chip[chips == "MTA-1_0"] <- "MTA 1.0"
platforms.affy$Chip[chips == "HG-U133A"] <- "HG-U133A"
platforms.affy$Chip[chips == "HG-U133A_2"] <- "HG-U133A 2.0"
platforms.affy$Chip[chips == "HT_HG-U133A"] <- "HT HG-U133A"
platforms.affy$Chip[chips == "HT_HG-U133_Plus_PM"] <- "HT HG-U133 Plus PM"
platforms.affy$Chip[chips == "Mouse430_2_Mm_REFSEQ"] <- "MG-430 2.0"
platforms.affy$Chip[chips == "Mouse4302_Mm_EntrezG"] <- "MG-430 2.0"
platforms.affy$Chip[chips == "Clariom_D_Human"] <- "Clariom D Human"
platforms.affy$Chip[chips == "Clariom_S_Human"] <- "Clariom S Human"
platforms.affy$Chip[chips == "Clariom_S_Human_HT"] <- "Clariom S Human HT"
platforms.affy$Chip[chips == "Clariom_D_Mouse"] <- "Clariom D Mouse"
platforms.affy$Chip[chips == "Clariom_S_Mouse"] <- "Clariom S Mouse"
platforms.affy$Chip[chips == "HuGene10stv1_Hs_ENSG"] <- "Human Gene 1.0 ST"
platforms.affy$Chip[chips == "HG-U133B"] <- "HG-U133B"
platforms.affy$Chip[chips == "HT_HG-U133B"] <- "HT HG-U133B"
platforms.affy$Chip[chips == "MG_U74Av2"] <- "MG-U74Av2"
platforms.affy$Chip[chips == "HTMG430PM_Mm_ENTREZG"] <- "HT MG-430 PM"
platforms.affy$Chip[chips == "Mouse4302_Mm_ENTREZG"] <- "MG-430 2.0"

platforms.affy$Chip[grep("U133 Plus 2.0|U133 Plus2", platforms.affy$Title)] <- "HG-U133 Plus 2.0"
platforms.affy$Chip[grep("Mouse Genome 430 2.0", platforms.affy$Title)] <- "MG-430 2.0"

write.table(platforms.affy, "~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/gse_platforms_summary_table.hs_mm_affymetrix.txt",
          quote = F, row.names = F, sep = "\t")


###### 

platforms.affy <- read.csv("~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/gse_platforms_summary_table.hs_mm_affymetrix.cleaned.csv",
                           header = T, stringsAsFactors = F)

platforms.affy <- platforms.affy %>% 
  filter(!is.na(Chip) & Chip != "Others") %>%
  arrange(desc(Series.Count))
  
View(platforms.affy)

t <- platforms.affy %>% group_by(Chip) %>%
  summarise(n.samples=sum(Samples.Count),
            n.series=sum(Series.Count)) %>%
  arrange(desc(n.series))

t
sum(t$n.series) # 28182
View(t)



######

samples <- read.csv("~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/gse_samples_summary_table.csv",
                    header = T, stringsAsFactors = F)
dim(samples)
which(duplicated(samples$Accession))

View(table(samples$Taxonomy))

samples.filtered <- samples %>% 
  filter(Taxonomy %in% c("Homo sapiens", "Mus musculus"),
         Platform %in% platforms.affy$Accession) %>% 
  separate_rows(Series, sep = ";")

dim(samples.filtered)
View(samples.filtered)

saveRDS(samples.filtered, 
        file = "~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/final_samples_list.hs_mm_affymetrix.RDS")


filter <- which(duplicated(paste0(samples.filtered$Series, "-", samples.filtered$Platform)))
filter

series.filtered <- samples.filtered[-filter,] %>% 
  select(Series, Platform, Taxonomy)
dim(series.filtered)
View(series.filtered)

idx <- match(series.filtered$Platform, platforms.affy$Accession)

series.filtered <- data.frame(series.filtered, Chip=platforms.affy$Chip[idx], Platform_Organism=platforms.affy$Taxonomy[idx])
series.filtered <- series.filtered %>% filter(Taxonomy == Platform_Organism)
dim(series.filtered)

accessions <- samples %>% 
  filter(Taxonomy %in% c("Homo sapiens", "Mus musculus"),
         Platform %in% platforms.affy$Accession,
         grepl(";", Series)) %>% 
  separate_rows(Series, sep = ";") %>%
  pull(Series) %>%
  unique()

View(accessions)


dest.dir <- "~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/microarray/series/"

scope <- "Series"
amount = 'brief'

super_series <- c()

for (accession in accessions) {
  
  soft.file <- file.path(dest.dir, paste0(accession, '.', scope, '.txt'))
  
  if (! file.exists(soft.file)) {
    downloadAccessionInfo(accession = accession, dest.dir = dest.dir, scope, amount)
  }
  
  txt <- readSOFTFileFun(soft.file, method = 'fread')
  
  if (sum(grepl("SuperSeries of", txt$V2))>0) {
    print(accession)
    super_series <- c(super_series, accession)
  }
  
}

super_series

series.filtered$Super_Series <- ifelse(series.filtered$Series %in% super_series, "Y", "N")

View(data.frame(accessions, ifelse(accessions %in% super_series, "Y", "N")))

sort(table(series.filtered$Chip), decreasing = T)

names(sort(table(series.filtered$Chip), decreasing = T))

chip_to_package <- c(
  "HG-U133 Plus 2.0"        = "pd.hgu133plus2.hs.gencodeg",
  "MG-430 2.0"              = "pd.mouse4302.mm.gencodeg",
  "Mouse Gene 1.0 ST"       = "pd.mogene10st.mm.gencodeg",
  "Human Gene 1.0 ST"       = "pd.hugene10st.hs.gencodeg",
  "Mouse Gene 2.0 ST"       = "pd.mogene20st.mm.gencodeg",
  "HG-U133A"                = "pd.hgu133a.hs.gencodeg",
  "Human Gene 2.0 ST"       = "pd.hugene20st.hs.gencodeg",
  "HTA 2.0"                 = "pd.hta20.hs.gencodeg",
  "Human Exon 1.0 ST"       = "pd.huex10st.hs.gencodeg",
  "HG-U133A 2.0"            = "pd.hgu133a2.hs.gencodeg",
  "Clariom S Human"         = "pd.clariomshuman.hs.gencodeg",
  "MG-U74Av2"               = "pd.mgu74av2.mm.gencodeg",
  "Human PrimeView"         = "pd.primeview.hs.gencodeg",
  "Clariom S Mouse"         = "pd.clariomsmouse.mm.gencodeg",
  "MG-430A 2.0"             = "pd.mouse430a2.mm.gencodeg",
  "MG-430A"                 = "pd.moe430a.mm.gencodeg",
  "Mouse Exon 1.0 ST"       = "pd.moex10st.mm.gencodeg",
  "MTA 1.0"                 = "pd.mta10.mm.gencodeg",
  "Clariom D Human"         = "pd.clariomdhuman.hs.gencodeg",
  "HG-U95Av2"               = "pd.hgu95av2.hs.gencodeg",
  "HT MG-430 PM"            = "pd.htmg430pm.mm.gencodeg",
  "HG-U219"                 = "pd.hgu219.hs.gencodeg",
  "Human Gene 2.1 ST"       = "pd.hugene21st.hs.gencodeg",
  "Mouse Gene 2.1 ST"       = "pd.mogene21st.mm.gencodeg",
  "Human Gene 1.1 ST"       = "pd.hugene11st.hs.gencodeg",
  "Mouse Gene 1.1 ST"       = "pd.mogene11st.mm.gencodeg",
  "HG-U133B"                = "pd.hgu133b.hs.gencodeg",
  "HT HG-U133 Plus PM"      = "pd.hthgu133pluspm.hs.gencodeg",
  "MG-430B"                 = "pd.moe430b.mm.gencodeg",
  "HG-Focus"                = "pd.hgfocus.hs.gencodeg",
  "HT HG-U133A"             = "pd.hthgu133a.hs.gencodeg",
  "MG-U74Bv2"               = "pd.mgu74bv2.mm.gencodeg",
  "MG-U74Cv2"               = "pd.mgu74cv2.mm.gencodeg",
  "HG-U95A"                 = "pd.hgu95a.hs.gencodeg",
  "Clariom S Human HT"      = "pd.clariomshumanht.hs.gencodeg",
  "Clariom S Mouse HT"      = "pd.clariomsmouseht.mm.gencodeg",
  "HG-U133 X3P"             = "pd.u133x3p.hs.gencodeg",
  "HT MG-430A"              = "pd.htmg430a.mm.gencodeg",
  "MG-U74A"                 = "pd.mgu74a.mm.gencodeg",
  "HG-U133A AofAv2"         = "pd.u133aaofav2.hs.gencodeg",
  "HG-U95B"                 = "pd.hgu95b.hs.gencodeg",
  "HG-U95C"                 = "pd.hgu95c.hs.gencodeg",
  "HG-U95D"                 = "pd.hgu95d.hs.gencodeg",
  "HG-U95E"                 = "pd.hgu95e.hs.gencodeg",
  "Clariom GO Screen Human" = "pd.goscreenhu.hs.gencodeg",
  "HT HG-U133B"             = "pd.hthgu133b.hs.gencodeg",
  "Human Almac Xcel"        = "pd.xcel.hs.gencodeg",
  "Clariom D Mouse"         = "pd.clariomdmouse.mm.gencodeg",
  "HT MG-430B"              = "pd.htmg430b.mm.gencodeg"
)

series.filtered$Package <- chip_to_package[series.filtered$Chip]

View(series.filtered)

saveRDS(series.filtered, 
        file = "~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/final_series_list.hs_mm_affymetrix.RDS")


#######

library(GEOquery)
library(httr)

# series_table <- read.csv("~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/gse_series_summary_table.csv",
#                          header = T, stringsAsFactors = F)
# 
# dim(series_table)

samples.filtered <- readRDS(file = "~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/final_samples_list.hs_mm_affymetrix.RDS")

series.filtered <- readRDS(file = "~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/final_series_list.hs_mm_affymetrix.RDS")
series.filtered <- series.filtered %>% filter(Super_Series == "N")

dim(series.filtered) # 25129

series <- paste0(series.filtered$Series, "-", series.filtered$Platform)
sum(paste0(samples.filtered$Series, "-", samples.filtered$Platform) %in% series) # 683738


stats <- samples.filtered %>% group_by(Series, Platform) %>%
  summarise(N=length(Accession))

idx <- match(paste0(series.filtered$Series, "-", series.filtered$Platform),
      paste0(stats$Series, "-", stats$Platform))

series.filtered$N_Samples <- stats$N[idx]

series.filtered <- series.filtered %>% arrange(N_Samples)
series.filtered


pd_pkgs <- unique(series.filtered$Package)
pd_pkgs

for(pkg in pd_pkgs) {
  library(pkg, character.only = TRUE)
}


base.dir <- "/Users/rli012/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/microarray/affymetrix"

#View(series.filtered)

# View(series.filtered[which(series.filtered$Series %in% series.filtered$Series[which(duplicated(series.filtered$Series))]),])

completed <- gsub("_Normalized_Data.rds", "", list.files(file.path(base.dir, "normalized")))
#completed

remaining <- which(!paste0(series.filtered$Series, "-", series.filtered$Platform) %in% completed)
remaining

which(series.filtered$Series[remaining] == "GSE117468")
which(remaining==24527)

series.filtered$Series[remaining[1665]]

# 23923 huge

for (i in remaining[-c(1:1888)]) {
  
  gse <- series.filtered$Series[i]
  anno <- series.filtered$Package[i]
  gpl <- series.filtered$Platform[i]
  
  taxonomy <- series.filtered$Taxonomy[i]
  
  message(sprintf("----- # %s: %s - %s -----", i, gse, gpl))
  print(sprintf("%s - %s", gpl, anno))
  
  if (file.exists(file.path(base.dir, "normalized", sprintf("%s-%s_Normalized_Data.rds", gse, gpl)))) {
    print ("Expression data has already been calculated !")
    next
  }
  
  gse_nnn <- gsub("(\\d{3})$", "nnn", gse)
  
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", gse_nnn, gse)
  response <- HEAD(url)
  
  if (status_code(response) != 200) {
    message("supp is not available!")
    next
  }
  
  tar_url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s_RAW.tar", gse_nnn, gse, gse)
  response <- HEAD(tar_url)
  
  if (status_code(response) != 200) {
    message("RAW is not available!")
    next
  }
  
  size_mb = as.numeric(headers(response)$`content-length`) / (1024^2)
  print(size_mb)
  
  if (size_mb > 2000) {
    message("Large file, skip!")
    next
  }
  
  gse.samples <- samples.filtered %>% filter(Series==gse, Platform==gpl) %>% pull(Accession)
  # print(sprintf("GSE samples - %s", length(gse.samples)))
  
  filePaths = getGEOSuppFiles(gse, baseDir = file.path(base.dir, "raw"), makeDirectory = FALSE, filter_regex = 'RAW')
  untar(file.path(base.dir, "raw", paste0(gse, '_RAW.tar')), exdir = file.path(base.dir, "raw", paste0(gse, '_RAW')))
  
  celFiles = list.celfiles(file.path(base.dir, "raw", paste0(gse, '_RAW')), full.names=T, listGzipped=T)
  
  # print(sprintf("CEL samples - %s", length(celFiles)))
  
  .samples <- unlist(lapply(basename(celFiles), function(x) strsplit(x, '_|\\.')[[1]][1]))
  
  idx <- which(.samples %in% gse.samples)
  celFiles <- celFiles[idx]
  
  # print(sprintf("CEL samples - %s", length(celFiles)))
  
  print(sprintf("Number of samples: %s", length(celFiles)))
  
  if (length(celFiles) != length(gse.samples) | length(celFiles) == 0) {
    print("Sample size is not consistent !")
    next
  }
  
  rawData = read.celfiles(celFiles, pkgname = anno, verbose = F)
  
  probesetData = oligo::rma(rawData)
  
  exprData = exprs(probesetData)
  
  rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '_|\\.')[[1]][1]))
  colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '_|\\.')[[1]][1]))
  
  if (taxonomy == "Mus musculus") {
    ensembl = "ENSMUSG"
  } else if (taxonomy == "Homo sapiens") {
    ensembl = "ENSG"
  } 
  
  filter <- which(!startsWith(rownames(exprData), ensembl))
  
  if (length(filter) > 0) {
    exprData <- exprData[-filter,]
  }
  
  if (sum(is.na(expr))>0) {
    next
  }
  
  dim(exprData)
  
  saveRDS(exprData, file=file.path(base.dir, "normalized", sprintf("%s-%s_Normalized_Data.rds", gse, gpl)))
  
  exprData <- NULL
  rawData <- NULL
  
}


# expr <- readRDS("data/microarray/affymetrix/normalized/GSE180556-GPL6246_Normalized_Data.rds")
View(expr)


### Illumina










# filter <- which(duplicated(paste0(samples.filtered$Series, "-", samples.filtered$Platform)))
# filter
# 
# View(samples.filtered[-filter,])
# 
# samples.filtered <- samples.filtered[-filter,]
# View(samples.filtered[samples.filtered$Series=="GSE153298",])
# 
# series.filtered1 <- samples.filtered %>% filter(grepl(";", Series)) %>%
#   separate_rows(Series, sep = ";") %>%
#   group_by(Series) %>%
#   filter(n() == 1) %>%
#   ungroup() %>% 
#   filter(Platform %in% platforms.affy$Accession)
#   
# series.filtered1
# 
# series.filtered2 <- samples.filtered %>% filter(!grepl(";", Series)) %>% 
#   filter(Platform %in% platforms.affy$Accession)
# 
# View(series.filtered2)
# 
# series.filtered <- bind_rows(series.filtered1, series.filtered2)
# 
# dim(series.filtered)
# 
# # samples.filtered <- samples.filtered[-filter,] %>%
# #   separate_rows(Series, sep = ";") #%>%
# #   group_by(Series) %>%
# #   filter(n() > 1) %>%
# #   ungroup() %>%
# #   filter(Platform %in% platforms.affy$Accession)
# 
# which(duplicated(series.filtered$Series))
# 
# 
# #samples.filtered <- samples %>% filter(Platform == gpl)
# #samples.filtered
# 
# View(series.filtered)
# View(samples)
# 
# idx <- match(series.filtered$Platform, platforms.affy$Accession)
# 
# series.filtered <- data.frame(series.filtered, Chip=platforms.affy$Chip[idx], Platform.Organism=platforms.affy$Taxonomy[idx])
# series.filtered
# 
# platforms.affy
# 
# series.filtered <- series.filtered %>% filter(Taxonomy == Platform.Organism)
# 
# # gpl <- platforms.affy %>% filter(Chip == "HG-U133 Plus 2.0") %>% 
# #   pull(Accession)
# # 
# # gpl <- platforms.affy$Accession[1]
# 
# 
# ### 
# dest.dir <- "~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/microarray/gpl"
# scope <- "Platform"
# amount = 'brief'
# 
# accession <- platforms.affy$Accession[1] ## GPL570
# accession
# 
# soft.file <- file.path(dest.dir, paste0(accession, '.', scope, '.txt'))
# 
# if (! file.exists(soft.file)) {
#   downloadAccessionInfo(accession = accession, dest.dir = dest.dir, scope, amount)
# }
# 
# txt <- readSOFTFileFun(soft.file, method = 'fread')
# gse_accessions <- txt %>% filter(V1=="!Platform_series_id") %>% pull(V2)
# 
# gse_accessions
# 
# 
# series <- series.filtered$Series[series.filtered$Platform==accession]
# series
# 
# gse_accessions <- gse_accessions[which(! gse_accessions %in% unlist(strsplit(series, ";")))]
# gse_accessions
# 
# dest.dir <- "~/Documents/BioinfoCamp/projects/cloud_platform/OmicsMine/backend/data/microarray/affymetrix/metadata/"
# 
# scope <- "Series"
# amount = 'brief'
# 
# super_series <- c()
# 
# for (accession in gse_accessions) {
#   
#   soft.file <- file.path(dest.dir, paste0(accession, '.', scope, '.txt'))
#   
#   if (! file.exists(soft.file)) {
#     downloadAccessionInfo(accession = accession, dest.dir = dest.dir, scope, amount)
#   }
#   
#   txt <- readSOFTFileFun(soft.file, method = 'fread')
#   
#   if (sum(grepl("SuperSeries of", txt$V2))>0) {
#     print(accession)
#     super_series <- c(super_series, accession)
#   }
#   
# }
# 
# super_series
# 
# gse_accessions[which(! gse_accessions %in% super_series)]

###














gse <- 'GSE12378'
gse <- 'GSE30521'
anno <- 'pd.huex10st.hs.gencodeg'

gse <- 'GSE2443'
anno <- 'pd.hgu133a.hs.gencodeg'

gse <- 'GSE26910' # 12 samples
gse <- 'GSE32448'
anno <- 'pd.hgu133plus2.hs.gencodeg'

gse <- 'GSE6956'
gse <- 'GSE7055'
anno <- 'pd.hgu133a2.hs.gencodeg'

filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE, filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)

rawData = read.celfiles(celFiles, pkgname = anno)

if (gse=='GSE26910') {
  rawData = read.celfiles(celFiles[1:12], pkgname = anno)
}

probesetData = oligo::rma(rawData)

exprData = exprs(probesetData)
exprData[1:5,1:5]

rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '_|\\.')[[1]][1]))
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '_|\\.')[[1]][1]))

filter <- which(!startsWith(rownames(exprData), 'ENSG'))
filter

if (length(filter) > 0) {
  exprData <- exprData[-filter,]
}

dim(exprData)

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression.RDS'))










### Illumina Arrays

### select most informative probe, MAX IQR

selectProbeFun <- function(expr) {
  expr$IQR <- apply(expr[,-which(colnames(expr)=='ID')], 1, IQR)
  
  expr <- expr %>% group_by(ID) %>% 
    filter(row_number() == which.max(IQR)) %>%
    column_to_rownames('ID')
  
  expr <- expr[,-which(colnames(expr)=='IQR')]
  
  return(expr)
  
}


library(dplyr)
library(tibble)

pcadb.anno <- readRDS('data/Homo_Sapiens_Gene_Annotation_ENSEMBL_HGNC_ENTRE.RDS')
table(pcadb.anno$chromosome)

# BiocManager::install('illuminaHumanv3.db')
library(illuminaHumanv4.db)
#columns(illuminaHumanv4.db)

library(illuminaHumanv3.db)
#columns(illuminaHumanv3.db)

IlluminaFun <- function(exprData, anno='illuminaHumanv4.db') {
  
  if (!exists('pcadb.anno')) {
    pcadb.anno <- readRDS('data/Homo_Sapiens_Gene_Annotation_ENSEMBL_HGNC_ENTRE.RDS')
  }
  
  if (anno=='illuminaHumanv4.db') {
    
    illumina.anno <- AnnotationDbi::select(illuminaHumanv4.db,
                                           keys = rownames(exprData),
                                           columns=c('ENSEMBL', 'GENETYPE', "SYMBOL","ENTREZID"), # "SYMBOL","ENTREZID",
                                           keytype="PROBEID")
    
    
  } else if (anno=='illuminaHumanv3.db') {
    
    illumina.anno <- AnnotationDbi::select(illuminaHumanv3.db,
                                           keys = rownames(exprData),
                                           columns=c('ENSEMBL', 'GENETYPE', "SYMBOL","ENTREZID"), # "SYMBOL","ENTREZID",
                                           keytype="PROBEID")
    
  }
  
  idx <- which(startsWith(illumina.anno$ENSEMBL, 'ENSG'))
  illumina.anno <- illumina.anno[idx,]
  
  filter <- which(!illumina.anno$ENSEMBL %in% pcadb.anno$ensembl_id)
  illumina.anno <- illumina.anno[-filter,]
  
  illumina.anno <- illumina.anno %>% group_by(PROBEID) %>% 
    filter(if (length(GENETYPE)>1) GENETYPE == 'protein-coding' else !is.na(GENETYPE))
  
  ### optional (if a probe mapped to multiple genes with some of them are novel proteins)
  
  illumina.anno$PROBE_ENTREZ_ID <- paste0(illumina.anno$PROBEID, '_', illumina.anno$ENTREZID)
  illumina.anno$PROBE_ENTREZ_ID
  
  dup.probe.entrez <- illumina.anno$PROBE_ENTREZ_ID[which(duplicated(illumina.anno$PROBE_ENTREZ_ID))]
  dup.probe.entrez
  
  idx <- which(illumina.anno$PROBE_ENTREZ_ID %in% dup.probe.entrez)
  idx
  
  illumina.anno$ENSEMBL_ENTREZ_ID <- paste0(illumina.anno$ENSEMBL, '_', illumina.anno$ENTREZID)
  illumina.anno$ENSEMBL_ENTREZ_ID
  
  i <- match(illumina.anno$ENSEMBL_ENTREZ_ID[idx], paste0(pcadb.anno$ensembl_id, '_', pcadb.anno$entrez_id))
  i <- i[which(!is.na(pcadb.anno$tax_id[i]))]
  
  filter <- idx[which(!illumina.anno$ENSEMBL_ENTREZ_ID[idx] %in% paste0(pcadb.anno$ensembl_id, '_', pcadb.anno$entrez_id)[i])]
  
  illumina.anno <- illumina.anno[-filter,]
  
  #### probe mapped to multiple genes
  
  probes <- illumina.anno$PROBEID[which(duplicated(illumina.anno$PROBEID))]
  filter <- which(illumina.anno$PROBEID %in% probes)
  illumina.anno <- illumina.anno[-filter,]
  illumina.anno
  
  
  ###
  exprData <- data.frame(exprData[illumina.anno$PROBEID,])
  exprData$ID <- illumina.anno$ENSEMBL
  
  idx <- which(!is.na(rowSums(exprData[,-ncol(exprData)])))
  exprData <- selectProbeFun(exprData[idx,])
  
  if (max(exprData) > 100) {
    exprData <- log2(exprData)
  }
  
  rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '_|\\.', )[[1]][1]))
  
  return (exprData)
  
}

gse <- 'GSE28680' # Illumina HumanHT-12 V4.0 Expression Beadchip
gse <- 'GSE29650' # Illumina HumanHT-12 V3.0 Expression Beadchip
gse <- 'GSE32571' # Illumina HumanHT-12 V3.0 Expression Beadchip

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
i.matrix <- 1
exprData <- exprs(seriesMatrix[[i.matrix]])
exprData

exprData <- IlluminaFun(exprData, anno = 'illuminaHumanv4.db')
exprData

exprData <- IlluminaFun(exprData, anno = 'illuminaHumanv3.db')
exprData

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression.RDS'))





