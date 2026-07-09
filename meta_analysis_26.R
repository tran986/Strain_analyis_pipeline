library(dplyr)
library(tidyverse)
library(readxl)
library(readr)
library(purrr)
library(dplyr)
library(tidyverse)
library(gdata)
library(devtools)
devtools::load_all("/Users/tran.986/Desktop/phylogenize/package/repermulize")
devtools::load_all("/Users/tran.986/Desktop/phylogenize/package/phylogenize")

md = read_excel("~/Desktop/meta_analysis_26/metadata_with_cond.xls",
                          sheet = 1) # use this as a standard (FROM THE PAPER)
#remove T1D:
md = md[md$Status != "T1D", ]
View(md)
#---------Chinese dataset:
md_CHNs = read_csv("~/Desktop/meta_analysis_26/metadata_CHNs.csv") #CHNs only
#merge info to obtain sequence list for wget downloads:
md_CHNs["Sample"]<-substring(md_CHNs$`Sample Name`, 5)
CHNs_md_final = md %>% left_join(md_CHNs, by = "Sample") %>%
  filter(!is.na(Run))
write.csv(CHNs_md_final, "~/Desktop/meta_analysis_26/metadata_CHNs_final.csv")

#---------MHN dataset1:--Study: ERP004605
md_mhn1 = read_tsv("~/Desktop/meta_analysis_26/filereport_MHN_1.tsv")

md_mhn1$Sample <- str_extract(md_mhn1$submitted_ftp, "[^/]+$") %>%
  str_extract("[^/]+$") %>%   # keep after last /
  str_extract("^[^_]+") %>% #keep before _
  str_extract("^[^-]+")  

mhn1_md_final = md_mhn1 %>%
  left_join(md, by = c("Sample")) %>%
  filter(!is.na(Status)) %>%
  select(Sample, `Country subset`, Status, fastq_ftp)
write.csv(mhn1_md_final, "~/Desktop/meta_analysis_26/mhn1_md_final.csv")

#---------MHN dataset2: - only have ND Ctrl: --Study: ERP003612
md_mhn2 = read_tsv("~/Desktop/meta_analysis_26/filereport_MHN_2.tsv")
md_mhn2$Sample <- str_extract(md_mhn2$submitted_ftp, "[^/]+$") %>%
  str_extract("[^/]+$") %>%   # keep after last /
  str_extract("(?<=MetaHIT-)[^_]+")

mhn2_md_final = md_mhn2 %>% 
  left_join(md, by = c("Sample")) %>% 
  filter(!is.na(Status)) %>%
  select(Sample, `Country subset`, Status, fastq_ftp)
write.csv(mhn2_md_final, "~/Desktop/meta_analysis_26/mhn2_md_final.csv")

#---------MHN dataset3: --Study: ERP002469 
md_mhn3 = read_tsv("~/Desktop/meta_analysis_26/filereport_MHN_3.tsv")
md_mhn3$Sample <- str_extract(md_mhn3$submitted_ftp, "[^/]+$") %>%
  str_extract("[^/]+$") %>%
  str_extract("^[^_]+_[^_]+")

mhn3_md_final = md_mhn3 %>% 
  left_join(md, by = c("Sample")) %>% 
  filter(!is.na(Status)) %>%
  select(Sample, `Country subset`, Status, fastq_ftp)
write.csv(mhn3_md_final, "~/Desktop/meta_analysis_26/mhn3_md_final.csv")


#--------MHN dataset4: --Study: ERP002061 -- all ND CTRL
md_mhn4 = read_tsv("~/Desktop/meta_analysis_26/filereport_MHN_4.tsv")
md_mhn4$Sample <- str_extract(md_mhn4$submitted_ftp, "[^/]+$") %>%
  str_extract("[^/]+$") %>%
  str_extract("^[^_]+") %>% #keep before _
  str_extract("^[^-]+")

mhn4_md_final = md_mhn4 %>% 
  left_join(md, by = c("Sample")) %>% 
  filter(!is.na(Status)) %>%
  select(Sample, `Country subset`, Status, fastq_ftp)
write.csv(mhn4_md_final, "~/Desktop/meta_analysis_26/mhn4_md_final.csv")


#--------
#Default Functions:
#1. function to get url to run on HPC:
#e.g., study_id = CHN, ERP004605, ERP003612, ERP002061, ERP002469
wget_urls_func <- function(md_final, study_id) {
  
  out_file <- paste0("~/Desktop/meta_analysis_26/fastq_url/", study_id, "_fastq_urls.txt")
  
  wget_urls <- c(
    paste0(
      "https://ftp.sra.ebi.ac.uk/vol1/fastq/",
      substr(md_final$Run, 1, 6), "/",
      md_final$Run, "/",
      md_final$Run,
      "_1.fastq.gz"
    ),
    paste0(
      "https://ftp.sra.ebi.ac.uk/vol1/fastq/",
      substr(md_final$Run, 1, 6), "/",
      md_final$Run, "/",
      md_final$Run,
      "_2.fastq.gz"
    )
  )
  
  writeLines(wget_urls, out_file)
}

#2. function to read into .bracken files (push this to ASC):
import_bracken=function(study_id) { 
  file = list.files(path=paste0(working_dir, "/abund/", study_id))
  data = map(file, ~ { read_tsv(file.path(paste0(working_dir, "/abund/", study_id), .), show_col_types=FALSE) })
  samples = gsub("\\.bracken","", file)
  data_newcol = map2(data, samples, ~ {mutate(.x, sample = .y)})
  data_tidy = bind_rows(data_newcol)
  data_w_count = tidyr::pivot_wider(data_tidy, 
                                    names_from = sample, 
                                    values_from = new_est_reads,
                                    id_cols = name,
                                    values_fill = 0)
  write.table(data_w_count, file = paste0(working_dir, "/phylogenize_out/", study_id, "/data_w_count.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  return(data_w_count)
}

#3. function to create phylogenize input:
extract_phyloz_metadata=function(import_bracken_out, 
                                 metadata,
                                 study_id, #study_id = CHN, MHN1, MHN3, 4, etc
                                 envs_compared = c("ND CTRL", "T2D metformin-")) { 
  #import_bracken_out = test
  #metadata = CHNs_md_final 
  #clean up metadata:
  if (study_id == "CHN") #then it is from CHN sample "SRR"
    
  {
    metadata_clean=metadata[,c("Run", "Status")] %>%
      rename("env"="Status",
             "sample"="Run") %>%
      mutate(dataset = study_id)
    
  } else { #all MHN studies with "ERR" samples
    metadata_clean = metadata %>% 
      mutate(sample = str_extract(fastq_ftp, "(?<=/)[^/]+(?=/[^/]+\\.fastq\\.gz)"),
             dataset = study_id) %>%
      select(Status, sample, dataset) %>%
      rename("env"="Status")
  }
  
  #filter 2 envs that will be compared:
  metadata_filter=metadata_clean %>% filter(env %in% envs_compared)
  
  write.table(metadata_filter, file = paste0(working_dir, "/phylogenize_out/", study_id, "/metadata_filter.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  return(metadata_filter)
  
}


#4. apply phylogenize on inputs: metadata and count
dataset_dir = "/fs/project/bradley.720/projects/phylogenize_v2/phylogenize/package/phylogenize/inst/extdata"
phylogenize_run=function(extract_phyloz_metadata_out,
                         import_bracken_out,
                         study_id,
                         ref_env) { #"ND CTRL or T2D metformin-
  
  res=phylogenize(
    output_file = paste0(working_dir, "/phylogenize_out/", study_id, "/phylogenize.html"),
    out_dir = paste0(working_dir, "/phylogenize_out/", study_id),
    data_dir = dataset_dir,
    db = "human-gut",
    taxon_level = "family",
    type_16S=F,
    which_phenotype = "abundance",
    diff_abund_method = "ancombc2",
    abundance_file = import_bracken_out,
    metadata_file = extract_phyloz_metadata,
    which_envir = ref_env, 
    env_column = "env",
    sample_column = "sample",
    ncl=parallelly::availableCores(),
    core_method = "permutrate-rlm")
  
  saveRDS(res, 
          paste0(working_dir, "/phylogenize_out/", study_id, "/phylogenize.rds"))
  return(res)
}

#5. Function to run the other pipeline (pre-phylogenize merging pipeline):
#from bracken output -> apply ancombc:
ancomRun = function(count_tbl, study_id, metadata){
  
  count = count_tbl %>% column_to_rownames(var = "name")
  metadata_tbl = extract_phyloz_metadata(import_bracken_out = count,
                                         metadata = metadata,
                                         study_id = study_id,
                                         envs_compared = c("ND CTRL", "T2D metformin-"))
  metadata_tbl = metadata_tbl %>% column_to_rownames(var = "sample")
  
  #filter out those samples in the metadata but not in count_tbl:
  count_tbl_clean = count %>% dplyr::select(rownames(metadata_tbl))
  all(colnames(count_tbl_clean) == rownames(metadata_tbl))
  tse <- TreeSummarizedExperiment(assays = S4Vectors::SimpleList(counts = as.matrix(count_tbl_clean)),
                                  colData = S4Vectors::DataFrame(metadata_tbl)
  )
  #apply ancombc2:
  print("fitting ANCOMBC")
  res_ancom=ANCOMBC::ancombc2(data = tse,
                              assay_name = "counts",
                              fix_formula = "env")
  
  return(list(res = res_ancom,
              count_input = count_tbl_clean,
              metadata_input = metadata_tbl))
  saveRDS(res_ancom,
          paste0(working_dir, "/pipeline2/ancom/", study_id, "_ancombc_res.rds"))
}

