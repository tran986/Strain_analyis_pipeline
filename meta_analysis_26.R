library(dplyr)
library(tidyverse)
library(readxl)
library(readr)
library(purrr)
library(dplyr)
library(tidyverse)

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

#---------MHN dataset2: - only have ND Ctrl: --Study: ERP003612
md_mhn2 = read_tsv("~/Desktop/meta_analysis_26/filereport_MHN_2.tsv")
View(md_mhn2)

md_mhn2$Sample <- str_extract(md_mhn2$submitted_ftp, "[^/]+$") %>%
  str_extract("[^/]+$") %>%   # keep after last /
  str_extract("(?<=MetaHIT-)[^_]+")

mhn2_md_final = md_mhn2 %>% 
  left_join(md, by = c("Sample")) %>% 
  filter(!is.na(Status)) %>%
  select(Sample, `Country subset`, Status, fastq_ftp)

#---------MHN dataset3: --Study: ERP002469 
md_mhn3 = read_tsv("~/Desktop/meta_analysis_26/filereport_MHN_3.tsv")
View(md_mhn3)

md_mhn3$Sample <- str_extract(md_mhn3$submitted_ftp, "[^/]+$") %>%
  str_extract("[^/]+$") %>%
  str_extract("^[^_]+_[^_]+")

mhn3_md_final = md_mhn3 %>% 
  left_join(md, by = c("Sample")) %>% 
  filter(!is.na(Status)) %>%
  select(Sample, `Country subset`, Status, fastq_ftp)

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
import_bracken=function(bracken_dir_path) {
  file = list.files(path=bracken_dir_path)
  data = map(file, ~ { read_tsv(file.path(bracken_dir_path, .), show_col_types=FALSE) })
  samples = gsub("\\.bracken","", file)
  data_newcol = map2(data, samples, ~ {mutate(.x, sample = .y)})
  data_tidy = bind_rows(data_newcol)
  data_w_count = pivot_wider(data_tidy, 
                             names_from = sample, 
                             values_from = new_est_reads,
                             id_cols = name,
                             values_fill = 0)
  return(data_w_count)
}

#--------apply default functions:
study_id_ls = c("CHN", "MHN_1_ERP004605", "MHN_2_ERP003612", "MHN_3_ERP002469", "MHN_4_ERP002061")
md_final_ls = list(CHNs_md_final,
                mhn1_md_final,
                mhn2_md_final,
                mhn3_md_final,
                mhn4_md_final)

purrr::map2(md_final_ls, study_id_ls, ~ wget_urls_func(md_final = .x,
                                                       study_id = .y))

test=import_bracken(bracken_dir_path="~/Desktop/meta_analysis_26/bracken_outputs/CHN")
