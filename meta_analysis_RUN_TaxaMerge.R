working_dir = "/Users/tran.986/Desktop/meta_analysis_26"
source(paste0(working_dir, "/meta_analysis_26.R"))

#-----------------------================================-----------Expand to more datasets:
#             HMP_2019_t2d --> HMP (not started)
#         KarlssonFH_2013 --> MH3 (done)
#                LiJ_2014 --> MH1 (done)
#      MetaCardis_2020_a --> MCA (started - in progress) 
#              QinJ_2012 --> CHN (done)
# SankaranarayananK_2015 --> SKK (not started)

#-----------------------================================ORDER OF PROCESSES:
#study 1 count_tbl + metadata --> ANCOMBC --> ashr 
#study 2 count_tbl + metadata --> ANCOMBC --> ashr
#...do for all of the studies --> merge (combined) at taxa signals --> Phylogenize2

#==================================================MERGE count table to run ANCOMBC:
#ource("/fs/project/bradley.720/projects/meta_analysis_26/meta_analysis_26.R")
#working_dir = "/fs/project/bradley.720/projects/meta_analysis_26"

#---merge all count table from all datasets:
#import_ls = c("CHN", "ERP004605_MH1", "ERP002469_MH3", "MCA", "SKK")
#import_bracken_ls = setNames(lapply(import_ls, function(i) {
#  read_tsv(paste0(working_dir, "/mergeGene/phylogenize_out/", i, "/data_w_count.tsv"))
#}), import_ls)

#count_tbl_merge = Reduce(function(x,y) merge(x, y, by = "name", all = T), import_bracken_ls)

#count_tbl_merge[-1]<- lapply(count_tbl_merge[-1], function(x){
#  x_num <- suppressWarnings(as.numeric(as.character(x)))
#  x_num[is.na(x_num) & !is.na(x)] <- 0
#  x_num
#})

#count_tbl_merge[is.na(count_tbl_merge)] <- 0
#print("any non-numeric value:")
#anyNA(count_tbl_merge[-1])


#write.table(count_tbl_merge, 
#            file = paste0(working_dir, "/mergeTaxa/ancom/combine_fix/count_tbl_merged.tab"),
#            sep = "\t",
#            row.names = F,
#            quote = F)

#--merge all metadata
metadata_ls=setNames(lapply(import_ls, function(i) {
  df = read_tsv(paste0(working_dir, "/mergeGene/phylogenize_out/", i, "/metadata_filter.tsv")) |> dplyr::mutate(dataset = i)
}), import_ls)

metadata_tbl_merge = bind_rows(metadata_ls) |> dplyr::select(sample, env, dataset)
write.table(metadata_tbl_merge, 
            file = paste0(working_dir, "/mergeTaxa/ancom/combine_fix/metadata_merged.tab"),
            sep = "\t", row.names = FALSE)

#--run ANCOMBC:
#ancomMerge = ancomRun(count_tbl = count_tbl_merge,
#	 study_id = NULL,
#	 metadata = metadata_tbl_merge)

ancomMerge = readRDS(paste0(working_dir, "/mergeTaxa/ancom/combine_fix/merged_ancom_res.rds"))

#--run ashr:
print("fitting ashr")
ancomMerge_res = ancomMerge$res
ashMerge = ashr::ash(betahat = ancomMerge_res$`lfc_envT2D metformin-`,
                     sebetahat = ancomMerge_res$`se_envT2D metformin-`) 

#rename ash res to become a "provided" file for phylogenize
ashMerge_res = ashMerge$res
ashMerge_res$taxon <- ancomMerge_res$taxon
ashMerge_res = ashMerge_res %>%
  dplyr::select(taxon, PosteriorMean, PosteriorSD) %>%
  dplyr::rename("estimate"="PosteriorMean",
                "stderr"="PosteriorSD")

write.table(ashMerge_res,
            paste0(working_dir, "/mergeTaxa/ash/ashMerge_res.tab"),            
            sep = "\t", row.names = FALSE)   #this is "provided" file 

#--run phylogenize2:
phylogenize_run(provided_file_path = paste0(working_dir, "/mergeTaxa/ash/ashMerge_res.tab"),
                study_id = NULL,
                phenotype = "provided",
                ref_env = "ND CTRL")







