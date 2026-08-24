#working_dir = "/fs/project/bradley.720/projects/meta_analysis_26"
working_dir = "/Users/tran.986/Desktop/meta_analysis_26"
source(paste0(working_dir, "/meta_analysis_26.R"))

#-----------------------================================-----------Expand to more datasets:
#             HMP_2019_t2d --> HMP (not started)
#         KarlssonFH_2013 --> MH3 (done)
#                LiJ_2014 --> MH1 (done)
#      MetaCardis_2020_a --> MCA (started - in progress) 
#              QinJ_2012 --> CHN (done)
# SankaranarayananK_2015 --> SKK (not started)

#---CHN only:
#phylogenize_full_func(study_id = "CHN",
#                      metadata_dir_path = paste0(working_dir,"/metadata/CHN_md_final.csv"),
#                      ref_env = "T2D metformin-",
#                      envs_compared = c("ND CTRL", "T2D metformin-"))	

#---ERP002469_MH3 only:
#test = import_bracken(study_id = "ERP002469_MH3")
#metadata = read_csv(paste0(working_dir,"/metadata/mhn3_md_final.csv"))
#test_md = extract_phyloz_metadata(import_bracken_out = test,
#				  metadata = metadata,
#				  study_id = "ERP002469_MH3",
#				  envs_compared = c("ND CTRL", "T2D metformin-"))
#test_md
#phylogenize_full_func(study_id = "ERP002469_MH3",
#		      metadata_dir_path = paste0(working_dir,"/metadata/mhn3_md_final.csv"),
#		      ref_env = "T2D metformin-",
#		      envs_compared = c("ND CTRL", "T2D metformin-"))

#---ERP004605_MH1 only:
#phylogenize_full_func(study_id = "ERP004605_MH1",
#		      metadata_dir_path = paste0(working_dir,"/metadata/mhn1_md_final.csv"),
#		      ref_env = "T2D metformin-",
#		      envs_compared = c("ND CTRL", "T2D metformin-"))

#---MCA only:
#metadata_MCA_ctrl = metadataRetrieve(study_id = "MCA")[["healthy"]]
#metadata_MCA_t2d = metadataRetrieve(study_id = "MCA")[["T2D"]]
#ctrl_MCA_seqid=ftpRetrieve(metadata_MCA_ctrl)
#t2d_MCA_seqid=ftpRetrieve(metadata_MCA_t2d)

#save fastq_df to run download.sh
#MCA_seqid = rbind(ctrl_MCA_seqid$fastq_df,
#                  t2d_MCA_seqid$fastq_df)

#write.table(MCA_seqid[c("fastq_ftp")], file = paste0(working_dir,"/fastq_url/MCA_fastq_url.txt"), sep = "\t", row.names = FALSE)
#process $sample_list to run extract_phyloz_metadata --> phylogenize_full_func:
#MCA_md_final=rbind(ctrl_MCA_seqid$sample_list["sample_id"] |> dplyr::mutate(env="ND CTRL"),
#                   t2d_MCA_seqid$sample_list["sample_id"] |> dplyr::mutate(env = "T2D metformin-")) |>
#             dplyr::rename("sample"="sample_id")

#write.csv(MCA_md_final, file = paste0(working_dir, "/metadata/MCA_md_final.csv"))
#phylogenize_full_func(study_id = "MCA",
#                      metadata_dir_path = paste0(working_dir, "/metadata/MCA_md_final.csv"),
#                      ref_env = "T2D metformin-",
#                      envs_compared = c("ND CTRL", "T2D metformin-"))

#---SKK only:
#metadata_SKK_ctrl = metadataRetrieve(study_id = "SKK")[["healthy"]]
#metadata_SKK_t2d = metadataRetrieve(study_id = "SKK")[["T2D"]]
#ctrl_SKK_seqid=ftpRetrieve(metadata_SKK_ctrl)
#t2d_SKK_seqid=ftpRetrieve(metadata_SKK_t2d)

#save fastq_df to run download.sh
#SKK_seqid = rbind(ctrl_SKK_seqid$fastq_df,
#                  t2d_SKK_seqid$fastq_df)

#write.table(SKK_seqid[c("fastq_ftp")], file = paste0(working_dir,"/fastq_url/SKK_fastq_url.txt"), sep = "\t", row.names = FALSE)

#SKK_md_final=rbind(ctrl_SKK_seqid$sample_list["sample_id"] |> dplyr::mutate(env="ND CTRL"),
#                   t2d_SKK_seqid$sample_list["sample_id"] |> dplyr::mutate(env = "T2D metformin-")) |>
#             dplyr::rename("sample"="sample_id")
#write.csv(SKK_md_final, paste0(working_dir, "/metadata/SKK_md_final.csv"))                   
#phylogenize_full_func(study_id = "SKK",
#                      metadata_dir_path = paste0(working_dir, "/metadata/SKK_md_final.csv"),
#                      ref_env = "T2D metformin-",
#                      envs_compared = c("ND CTRL", "T2D metformin-"))

#---MetaCardis -- Molinaro et al: Imidazole propionate is increased in diabetes and associated with dietary patterns and altered microbial ecology
metadata_MCA_ctrl = metadataRetrieve(study_id = "MCA")[["healthy"]]
metadata_MCA_t2d = metadataRetrieve(study_id = "MCA")[["T2D"]]
ctrl_MCA_seqid=ftpRetrieve(metadata_MCA_ctrl)
t2d_MCA_seqid=ftpRetrieve(metadata_MCA_t2d)


#save fastq_df to run download.sh
MCA_seqid = rbind(ctrl_MCA_seqid$fastq_df,
                  t2d_MCA_seqid$fastq_df)

#write.table(MCA_seqid[c("fastq_ftp")], file = paste0(working_dir,"/fastq_url/MCA_fastq_url.txt"), sep = "\t", row.names = FALSE)
#process $sample_list to run extract_phyloz_metadata --> phylogenize_full_func:
MCA_md_final=rbind(ctrl_MCA_seqid$sample_list["sample_id"] |> dplyr::mutate(env="ND CTRL"),
                   t2d_MCA_seqid$sample_list["sample_id"] |> dplyr::mutate(env = "T2D metformin-")) |>
  dplyr::rename("sample"="sample_id")

write.csv(MCA_md_final, file = paste0(working_dir, "/metadata/MCA_md_final.csv"))
phylogenize_full_func(study_id = "MCA",
                      metadata_dir_path = paste0(working_dir, "/metadata/MCA_md_final.csv"),
                      ref_env = "T2D metformin-",
                      envs_compared = c("ND CTRL", "T2D metformin-"))


#---SankaranarayananK_2015: Gut Microbiome Diversity among Cheyenne and Arapaho Individuals from Western Oklahoma
metadata_SKK_ctrl = metadataRetrieve(study_id = "SKK")[["healthy"]]
metadata_SKK_t2d = metadataRetrieve(study_id = "SKK")[["T2D"]]
ctrl_SKK_seqid=ftpRetrieve(metadata_SKK_ctrl)
t2d_SKK_seqid=ftpRetrieve(metadata_SKK_t2d)

#save fastq_df to run download.sh
SKK_seqid = rbind(ctrl_SKK_seqid$fastq_df,
                  t2d_SKK_seqid$fastq_df)

#write.table(SKK_seqid[c("fastq_ftp")], file = paste0(working_dir,"/fastq_url/SKK_fastq_url.txt"), sep = "\t", row.names = FALSE)

SKK_md_final=rbind(ctrl_SKK_seqid$sample_list["sample_id"] |> dplyr::mutate(env="ND CTRL"),
                   t2d_SKK_seqid$sample_list["sample_id"] |> dplyr::mutate(env = "T2D metformin-")) 
write.csv(SKK_md_final, paste0(working_dir, "/metadata/SKK_md_final.csv"))                   
#phylogenize_full_func(study_id = "SKK",
#                      metadata_dir_path = paste0(working_dir, "/metadata/SKK_md_final.csv"),
#                      ref_env = "T2D metformin-",
#                      envs_compared = c("ND CTRL", "T2D metformin-"))

#------------------------------ANALYSIS OF FINAL BIOLOGICAL HITS:
core_out=readRDS(paste0(working_dir, "/core_output.rds"))
annotation_df=core_out$list_pheno$pz.db$gene.to.fxn
gene_pres=core_out$list_pheno$pz.db$gene.presence
tree=core_out$list_pheno$pz.db$trees

#----------pipeline 2 MERGE AT GENE LEVEL: (run each dataset separately) results import:
#calculate for Faith's Diversity:
#--CHN:
CHN_all_res=readRes_n_pdCal(study_id = "CHN")
CHN_anno=left_join(CHN_all_res, annotation_df, by = "gene")


#ERP004605_MH1:
MH1_all_res=readRes_n_pdCal(study_id = "MH1")
#MH1_anno=left_join(MH1_all_res, annnotation_df, by = "gene")


#ERP002469_MH3:
MH3_all_res=readRes_n_pdCal(study_id = "MH3")
#MH3_anno=left_join(MH3_all_res, annnotation_df, by = "gene")

#MCA:
MCA_all_res = readRes_n_pdCal(study_id = "MCA")

#SKK:
SKK_all_res = readRes_n_pdCal(study_id = "SKK")

#======================================= Pipeline 2A
#-------step 1: Recover SE from estimate and p-values - do it per gene/per taxa for each study:
study_id_ls = c("CHN", "MH1", "MH3", "MCA", "SKK") #list of md_final csv, add MCK and SKK if needed

#extract DF (# of case and ctrl subjects) for each study
extractDF_study_list=extractDF_study_func(study_id_ls = study_id_ls)

#input: CHN_all_res, MH1_all_res, MH3_all_res, SKK_all_res, MCA_all_res + their id: "MCA", "SKK", "MH1", "MH3", etc
all_res_list = list(CHN = CHN_all_res,
     MH1 = MH1_all_res,
     MH3 = MH3_all_res,
     MCA = MCA_all_res,
     SKK = SKK_all_res)

all_variance = lapply(names(all_res_list), function(study_name) {
  
  #convert p-value to t-statistics from each study df and p-value:
  all_res_list[[study_name]]$t.statistic = tstatCal(
    study_id_res = all_res_list[[study_name]],
    study_id = study_name,
    extractDF_study_id = extractDF_study_list[[study_name]]
  )
  #back-calculate se from effect size and t-statistic
  all_res_list[[study_name]]$se = se_recCal(
    study_id_res = all_res_list[[study_name]],
    tstatCal_out = all_res_list[[study_name]]$t.statistic
  )
  #calculate variance for each study from SE
  all_res_list[[study_name]]$variance = varianceCal(
    se_recCal_out = all_res_list[[study_name]]$se
  )
  
  #add study_name to easy combine and group by later:
  all_res_list[[study_name]]$study_name <- study_name
  
  #add DF for each study:
  all_res_list[[study_name]]$study_DF <-extractDF_study_list[[study_name]]
  
  all_res_list[[study_name]]
  
})
names(all_variance) = names(all_res_list)

#------- step 2: Plug variance, and study DF into Satterthwaite formula:
#combined all_variance with all studies:
combined_variance=bind_rows(all_variance)
satt_df=Satt_DFCal(combined_variance) 

#------- step 3: Compute pooled estimate and pooled SE -> t-statistics -> p-value
pooled_res = satt_df |> mutate(t_pooled = effsize_pooled / se_pooled,
         p_value_satt = 2 * pt(abs(t_pooled), df = Satt_DF, lower.tail = FALSE)) 

#------- step 4: multiple hypothesis correction: - use either BH or qvalue storey method (less conservative):
pooled_res$q_value_satt <-  p.adjust(pooled_res$p_value_satt, method = "BH")

#==================================================== Inspection of why we have so little overlapping genes signals:
#import the .rds from phylogenize2 for each dataset:
merge_core = lapply(study_id_ls, function(i) 
  readRDS(paste0(working_dir, "/core_output_", i,"_mergeGene.rds")))

phenotype = lapply(merge_core, function(i) {
  enframe(i$list_pheno$phenotype_results$phenotype, 
          name = "taxon",
          value = "phenotype values") |>
    left_join(i$list_pheno$pz.db$taxonomy,
              by = c("taxon"="cluster"))
})
  
#step 1: Determine which genes each study "has access to"
detected_taxa = setNames(lapply(phenotype, function(i) {
  i$taxon
}), study_id_ls)
  
#subset gene matrix of Corio and Lachno:
gene_pres_Corio = gene_pres$Coriobacteriaceae
gene_pres_Lachno = gene_pres$Lachnospiraceae

gene_sets_per_study_Corio = setNames(lapply(study_id_ls, function(s) {
  taxa_in_study = detected_taxa[[s]]
  
  # subset to only the columns (taxa) detected in this study
  taxa_in_study_present = intersect(taxa_in_study, colnames(gene_pres_Corio))
  sub_matrix = gene_pres_Corio[, taxa_in_study_present, drop = FALSE]
  
  # a gene is "available" in this study if ANY of its detected taxa carry it
  gene_present = rownames(sub_matrix)[rowSums(sub_matrix > 0) > 0]
  gene_present
}), study_id_ls)

gene_sets_per_study_Lachno = setNames(lapply(study_id_ls, function(s) {
  taxa_in_study = detected_taxa[[s]]
  
  # subset to only the columns (taxa) detected in this study
  taxa_in_study_present = intersect(taxa_in_study, colnames(gene_pres_Lachno))
  sub_matrix = gene_pres_Lachno[, taxa_in_study_present, drop = FALSE]
  
  # a gene is "available" in this study if ANY of its detected taxa carry it
  gene_present = rownames(sub_matrix)[rowSums(sub_matrix > 0) > 0]
  gene_present
}), study_id_ls)

#step 2: Build the gene × study presence table (same as before)
all_genes_Corio = rownames(gene_pres_Corio)  # use full gene list as reference, not just union
all_genes_Lachno = rownames(gene_pres_Lachno)

#Corio:
gene_presence_Corio = data.frame(
  gene = all_genes_Corio,
  CHN = all_genes_Corio %in% gene_sets_per_study_Corio$CHN,
  MH1 = all_genes_Corio %in% gene_sets_per_study_Corio$MH1,
  MH3 = all_genes_Corio %in% gene_sets_per_study_Corio$MH3,
  SKK = all_genes_Corio %in% gene_sets_per_study_Corio$SKK,
  MCA = all_genes_Corio %in% gene_sets_per_study_Corio$MCA
)

#Lachno
gene_presence_Lachno = data.frame(
  gene = all_genes_Lachno,
  CHN = all_genes_Lachno %in% gene_sets_per_study_Lachno$CHN,
  MH1 = all_genes_Lachno %in% gene_sets_per_study_Lachno$MH1,
  MH3 = all_genes_Lachno %in% gene_sets_per_study_Lachno$MH3,
  SKK = all_genes_Lachno %in% gene_sets_per_study_Lachno$SKK,
  MCA = all_genes_Lachno %in% gene_sets_per_study_Lachno$MCA
)

#step 3: Summarize
#Corio
genes_in_all_Corio = gene_presence_Corio %>% filter(CHN & MH1 & MH3 & SKK & MCA) %>% pull(gene) #genes that are present across all studies
genes_dropped_somewhere_Corio = gene_presence_Corio %>% filter(!(CHN & MH1 & MH3 & SKK & MCA)) %>% pull(gene)

gene_presence_Corio %>%
  summarise(
    n_all_5 = sum(CHN & MH1 & MH3 & SKK & MCA),
    n_missing_at_least_one = sum(!(CHN & MH1 & MH3 & SKK & MCA)),
    n_only_CHN = sum(CHN & !MH1 & !MH3 & !SKK & !MCA),
    n_only_MH1 = sum(!CHN & MH1 & !MH3 & !SKK & !MCA),
    n_only_MH3 = sum(!CHN & !MH1 & MH3 & !SKK & !MCA),
    n_only_SKK = sum(!CHN & !MH1 & !MH3 & SKK & !MCA),
    n_only_MCA = sum(!CHN & !MH1 & !MH3 & !SKK & MCA)
  )

upset_input_Corio = gene_presence_Corio %>%
  mutate(across(c(CHN, MH1, MH3, SKK, MCA), as.integer)) %>%
  select(CHN, MH1, MH3, SKK, MCA)

upset(upset_input_Corio, sets = c("CHN", "MH1", "MH3", "SKK", "MCA"), order.by = "freq")

#Lachno
genes_in_all_Lachno = gene_presence_Lachno %>% filter(CHN & MH1 & MH3 & SKK & MCA) %>% pull(gene) #genes that are present across all studies
genes_dropped_somewhere_Lachno = gene_presence_Lachno %>% filter(!(CHN & MH1 & MH3 & SKK & MCA)) %>% pull(gene)

gene_presence_Lachno %>%
  summarise(
    n_all_5 = sum(CHN & MH1 & MH3 & SKK & MCA),
    n_missing_at_least_one = sum(!(CHN & MH1 & MH3 & SKK & MCA)),
    n_only_CHN = sum(CHN & !MH1 & !MH3 & !SKK & !MCA),
    n_only_MH1 = sum(!CHN & MH1 & !MH3 & !SKK & !MCA),
    n_only_MH3 = sum(!CHN & !MH1 & MH3 & !SKK & !MCA),
    n_only_SKK = sum(!CHN & !MH1 & !MH3 & SKK & !MCA),
    n_only_MCA = sum(!CHN & !MH1 & !MH3 & !SKK & MCA)
  )

upset_input_Lachno = gene_presence_Lachno %>%
  mutate(across(c(CHN, MH1, MH3, SKK, MCA), as.integer)) %>%
  select(CHN, MH1, MH3, SKK, MCA)

upset(upset_input_Lachno, sets = c("CHN", "MH1", "MH3", "SKK", "MCA"), order.by = "freq")



