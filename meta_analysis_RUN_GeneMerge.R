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

