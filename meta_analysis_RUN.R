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


#---pipeline 1: merge at gene-level (=)
#---CHN only: (Qin et al, 2012) - a metagenome wide association study...
#phylogenize_full_func(study_id = "CHN",
#                      metadata_dir_path = paste0(working_dir,"/metadata/CHN_md_final.csv"),
#                      ref_env = "T2D metformin-",
#                      envs_compared = c("ND CTRL", "T2D metformin-"))	

#---ERP002469_MH3 only: (karlssonFH et al, 2013)
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

#---ERP004605_MH1 only: "(Li et al, 2014) An integrated catalog of reference genes in the human gut microbiome" 
#phylogenize_full_func(study_id = "ERP004605_MH1",
#		      metadata_dir_path = paste0(working_dir,"/metadata/mhn1_md_final.csv"),
#		      ref_env = "T2D metformin-",
#		      envs_compared = c("ND CTRL", "T2D metformin-"))

#-----pipeline2: merge pre-phylogenize:
#--step 1: run ancombc2
#import_bracken_CHN=read_tsv(paste0(working_dir, "/phylogenize_out/CHN/data_w_count.tsv"))
#import_bracken_MH1=read_tsv(paste0(working_dir, "/phylogenize_out/ERP004605_MH1/data_w_count.tsv"))
#import_bracken_MH3=read_tsv(paste0(working_dir, "/phylogenize_out/ERP002469_MH3/data_w_count.tsv"))

#count_tbl_ls = list(#CHN = import_bracken_CHN,
#     ERP004605_MH1 = import_bracken_MH1,
#     ERP002469_MH3 = import_bracken_MH3)

#metadata_tbl_ls = list(#CHN = read_csv(paste0(working_dir,"/metadata/CHN_md_final.csv")),
#ERP004605_MH1 =  read_csv(paste0(working_dir,"/metadata/mhn1_md_final.csv")),
#		      ERP002469_MH3 = read_csv(paste0(working_dir,"/metadata/mhn3_md_final.csv")))

study_id_ls = c("CHN", "ERP004605_MH1", 
                "ERP002469_MH3")


#run through all of the items in the list:
#all_ancom_res <- lapply(seq_along(study_id_ls), function(i) {
#	cat("\nRunning", study_id_ls[i], "\n")
#       ancomRun(count_tbl = count_tbl_ls[[i]],
#	study_id = study_id_ls[i],
#	metadata = metadata_tbl_ls[[i]])
#})

#names(all_ancom_res) <- study_id_ls

#saveRDS(all_ancom_res,
#	paste0(working_dir, "/pipeline2/ancom/all_ancom_res.rds"))

#---step 2: combine ancombc w fixed effect?
#ancom_out_ls = lapply(study_id_ls, function(id) 
#       readRDS(paste0(working_dir, "/pipeline2/ancom/", id, "_ancombc_res.rds"))
#       )
#--calculate weight per taxon for each study_id:
#weight_per_taxon=weight_per_taxonCal(ancomRun_output_ls = ancom_out_ls,
#		    study_id_ls = study_id_ls)

#--calculate the sum of weights per taxon across studies:
#sigma_weight_per_taxon=sigma_weightCal(weight_per_taxon_ls = weight_per_taxon)
#write.csv(sigma_weight_per_taxon, paste0(working_dir, "/pipeline2/ancom/combine/sigma_weight.csv"))

#--calculate the SE pooled:
#se_pooled=SE_poolCal(sigma_weightCal_out = sigma_weight_per_taxon)
#write.csv(se_pooled, paste0(working_dir, "/pipeline2/ancom/combine/se_pooled.csv"))

#--for each study_id, 
#--calculate the product of weight (per_taxon) * effect_size (of that taxon) 
#sum_prod_weight_mu=sum_prod_weightCal(weight_per_taxon_ls = weight_per_taxon,
#				      study_id_ls = study_id_ls)
#write.csv(sum_prod_weight_mu, paste0(working_dir, "/pipeline2/ancom/combine/sum_prod_weight_mu.csv"))

#--calculate for pooled effect size
#eff_size_pooled = mu_poolCal(sum_prod_weightCal_out = sum_prod_weight_mu, 
#			     sigma_weightCal_out = sigma_weight_per_taxon)
#write.csv(eff_size_pooled, paste0(working_dir, "/pipeline2/ancom/combine/eff_size_pooled.csv"))

#---step 3a: run ashr only once on the combined effect 
#ash_res_pooled=ashRun(mu_poolCal_out = eff_size_pooled,
#       SE_poolCal_out = se_pooled)
#write.csv(ash_res_pooled, paste0(working_dir, "/pipeline2/ash/ash_res_pooled.csv"))

#---step 3b: run phylogenize repermulize()
#preping for provided phenotype file:
#se_poolCal_out = read.csv(paste0(working_dir,"/pipeline2/ancom/combine/se_pooled.csv"))
#mu_poolCal_out = read.csv(paste0(working_dir, "/pipeline2/ancom/combine/eff_size_pooled.csv"))

#provided_file=inner_join(mu_poolCal_out, se_poolCal_out, 
#                         by = "taxon")  %>%
#  dplyr::filter(n_studies == 3) %>%
#  dplyr::select(taxon, se_pooled, mu_pooled) %>% 
#  dplyr::rename("estimate"="mu_pooled",
#                "stderr"="se_pooled")

#write.table(provided_file, 
#            paste0(working_dir, "/pipeline2/phylogenize/input/phenotype_file.tsv"),
#            sep = "\t",
#            row.names = F, 
#            quote = F)

#phylogenize_run(provided_file = paste0(working_dir,
#                                       "/pipeline2/phylogenize/input/phenotype_file.tsv"),
#                study_id =NULL,
#                phenotype = "provided",
#                ref_env = c("T2D metformin-"))

#---------------------------------------------pipeline2: merge pre-phylogenize:
#--step 1: run ancombc2
#import_bracken_CHN=read_tsv(paste0(working_dir, "/pipeline1/phylogenize_out/CHN/data_w_count.tsv"))
#import_bracken_MH1=read_tsv(paste0(working_dir, "/pipeline1/phylogenize_out/ERP004605_MH1/data_w_count.tsv"))
#import_bracken_MH3=read_tsv(paste0(working_dir, "/pipeline1/phylogenize_out/ERP002469_MH3/data_w_count.tsv"))
#import_bracken_MCA=read_tsv(paste0(working_dir, "/pipeline1/phylogenize_out/MCA/data_w_count.tsv"))
#import_bracken_SKK = read_tsv(paste0(working_dir, "/pipeline1/phylogenize_out/SKK/data_w_count.tsv"))


#count_tbl_ls = list(MCA = import_bracken_MCA,
#		    SKK = import_bracken_SKK)
#CHN = import_bracken_CHN,
#     ERP004605_MH1 = import_bracken_MH1,
#     ERP002469_MH3 = import_bracken_MH3)

#metadata_tbl_ls = list(MCA = read_csv(paste0(working_dir, "/metadata/MCA_md_final.csv")),
#		       SKK = read_csv(paste0(working_dir, "/metadata/SKK_md_final.csv")))


#CHN = read_csv(paste0(working_dir,"/metadata/CHN_md_final.csv")),
#ERP004605_MH1 =  read_csv(paste0(working_dir,"/metadata/mhn1_md_final.csv")),
#ERP002469_MH3 = read_csv(paste0(working_dir,"/metadata/mhn3_md_final.csv")))

study_id_ls = c("MCA", "SKK", "CHN", "ERP004605_MH1","ERP002469_MH3")


#run through all of the items in the list:
#all_ancom_res <- lapply(seq_along(study_id_ls), function(i) {
#	cat("\nRunning", study_id_ls[i], "\n")
#        ancomRun(count_tbl = count_tbl_ls[[i]],
#	study_id = study_id_ls[i],
#	metadata = metadata_tbl_ls[[i]])
#})

#names(all_ancom_res) <- study_id_ls

#saveRDS(all_ancom_res,
#	paste0(working_dir, "/pipeline2/ancom/all_ancom_res_MCA_SKK.rds"))

#---step 2: combine ancombc w fixed effect?
#ancom_out_ls = lapply(study_id_ls, function(id) 
#       readRDS(paste0(working_dir, "/pipeline2/ancom/", id, "_ancombc_res.rds"))
#       )
#--calculate weight per taxon for each study_id:
#weight_per_taxon=weight_per_taxonCal(ancomRun_output_ls = ancom_out_ls,
#		    study_id_ls = study_id_ls)

#--calculate the sum of weights per taxon across studies:
#sigma_weight_per_taxon=sigma_weightCal(weight_per_taxon_ls = weight_per_taxon)
#write.csv(sigma_weight_per_taxon, paste0(working_dir, "/pipeline2/ancom/combine/sigma_weight_5studies.csv"))

#--calculate the SE pooled:
#se_pooled=SE_poolCal(sigma_weightCal_out = sigma_weight_per_taxon)
#write.csv(se_pooled, paste0(working_dir, "/pipeline2/ancom/combine/se_pooled_5studies.csv"))

#--for each study_id, 
#--calculate the product of weight (per_taxon) * effect_size (of that taxon) 
#sum_prod_weight_mu=sum_prod_weightCal(weight_per_taxon_ls = weight_per_taxon,
#				      study_id_ls = study_id_ls)
#write.csv(sum_prod_weight_mu, paste0(working_dir, "/pipeline2/ancom/combine/sum_prod_weight_mu_5studies.csv"))

#--calculate for pooled effect size
#eff_size_pooled = mu_poolCal(sum_prod_weightCal_out = sum_prod_weight_mu, 
#			     sigma_weightCal_out = sigma_weight_per_taxon)
#write.csv(eff_size_pooled, paste0(working_dir, "/pipeline2/ancom/combine/eff_size_pooled_5studies.csv"))

#---step 3: run ashr only once on the combined effect
#ash_res_pooled=ashRun(mu_poolCal_out = eff_size_pooled,
#       SE_poolCal_out = se_pooled)
#write.csv(ash_res_pooled, paste0(working_dir, "/pipeline2/ash/ash_res_pooled_5studies.csv"))

#---step 4: run phylogenize repermulize()
#preping for provided phenotype file:
#se_poolCal_out = read.csv(paste0(working_dir,"/pipeline2/ancom/combine/se_pooled_5studies.csv"))
#mu_poolCal_out = read.csv(paste0(working_dir, "/pipeline2/ancom/combine/eff_size_pooled_5studies.csv"))

#provided_file=inner_join(mu_poolCal_out, se_poolCal_out, 
#                         by = "taxon")  %>%
#  dplyr::filter(n_studies == 5) %>% #fix "n_studies" if more datasets are added... 
#  dplyr::select(taxon, se_pooled, mu_pooled) %>% 
#  dplyr::rename("estimate"="mu_pooled",
#                "stderr"="se_pooled")

#write.table(provided_file, 
#            paste0(working_dir, "/pipeline2/phylogenize/input/phenotype_file_5studies.tsv"),
#            sep = "\t",
#            row.names = F, 
#            quote = F)

#run phylogenize "provided" phenotype:
#phylogenize_run(provided_file_path = paste0(working_dir,"/pipeline2/phylogenize/input/phenotype_file_5studies.tsv"),
#                study_id = NULL,
#                phenotype = "provided",
#                ref_env = "T2D metformin-")

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

#----------pipeline 1 (run each dataset separately) results import:
#calculate for Faith's Diversity:
#--CHN:
CHN_all_res=readRes_n_pdCal(study_id = "CHN")
CHN_anno=left_join(CHN_all_res, annnotation_df, by = "gene")


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

#======================================= Pipeline 1A
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

#------- step 4: multiple hypothesis correction:
pooled_res$q_value_satt <-  p.adjust(pooled_res$p_value_satt, method = "BH")
pooled_res |> View()

#------Make a figure for those genes in those taxa that are from >1 studies:
genes_to_plot = c("GUT_GENOME029728_01359","GUT_GENOME076986_01167", "UniRef50_A0A060CPZ1", "UniRef50_A0A417LEQ2", "UniRef50_UPI00102F5044")
overlap_genes_study = combined_variance |> filter(gene %in% genes_to_plot) |>
  mutate(se_val = sqrt(variance),
         ci_lower = effect.size - 1.96 * se_val, #ci have to be computed from original effect.size for each gene/study
         ci_upper = effect.size + 1.96 * se_val)

overlap_genes_pooled = pooled_res |>
  filter(gene %in% genes_to_plot) |>
  mutate(
    study_name = "Pooled",
    se_val = se_pooled,
    ci_lower = effsize_pooled - 1.96 * se_pooled,
    ci_upper = effsize_pooled + 1.96 * se_pooled,
    effect.size = effsize_pooled
  ) %>%
  select(study_name, taxon, gene, effect.size, se_val, ci_lower, ci_upper)

forest_data_full =  bind_rows(overlap_genes_study[, c("study_name", "taxon", "gene", "effect.size", "se_val", "ci_lower", "ci_upper")],
                              overlap_genes_pooled) 
forest_data_full$study_name<-factor(forest_data_full$study_name, levels = c("Pooled", "CHN", "MH1", "MH3", "MCA", "SKK"))

test_pipeline1A=ggplot(forest_data_full, aes(x = effect.size, y = study_name, color = study_name == "Pooled")) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ gene + taxon, scales = "free_x", ncol = 1) +   # <-- ncol = 1 forces single column
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), guide = "none") +
  labs(x = "Effect size (95% CI)", y = NULL, title = "Effect size comparison across studies") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold",
                                  size = 14, 
                                  hjust = 0.5))

ggsave(plot = test_pipeline1A,
       filename= "~/Desktop/test_1A.png",
       width = 5,
       height = 14,
       units = "in",
       dpi = 600)

#======================================= Pipeline 2 FIX USING RANDOM EFFECT INSTEAD - locally run
study_ids = c("CHN", "ERP004605_MH1", "ERP002469_MH3", "MCA", "SKK")

#read into list of ancombc results + ashr 
ancombc_list = setNames(lapply(study_ids, function(id) {

  res = readRDS(paste0(working_dir, "/ancom/", id, "_ancombc_res.rds"))
  #run ancombc
  extracted = ancomExtract(ancomRun_output = res)
  #run ashr:
  ash_fit = ash(betahat = extracted$mu, sebetahat = extracted$se)
  
  extracted$lfc_shrunk = ash_fit$result$PosteriorMean
  extracted$se_shrunk  = ash_fit$result$PosteriorSD
  
  extracted$study<-id
  extracted %>% select(taxon, lfc_shrunk, se_shrunk, study)
  
}), study_ids)

combined_ancombc = bind_rows(ancombc_list)

#Run random-effects meta-analysis per taxon
taxon_meta_results = combined_ancombc %>%
  group_by(taxon) %>%
  group_modify(~ {
    # need at least 2 studies to run meta-analysis
    if (nrow(.x) < 2) {
      return(tibble(
        theta_pooled = NA, se_pooled = NA, tau2 = NA,
        I2 = NA, p_value = NA, n_studies = nrow(.x)
      ))
    }
    
    fit = tryCatch(
      rma(yi = .x$lfc_shrunk, sei = .x$se_shrunk, method = "REML"),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      return(tibble(
        theta_pooled = NA, se_pooled = NA, tau2 = NA,
        I2 = NA, p_value = NA, n_studies = nrow(.x)
      ))
    }
    
    tibble(
      theta_pooled = as.numeric(fit$b), #random-effects pooled log fold-change for that taxon
      se_pooled    = fit$se, #SE of the pooled estimate (already reflects both within-study and between-study variance)
      tau2         = fit$tau2, #between-study variance
      I2           = fit$I2, #% of total variability due to heterogeneity rather than sampling error
      p_value      = fit$pval,
      n_studies    = nrow(.x)
    )
  }) %>%
  ungroup() 

write.csv(taxon_meta_results, paste0(working_dir, "/random_effect_out.csv"))


#======================================================================DRAFT=======================
