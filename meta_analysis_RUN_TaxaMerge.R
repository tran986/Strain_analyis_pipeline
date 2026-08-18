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

#============================================================#======================================= Pipeline 1 merge at taxa level:
#--step 1: run ancombc2
#import_bracken_CHN=read_tsv(paste0(working_dir, "/phylogenize_out/CHN/data_w_count.tsv"))
#import_bracken_MH1=read_tsv(paste0(working_dir, "/phylogenize_out/ERP004605_MH1/data_w_count.tsv"))
#import_bracken_MH3=read_tsv(paste0(working_dir, "/phylogenize_out/ERP002469_MH3/data_w_count.tsv"))
#import_bracken_MCA=read_tsv(paste0(working_dir, "/pipeline1/phylogenize_out/MCA/data_w_count.tsv"))
#import_bracken_SKK = read_tsv(paste0(working_dir, "/pipeline1/phylogenize_out/SKK/data_w_count.tsv"))

#count_tbl_ls = list(MCA = import_bracken_MCA,
#		    SKK = import_bracken_SKK,
#		    CHN = import_bracken_CHN,
#	       	    ERP004605_MH1 = import_bracken_MH1,
#	       	    ERP002469_MH3 = import_bracken_MH3)

#metadata_tbl_ls = list(MCA = read_csv(paste0(working_dir, "/metadata/MCA_md_final.csv")),
#		       SKK = read_csv(paste0(working_dir, "/metadata/SKK_md_final.csv")),
#		       CHN = read_csv(paste0(working_dir,"/metadata/CHN_md_final.csv")),
#		       ERP004605_MH1 =  read_csv(paste0(working_dir,"/metadata/mhn1_md_final.csv")),
#		       ERP002469_MH3 = read_csv(paste0(working_dir,"/metadata/mhn3_md_final.csv")))

study_id_ls = c("CHN", "ERP004605_MH1", 
                "ERP002469_MH3")


#-------------------run through all of the items in the list:
#all_ancom_res <- lapply(seq_along(study_id_ls), function(i) {
#	cat("\nRunning", study_id_ls[i], "\n")
#       ancomRun(count_tbl = count_tbl_ls[[i]],
#	study_id = study_id_ls[i],
#	metadata = metadata_tbl_ls[[i]])
#})

#names(all_ancom_res) <- study_id_ls

#saveRDS(all_ancom_res,
#	paste0(working_dir, "/pipeline2/ancom/all_ancom_res.rds"))

#---step 2:--------------------combine ancombc w fixed effect?
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

#use ashr result as a provided_file instead:
#ash_res_pooled = read.csv(paste0(working_dir, "/pipeline2/ash/ash_res_pooled_5studies.csv"))

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













#create a provided_file where n_studies >= 4
provided_file = taxon_meta_results %>%
  dplyr::filter(n_studies == 5) |>
  dplyr::select(taxon, se_pooled, theta_pooled) |>
  dplyr::rename("estimate"="theta_pooled",
                "stderr"="se_pooled") %>% 
  drop_na(stderr, estimate)

write.table(provided_file, 
            paste0(working_dir, "/provided_file_REML_w_ash.tsv"),
            sep = "\t",
            row.names = F, 
            quote = F)

#phylogenize is run on HPC:...

#read into phylogenize results of merged 5 studies, w ash, REML:
all_res_reml5= readRes_n_pdCal(study_id = "ash-5studies-REML") |> #506 hits in total
  mutate(neglog10q = -log10(q.value),
         sig_status = case_when(q.value < 0.05 & effect.size > 0 ~ "Up, significant",
                                q.value < 0.05 & effect.size < 0 ~ "Down, significant"),
         plot_color = ifelse(neglog10q >= -log10(0.05), 
                             as.character(taxon),
                             "not significant"))

pipeline2_plot=ggplot(all_res_reml5, aes(x = effect.size, y = neglog10q, color = taxon)) +
  geom_point(alpha = 0.65, size = 1.9) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(x = "Effect size", y = expression(-log[10](q)), title = "Pipeline 2: gene associations by taxon across 5 studies (n=506)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave(plot = pipeline2_plot,
       filename = paste0(working_dir, "/test_pipeline2_plot.png"),
       width = 9,
       height = 5,
       units = "in",
       dpi = 600)

study_ids = c("CHN", "ERP004605_MH1", "ERP002469_MH3", "MCA", "SKK")


