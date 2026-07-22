#working_dir = "/fs/project/bradley.720/projects/meta_analysis_26"
working_dir = "/Users/tran.986/Desktop/meta_analysis_26"
source(paste0(working_dir, "/meta_analysis_26.R"))

#---pipeline 1: merge post-phylogenize
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
MH1_anno=left_join(MH1_all_res, annnotation_df, by = "gene")


#ERP002469_MH3:
MH3_all_res=readRes_n_pdCal(study_id = "MH3")
MH3_anno=left_join(MH3_all_res, annnotation_df, by = "gene")


#-------FOR TAXA LEVEL
#Make an upset plot (intersecting taxa identity)
#positive effect size:
anno_list <- list(
  CHN = CHN_anno,
  MH1 = MH1_anno,
  MH3 = MH3_anno
  # add more datasets as you have them
)

upset_plotMake(taxa_list=anno_list,
               direction_es = "negative")

upset_plotMake(taxa_list=anno_list,
               direction_es = "positive")

#Make a heatmap (to pair with the upset plot)
taxaLevel_plotMake(taxa_list=anno_list,
                   direction_es = "positive",
                   plot_type = "heatmap")

taxaLevel_plotMake(taxa_list=anno_list,
                   direction_es = "negative",
                   plot_type = "heatmap")


#-------FOR GENE LEVEL that consensus in each species
res_list=list("CHN"=CHN_all_res,
     "MH1"=MH1_all_res,
     "MH3"=MH3_all_res)

consensus_geneFind(all_res_ls=res_list,
                   direction_es = "positive") 

consensus_geneFind(all_res_ls=res_list,
                   direction_es = "negative")

#----------pipeline 2:
merge_res=read.csv(paste0(working_dir,"/all-results-merged-3ds.csv")) 

#---TAXON-GENE LEVEL:
Pre_mergePlot(mergeRes = merge_res)


#----------pipeline1 AND pipeline 2:
#how much of each study_id's all-result 
merge_res
MH1_all_res
CHN_all_res
MH3_all_res


#---Expand to more datasets:
library(dplyr)
data("sampleMetadata")   # explicitly loads it into your environment
sampleMetadata |> filter(study_condition == "T2D") |> View()





