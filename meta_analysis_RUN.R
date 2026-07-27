#working_dir = "/fs/project/bradley.720/projects/meta_analysis_26"
working_dir = "/Users/tran.986/Desktop/meta_analysis_26"
source(paste0(working_dir, "/meta_analysis_26.R"))

#---pipeline 1: merge post-phylogenize
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
#count_feature shows how many study_id that gene/taxa shows up:
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
#function returns 2 heatmaps:
#1 is Taxa-level:
#how much of each study_id's all-result (post-Phylogenize merging) 
#taxa overlaps with all-results of pre-Phylogenize merging:

#2 is GENE LEVEL:
#Identify taxa AND genes that overlaps bw
# post-Phylogenize 
# and pre-Phylogenize merging-> compare their ES:

preVpost_HeatmapMake(direction_es="positive",
                     merge_resPre = merge_res,
                     res_listPost = res_list)


preVpost_HeatmapMake(direction_es = "negative",
                     merge_resPre = merge_res,
                     res_listPost = res_list)


#-----------------------================================-----------Expand to more datasets:
#             HMP_2019_t2d --> HMP (not started)
#         KarlssonFH_2013 --> MH3 (done)
#                LiJ_2014 --> MH1 (done)
#      MetaCardis_2020_a --> MCA (started - in progress) 
#              QinJ_2012 --> CHN (done)
# SankaranarayananK_2015 --> SKK (not started)

#---MetaCardis -- Molinaro et al: Imidazole propionate is increased in diabetes and associated with dietary patterns and altered microbial ecology
metadata_MCA_ctrl = metadataRetrieve(study_id = "MCA")[["healthy"]]
metadata_MCA_t2d = metadataRetrieve(study_id = "MCA")[["T2D"]]

#---SankaranarayananK_2015: Gut Microbiome Diversity among Cheyenne and Arapaho Individuals from Western Oklahoma
metadata_SKK_ctrl = metadataRetrieve(study_id = "SKK")[["healthy"]]
metadata_SKK_t2d = metadataRetrieve(study_id = "SKK")[["T2D"]]

#run seq_Retrieve to obtain the URL --> URL will download an API -->
#API will return fastq_ftp on HPC:

#---MetaCardis:-this takes extreme long- DONOT RERUN - just read it locally"
ctrl_MCA_seqid=ftpRetrieve(metadata_MCA_ctrl)
t2d_MCA_seqid=ftpRetrieve(metadata_MCA_t2d)
        
MCA_seqid = rbind(ctrl_MCA_seqid,
                  t2d_MCA_seqid)

#write.table(MCA_seqid, file = paste0(working_dir,"/fastq_url/MCA_fastq_url.txt"), sep = "\t", row.names = FALSE)


#---SankaranarayananK_2015:
ctrl_SKK_seqid=ftpRetrieve(metadata_SKK_ctrl)
t2d_SKK_seqid=ftpRetrieve(metadata_SKK_t2d)

#SKK_seqid = rbind(ctrl_SKK_seqid,
#                  t2d_SKK_seqid)

#write.table(SKK_seqid, file = paste0(working_dir,"/fastq_url/SKK_fastq_url.txt"), sep = "\t", row.names = FALSE)







 


