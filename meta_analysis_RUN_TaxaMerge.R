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
#source("/fs/project/bradley.720/projects/meta_analysis_26/meta_analysis_26.R")
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
#---merge all count table from all studies:
#import_ls = c("CHN", "MH1", "MH3", "MCA_fixed_effect") #, "SKK")

#read into bracken + save into fix effect dir
#print("reading into bracken dirs...")
#import_bracken(study_id = "CHN",
#	       count_tbl_outdir = paste0(working_dir, "/fix_effect/count_tbl_fix_eff/CHN"))

#lapply(import_ls, function(id) {
#   import_bracken(study_id = id,
#                  count_tbl_outdir = paste0(working_dir, "/fix_effect/count_tbl_fix_eff/", id))
#})

#print("put them into list...")
#import_bracken_ls = setNames(lapply(import_ls, function(i) {
#  read_tsv(paste0(working_dir, "/fix_effect/count_tbl_fix_eff/", i, "/data_w_count.tsv"))
#}), import_ls)

#print("merging count tbl ...")
#count_tbl_merge = Reduce(function(x,y) merge(x, y, by = "name", all = T), import_bracken_ls)

#count_tbl_merge[-1]<- lapply(count_tbl_merge[-1], function(x){
#x_num <- suppressWarnings(as.numeric(as.character(x)))
#    x_num[is.na(x_num) & !is.na(x)] <- 0
#    x_num
#})

#count_tbl_merge[is.na(count_tbl_merge)] <- 0
#print("any non-numeric value:")
#anyNA(count_tbl_merge[-1])

#write.table(count_tbl_merge, 
#	  file = paste0(working_dir, "/fix_effect/count_tbl_fix_eff/count_tbl_merged.tab"),
#	  sep = "\t",
#	  row.names = F,
#	  quote = F)

count_tbl_merge = read.delim(
  paste0(working_dir, "/fix_effect/count_tbl_fix_eff/count_tbl_merged.tab"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(count_tbl_merge) <- NULL

#--merge all metadata
#metadata_ls=setNames(lapply(import_ls, function(i) {
# df = read.csv(paste0(working_dir, "/fix_effect/metadata_fix_eff/", i, "/metadata_filter.csv")) 
#}), import_ls)

print("merging metadata")
#metadata_tbl_merge = bind_rows(metadata_ls) |> dplyr::select(sample, disease, treatment, dataset)
#write.table(metadata_tbl_merge, 
#	    file = paste0(working_dir, "/fix_effect/metadata_fix_eff/metadata_merged.tab"),
#	    sep = "\t", row.names = FALSE)

metadata_tbl_merge <- read.delim(
  paste0(working_dir, "/fix_effect/metadata_fix_eff/metadata_merged.tab"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(metadata_tbl_merge) <- NULL

#--run ANCOMBC:
#print("fitting ANCOMBC...")
#ancomMerge = ancomRun(count_tbl = count_tbl_merge,
#	 study_id = NULL,
#	 metadata = metadata_tbl_merge,
#	 fixed_effect_used = T)

#ancomMerge = readRDS(paste0(working_dir, "/fix_effect/ancom/merged_ancom_res.rds"))

#--run ashr:
#ancomMerge_res = ancomMerge$res

#print("fitting ashr for T2D vs non-T2D")
#ashMerge_t2d = ashr::ash(betahat = ancomMerge_res$lfc_diseaseT2D, sebetahat = ancomMerge_res$se_diseaseT2D) 

#print("fitting ashr for metformin vs non-metformin")
#ashMerge_met = ashr::ash(betahat = ancomMerge_res$lfc_treatmentyes, sebetahat = ancomMerge_res$se_treatmentyes)

#print("rename ash res to become a 'provided' file for phylogenize")
#ashMerge_t2d_res = ashMerge_t2d$res
#ashMerge_met_res = ashMerge_met$res

#t2d:
#ashMerge_t2d_res$taxon <- ancomMerge_res$taxon
#ashMerge_t2d_res = ashMerge_t2d_res %>%
#	dplyr::select(taxon, PosteriorMean, PosteriorSD) %>%
#	dplyr::rename("estimate"="PosteriorMean",
#		      "stderr"="PosteriorSD")

#write.table(ashMerge_t2d_res,
#	    paste0(working_dir, "/fix_effect/ashr/ashMerge_t2d_res.tab"),            
#	    sep = "\t", row.names = FALSE)   #this is "provided" file 

#met:
#ashMerge_met_res$taxon<- ancomMerge_res$taxon
#ashMerge_met_res = ashMerge_met_res %>%
#	dplyr::select(taxon, PosteriorMean, PosteriorSD) %>%
#	dplyr::rename("estimate"="PosteriorMean",
#		      "stderr"="PosteriorSD")

#write.table(ashMerge_met_res,
#	    paste0(working_dir, "/fix_effect/ashr/ashMerge_met_res.tab"),
#	    sep = "\t", row.names = FALSE)   #this is "provided" file 

#--run phylogenize2:
#for t2d:
#phylogenize_run(provided_file_path = paste0(working_dir, "/fix_effect/ashr/ashMerge_t2d_res.tab"),
#		study_id = NULL,
#		phenotype = "provided",
#		ref_env = "T2D",
#	        env_col = "disease",	
#	       	out_dir = paste0(working_dir, "/fix_effect/phylogenize/t2d"))

#for met:
#phylogenize_run(provided_file_path = paste0(working_dir,"/fix_effect/ashr/ashMerge_met_res.tab"),
#		study_id = NULL,
#		phenotype = "provided",
#		ref_env = "yes",
#		env_col = "treatment",
#		out_dir = paste0(working_dir, "/fix_effect/phylogenize/met"))


#-- look at the top 5 most used drugs: -- for MCA: (SKK study does not include that information in their study)
treatment_counts <- sampleMetadata[sampleMetadata$study_name == "MetaCardis_2020_a", ] |>
  select(treatment) |>
  separate_longer_delim(treatment, delim = ";") |>
  filter(!is.na(treatment), treatment != "") |>
  count(treatment, name = "n_people", sort = TRUE) |>
  arrange(-n_people)

View(sampleMetadata[sampleMetadata$study_name == "MetaCardis_2020_a", ] )

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

#=========================After Phylogenize:
#---read into phylogenize result:
RandMerge_phyloz_out = read.csv(paste0(working_dir, "/all-results-randTaxaMerge.csv"))
fix_eff_phyloz_out_t2d = read.csv(paste0(working_dir, "/all-results-t2d-fix-eff.csv"))
fix_eff_phyloz_out_met = read.csv(paste0(working_dir, "/all-results-met-fix-eff.csv"))

#run function 24 to make a figure:
figMake_phyloz(phyloz_out = fix_eff_phyloz_out_t2d)
figMake_phyloz(phyloz_out = fix_eff_phyloz_out_met)

#add PD calculation:
core_outRandMerge = readRDS(paste0(working_dir, "/core_output/core_output_randTaxaMerge.rds"))
core_out_fix_eff_t2d = readRDS(paste0(working_dir, "/core_output/core_output-fix-eff-t2d.rds"))
core_out_fix_eff_met = readRDS(paste0(working_dir, "/core_output/core_output-fix-eff-met.rds"))

pdAdd(core_out=core_outRandMerge,
      phyloz_out = RandMerge_phyloz_out)

pdAdd(core_out=core_out_fix_eff_t2d,
      phyloz_out = fix_eff_phyloz_out_t2d)

pdAdd(core_out=core_out_fix_eff_met,
      phyloz_out = fix_eff_phyloz_out_met)

#================================Validating Corio and Lachno hits by different means:
#-----------------------option 1:Applies Aldex3 on the same data:
library(ALDEx3)
#import X (count tbl) and Y (metadata):
#count_tbl_merge = read.delim(paste0(working_dir, "/mergeTaxa/ancom/combine_fix/count_tbl_merged.tab")) 
#metadata_tbl_merge = read.delim(paste0(working_dir, "/mergeTaxa/ancom/combine_fix/metadata_merged.tab")) |> column_to_rownames(var = "sample")

#if fix effect of treatment is used:
metadata_tbl_merge = read.delim(paste0(working_dir, "/fix_effect/metadata_fix_eff/metadata_merged.tab"))
count_tbl_merge = read.delim(paste0(working_dir, "/fix_effect/count_tbl_fix_eff/count_tbl_merged.tab"))


#print("fitting Aldex3 with only Disease only...")
#saveRDS(aldex_fit_diseaseOnly, paste0(working_dir, "/fix_effect/aldex/aldex3_fit_diseaseOnly.rds"))

print("fitting Aldex3 with disease & treatment as fixed effects ...")
#aldex_fit_fix_effs = aldexRun(count_tbl = count_tbl_merge,
#			      metadata_tbl = metadata_tbl_merge,
#			      fix_effect_used = TRUE)

#saveRDS(aldex_fit_fix_effs, paste0(working_dir, "/fix_effect/aldex/aldex3_fit_fix_effs.rds"))

aldex_fit_fix_effs = readRDS(paste0(working_dir, "/aldex3_fit_fix_effs.rds"))
#use this aldex3 fit for ashr -> phylogenize2 
#print("saving aldex3 fit res for Phylogenize...")
#aldex_fit_diseaseOnly_res = summary(aldex_fit_diseaseOnly)
print("processing aldex3 output-mean")
aldex_fix_eff_res = aldexExtract(aldex_output = aldex_fit_fix_effs)

#combining mean and se into 1 df -- diseaseT2D:
diseaseT2D = dplyr::inner_join(data.frame(taxon = aldex_fix_eff_res$estimate$taxon, 
                                          estimate = aldex_fix_eff_res$estimate$diseaseT2D),
                               data.frame(taxon = aldex_fix_eff_res$std.error$taxon,
                                          std.error = aldex_fix_eff_res$std.error$diseaseT2D),
                               by = "taxon"
)


print("applied ashr on aldex3 output...") #saved into "cross_check" dir under mergeTaxa
aldex_ash_fix_effT2D = ashr::ash(betahat = diseaseT2D$estimate,
                                sebetahat = diseaseT2D$std.error)

aldex_ash_fix_eff_T2D_res = aldex_ash_fix_effT2D$res
aldex_ash_fix_eff_T2D_res$taxon <- diseaseT2D$taxon

aldex_ash_fix_eff_res = aldex_ash_fix_eff_T2D_res$taxon %>%
  dplyr::select(taxon, PosteriorMean, PosteriorSD) %>%
  dplyr::rename("estimate"="PosteriorMean",
                "stderr"="PosteriorSD")

#save "provided" phenotype file:
write.table(aldex_ash_fix_eff_res,
            paste0(working_dir, "/fix_effect/aldex/aldex3_fix_eff_ash_res.tab"),            
            sep = "\t", row.names = FALSE)

#phylogenize_run(provided_file_path = paste0(working_dir, "/fix_effect/aldex/aldex_ash_res.tab"),
#		study_id = NULL,
#		phenotype = "provided",
#		ref_env = "ND CTRL",
#		out_dir = paste0(working_dir, "/phylogenze_aldex"),
#		output_file = paste0(working_dir, "/mergeTaxa/cross_check/phylogenize_aldex/phylogenize.html"))	    

#-----------------
#add PD to Aldex3 run as well:
aldex_core_output = readRDS(paste0(working_dir, "/core_output_aldex3_taxaMerge_crosscheck.rds"))
aldex_phyloz_PD=gene_dist_func(phyloz_output = aldex_phyloz_out,
               core_out = aldex_core_output)
aldex_phyloz_PD_func = left_join(aldex_phyloz_PD, 
                                 aldex_core_output$list_pheno$pz.db$gene.to.fxn,
                                 by = "gene") |> arrange(-PD) |> View()

taxa_overlap_aldex_ancom = intersect(aldex_phyloz_out[aldex_phyloz_out$q.value < 0.05, ]$taxon,
                                   RandMerge_phyloz_out[RandMerge_phyloz_out$q.value < 0.05, ]$taxon)

#---------Aldex3 when run with fix effect:
aldex_core_fix_eff_t2d = readRDS(paste0(working_dir, "/core_output/core_output-aldex-fix-eff-t2d.rds"))
aldex_core_fix_eff_met = readRDS(paste0(working_dir, "/core_output/core_output-aldex-fix-eff-met.rds"))

aldex_fix_eff_phyloz_t2d = read.csv(paste0(working_dir, "/all-results-aldex-fix-eff-t2d.csv"))
aldex_fix_eff_phyloz_met = read.csv(paste0(working_dir, "/all-results-aldex-fix-eff-met.csv"))

#create fig for phylogenize output:
figMake_phyloz(aldex_fix_eff_phyloz_t2d)
figMake_phyloz(aldex_fix_eff_phyloz_met)

#add PD:
pdAdd(core_out = aldex_core_fix_eff_t2d,
      phyloz_out = aldex_fix_eff_phyloz_t2d)

pdAdd(core_out = aldex_core_fix_eff_met,
      phyloz_out = aldex_fix_eff_phyloz_met)

#-----------------------option 2:Making sure what we see is due to technical (length) of the seq:
study_id_ls = c("CHN", "ERP002469_MH3", "ERP004605_MH1", "MCA", "SKK")
import_bracken_truncate_ls = setNames(lapply(study_id_ls, function(x) 
  {import_bracken_truncate(study_id = x)}), study_id_ls)

taxa_counts <- imap_dfr(
  import_bracken_truncate_ls,
  function(study_data, study_name) {
    
    map_dfr(
      names(study_data),
      function(len) {
        
        df <- study_data[[len]]
        
        tibble(
          study = study_name,
          length = as.numeric(len),
          sample = names(df)[-1],
          n_taxa = colSums(df[-1] > 0))
      }
    )
  }
)

#make a plot:
technical_check_plot = ggplot(
  taxa_counts,
  aes(x = length, y = n_taxa, group = sample)
) +
  geom_line(alpha = 0.2) +
  geom_point(alpha = 0.4) +
  geom_smooth(
    aes(group = study),
    method = "lm",
    se = TRUE
  ) +
  facet_wrap(~ study) +
  labs(
    x = "Sequence length (bp)",
    y = "Number of detected taxa"
  ) +
  theme_classic()

#apply correlation len vs count of taxa:
cross_check = setNames(lapply(study_id_ls, function(id) {
  cross_check_func(study_id = id,
                   import_bracken_truncate_ls = import_bracken_truncate_ls)
}), study_id_ls)


#CHN:-- all good
cross_check$CHN$lm_res
hist(cross_check$CHN$glm_res)

#MH3:-- all good
cross_check$ERP002469_MH3$lm_res
hist(cross_check$ERP002469_MH3$glm_res)

#MH1: --all good
cross_check$ERP004605_MH1$lm_res
hist(cross_check$ERP004605_MH1$glm_res)

#MCA:-weird lmao
cross_check$MCA$lm_res
hist(cross_check$MCA$glm_res)

#SKK:
cross_check$SKK$lm_res
hist(cross_check$SKK$glm_res) #this is ok

#-----------#-----------#-----------------------option 3:


