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

#---read into phylogenize result:
RandMerge_phyloz_out = read.csv(paste0(working_dir, "/all-results-randTaxaMerge.csv"))

RandMerge_results = RandMerge_phyloz_out %>%
  mutate(
    neglog10q = -log10(q.value),
    sig_status = case_when(
      q.value < 0.05 & effect.size > 0 ~ "Up, significant",
      q.value < 0.05 & effect.size < 0 ~ "Down, significant",
      TRUE ~ "Not significant"),
    plot_color = if_else(neglog10q >= -log10(0.05), as.character(taxon), "Not significant")
  )

taxon_levels = sort(unique(RandMerge_results$plot_color[RandMerge_results$plot_color != "Not significant"]))
taxon_colors = setNames(hue_pal()(length(taxon_levels)), taxon_levels)
color_values = c(taxon_colors, "Not significant" = "grey80")

RandMerge_results$plot_color = factor(
  RandMerge_results$plot_color,
  levels = c("Not significant", taxon_levels)
)

ggplot(RandMerge_results %>% arrange(plot_color), aes(x = effect.size, y = neglog10q, color = plot_color)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = color_values, name = "Taxon") +
  labs(x = "Effect size", y = expression(-log[10](q)), title = "Merge at GENE: 377 genes significant ") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

#testing/Drafting
RandMerge_results[RandMerge_results$sig_status == "Up, significant",] |> arrange(-effect.size) |> left_join(annotation_df, 
                                                                                                            by = "gene") |> View()

#================================Validating Corio and Lachno hits by different means:
#-----------------------option 1:Applies Aldex3 on the same data:
library(ALDEx3)
#import X (count tbl) and Y (metadata):
count_tbl_merge = read.delim(paste0(working_dir, "/count_tbl_merged.tab")) 
metadata_tbl_merge = read.delim(paste0(working_dir, "/metadata_merged.tab")) |> column_to_rownames(var = "sample")

#filter count tbl at 75%
count_tbl_merge = count_tbl_merge %>%
  dplyr::select(name, rownames(metadata_tbl_merge)) %>%  #where count_tbl = merged count_tbl
  column_to_rownames(var = "name")

keep_names=row.names(count_tbl_merge[((rowSums(count_tbl_merge==0))/ncol(count_tbl_merge))<=0.75,])
other_names=colSums(count_tbl_merge[((rowSums(count_tbl_merge==0))/ncol(count_tbl_merge))>0.75,])
count_tbl_merge <- count_tbl_merge[keep_names,]
count_tbl_merge_aldex <- rbind(count_tbl_merge, other_names) 

#clean up metadata:
metadata_tbl_merge$env <- factor(metadata_tbl_merge$env,
                                 levels=c("ND CTRL", "T2D metformin-"))

#fit into Aldex3 allowing uncertainty around the CLR-implied scale differences:
aldex_fit <- aldex(count_tbl_merge_aldex,
                   ~env,
                   metadata_tbl_merge,
                   nsample=1554,
                   scale=clr.sm,  # CLR assumption
                   gamma=1) 


#=======================================draft:
#count
data(gut_crohns_data)
Y <- gut_crohns_data$counts
keep_names <- row.names(Y[((rowSums(Y==0))/ncol(Y))<=0.75,])
other <- colSums(Y[((rowSums(Y==0))/ncol(Y))>0.75,])
Y <- Y[keep_names,]
Y <- rbind(Y, other)

#metadata:
X <- gut_crohns_data$metadata
X$Health.status <- factor(X$Health.status,
                          levels=c("Control", "CD"))

ncol(Y)
# Allow uncertainty around the CLR-implied scale differences
aldex.gut.raw <- aldex(Y,
                       ~Health.status,
                       X,
                       nsample=2000,
                       scale=clr.sm,  # CLR assumption
                       gamma=1)  

aldex.gut.summary <- summary(aldex.gut.raw)
