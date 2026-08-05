library(TreeSummarizedExperiment)
library(dplyr)
library(tidyverse)
library(readxl)
library(httr2)
library(readr)
library(purrr)
library(dplyr)
library(tidyverse)
library(devtools)
library(ape)
library(picante)
library(UpSetR)
library(circlize)
library(curatedMetagenomicData)
#devtools::load_all("/fs/project/bradley.720/projects/phylogenize_v2/phylogenize_repermulize/package/repermulize")
devtools::load_all("/Users/tran.986/Desktop/phylogenize/package/phylogenize")


#--------
#Default Functions:
#1. function to get url to run on HPC:
#e.g., study_id = CHN, ERP004605, ERP003612, ERP002061, ERP002469
wget_urls_func <- function(md_final, study_id) {
  
  out_file <- paste0(working_dir, "/fastq_url/", study_id, "_fastq_urls.txt")
  
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

#function that run full function 2, 3 and 4:
phylogenize_full_func=function(study_id, metadata_dir_path, envs_compared, ref_env) {
  #import bracken:
  bracken_out = import_bracken(study_id = study_id)
  
  #extract phyloz_metadata:
  metadata = read_csv(metadata_dir_path)
  ext_md = extract_phyloz_metadata(import_bracken_out = bracken_out,
                                   metadata = metadata,
                                   study_id = study_id,
                                   envs_compared = envs_compared)
  
  #run phylogenize to get .rds:
  phylogenize_run(#extract_phyloz_metadata_out = ext_md,
    #import_bracken_out = bracken_out,
    study_id = study_id,
    ref_env = ref_env,
    phenotype = "abundance")
  
}


#2. function to read into .bracken files (push this to ASC):
import_bracken=function(study_id) { 
  file = list.files(path=paste0(working_dir, "/abund/", study_id))
  data = map(file, ~ { read_tsv(file.path(paste0(working_dir, "/abund/", study_id), .), show_col_types=FALSE) })
  samples = gsub("\\.bracken","", file)
  data_newcol = map2(data, samples, ~ {mutate(.x, sample = .y)})
  data_tidy = bind_rows(data_newcol)
  data_w_count = tidyr::pivot_wider(data_tidy, 
                                    names_from = sample, 
                                    values_from = new_est_reads,
                                    id_cols = name,
                                    values_fill = 0)
  write.table(data_w_count, file = paste0(working_dir, "/phylogenize_out/", study_id, "/data_w_count.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  return(data_w_count)
}

#3. function to create phylogenize input:
extract_phyloz_metadata=function(import_bracken_out, 
                                 metadata,
                                 study_id, #study_id = CHN, MHN1, MHN3, 4, SKK, MCA, etc..
                                 envs_compared = c("ND CTRL", "T2D metformin-")) { 
  #import_bracken_out = test
  #metadata = CHNs_md_final 
  #clean up metadata:
  if (study_id == "CHN") #then it is from CHN sample "SRR"
    
  { #CHN dataset:
    metadata_clean=metadata[,c("Run", "Status")] %>%
      dplyr::rename("env"="Status",
                    "sample"="Run") %>%
      dplyr::mutate(dataset = study_id)
    
  } else { #all MHN studies with "ERR" samples
    if (study_id %in% c("ERP004605_MH1","ERP002468_MH3")) {
      metadata_clean = metadata %>% 
        dplyr::mutate(sample = str_extract(fastq_ftp, "(?<=/)[^/]+(?=/[^/]+\\.fastq\\.gz)"),
                      dataset = study_id) %>%
        dplyr::select(Status, sample, dataset) %>%
        dplyr::rename("env"="Status")
    } else { #SKK and MCA studies from curatedmetagenomics:
      metadata_clean = metadata |>
        dplyr::mutate(dataset = study_id) 
    }
  }
  
  #filter 2 envs that will be compared:
  metadata_filter=metadata_clean %>% filter(env %in% envs_compared)
  
  write.table(metadata_filter, file = paste0(working_dir, "/pipeline1/phylogenize_out/", study_id, "/metadata_filter.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  return(metadata_filter)
  
}


#4. apply phylogenize on inputs: metadata and count
#dataset_dir = "/fs/project/bradley.720/projects/phylogenize_v2/phylogenize/package/phylogenize/inst/extdata"
phylogenize_run=function(provided_file_path = NULL,
                         #extract_phyloz_metadata_out,
                         #import_bracken_out,
                         study_id = NULL,
                         phenotype, #"abundance" or "provided"
                         ref_env) { #"ND CTRL or T2D metformin-
  if (phenotype == "abundance") {
    res=phylogenize(
      output_file = paste0(working_dir, "/pipeline1/phylogenize_out/", study_id, "/phylogenize.html"),
      out_dir = paste0(working_dir, "/pipeline1/phylogenize_out/", study_id),
      data_dir = dataset_dir,
      db = "human-gut",
      taxon_level = "family",
      type_16S=F,
      which_phenotype = phenotype,
      diff_abund_method = "ancombc2",
      abundance_file = paste0(working_dir, "/pipeline1/phylogenize_out/", study_id, "/data_w_count.tsv"),
      metadata_file = paste0(working_dir, "/pipeline1/phylogenize_out/", study_id, "/metadata_filter.tsv"),
      which_envir = ref_env, 
      env_column = "env",
      sample_column = "sample",
      ncl=parallelly::availableCores())
    
  } else { #provided phenotype
    res = phylogenize(
      output_file = paste0(working_dir, "/pipeline2/phylogenize/out/phylogenize.html"),
      out_dir = paste0(working_dir, "/pipeline2/phylogenize/out"),
      data_dir = dataset_dir,
      db = "human-gut",
      taxon_level = "family",
      which_phenotype = "provided",
      diff_abund_method = "ancombc2",
      phenotype_file = provided_file_path,
      abundance_file = paste0(working_dir, "/phylogenize_out/CHN/data_w_count.tsv"),
      metadata_file = paste0(working_dir,"/phylogenize_out/CHN/metadata_filter.tsv"),
      which_envir = ref_env,
      env_column = "env",
      sample_column = "sample",
      ncl = parallelly::availableCores())
    
  }
  #saveRDS(res, paste0(working_dir, "/phylogenize_out/", study_id, "/phylogenize.rds"))
  return(res)
}

#PIPELINE2: function for pre-phylogenize merging:
#5. function to run ancombc inputs separately:
ancomRun = function(count_tbl, study_id, metadata){
  
  count = count_tbl %>% column_to_rownames(var = "name")
  metadata_tbl = extract_phyloz_metadata(import_bracken_out = count,
                                         metadata = metadata,
                                         study_id = study_id,
                                         envs_compared = c("ND CTRL", "T2D metformin-"))
  metadata_tbl = metadata_tbl %>% column_to_rownames(var = "sample") 
  
  #filter out those samples in the metadata but not in count_tbl:
  count_tbl_clean = count %>% dplyr::select(rownames(metadata_tbl))
  tse <- TreeSummarizedExperiment(assays = S4Vectors::SimpleList(counts = as.matrix(count_tbl_clean)),
                                  colData = S4Vectors::DataFrame(metadata_tbl)
  )
  #apply ancombc2:
  print("fitting ANCOMBC")
  res_ancom=ANCOMBC::ancombc2(data = tse,
                              assay_name = "counts",
                              fix_formula = "env")
  
  saveRDS(res_ancom,
          paste0(working_dir, "/pipeline2/ancom/", study_id, "_ancombc_res.rds"))
  
  return(list(res = res_ancom,
              count_input = count_tbl_clean,
              metadata_input = metadata_tbl))
}



#6. function to calculate combined mu and combined se:
#input is a list of ashr output of all study_id:
#taking both ancombc output and ashr outputs (parameter: ash_used = T/F)
#formula to compute for combined effect size: -for 1 study_id first
weightCal=function(SE) { #input: se of 1 study_id only 
  weight_per_study = 1/(SE**2)
  return(weight_per_study)  
}

#a function that loops through ancom_res of each study_id and extract for mu and se:
ancomExtract=function(ancomRun_output, 
                      envs_compared = c("ND CTRL",
                                        "T2D metformin-")) { #of 1 study_id only
  if (identical(envs_compared, c("ND CTRL", "T2D metformin-"))) {
    mu_se_df=data.frame(taxon = ancomRun_output$res$taxon,
                        mu=ancomRun_output$res$`lfc_envT2D metformin-`,
                        se=ancomRun_output$res$`se_envT2D metformin-`)}
  mu_se_df
}


weight_per_taxonCal=function(ancomRun_output_ls, study_id_ls){
  weight_per_taxon=lapply(ancomRun_output_ls, function(x) { #per study here:
    mu_se_df=ancomExtract(ancomRun_output = x,
                          envs_compared = c("ND CTRL", "T2D metformin-"))
    weight_per_taxon=data.frame(taxon=mu_se_df$taxon,
                                weight=weightCal(SE = mu_se_df$se),
                                mu = mu_se_df$mu)
    weight_per_taxon #return weight (per_taxon) for each study_id in a list
  })
  names(weight_per_taxon)<-study_id_ls
  return(weight_per_taxon) #should be a ls of df with 2 columns - taxon and each of their weight
  
}

sigma_weightCal <- function(weight_per_taxon_ls) {
  # combine list of dfs into one long df, tagging study
  combined <- purrr::imap_dfr(weight_per_taxon_ls, ~ dplyr::mutate(.x, study = .y))
  
  # sum weights per taxon across studies (NA-safe: taxa absent in a study just don't contribute)
  combine_df = combined %>% group_by(taxon) %>% summarise(sum_weight = sum(weight, na.rm = TRUE), 
                                                          n_studies  = n(),  # how many studies contributed to this taxon
                                                          .groups = "drop")
  combine_df
}


#loop throgh all studyId and apply weightCal--> to get sigma_weightCal_out
SE_poolCal = function(sigma_weightCal_out) { #input: sum of all weigtCal_out
  sigma_weightCal_out["se_pooled"] <- sqrt(1/sigma_weightCal_out$sum_weight)
  return(sigma_weightCal_out)
}

sum_prod_weightCal = function(weight_per_taxon_ls,
                              study_id_ls) 
{
  prod_ls=lapply(weight_per_taxon_ls, function(i) {
    i["prod_weight_mu"]<- i$mu * i$weight
    i 
  })				      
  names(prod_ls)<-study_id_ls
  combined <- purrr::imap_dfr(prod_ls, ~ dplyr::mutate(.x, study = .y))
  combined %>% 
    dplyr::group_by(taxon) %>%
    dplyr::summarise(
      sum_prod_weight_mu = sum(prod_weight_mu, na.rm = T),
      n_studies = n(),
      .groups = "drop")
}

mu_poolCal = function(sum_prod_weightCal_out,
                      sigma_weightCal_out) { #in = sum of all w (from weightCal)
  df_comb = inner_join(sum_prod_weightCal_out, sigma_weightCal_out,
                       by = "taxon")
  df_comb["mu_pooled"] = df_comb$sum_prod_weight_mu / df_comb$sum_weight
  df_comb=df_comb %>% dplyr::select(taxon, 
                                    sum_prod_weight_mu,
                                    n_studies.x,
                                    sum_weight,
                                    mu_pooled
  )
  return(df_comb)
}


#7. function to run ashr:- clean taxon name:
#works only for 1 study_id: all_ancom_res[[1]]
ashRun = function(mu_poolCal_out, SE_poolCal_out) {
  #combine 2 inputs so the taxon matched:
  interDf=inner_join(mu_poolCal_out, SE_poolCal_out, by = "taxon") %>%
    dplyr::select(taxon, n_studies, se_pooled, mu_pooled)
  bh = interDf$mu_pooled
  sebh = interDf$se_pooled
  
  ash_run = ashr::ash(betahat = bh,
                      sebetahat = sebh)
  ash_res = ash_run$result
  
  #attach taxon name from ancom to ash res:
  res_tbl = cbind(interDf$taxon,
                  ash_res)
  #rename taxon column:
  colnames(res_tbl)[1]<-"taxon"
  return(res_tbl)
  
}

#==============ANALYSIS OF FINAL BIOLOGICAL HITS:
#-----pipeline1:
#8. function to read into all-result.csv -> filter q.value<0.05 -> calculate for PD: ->return all info:
readRes_n_pdCal=function(study_id) { #study_id = CHN, MH1, MH3
  
  res=read.csv(paste0(working_dir, "/all-results-", study_id, ".csv")) #change if working on HPC
  sig=res[res$q.value < 0.05, ]
  sig_ls=split(sig, sig$taxon)
  length(sig_ls)
  
  #for each family, calculate PD for each gene:
  pd=lapply(seq_along(sig_ls), function(i){
    
    family=names(sig_ls)[i] #obtain the name of the family
    
    gene_family=gene_pres[[family]] #only take gene-pres of that family i
    genePres=gene_family[rownames(gene_family) %in% sig_ls[[i]]$gene, ,drop=F] #for each family, only take the sig. genes
    pres=as.matrix(genePres)
    
    #filter out trees: tip_label = species in the final presence gene matrix:
    tree_family = tree[[family]]
    tree_fin=keep.tip(tree_family, colnames(pres))
    
    #pd calculated:
    pd_df=pd(pres, tree_family)
    pd_df 
    
  }) 
  
  #merge pd info with all_res:
  res_final=bind_rows(pd) %>% 
    rownames_to_column(var="gene") %>%
    inner_join(sig, by = "gene") %>%
    arrange(-PD)
  return(res_final)
}

#9. Make a function that takes list of all res across all study_id
# -> and outputs an upset plot ("negative" or "positive" direction) Or
# an heatmap plot
taxaLevel_plotMake <- function(taxa_list, direction_es, plot_type) {

  taxa_list_filt <- lapply(taxa_list, function(f) {
    if (direction_es == "positive"){
      f[f$effect.size > 0, ]
    } else {
      f[f$effect.size < 0, ]
    }
  })
  if (plot_type == "upset") {
  long_df <- bind_rows(
    lapply(names(taxa_list_filt), function(ds) {
      data.frame(taxon = unique(taxa_list_filt[[ds]]$taxon), dataset = ds)
    })
  )
  
  binary_matrix <- long_df |>
    mutate(present = 1) |>
    pivot_wider(names_from = dataset, values_from = present, values_fill = 0)
  
  upset(
    as.data.frame(binary_matrix[, c("CHN","MH1","MH3")]),  # numeric 0/1 columns only
    sets = c("CHN","MH1","MH3"),
    order.by = "freq",
    mainbar.y.label = paste0("Intersection Size - ",direction_es),
    sets.x.label = "Overlapping Taxa",
    text.scale = 1.75
  )}
  else {
    heat_df <- taxa_list_filt |>
      lapply(\(a) {
        a |>
          count(taxon, name = "taxa_hits_count")}) |>
      imap(\(df, dataset) {
        df$dataset <- dataset
        df
      }) |>
      bind_rows() |>
      pivot_wider(
        names_from = dataset,
        values_from = taxa_hits_count,
        values_fill = 0
      )
    
    heat_mat = heat_df |> column_to_rownames(var="taxon") |>
      as.matrix()
    if (direction_es == "positive") {
      col_fun = circlize::colorRamp2(c(0, max(heat_mat) / 2, max(heat_mat)),
                                     c("white", "#9a0d1b", "#400000"))
    } else {
      col_fun <- circlize::colorRamp2(c(0, max(heat_mat) / 2, max(heat_mat)),
                                      c("white","#03396c", "#131e3a"))
    }
    #reordering block of heatmap so higher values go on top
    heat_mat <- heat_mat[
      order(apply(heat_mat, 1, max), decreasing = TRUE),
    ]
    
    ComplexHeatmap::Heatmap(
      heat_mat,
      name = "# significant genes",
      col = col_fun,
      cluster_rows = FALSE
    )
  }
  
}

#10.a function that takes in the list of all_res dfs (e.g CHN_all_res outputed from readRes_n_pdCal func)
#and outputs genes in taxa that are consistent across study_ids
#count_feature column = number of dataset that gene shows up
consensus_geneFind=function(all_res_ls, direction_es) {
  feature_df=lapply(all_res_ls, function(s) {
    s$"feature"<- paste0(s$taxon, "__", s$gene)
    s
  }) |>
    imap(\(df, dataset) {
      df$dataset <- dataset
      df
    }) 
  
  if (direction_es == "positive") {
    feature_df = bind_rows(feature_df) |>
      filter(effect.size > 0) 
  } else {
    feature_df = bind_rows(feature_df) |>
      filter(effect.size < 0) 
  }
  feature_df |> group_by(feature) |>
    summarize(count_feature = n()) %>% 
    arrange(-count_feature) |> 
    separate(feature, into = c("taxon", "gene"), sep = "__") |>
    left_join(annotation_df, by=c("gene")) 
}

#11.A function that inputs all-results.csv from pre-Phylogenize merging
#and outputs bar plots of # of sig.genes across taxa:

Pre_mergePlot=function(mergeRes) #all-results from pre-Phylogenize merging:
{
  mergeRes = merge_res
  pos_merge=mergeRes |>
    dplyr::filter(q.value < 0.05, effect.size > 0) |>
    group_by(taxon) |>
    summarize(count_taxon = n(), .groups = "drop") |>
    arrange(desc(count_taxon)) |>
    mutate(es = "positiveES")
  
  neg_merge=mergeRes |>
    dplyr::filter(q.value < 0.05, effect.size < 0) |>
    group_by(taxon) |>
    summarize(count_taxon = n(), .groups = "drop") |>
    arrange(desc(count_taxon)) |>
    mutate(es = "negativeES")
  
  merge_df=rbind(pos_merge,neg_merge) 
  
  merge_fix = pos_merge |> rbind(
    data.frame(
      taxon=anti_join(merge_df, pos_merge, by = "taxon")$taxon,
      count_taxon = 0,
      es = "positiveES"
    )) |> rbind(
      data.frame(
        taxon=anti_join(merge_df, neg_merge, by = "taxon")$taxon,
        count_taxon = 0,
        es = "negativeES"
      ) |> rbind(neg_merge))
  
  merge_fix |> ggplot(aes(x = taxon, y = count_taxon, fill = es)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    theme_bw() +
    labs(
      x = "Taxon",
      y = "Count of significant genes (pre-Phylogenize merge)"
    ) + theme(
      axis.title.x = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")) + 
    coord_flip() +
    scale_fill_manual(
      values = c("negativeES" = "#131e3a", "positiveES" = "darkred"))
}

#-----------PIPELINE 1 VS PIPELINE 2:
#11 OPTIONAL.a function how much of each study_id's all-result (post-Phylogenize merging) 
#taxa overlaps with all-results of pre-Phylogenize merging:
pre_v_postFind_eachStudy <- function(merge_resPre, res_listPost, direction_es) {
  
  ## PRE-PHYLOGENIZE
  if (direction_es == "positive") {
    gene_preMerge <- merge_resPre %>%
      filter(effect.size > 0, q.value < 0.05)
  } else {
    gene_preMerge <- merge_resPre %>%
      filter(effect.size < 0, q.value < 0.05)
  }
  
  taxa_preMerge <- gene_preMerge %>%
    count(taxon, name = "taxon_count_pre")
  
  ## POST-PHYLOGENIZE
  taxa_post <- imap(res_listPost, function(df, dataset) {
    
    col_name <- paste0("taxon_count_post_", dataset)
    
    if (direction_es == "positive") {
      df <- df %>% filter(effect.size > 0)
    } else {
      df <- df %>% filter(effect.size < 0)
    }
    
    df %>%
      count(taxon, name = col_name)
    
  })
  
  ## Merge all post datasets by taxon
  taxa_postMerge <- reduce(
    taxa_post,
    full_join,
    by = "taxon"
  )
  
  ## Merge with pre counts
  result <- taxa_preMerge %>%
    full_join(taxa_postMerge, by = "taxon") %>%
    mutate(
      across(where(is.numeric), ~ replace_na(.x, 0))
    )
  
  return(result)
}

#12 OPTIONAL. a function to make 2 heatmaps from function 11:
#1 is at taxa level: which each study_id (post Merge) taxa overlapping w that to pre-Merge
#2 is at gene level: same as above but es comapred for those shared gene/taxa.
preVpost_HeatmapMake = function(direction_es, merge_resPre, res_listPost) {
  
  # TAXA LEVEL:
  #building a df for TAXA heatmap
  if (direction_es == "positive") {
    preVpost=pre_v_postFind_eachStudy(merge_resPre = merge_resPre,
                                          res_listPost = res_listPost,
                                          direction_es = "positive") |>
      column_to_rownames("taxon") |>
      as.matrix() 
    
    #set up color for making heatmap afterward
    col_setup=circlize::colorRamp2(c(0, max(preVpost) / 2, max(preVpost)),
                         c("white", "#9a0d1b", "#400000"))
    
  } else {
    preVpost=pre_v_postFind_eachStudy(merge_resPre = merge_resPre,
                                          res_listPost = res_listPost,
                                          direction_es = "negative") |>
      column_to_rownames("taxon") |>
      as.matrix() 
    
    col_setup=circlize::colorRamp2(c(0, max(preVpost) / 2, max(preVpost)),
                                   c("white","#03396c", "#131e3a"))
  }
  
  #making heatmap:
  col_group <- data.frame(
    Pipeline = c("Pre-Phylogenize Merging", rep("Post-Phylogenize Merging", length(res_listPost)))
  )
  
  rownames(col_group) <- colnames(preVpost)
  colnames(preVpost) <- c("Pre", names(res_listPost))
  
  ha <- HeatmapAnnotation(
    Pipeline = col_group$Pipeline,
    col = list(
      Pipeline = c(
        "Pre-Phylogenize Merging"  = "#e69f00",   # orange
        "Post-Phylogenize Merging" = "#009e73"    # teal
      )
    )
  )
  
  taxa_preVpost_hm=ComplexHeatmap::Heatmap(matrix = preVpost, 
                          cluster_columns = F,
                          top_annotation = ha,
                          name = paste0("# of significant genes\n", direction_es),
                          col = col_setup)
  
  # GENE LEVEL:
  genePrevPost=lapply(seq_along(res_listPost), function(i) {
      if (direction_es == "positive") {
      test_gene = inner_join(merge_resPre[merge_resPre$q.value<0.05 & merge_resPre$effect.size > 0,],
                                res_listPost[[i]][res_listPost[[i]]$effect.size > 0,],
                                by = c("gene", "taxon")) |> 
        dplyr::select(taxon, gene, effect.size.x, effect.size.y) |>
        dplyr::rename("effect.size_pre"="effect.size.x",
                      "effect.size"="effect.size.y") |>
        left_join(annotation_df, by = "gene") |>
        mutate(study_id = names(res_listPost)[i])
      } else {
        test_gene = inner_join(merge_resPre[merge_resPre$q.value<0.05 & merge_resPre$effect.size < 0,],
                                  res_listPost[[i]][res_listPost[[i]]$effect.size < 0,],
                                  by = c("gene", "taxon")) |> 
          dplyr::select(taxon, gene, effect.size.x, effect.size.y) |>
          dplyr::rename("effect.size_pre"="effect.size.x",
                        "effect.size"="effect.size.y") |>
          left_join(annotation_df, by = "gene") |>
          mutate(study_id = names(res_listPost)[i])
      }
      
      if (nrow(test_gene) > 0) {
        test_gene = test_gene |> select(-any_of(c("accession","function")))
        ## Matrix of effect sizes
        heat_mat <- test_gene |>
          dplyr::select(effect.size_pre, effect.size) |>
          mutate(across(everything(), as.numeric)) |>
          as.matrix()
        
        rownames(heat_mat) <- test_gene$gene
        colnames(heat_mat) <- c("Pre", "Post")
        
        ## Row annotations
        taxa <- unique(test_gene$taxon)
        
        cols <- RColorBrewer::brewer.pal(max(3, length(taxa)), "Set2")[seq_along(taxa)]
        
        taxa_colors <- setNames(cols, taxa)
        row_ha <- rowAnnotation(
          Study = test_gene$study_id,
          Taxon = test_gene$taxon,
          col = list(
            Study = c(
              CHN = "#66C2A5",
              MH1 = "#FC8D62",
              MH3 = "#8DA0CB"
            ),
            Taxon = taxa_colors
          ),
          annotation_name_side = "top"
        )
        
        ## Color scale
        lim <- max(abs(heat_mat))
        if(direction_es=="positive"){
          col_fun = colorRamp2(
            c(0, lim/2, lim),
            c("white","#FB9A99","#B2182B")
          )
        }else{
          col_fun = colorRamp2(
            c(-lim,-lim/2,0),
            c("#131e3a","#03396c","white")
          )
        }
        
        ## Heatmap
        gene_preVpost_hm=Heatmap(
          heat_mat,
          name = "Effect size",
          left_annotation = row_ha,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_names_side = "right",
          column_names_side = "bottom",
          col = col_fun)
      } #close of if
    else {
      return(NULL)
    }
      
    return(gene_preVpost_hm)
    })
  
  return(list(taxa_level_heatmap = taxa_preVpost_hm,
              gene_level_heatmap = genePrevPost))
}
    
#------Expand to other datasets:
# explicitly loads it into your environment:
#13. function to retrieve metadata info for each T2D studies found in curatedMetagenomicData:
#output a list of metadata control and metadata t2d:
metadataRetrieve=function(study_id, confounder_para = c("metformin")) #
{

  study_code=data.frame(
    name_CM = c("SankaranarayananK_2015","MetaCardis_2020_a","HMP_2019_t2d"),
    id = c("SKK","MCA","HMP"))
  
  std_name = study_code[study_code$id == study_id,]$name_CM
  metadata = lapply(c("healthy","T2D"), function(x) #create metadata where [[1]] is for control, [[2]] is for T2D.
  {
    subj = sampleMetadata |> 
      filter(study_name == std_name,
             disease == x) |>
      dplyr::select(study_name, antibiotics_current_use, sample_id,
                    disease, NCBI_accession, treatment)  
    subj
  })
  names(metadata)<-c("healthy","T2D")
  
  #----filter out those used antibiotics for control:
  #----filter out those used metformin:
  metadata = lapply(metadata, function(t) 
    t |> filter(!grepl(confounder_para, treatment),
                antibiotics_current_use != "yes"))
  
  return(metadata)
}

#14. a function to get sample ID --> obtain URL --> use that URL to retrieve the fastq_ftp:
ftpRetrieve <- function(metadataRetrieve_out) {
  #metadataRetrieve_out = metadata_MCA_t2d
  sampleDf <- metadataRetrieve_out |>
    filter(!is.na(NCBI_accession)) |>
    pull(NCBI_accession) |>
    strsplit(";") |>
    unlist() |>
    trimws() |>
    data.frame(sample_id = _)
  
  sampleDf <- sampleDf |>
    filter(sample_id != "")
  
  sampleDf$URL <- paste0(
    "https://www.ebi.ac.uk/ena/portal/api/filereport?",
    "accession=", sampleDf$sample_id,
    "&result=read_run",
    "&fields=fastq_ftp"
  )
  
  colnames(sampleDf)<-c("sample_id", "URL")
  
  #map URL with fastp_ftp:
  fastq_ftp <- lapply(seq_len(nrow(sampleDf)), function(i){
    url = sampleDf[i, ]$URL
    txt = httr::content(httr::GET(url),
                  as = "text", 
                  encoding = "UTF-8")
    fastq_ftp = strsplit(txt, "\n")[[1]][2] |> strsplit("\t") |> (\(y) y[[1]][2])()
  }) 
  
  fastq_vec <- unlist(fastq_ftp)  # character vector, one string per row, "url1;url2"
  
  fastq_df <- data.frame("fastq_ftp"=fastq_vec) %>%
    cbind(sampleDf) %>%
    separate_rows("fastq_ftp", sep = ";") 
  
  return(list(sample_list = sampleDf,
              fastq_df = fastq_df))
}

#--------------PIPELINE 1A: using t-statistics and Satterthwaite for DF estimation:
#15. extract DF (# of case and ctrl subjects) for each study
extractDF_study_func = function(study_id_ls)  {
  setNames(lapply(study_id_ls, function(id){
    if (id %in% c("CHN", "MH1", "MH3")) {
      metadata = read.csv(paste0(working_dir, "/", id, "_md_final.csv"))
      df_study = length(unique(metadata$Sample))
    } else { #for SKK and CHN: read into outputs from metadataRetrieve() outputs:
      metadata = bind_rows(metadataRetrieve(study_id = id))
      df_study = length(unique(metadata$sample_id))
    }
  }), study_id_ls)
}

#16. convert p-value to t-statistics from each study df and p-value:
#input: CHN_all_res, MH1_all_res, MH3_all_res, SKK_all_res, MCA_all_res + their id: "MCA", "SKK", "MH1", "MH3", etc
tstatCal = function(study_id_res, study_id, extractDF_study_id) {
  #study_id_res = CHN_all_res
  #study_id = "CHN"
  #extractDF_study_out= extractDF_study_out[[study_id]]
  p_capped = pmin(pmax(study_id_res$p.value, 1e-300), 1 - 1e-16) #cap p-values to avoid Inf/0 in t-statistic and SE
  qt(p_capped/2, 
     df = extractDF_study_id - 2,
     lower.tail = F) #n_case + n_ctrl - 2
}

#17. back-calculate se from effect size and t-statistic
se_recCal=function(study_id_res, tstatCal_out) {
  abs(study_id_res$effect.size) / tstatCal_out
}

#18. get variance from se recovered:
varianceCal = function(se_recCal_out) {
  se_recCal_out^2
}

#19. functions to compute Satterthwaite DF:
Satt_DFCal = function(varianceCal_out_combined) #input: list of all_res_with
{ varianceCal_out_combined |> group_by(taxon, gene) |>
  summarize(sum_var = sum(variance, na.rm = T),
            sum_v2_df = sum(variance^2/study_DF, na.rm = T),
            Satt_DF = sum_var / sum_v2_df,
            weight_sum = sum(1 / variance, na.rm = TRUE),
            effsize_pooled = sum(effect.size / variance, na.rm = TRUE) / weight_sum,
            se_pooled = sqrt(1 / weight_sum),
            n_studies  = n(),
            .groups = "drop")
}
  


