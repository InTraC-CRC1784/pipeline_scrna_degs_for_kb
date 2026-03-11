#!/usr/bin/env Rscript
library(edgeR)
library(dplyr)
library(magrittr)
library(tidyr) 

###-------- Arg Parser --------
parse_args <- function(args) {
  res <- list()
  for (a in args) {
    if (grepl("^--", a)) {
      kv <- sub("^--", "", a)
      key <- sub("=.*", "", kv)
      value <- sub("^[^=]*=", "", kv)
      res[[key]] <- value
    }
  }
  res
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
cli_args <- parse_args(args)

# --- Configuration ---
RESULTS_DIR                 <- cli_args$RESULTS_DIR %||% "EDGER_RESULT"
pb_whole_path               <- cli_args$pb_whole_path %||% "pseudobulk_whole.csv"
pb_filtering_path           <- cli_args$pb_filtering_path %||% "pseudobulk_filtering.csv"
cell_state_translation_path <- cli_args$cell_state_translation_path %||% "cell_state_translation_table.csv"
cell_state_number_path      <- cli_args$cell_state_number_path %||% "cell_state_number.csv"
ontology_path               <- cli_args$ontology_path %||% NULL
gene_mapping_path           <- cli_args$gene_mapping_path %||% NULL
threshold                   <- as.numeric(cli_args$threshold %||% 0.0125)
output_file_name            <- cli_args$output_file_name %||% "edgeR_results.csv"
Cell_Type_Col               <- cli_args$Cell_Type_Col %||% "cell_origin"
Cell_State_Col              <- cli_args$Cell_State_Col %||% "cell_type"
CELL_LEVEL                  <- cli_args$CELL_LEVEL %||% "celltype"
REGION                      <- cli_args$REGION %||% NA
Condition_Col               <- cli_args$Condition_Col %||% "condition"
Sample_Information_Col      <- cli_args$Sample_Information_Col %||% NA
comparison_normal_value     <- cli_args$comparison_normal_value %||% "control"
min_cells_per_state         <- as.integer(cli_args$min_cells_per_state %||% 5)
fdr_threshold               <- as.numeric(cli_args$fdr_threshold  %||% 1)

# --- Metadata Arguments ---
species_arg    <- cli_args$species    %||% NA
study_name_arg <- cli_args$study_name %||% NA
study_id_arg   <- cli_args$study_id   %||% NA
study_year_arg <- cli_args$study_year %||% NA
disease_id_arg <- cli_args$disease_id %||% NA
disease_arg    <- cli_args$disease    %||% NA
tissue_arg     <- cli_args$tissue     %||% NA
tissue_id_arg  <- cli_args$tissue_id  %||% NA

# filename for filtered output
ext <- strsplit(output_file_name, ".", fixed=T)[[1]][-1]
output_filtered_file_name <- gsub(paste0("\\.", ext, "$"), paste0("_filtered.", ext), output_file_name, perl=TRUE)

# --- 1. Data Loading ---
cat("Loading data\n")
cat(paste0("Pseudobulk: ", pb_whole_path, "\n"))
x <- read.csv(pb_whole_path, row.names = 1, check.names = FALSE)
cat(paste0("Gene filtering list: ", pb_filtering_path, "\n"))
genes_tofilter <- read.csv(pb_filtering_path, row.names = 1, check.names = FALSE)
cat(paste0("cell state to cell type mapping: ", cell_state_translation_path, "\n"))
CELLTYPE_STATE <- read.csv(cell_state_translation_path, check.names = FALSE)
colnames(CELLTYPE_STATE) <- c(Cell_State_Col, Cell_Type_Col)

CELLTYPE_FILTER <- read.csv(cell_state_number_path, check.names = FALSE)
CELLTYPE_FILTER[, -1] <- lapply(CELLTYPE_FILTER[, -1, drop=FALSE], function(i) i >= min_cells_per_state)

l <- strsplit(colnames(x), "__")
meta.data <- as.data.frame(do.call(rbind, l))
rownames(meta.data) <- colnames(x)
colnames(meta.data) <- c(Condition_Col, Cell_State_Col, Sample_Information_Col)

GENOTYPES <- unique(meta.data[[Condition_Col]][meta.data[[Condition_Col]] != comparison_normal_value])
results_list <- list()

# --- 2. Main Analysis Loop ---
cat("Running edgeR\n")
for (GENOTYPE in as.character(GENOTYPES)) {
  current_states <- unique(meta.data[[Cell_State_Col]])
  for (STATE_ID in as.character(current_states)) {
    CELLTYPE <- CELLTYPE_STATE[[Cell_Type_Col]][CELLTYPE_STATE[[Cell_State_Col]] == STATE_ID][1]
    if (is.na(CELLTYPE)) next
    
    samples_mask <- meta.data[[Condition_Col]] %in% c(comparison_normal_value, GENOTYPE) &
                    meta.data[[Cell_State_Col]] == STATE_ID
    
    meta_sub <- meta.data[samples_mask, , drop = FALSE]
    state_filter_row <- CELLTYPE_FILTER[CELLTYPE_FILTER[[1]] == STATE_ID, ]
    if (nrow(state_filter_row) == 0) next
    
    PATIENTS_TOKEEP <- colnames(state_filter_row)[which(as.logical(state_filter_row))]
    clean_meta_ids <- gsub("[._-]", "", as.character(meta_sub[[Sample_Information_Col]]))
    clean_keep_ids <- gsub("[._-]", "", PATIENTS_TOKEEP)
    
    final_mask <- clean_meta_ids %in% clean_keep_ids
    meta_sub <- meta_sub[final_mask, , drop = FALSE]
    x_sub <- x[, rownames(meta_sub), drop = FALSE]
    
    cond_counts <- table(as.character(meta_sub[[Condition_Col]]))
    
    if (length(cond_counts) == 2 && all(cond_counts >= 3)) {
      meta_sub[[Condition_Col]] <- factor(meta_sub[[Condition_Col]], levels = c(comparison_normal_value, GENOTYPE))
      dge <- DGEList(counts = x_sub, group = meta_sub[[Condition_Col]])
      design <- model.matrix(as.formula(paste0("~", Condition_Col)), data = meta_sub)
      dge <- estimateDisp(dge, design)
      fit <- glmQLFit(dge, design)
      qlf <- glmQLFTest(fit)
      
      # --- Calculation of Means per Group ---
      pat_pattern <- paste0("__", STATE_ID, "__")
      sub_filter  <- genes_tofilter[, grep(pat_pattern, colnames(genes_tofilter)), drop=FALSE]
      mean_genotype <- rowMeans(sub_filter[, grepl(paste0("^", GENOTYPE, "__"), colnames(sub_filter)), drop=FALSE])
      mean_control  <- rowMeans(sub_filter[, grepl(paste0("^", comparison_normal_value, "__"), colnames(sub_filter)), drop=FALSE])
      
      mean_df <- data.frame(
        Gene = rownames(genes_tofilter),
        mean_exp_control  = mean_control,
        mean_exp_genotype = mean_genotype,
        stringsAsFactors = FALSE
      )

      keep_genes    <- rownames(genes_tofilter)[mean_control > threshold | mean_genotype > threshold]
      
      tt_all  <- topTags(qlf, n = Inf)$table %>% mutate(Gene = rownames(.))
      tt_filt <- topTags(qlf[intersect(keep_genes, rownames(qlf)), ], n = Inf)$table %>% 
                 mutate(Gene = rownames(.)) %>% select(Gene, FDR_onlyhigh = FDR)
      
      # --- Merging Stats, FDR, and Means ---
      tt_merged <- tt_all %>%
        left_join(tt_filt, by = "Gene") %>%
        left_join(mean_df, by = "Gene") %>%
        mutate(
          low_expression = ifelse(is.na(FDR_onlyhigh), "T", "F"),
          FDR_onlyhigh   = ifelse(is.na(FDR_onlyhigh), 1, FDR_onlyhigh),
          cell_type_name_raw = STATE_ID,
          cell_heart     = CELLTYPE,
          comparison     = paste0(comparison_normal_value, "_", GENOTYPE),
          Observations_genotype = as.numeric(cond_counts[GENOTYPE]),
          Observations_reference = as.numeric(cond_counts[[comparison_normal_value]])
        )
      results_list[[paste(GENOTYPE, STATE_ID, sep="_")]] <- tt_merged
    }
  }
}

# --- 3. Final Formatting & Save ---
cat("formatting results\n")
if (length(results_list) > 0) {
  final_df <- do.call(rbind, results_list)
  rownames(final_df) <- NULL
  
  # Assign Metadata
  final_df$species     <- species_arg
  final_df$study_name  <- study_name_arg
  final_df$dataset_id  <- study_id_arg
  final_df$study_year  <- study_year_arg
  final_df$disease_id  <- disease_id_arg
  final_df$disease     <- disease_arg
  final_df$test        <- "edgeR"
  final_df$Region      <- REGION
  final_df$tissue      <- tissue_arg
  final_df$tissue_id   <- tissue_id_arg

  # Standard Rename block
  final_df <- final_df %>%
    rename(
      gene           = Gene,
      cell_type_name = cell_type_name_raw,
      log2fc         = logFC,
      p_value        = PValue,
      adj_p          = FDR_onlyhigh,
      comparison_id  = comparison
    )

  final_df <- final_df %>%
    separate(comparison_id, into = c("control", "case"), sep = "_", remove = FALSE)

  # --- Ontology Mapping ---
  if (!is.null(ontology_path) && file.exists(ontology_path)) {
    ontology_df <- read.csv(ontology_path, stringsAsFactors = FALSE)
    ontology_df <- ontology_df %>%
      rename(cell_type_id = ontology_id, cell_type = ontology_name)
    
    if ("cell_type_name" %in% colnames(ontology_df)) {
      final_df <- final_df %>% left_join(ontology_df, by = "cell_type_name", relationship = "many-to-many")
    }
  }

  # --- Gene ID Mapping ---
  
  if (!is.null(gene_mapping_path) && file.exists(gene_mapping_path)) {
    gene_map <- read.csv(gene_mapping_path, stringsAsFactors = FALSE)
    if ("gene" %in% colnames(gene_map)) {
      # Use many-to-many to allow one gene symbol to match multiple IDs/rows as expected
      final_df <- final_df %>% left_join(gene_map, by = "gene", relationship = "many-to-many")
    }
  }

  if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
  write.csv(final_df, file.path(RESULTS_DIR, output_file_name), row.names = FALSE)
  message("SUCCESS: Results saved to ", RESULTS_DIR)
    
  # --- filtering by FDR and saving ---
  filtered_df <- final_df %>% filter(adj_p <= fdr_threshold)
  write.csv(filtered_df, file.path(RESULTS_DIR, output_filtered_file_name), row.names = FALSE)
    
} else {
  message("ERROR: No results generated.")
}