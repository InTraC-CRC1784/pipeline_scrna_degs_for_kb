library(edgeR)
library(dplyr)
library(magrittr)
###-------- Arg Parser
args <- commandArgs(trailingOnly = TRUE)

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

cli_args <- parse_args(args)

if (!is.null(cli_args$RESULTS_DIR)) {
  RESULTS_DIR <- cli_args$RESULTS_DIR
}

if (!is.null(cli_args$pb_whole_path)) {
  pb_whole_path <- cli_args$pb_whole_path
}

if (!is.null(cli_args$pb_filtering_path)) {
  pb_filtering_path <- cli_args$pb_filtering_path
}

if (!is.null(cli_args$cell_state_translation_path)) {
  cell_state_translation_path <- cli_args$cell_state_translation_path
}

if (!is.null(cli_args$cell_state_number_path)) {
  cell_state_number_path <- cli_args$cell_state_number_path
}

if (!is.null(cli_args$threshold)) {
  threshold <- as.numeric(cli_args$threshold)
}

if (!is.null(cli_args$output_file_name)) {
  output_file_name <- cli_args$output_file_name
  ext <- strsplit(output_file_name, ".", fixed=T)[[1]][-1]
  output_filtered_file_name <- gsub(paste0("\\.", ext, "$"), paste0("_filtered.", ext), output_file_name, perl=TRUE)
}

# Cell_Type_Col
if(!is.null(cli_args$Cell_Type_Col)){
  Cell_Type_Col<- cli_args$Cell_Type_Col
}

# Cell_State_Col
if(!is.null(cli_args$Cell_State_Col)){
  Cell_State_Col <- cli_args$Cell_State_Col
}
#CELL_LEVEL
if(!is.null(cli_args$CELL_LEVEL)){
  CELL_LEVEL<- cli_args$CELL_LEVEL
}
# REGION
if(!is.null(cli_args$REGION)){
  REGION<- cli_args$REGION
}
# Condition_Col
if(!is.null(cli_args$Condition_Col)){
  Condition_Col<- cli_args$Condition_Col
}
# Sample_Information_Col
if(!is.null(cli_args$Sample_Information_Col)){
  Sample_Information_Col<- cli_args$Sample_Information_Col
}
# comparison_normal_value
if(!is.null(cli_args$comparison_normal_value)){
  comparison_normal_value<- cli_args$comparison_normal_value
}
# species
if(!is.null(cli_args$species)){
  species <- cli_args$species
}
# year
if(!is.null(cli_args$year)){
  year<- cli_args$year
}
# paper
if(!is.null(cli_args$paper)){
  paper<- cli_args$paper
}

if(!is.null(cli_args$min_cells_per_state)){
  min_cells_per_state <- as.integer(cli_args$min_cells_per_state)
}

if(!is.null(cli_args$fdr_threshold)){
  fdr_threshold <- as.numeric(cli_args$fdr_threshold)
} else {
  fdr_threshold <- 1
}

###-------- 


# --- 1. Data Loading & Preprocessing ---
x <- read.csv(pb_whole_path, row.names = 1)
genes_tofilter <- read.csv(pb_filtering_path)
CELLTYPE_STATE <- read.csv(cell_state_translation_path)
colnames(CELLTYPE_STATE) <- c(Cell_State_Col, Cell_Type_Col)

CELLTYPE_FILTER <- read.csv(cell_state_number_path, check.names = FALSE)
# Convert counts to boolean (TRUE if >= 5)

CELLTYPE_FILTER[, -1] <- lapply(CELLTYPE_FILTER[, -1], function(i) i >= min_cells_per_state)

# Prepare metadata from column names (Format: Condition__CellState__SampleID)
l <- strsplit(colnames(x), "__")
meta.data <- as.data.frame(do.call(rbind, l))
rownames(meta.data) <- colnames(x)
colnames(meta.data) <- c(Condition_Col, Cell_State_Col, Sample_Information_Col)

GENOTYPES <- unique(meta.data[[Condition_Col]][meta.data[[Condition_Col]] != comparison_normal_value])
results_list <- list() # Store results here instead of rbind in a loop

# --- 2. Main Analysis Loop ---
for (GENOTYPE in as.character(GENOTYPES)) {
  message("\n>>> Processing sample category: ", GENOTYPE)
  
  # Get unique states for this specific genotype/control set
  current_states <- unique(meta.data[[Cell_State_Col]])
  
  for (STATE_ID in as.character(current_states)) {
    message("Checking State: ", STATE_ID)

    # Define Column Names for filtering mean calculation
    CONTROL_LABEL  <- paste0(comparison_normal_value, STATE_ID)
    GENOTYPE_LABEL <- paste0(GENOTYPE, "_", STATE_ID)

    # Get Cell Type Mapping
    CELLTYPE <- CELLTYPE_STATE[[Cell_Type_Col]][CELLTYPE_STATE[[Cell_State_Col]] == STATE_ID][1]
    if (is.na(CELLTYPE)) {
      message("Skipping: No matching cell type for ", STATE_ID)
      next
    }

    # Subset Data for SR and current GENOTYPE
    samples_to_keep_mask <- meta.data[[Condition_Col]] %in% c(comparison_normal_value, GENOTYPE) & meta.data[[Cell_State_Col]] == STATE_ID
    x_sub <- x[, samples_to_keep_mask, drop = FALSE]
    meta_sub <- meta.data[samples_to_keep_mask, , drop = FALSE]

    # Filter based on PATIENTS_TOKEEP (cell_state_number.csv)
    state_filter_row <- CELLTYPE_FILTER[CELLTYPE_FILTER[[1]] == STATE_ID, ]
    if (nrow(state_filter_row) == 0) next

    patients_log <- unlist(state_filter_row[1, -1], use.names = TRUE)  # ohne erste Spalte (STATE_ID)
    PATIENTS_TOKEEP <- names(which(patients_log))
    
    
    final_sample_mask <- meta_sub[[Sample_Information_Col]] %in% PATIENTS_TOKEEP
    meta_sub[[Sample_Information_Col]] <- sub("\\..*$", "", meta_sub[[Sample_Information_Col]])
    PATIENTS_TOKEEP <- sub("\\..*$", "", PATIENTS_TOKEEP)
    
    x_sub <- x_sub[, final_sample_mask, drop = FALSE]
    meta_sub <- meta_sub[final_sample_mask, , drop = FALSE]
   

    # Validate Sample Size (Need 2 groups and at least 3 reps per group)
    cond_counts <- table(as.character(meta_sub[[Condition_Col]]))
    if (length(cond_counts) == 2 && all(cond_counts >= 3)) {
      
      message("Running edgeR for: ", STATE_ID)
      
      # Setup edgeR Object
      meta_sub[[Condition_Col]] <- factor(meta_sub[[Condition_Col]], levels = c(comparison_normal_value, GENOTYPE))
      dge <- DGEList(counts = x_sub, group = meta_sub[[Condition_Col]])
      
      # Advanced Filtering using genes_tofilter
      # Match columns using regex: Condition__State__Patient
      col_pattern_genotype <- paste0("^", GENOTYPE, "__", STATE_ID, "__")
      col_pattern_control  <- paste0("^", comparison_normal_value, "__", STATE_ID, "__")
      
      cols_genotype <- grep(col_pattern_genotype, colnames(genes_tofilter))
      cols_control  <- grep(col_pattern_control, colnames(genes_tofilter))
      
      # Only keep columns belonging to the PATIENTS_TOKEEP list
      cols_genotype <- cols_genotype[gsub(".*__.*__", "", colnames(genes_tofilter)[cols_genotype]) %in% PATIENTS_TOKEEP]
      cols_control  <- cols_control[gsub(".*__.*__", "", colnames(genes_tofilter)[cols_control]) %in% PATIENTS_TOKEEP]

      # Calculate means for filtering
      mean_genotype <- rowMeans(genes_tofilter[, cols_genotype, drop = FALSE])
      mean_control  <- rowMeans(genes_tofilter[, cols_control, drop = FALSE])
      
      # Expression threshold logic
      keep_genes <- which(mean_control > threshold | mean_genotype > threshold)     # threshold anpassen?
      
      if (length(keep_genes) == 0) {
        message("No genes passed threshold.")
        next
      }

      # edgeR Statistics
      design <- model.matrix( as.formula(paste0("~ ", Condition_Col)),
                              data = meta_sub)
      dge <- estimateDisp(dge, design)
      fit <- glmQLFit(dge, design)
      qlf <- glmQLFTest(fit)
   
      # Extract Results
      tt_all <- topTags(qlf, n = Inf)$table
      tt_all$Gene <- rownames(tt_all)
      
      # Filtered tags (for specific FDR calculation if needed)
      tt_filt <- topTags(qlf[keep_genes, ], n = Inf)$table
      tt_filt$Gene <- rownames(tt_filt)
      colnames(tt_filt)[colnames(tt_filt) == "FDR"] <- "FDR_onlyhigh"
        
      mean_df <- data.frame(
        Gene = rownames(genes_tofilter),
        mean_exp_control  = mean_control,
        mean_exp_genotype = mean_genotype,
        stringsAsFactors = FALSE
      )

      # Merge Stats
      tt_merged <- tt_all %>%
        left_join(tt_filt[, c("Gene", "FDR_onlyhigh")], by = "Gene") %>%
        left_join(mean_df, by = "Gene") %>%
        mutate(
          low_expression = ifelse(is.na(FDR_onlyhigh), "T", "F"),
          FDR_onlyhigh   = ifelse(is.na(FDR_onlyhigh), 1, FDR_onlyhigh),
          FDR_plot       = -log10(FDR_onlyhigh),
          Region            = REGION,
          annotation_level  = CELL_LEVEL,
          cell_type         = STATE_ID,
          cell_heart        = CELLTYPE,
          comparison        = paste0(comparison_normal_value,"_", GENOTYPE),
          Observations_genotype = as.numeric(cond_counts[GENOTYPE]),
          Observations_reference = as.numeric(cond_counts[[comparison_normal_value]])
        )

      results_list[[paste(GENOTYPE, STATE_ID, sep="_")]] <- tt_merged
      message("Finished: ", STATE_ID)
      
    } else {
      message("Skipping ", STATE_ID, ": Insufficient replicates.")
    }
  }
}

# --- 3. Final Compilation & Save ---
if (length(results_list) > 0) {
  final_df <- do.call(rbind, results_list)
  rownames(final_df) <- NULL
  
  # Rename for clarity as requested
  colnames(final_df)[colnames(final_df) == "FDR"] <- "FDR_all"

    
    #added these three variables, needed to make canonical table in the end
  final_df$species <- species 
  final_df$paper   <- paper
  final_df$year    <- year
  final_df$test <- "edgeR"

  if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
  write.csv(final_df, file.path(RESULTS_DIR, output_file_name), row.names = FALSE)

  # --- 4. filtering by FDR and saving ---
  filtered_df <- final_df %>% filter(FDR_onlyhigh <= fdr_threshold)
  write.csv(filtered_df, file.path(RESULTS_DIR, output_filtered_file_name), row.names = FALSE)

  message("SUCCESS: Results saved to ", RESULTS_DIR)

} else {
  message("ERROR: No results generated.")
}





