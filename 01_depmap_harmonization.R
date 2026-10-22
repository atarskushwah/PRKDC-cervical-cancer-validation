# scripts/01_depmap_harmonization.R
library(tidyverse)
library(data.table)

# 1. Indestructible Base R Loader
load_safe <- function(path, name) {
  df <- fread(path, data.table = FALSE, check.names = TRUE)
  
  # Safely assign the primary ID without creating duplicates
  if ("DepMap_ID" %in% names(df)) {
    names(df)[names(df) == "DepMap_ID"] <- "ModelID"
  } else if (!("ModelID" %in% names(df))) {
    names(df)[1] <- "ModelID"
  }
  
  # If any duplicate 'ModelID' columns were generated, destroy them
  id_cols <- which(names(df) == "ModelID")
  if (length(id_cols) > 1) {
    df <- df[, -id_cols[-1]]
  }
  
  cat(sprintf("[%s] Loaded %d cell lines\n", name, nrow(df)))
  return(df)
}

models <- load_safe("data/raw/Model.csv", "Models")

# 2. CRISPR (Outcome)
crispr <- load_safe("data/raw/CRISPRGeneEffect.csv", "CRISPR")
prkdc_idx <- grep("PRKDC|5591", names(crispr), ignore.case = TRUE)[1]

crispr_out <- crispr[, c(1, prkdc_idx)]
names(crispr_out)[2] <- "PRKDC_effect"

# Strip ALL hidden characters/spaces before forcing to numeric
crispr_out$PRKDC_effect <- as.numeric(gsub("[^0-9.-]", "", as.character(crispr_out$PRKDC_effect)))
crispr_out <- crispr_out[!is.na(crispr_out$PRKDC_effect), ]
cat("[CRISPR] Valid PRKDC scores retained:", nrow(crispr_out), "\n")

# 3. Expression 
expr <- load_safe("data/raw/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", "Expression")
prkdc_exp_idx <- grep("PRKDC|5591", names(expr), ignore.case = TRUE)[1]
names(expr)[prkdc_exp_idx] <- "PRKDC_expression"

# 4. Copy Number
cn <- load_safe("data/raw/OmicsCNGeneWGS.csv", "Copy Number")
prkdc_cn_idx <- grep("PRKDC|5591", names(cn), ignore.case = TRUE)[1]
names(cn)[prkdc_cn_idx] <- "PRKDC_CN"

# 5. Mutations
mut_raw <- load_safe("data/raw/OmicsSomaticMutations.csv", "Mutations")
var_col <- grep("Variant|Classification|Consequence", names(mut_raw), ignore.case = TRUE, value = TRUE)[1]
gene_col <- grep("Hugo|Gene", names(mut_raw), ignore.case = TRUE, value = TRUE)[1]

mut_sub <- mut_raw[, c("ModelID", gene_col, var_col)]
names(mut_sub) <- c("ModelID", "HugoSymbol", "VariantInfo")

is_damaging <- grepl("frame|nonsense|splice|stop", mut_sub$VariantInfo, ignore.case = TRUE)
mut_sub <- mut_sub[is_damaging, ]

if (nrow(mut_sub) > 0) {
  mut_sub$Mutated <- 1
  mut_sub <- mut_sub[!duplicated(paste(mut_sub$ModelID, mut_sub$HugoSymbol)), ]
  mutations <- pivot_wider(mut_sub[, c("ModelID", "HugoSymbol", "Mutated")],
                           names_from = HugoSymbol,
                           values_from = Mutated,
                           values_fill = list(Mutated = 0))
} else {
  mutations <- data.frame(ModelID = character())
}
cat("[Mutations] Damaging mutations parsed:", nrow(mutations), "\n")

# 6. Harmonize
cat("\n--- Starting Harmonization ---\n")
depmap_cohort <- as_tibble(crispr_out) %>%
  left_join(as_tibble(models), by = "ModelID") %>%
  left_join(as_tibble(expr), by = "ModelID") %>%
  left_join(as_tibble(cn), by = "ModelID", suffix = c("_exp", "_cn")) %>%
  left_join(as_tibble(mutations), by = "ModelID", suffix = c("", "_mut"))

cat("Harmonization complete. Final cell lines ready for modeling:", nrow(depmap_cohort), "\n")
write_csv(depmap_cohort, "data/processed/01_harmonized_depmap_discovery.csv")
cat("Saved to data/processed/01_harmonized_depmap_discovery.csv\n")