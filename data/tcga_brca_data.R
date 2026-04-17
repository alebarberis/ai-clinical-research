# Script used to create the dataset for the laboratory on Clustering.

# 1. Setup ---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("TCGAbiolinks", "SummarizedExperiment", "dplyr", "matrixStats")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p)
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(matrixStats)
library(stringr)

# 2. Query RNA-seq data --------------------------------------------------------
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query = query, directory = file.path("data", "GDCdata"))

rnaseq <- GDCprepare(
  query = query, 
  directory = file.path("data", "GDCdata")
)
# saveRDS(object = rnaseq, file = file.path("data", "tcga_brca_largeRangedSE.rds"))

# 3. Select PRIMARY TUMOURS ----------------------------------------------------
sample_types <- colData(rnaseq)$shortLetterCode

# "TP" = Primary Tumor
keep_samples <- sample_types == "TP"
rnaseq <- rnaseq[, keep_samples]

# 3. Expression (TPM) ----------------------------------------------------------
expr <- assay(rnaseq, "tpm_unstrand")
expr <- log2(expr + 1)

# 4. Gene IDs: ENSEMBL -> clean + HGNC -----------------------------------------

# Get gene info
gene_info <- rowData(rnaseq) %>% as.data.frame()

# Keep protein coding
keep_genes <- gene_info$gene_type == "protein_coding"

expr <- expr[keep_genes, ]
gene_info <- gene_info[keep_genes, ]

# Remove ENSG versions (es: ENSG000001.1 -> ENSG000001)
gene_info$ensembl_id <- str_replace(gene_info$gene_id, "\\..*", "")

# Gene symbol
gene_info$symbol <- gene_info$gene_name

# Fallback: if symbol is missing, use ensembl_id
gene_info$symbol[is.na(gene_info$symbol) | gene_info$symbol == ""] <- 
  gene_info$ensembl_id[is.na(gene_info$symbol) | gene_info$symbol == ""]

# 5. Collapse duplicated (same gene symbol) ------------------------------------
expr_df <- as.data.frame(expr)
expr_df$symbol <- gene_info$symbol

expr_collapsed <- expr_df %>%
  group_by(symbol) %>%
  summarise(across(where(is.numeric), mean))

expr_collapsed <- as.data.frame(expr_collapsed)
rownames(expr_collapsed) <- expr_collapsed$symbol
expr_collapsed$symbol <- NULL

expr <- as.matrix(expr_collapsed)

# 6. Filtering low expression --------------------------------------------------
keep <- rowMeans(expr) > 1
expr <- expr[keep, ]

# 7. Variable genes ------------------------------------------------------------
vars <- rowVars(expr)
top_genes <- order(vars, decreasing = TRUE)[1:4000]
expr <- expr[top_genes, ]

# 8. Trasposizione (samples x genes) -------------------------------------------
expr_final <- as.data.frame(t(expr))
expr_final$sample_id <- rownames(expr_final)

# 9. Clinical data -------------------------------------------------------------
clinical <- GDCquery_clinic("TCGA-BRCA", "clinical")

clinical$patient_id <- substr(clinical$submitter_id, 1, 12)

# Mapping sample → patient
samples <- colnames(expr)
patients <- substr(samples, 1, 12)

sample_map <- data.frame(
  sample_id = samples,
  patient_id = patients
)

metadata <- sample_map %>%
  left_join(clinical, by = "patient_id")

# 10. Order --------------------------------------------------------------------
metadata <- metadata %>%
  filter(sample_id %in% expr_final$sample_id)

expr_final <- expr_final %>%
  filter(sample_id %in% metadata$sample_id)

# Same order
metadata <- metadata[match(expr_final$sample_id, metadata$sample_id), ]

# Select variables
metadata_sel <- metadata %>%
  dplyr::select(
    sample_id,
    patient_id,
    gender,
    age_at_diagnosis,
    stage = ajcc_pathologic_stage,
    vital_status,
    days_to_death,
    days_to_last_follow_up,
    race
  )

metadata_sel <- metadata_sel %>%
  mutate(age_years = round(age_at_diagnosis / 365.25, 1))

coldata <- colData(rnaseq) %>% as.data.frame()

clinical_sel <- coldata %>%
  dplyr::select(
    sample_id = barcode,
    subtype = paper_BRCA_Subtype_PAM50
  )

sample_metadata <- metadata_sel %>%
  left_join(clinical_sel, by = "sample_id")

# Same order
sample_metadata <- sample_metadata[match(expr_final$sample_id, sample_metadata$sample_id), ]

library(RTCGA.clinical)

data(BRCA.clinical)

rtcga <- BRCA.clinical
rtcga <- rtcga %>% select(
  patient.bcr_patient_barcode, 
  patient.breast_carcinoma_estrogen_receptor_status, 
  patient.breast_carcinoma_progesterone_receptor_status
)
rtcga <- rtcga %>% rename(
  patient_id = patient.bcr_patient_barcode,
  ER_status = patient.breast_carcinoma_estrogen_receptor_status,
  PR_status = patient.breast_carcinoma_progesterone_receptor_status) %>%
  filter(ER_status %in% c("positive", "negative")) %>%
  mutate(ER_status = factor(ifelse(ER_status=="positive","ER+","ER-"))) %>%
  filter(PR_status %in% c("positive", "negative")) %>%
  mutate(PR_status = factor(ifelse(PR_status=="positive","PR+","PR-")))
rtcga$patient_id <- toupper(rtcga$patient_id)

sample_metadata <- sample_metadata %>%
  left_join(rtcga, by = "patient_id")

# Filter by gender
sample_metadata <- sample_metadata %>%
  filter(gender %in% c("female"))

# Setup survival time
sample_metadata <- sample_metadata %>%
  mutate(
    vital_status = tolower(vital_status),
    
    event = ifelse(vital_status == "dead", 1, 0),
    
    time = ifelse(
      event == 1,
      as.numeric(days_to_death),
      as.numeric(days_to_last_follow_up)
    )
  )

tcga_brca_sample_metadata <- sample_metadata %>%
  filter(!is.na(time) & time > 0) %>%
  filter(!is.na(ER_status)) %>%
  filter(!is.na(PR_status)) %>%
  filter(!is.na(stage))

tcga_brca_log2_tpm <- expr_final %>%
  filter(sample_id %in% tcga_brca_sample_metadata$sample_id)

# Same order
tcga_brca_sample_metadata <- tcga_brca_sample_metadata[match(tcga_brca_log2_tpm$sample_id, tcga_brca_sample_metadata$sample_id), ]
tcga_brca_sample_metadata <- unique(tcga_brca_sample_metadata)

# Check duplicates
if(any(duplicated(tcga_brca_sample_metadata$patient_id))){
  # Keep vial 'A'
  tcga_brca_sample_metadata <- tcga_brca_sample_metadata %>%
    mutate(
      vial = substr(sample_id, 16, 16)  # estrae "01A", "01B"
    ) 
  tmp <- tcga_brca_sample_metadata[which(duplicated(tcga_brca_sample_metadata$patient_id)),]
  rmid <- tmp[which(tmp$vial != 'A'),"sample_id"]
  tcga_brca_sample_metadata <- tcga_brca_sample_metadata %>%
    filter(!(sample_id %in% rmid))
  
  # Keep first
  tcga_brca_sample_metadata <- tcga_brca_sample_metadata %>%
    group_by(patient_id) %>%
    slice(1) %>%
    ungroup()
  
  tcga_brca_sample_metadata$vial <- NULL
}

tcga_brca_log2_tpm <- expr_final %>%
  filter(sample_id %in% tcga_brca_sample_metadata$sample_id)

# Same order
tcga_brca_sample_metadata <- tcga_brca_sample_metadata[match(tcga_brca_log2_tpm$sample_id, tcga_brca_sample_metadata$sample_id), ]

# Remove sample id col
tcga_brca_log2_tpm$sample_id <- NULL
tcga_brca_log2_tpm <- as.matrix(tcga_brca_log2_tpm)

# Assign names
rownames(tcga_brca_log2_tpm) <- tcga_brca_sample_metadata$patient_id
rownames(tcga_brca_sample_metadata) <- tcga_brca_sample_metadata$patient_id

# Remove cols
# tcga_brca_sample_metadata$sample_id <- NULL
# tcga_brca_sample_metadata$race <- NULL
# tcga_brca_sample_metadata$days_to_death <- NULL
# tcga_brca_sample_metadata$days_to_last_follow_up <- NULL

# 10. Reduce size --------------------------------------------------------------
# The data above is around 14MB. It is better to reduce the size so we can store
# it on GitHub. For example, we can round the numeric values to 5 digits.

# 11. Round the values ---------------------------------------------------------
# Rounding
tcga_brca_log2_tpm <- round(x = tcga_brca_log2_tpm, 5)

# 12. Save ---------------------------------------------------------------------
# write.csv(tcga_brca_log2_tpm, "expression_matrix.csv", row.names = FALSE)
# write.csv(tcga_brca_sample_metadata, "sample_metadata.csv", row.names = FALSE)

save(tcga_brca_log2_tpm, tcga_brca_sample_metadata, file = file.path("data","tcga_brca_mini.rda"), compress = "xz")

# saveRDS(
#   list(expr = tcga_brca_log2_tpm, meta = tcga_brca_sample_metadata),
#   file = file.path("data","tcga_brca_mini.rds"),
#   compress = "xz"
# )
# library(arrow)
# write_parquet(x = as.data.frame(tcga_brca_log2_tpm), sink = file.path("data","tcga_brca_mini.parquet"))
cat("Done!\n")
cat("Expression matrix:", dim(tcga_brca_log2_tpm), "\n")
cat("Metadata:", dim(tcga_brca_sample_metadata), "\n")
