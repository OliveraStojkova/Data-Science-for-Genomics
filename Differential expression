# Investigate how the drug affects cellular mRNA levels in general, and whether the expression of key groups of genes are affected. 

# 4 files are supplied:
# `studyDesign.tsv`: File describing treatment of the 10 samples included in the study.
# `countMatrix.tsv`: Number of RNA-Seq reads mapping to each of the genes.
# `normalizedMatrix.tsv`: Normalized expression to each of the genes.
# `diabetesGene.tsvs`: Collection of genes known to be involved in type 2 diabetes.

### Load data
# study design
study_design <- read_tsv("data_q3/studyDesign.tsv", show_col_types = FALSE)
samples_study_design <- study_design$Sample

# count matrix
count_m <- read_table("data_q3/countMatrix.tsv", col_names = c("Gene", paste0("Sample", 1:10)), show_col_types = FALSE) # add column names

count_matrix<-count_m[-1, ]%>%                        # remove first row 
  mutate(across(starts_with("Sample"), as.numeric))   # convert the last column from character to numeric
   

# check whether the samples from the study design and count matrix are the same 
samples_count_m <- colnames(count_matrix)[-1]
samples_sd_m <- all(samples_study_design %in% samples_count_m) & all(samples_count_m %in% samples_study_design)

# normalized matrix
normalized_m <- read_table("data_q3/normalizedMatrix.tsv", col_names = c("Gene", paste0("Sample", 1:10)), show_col_types = FALSE)

normalized_matrix<-normalized_m[-1, ]%>%
    mutate(across(starts_with("Sample"), as.numeric))

# check whether the samples from the study design and normalized matrix are the same 
samples_count_nm <- colnames(normalized_matrix)[-1]
samples_sd_nm <- all(samples_study_design %in% samples_count_nm) & all(samples_count_nm %in% samples_study_design)

## perform a hierarchical clustering on the samples using euclidian distance and average linkage, and plot the resulting tree
normalized_matrix_no_gene <- normalized_matrix %>% select(-Gene) # remove the gene column

# calculate distance matrix from the normalized matrix
normalized_matrix_t <- as.matrix(t(normalized_matrix_no_gene))
distance_matrix <- dist(normalized_matrix_t)

# perform hierarchical clustering with hclust, method = "avergae"
hc_clustering <- hclust(distance_matrix, method = "average")

# plot
plot(hc_clustering, main="Hierarchical Clustering of Samples", 
     xlab = "Samples",
     ylab = "Distance")

### Differential expression
# remove the gene column of the outlier (Sample 5)
count_matrix_no_outliers <- count_matrix %>%
  select(-c(Gene, Sample5))%>%
  as.matrix()

row.names(count_matrix_no_outliers) <- count_matrix$Gene

# remove sample5 column
study_design_1 <- study_design%>%
  filter(!(Sample %in% c("Sample5")))

# save data as DESeqDatSet-object
dds <- DESeqDataSetFromMatrix(countData = count_matrix_no_outliers, # Count matrix without outliers
  colData = study_design_1, #  design matrix: study design
  design = ~ Condition    #  column that 'splits' the data - Condition: Ctrl and Trt
) 

# perform DESeq2
dds <- DESeq(dds)

# save results
results_dds <- results(dds, 
contrast=c("Condition", "Trt", "Ctrl"), # Comparison (LFC > 0.25, means higher expression in treatment)
lfcThreshold=0.25,                      # adjusted logFC cutoff
alpha=0.05)                             # adjusted FDR cutoff)

summary(results_dds)

# save DESeq2 results in a tibble
results_dds_tible <- results(dds, 
                             contrast=c("Condition", "Trt", "Ctrl"), 
                             lfcThreshold=0.25,
                             alpha=0.05,
                             tidy = TRUE) %>% 
  as_tibble%>%
  rename(row = "GeneID")

# MA plot
ggplot(results_dds_tible, 
       aes(x=baseMean, 
           y=log2FoldChange,
           col = padj < 0.05)) + 
  geom_point(size = 0.5, alpha=0.7) + 
  scale_x_log10() + 
  geom_hline(yintercept = 0, alpha = 0.75)+
  theme_bw()

# find up and down regulated genes
upregulated_genes <- results_dds_tible[which(results_dds$log2FoldChange > 0.25 & results_dds$padj < 0.05), ]
downregulated_genes <- results_dds_tible[which(results_dds$log2FoldChange < -0.25 & results_dds$padj < 0.05), ]

# count up and down regulated genes
no_upregulated_genes <- nrow(upregulated_genes)
no_downregulated_genes <- nrow(downregulated_genes)
 
# Sort the DE statistics table  that you get from DESeq2 to report the top 5 genes sorted by
# 1. positive logFC (highest on top)
# 2. negative logFC (lowest on top)
# only looking at significantly differentially expressed genes*

# top 5 significant genes with positive logFC
highest5_logFC <- results_dds_tible%>%
  filter(padj < 0.05)%>%
  arrange(desc(log2FoldChange))%>%
  head(5)

# top 5 significant genes with negative logFC
lowest5_logFC <- results_dds_tible%>%
  filter(padj < 0.05)%>%
  arrange(log2FoldChange)%>%
  head(5)

