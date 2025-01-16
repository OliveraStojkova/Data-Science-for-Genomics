# The gtex project has made RNA-seq experiments for a large number of human tissues and individuals. 
# The gtex_data.txt file is a summary of these data sets: it is the median TPM expression for each gene in each tissue. 

#  What 10 genes have the highest summed expression across tissues? 
# Load dataset
gtex_data <- read_tsv("data_q1/gtex_data.txt", show_col_types = FALSE) # 44 219 rows (genes) and 57 columns (54 are median TPM expression for each gene in each tissue)

# find the highest 10 expressed genes by sorting the last column (total) in descending order 
highest10_exp <- gtex_data%>%
  arrange(desc(total))%>%
  head(10)

summed_exp_tibble <- highest10_exp%>%
  select(Name, Description, total)

print(summed_exp_tibble)

# Which 10 genes have the highest coefficient of variation?
# CV = SD/mean
# Add columns for sd, mean and CV to the existing gtex_data tibble

highest10_genes_cv <- gtex_data%>%
  rowwise()%>% # ensures that sd, mean and cv are applies to each row individually
  mutate(
    sd_exp = sd(c_across(3:(ncol(gtex_data)-1))), #only include the columns that have tissue expression data
    mean_exp = mean(c_across(3:(ncol(gtex_data)-1))),
    coef_variation_percent = ((sd_exp/mean_exp)*100)
    ) %>%                                         
  ungroup()%>%         # reverts the rowwise function, back to regular tibble 
  dplyr::arrange(desc(coef_variation_percent))%>%
  slice_head(n=10)

cv_exp_tibble <- highest10_genes_cv%>%
  select(Name, Description, coef_variation_percent)

print(cv_exp_tibble)

# Using the gost function in the gprofiler2 package and an appropriate background, what are the 20 most significantly enriched GO terms for the 100 genes 
# that have the highest summed expression across tissues? Use FDR as correction method, and only report terms where FDR<0.05. 
# Show this as a readable bar plot which shows significance, as -log10 (FDR) and GO term names.  

highest100_summed_exp <- gtex_data%>%
  arrange(desc(total))%>%
  head(100)

highest100_genes <- highest100_summed_exp$Description

highest100_gost <- gost(query = highest100_genes,
     organism = "hsapiens", 
     correction_method = "fdr", 
     significant = TRUE,          # keep only significant value 
     domain_scope = "known",      # all known genes 
     sources = "GO",              # only keep GO terms
     user_threshold = 0.05)       # set significance level

most_enriched20 <- highest100_gost$result%>%
  arrange(p_value)%>%
  head(20)

go_terms <- most_enriched20$term_name

plot_most_enriched20 <- most_enriched20 %>%
  mutate(log10_FDR = -log10(p_value)) %>%   # add a column for -log10(FDR)
  arrange(desc(log10_FDR))

ggplot(plot_most_enriched20, 
       aes(x = reorder(term_name, log10_FDR), y = log10_FDR)) + # orders GO terms, based on -log10FDR
  geom_bar(stat = "identity", fill = cividis(1)) +
  coord_flip()+
  labs(title = "20 Most Enriched GO Terms",
    x = "GO Terms",
    y = "-log10(FDR)") +
  theme_bw()

# What are the 20 most significantly enriched GO terms for the 100 genes that have the  highest coefficient of variation? 
# Use FDR as correction method, and only report terms where FDR<0.05. Show this as a readable bar plot which shows significance, as -log10 (FDR) and GO term names.
highest100_varied_exp <- gtex_data%>%
  rowwise()%>% # ensures that sd, mean and cv are applies to each row individually
  mutate(
    sd_exp = sd(c_across(3:(ncol(gtex_data)-1))), #only include the columns that have tissue expression data
    mean_exp = mean(c_across(3:(ncol(gtex_data)-1))),
    coef_variation_percent = ((sd_exp/mean_exp)*100)
    ) %>%                                         
  ungroup()%>%         # reverts the rowwise function, back to regular tibble 
  dplyr::arrange(desc(coef_variation_percent))%>%
  slice_head(n=100)

highest100_varied_genes <- highest100_varied_exp$Description

highest100_varied_gost <- gost(query = highest100_varied_genes,
     organism = "hsapiens", 
     correction_method = "fdr",           
     domain_scope = "annotated",     # all annotated genes 
     sources = "GO",              # only keep GO terms
     significant = TRUE,
     user_threshold = 0.05)       # set significance level

most_enriched20_varied <- highest100_varied_gost$result%>%
  arrange(p_value)%>%
  head(20)

plot_most_enriched20_varied <- most_enriched20_varied %>%
  mutate(log10_FDR = -log10(p_value)) %>%   # add a column for -log10(FDR)
  arrange(desc(log10_FDR))

ggplot(plot_most_enriched20_varied, 
       aes(x = reorder(term_name, log10_FDR), y = log10_FDR)) + # orders GO terms, based on -log10FDR
  geom_bar(stat = "identity", fill = cividis(1)) +
  coord_flip()+
  labs(title = "20 Most Enriched GO Terms",
    x = "GO Terms",
    y = "-log10(FDR)") +
  theme_bw() 
