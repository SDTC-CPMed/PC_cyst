#####This code is from PC cyst project, where I use 132 genes and TCGA DEGs to find the comon genes and make the heatmaps
##written by firoj.mahmud@ki.se
library(pheatmap) # For plotting heatmap
library(readxl)  # For reading Excel files
library(tidyverse)  # For data manipulation and plotting
library(ggplot2)  # For creating box plots


setwd("/data/sharedData/Firoj/BACKup_Documents/EarlyStagePANCREATICcancer/PCcystMANUSCRIPT_FILES_UPDATED_29jan_2024/Validation_PROTEINdata_NatureCancer")
      
validation_gwas_genes <- c("ACTB","ACTR3","ACTG1","ARPC1B","ARPC2","ARPC3","ARPC5","ARPC4")
file_proteomics_valid <- "/data/sharedData/Firoj/BACKup_Documents/EarlyStagePANCREATICcancer/PCcystMANUSCRIPT_FILES_UPDATED_29jan_2024/Validation_PROTEINdata_NatureCancer/20221227_PDAC_PRO_exp.xlsx"
data_prt_valid <- read_excel(file_proteomics_valid)
#####NOW I have to search the proteins to be validated in the data sets#############################

data_long <- data_prt_valid %>%
  pivot_longer(
    cols = -...1, 
    names_to = "sample", 
    values_to = "expression"
  ) %>%
  mutate(
    protein = ...1,
    condition = if_else(str_detect(sample, "_N_PRO"), "Normal", "Tumor"),
    sample_id = str_remove(sample, "_[NT]_PRO")
  ) %>%
  select(-...1, -sample) %>%
  filter(protein %in% c("ACTB","ACTR3","ACTG1","ARPC1B","ARPC2","ARPC3","ARPC5","ARPC4")) ###this are the validated proteins

ggplot(data_long, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot() +
  facet_wrap(~protein, scales = "free_y") +
  theme_minimal() +
  labs(title = "Protein Expression in Normal vs. Tumor Samples", y = "Expression Level") +
  scale_fill_brewer(palette = "Set1")
dev.off()
###################################################calculate the significance and print in the box plots#
results <- data_long %>%
  group_by(protein) %>%
  summarise(p_value = t.test(expression[condition == "Normal"], expression[condition == "Tumor"], paired = FALSE)$p.value)

print(results)


significance_stars <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns") # Not significant
  }
}

results$significance <- sapply(results$p_value, significance_stars)
data_long <- left_join(data_long, results, by = "protein")
head(data_long)

 ggplot(data_long, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot() +
  facet_wrap(~protein, scales = "free_y") +
  theme_minimal() +
  labs(title = "Protein Expression in Normal vs. Tumor Samples", y = "Expression Level") +
  scale_fill_brewer(palette = "Set1") +
  geom_text(data = distinct(data_long, protein, .keep_all = TRUE),
            aes(label = significance, x = 1.5, y = Inf), 
            position = position_nudge(y = -0.05), hjust = -0.1, vjust = 1, check_overlap = TRUE)


##############Make a bubble plot with significance and fold change


################################### TCGA AND CPTAC DATA###########

TCGA_CPTAC=read.csv("/data/sharedData/Firoj/BACKup_Documents/EarlyStagePANCREATICcancer/NEW_clean_RESULT_MANUSCRIPT_PCcyst/TCGA_CPTAC_Concordants_filtered.csv")


# Filter out specific genes
genes_of_interest <- c("ACTB","ACTR3","ACTG1","ARPC1B","ARPC2","ARPC3","ARPC5","ARPC4")
filtered_data <- TCGA_CPTAC %>% 
  filter(gene %in% genes_of_interest)
filtered_data
# Divide the data into two subsets for CPTAC and TCGA
cptac_data <- filtered_data %>%
  select(gene, logFC = logFC_CPTAC, pValue = pValue_CPTAC)
head(cptac_data)
tcga_data <- filtered_data %>%
  select(gene, logFC = logFC_TCGA, pValue = Padj_TCGA)
head(tcga_data)
#############################################PLOT BASED ON PVALUE AS BUBBLE SIZE#########
significance_color <- function(p) {
  ifelse(p < 0.001, "High", ifelse(p < 0.01, "Medium", ifelse(p < 0.05, "Low", "Not Significant")))
}

# Create bubble plots for CPTAC
P2=ggplot(cptac_data, aes(x = reorder(gene, logFC), y = logFC, size = -log10(pValue), color = significance_color(pValue))) +
  geom_point(alpha = 0.6) +
  scale_size_continuous(range = c(3, 12), name = "-log10(p-value)") +
  scale_color_manual(values = c("High" = "red", "Medium" = "orange", "Low" = "blue", "Not Significant" = "grey")) +
  labs(title = "CPTAC Dataset: Protein Expression Fold Change",
       x = "Protein",
       y = "Log Fold Change") +
  theme_minimal() +
  theme(legend.position = "right")

# Create bubble plots for TCGA 
P1=ggplot(tcga_data, aes(x = reorder(gene, logFC), y = logFC, size = -log10(pValue), color = significance_color(pValue))) +
  geom_point(alpha = 0.6) +
  scale_size_continuous(range = c(3, 12), name = "-log10(p-value)") +
  scale_color_manual(values = c("High" = "red", "Medium" = "orange", "Low" = "blue", "Not Significant" = "grey")) +
  labs(title = "TCGA Dataset: Gene Expression Fold Change",
       x = "Gene",
       y = "Log Fold Change") +
  theme_minimal() +
  theme(legend.position = "right")

