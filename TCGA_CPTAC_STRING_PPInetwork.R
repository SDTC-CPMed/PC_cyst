
######################################This code will take Differential expressed genes and differential expressed protein list from TCGA and CPTAC find the common regulatory pattern and based on the selected gene list further create a protein-protein interection network based on STRING db and further find the sub cluster with GWAS enriched genes
##Codes written by firoj mahmud, post doctoral researcher, medical digital twin group, karolinska institute, for further email to firoj.mahmud@ki.se

  library(dplyr)
  library(ggplot2)
  library(stats)
  library(ggraph)
  library(igraph)
  library(patchwork)
  library(graphlayouts)
  library(ggforce)
  library(STRINGdb)
  

  setwd("/home/firoj/Documents/EarlyStagePANCREATICcancer/TCGA_CPTAC_STRINGs")
  
  
  data_CPTAC<- read.table("/home/firoj/Documents/EarlyStagePANCREATICcancer/TCGAdataSets_DEseq/all-DE.txt", header=TRUE, sep="\t")
  data_TCGA <- read.table("/home/firoj/Documents/EarlyStagePANCREATICcancer/TCGAdataSets_DEseq/DEG_PAAD_LIMMA.txt", header=TRUE, sep="\t")
  
  head(data_CPTAC)
  head(data_TCGA)

  
  ###Now i want to find similar up or down regulated genes in both data sets 
  # filtering based on the absolute value of log fold change > 0.25
  filtered_CPTAC <- subset(data_CPTAC, abs(logFC_TvsN) > 0.25)
  filtered_TCGA <- subset(data_TCGA, abs(Log2.Fold.Change.) > 0.25)
  
  filtered_TCGA$gene <- filtered_TCGA$Gene.Symbol
  
  merged_data <- merge(filtered_CPTAC, filtered_TCGA, by = "gene")
  write.csv(merged_data, file="/home/firoj/Documents/EarlyStagePANCREATICcancer/TCGA_CPTAC_STRINGs/TCGA_CPTAC_Concordants.csv")
  
  same_sign_genes <- merged_data[(merged_data$logFC_TvsN * merged_data$Log2.Fold.Change.) > 0, ]
  final_data <- same_sign_genes[, c("gene", "logFC_TvsN", "pValue", "p.adj", "Log2.Fold.Change.", "adjp")]
 
  colnames(final_data) <- c("gene", "logFC_CPTAC", "pValue_CPTAC", "Padj_CPTAC", "logFC_TCGA", "Padj_TCGA") ##changed the colum as accordingly
  head(final_data)
  write.csv(final_data, file="/home/firoj/Documents/EarlyStagePANCREATICcancer/TCGA_CPTAC_STRINGs/TCGA_CPTAC_Concordants_filtered.csv")
  
  
  # Create a scatter plot with logFC_CPTAC on the x-axis and logFC_TCGA on the y-axis
  ggplot(final_data, aes(x = logFC_CPTAC, y = logFC_TCGA)) +
    geom_point() +  # Add points
    theme_minimal() +  # Use a minimal theme
    labs(x = "Log Fold Change in CPTAC", 
         y = "Log Fold Change in TCGA", 
         title = "Comparison of Log Fold Changes between CPTAC and TCGA") +
    geom_smooth(method = "lm", color = "blue", se = FALSE) + # Add a linear regression line without the shaded confidence interval
    theme(plot.title = element_text(hjust = 0.5)) # Center the plot title
  
  
  
  GWASgenes=read.table("/home/firoj/Documents/EarlyStagePANCREATICcancer/BOXplots_PCcystUKBIOBANK/genes_top_p_cyst_vs_HC.txt",header=TRUE, sep="\t")
  head(GWASgenes)
  GWAS_id=GWASgenes$V1
  gene_names=GWASgenes$V1
  
  
  # Check if any GWAS genes are found in the final_data dataset
  matching_genes <- final_data$gene[final_data$gene %in% GWASgenes$gene]
  


   # Subset for GWAS genes
  gw_genes <- subset(final_data, gene %in% GWASgenes$gene)

  final_data$highlight <- ifelse(final_data$gene %in% GWASgenes$gene, "GWAS Gene", "Other Gene")
  #pdf("scatter_plot_TCGA-CPTAC.pdf", width = 11, height = 8.5)
  #tiff("scatter_plot_TCGA-CPTAC.tiff", width = 11, height = 8.5, units = 'in', res = 300)
  ggplot(final_data, aes(x = logFC_CPTAC, y = logFC_TCGA, color = highlight, size = highlight)) +
    geom_point(alpha = 0.7) +  # Set transparency to see overlapping points
    scale_size_manual(values = c("GWAS Gene" = 3, "Other Gene" = 1)) + # Adjust point sizes
    scale_color_manual(values = c("GWAS Gene" = "red", "Other Gene" = "black")) + # Adjust point colors
    theme_minimal() +
    labs(x = "Log Fold Change in CPTAC", 
         y = "Log Fold Change in TCGA", 
         title = "Comparison of Log Fold Changes between CPTAC and TCGA",
         color = "Gene Type",
         size = "Gene Type") +
    geom_smooth(method = "lm", color = "blue", se = FALSE, linetype = "dashed") +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color = guide_legend(override.aes = list(size = 6)), # Adjust legend point size
           size = guide_legend(override.aes = list(size = 6)))
  
########################################################### Now create the PPI network based on the saved protein-gene list that are concordantly expressed
  
 
  
  string_db <- STRINGdb$new(version="11.0", species=9606, score_threshold=400)
  
  mapped_genes <- string_db$map(final_data, "gene", removeUnmappedRows = TRUE)
  string_interaction <- string_db$get_interactions(mapped_genes$STRING_id)
  
  write.table(interactions, "/data/sharedData/Firoj/BACKup_Documents/GWAS_gene_STRING/PPI_TCGA_CPTAC_data.tsv", sep="\t", row.names=FALSE, quote=FALSE)
  
  string_interaction= read.delim("/data/sharedData/Firoj/BACKup_Documents/GWAS_gene_STRING/PPI_TCGA_CPTAC_data.tsv")
  #Here I set up the string interection which is more than 0.95
  filtered_interaction <- string_interaction %>%
    filter(combined_score > 0.95)###string interection cut-off set to 0.95

  #Here I set up the string interection which is more than 0.95
  filtered_interaction <- string_interaction %>%
    filter(combined_score > 0.6)###string interection cut-off set to 0.95

  
  GWASgenes=read.table("/home/firoj/Documents/GWAS_gene_STRING/genes_1330_dermatologic_diseases_send.txt",header=TRUE, sep="\t")

  GWASgenes_vector <- as.vector(GWASgenes$gene)  # Adjust the column name if necessary
  
  g <- graph_from_data_frame(d = filtered_interaction[, c("node1", "node2")], 
                             directed = FALSE)
  E(g)$combined_score <- filtered_interaction$combined_score
  E(g)$coexpression <- filtered_interaction$coexpression
  node_degree <- degree(g)
  threshold <- quantile(node_degree, 0.9)
  hub_nodes <- names(which(node_degree > threshold))
  
  
 
  
  ##############isolate subcluster 
  
  centrality_scores <- degree(g, mode = "all")
  #centrality_scores <- degree(g, mode = "all")
  V(g)$centrality <- centrality_scores
  # Identify top central genes
  top_central_genes <- names(sort(centrality_scores, decreasing = TRUE)[1:50])
  
  # This includes top genes and any nodes directly connected to them, here i wanted to extract genes that are connected with each other by their first neighbour
  subgraph_nodes <- unique(c(top_central_genes, unlist(lapply(top_central_genes, function(gene) neighbors(g, gene)$name))))
  subgraph_g <- induced_subgraph(g, subgraph_nodes)
  ifelse(names(V(g)) %in% GWASgenes_vector, "green", "blue"))

  
    ggraph(subgraph_g, layout = "fr") + 
    geom_edge_link(color = "grey80", edge_alpha = 0.5) +
    geom_node_point(aes(color = name %in% top_central_genes), size = centrality) +
    scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
    geom_node_text(aes(label = ifelse(name %in% top_central_genes, name, "")), repel = TRUE, size = 3, color = "black", fontface = "bold") +
    theme_graph() +
    labs(title = "Isolated Cluster of Top Central Genes", subtitle = "With Immediate Neighbors")
  

  
  