library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplotify)
library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(Cairo)
library(data.table)
library(tibble)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(lme4)
library(lattice)
library(lme4qtl)
library(stats)
library(lmtest)
library(reshape2)
library(dendextend)
library(ggpubr)
library(ape)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)

load("HGT_pipeline.RData")

#### 1 HGT information (Basic identity/ length/ Rates/ Number/ After filter compare) -- Figure1 ####

#HGT_info <- read.table("HGT_Clusters_info_filter_conseve_Final.csv",sep=",",header=T,row.names = 1)
median(HGT_info$Length)
sd(HGT_info$Length)

median(HGT_info$Identity)
sd(HGT_info$Identity)

## 1.1 Figure 1a Phylogenetic tree of MAGs
MAGs <- c(HGT_info$MAG.1,HGT_info$MAG.2)
MAGs <- MAGs %>% unique()
#node_info <- read.csv("HGT_network/HGT_node_info_full.csv",header=T)
MAG_hgt <- data.frame(node = MAGs, hgt = "1")
node_info_hgt <- left_join(node_info, MAG_hgt, by = c("genome" = "node"))
node_info_hgt$hgt[is.na(node_info_hgt$hgt)] <- 0
#tree_mags <- read.tree("Tree_HGT_MAGs/concatenated_gsi(81).nwk")
tree_labels <- tree_mags$tip.label
labels_to_keep <- tree_labels[tree_labels %in% node_info_hgt$genome1]
filtered_tree_mags <- drop.tip(tree_mags, setdiff(tree_labels, labels_to_keep))

phylum_genomes_list <- node_info_hgt %>%
  group_by(Phylum) %>%
  summarise(genomes = list(genome1)) %>%
  deframe()

tree_mags_ano <- groupOTU(filtered_tree_mags,phylum_genomes_list)

ALl_tree_plot <- ggtree(tree_mags_ano,aes(color=group),
                        layout="fan",
                        linewidth=0.6,
                        show.legend=F) +
  geom_point(mapping=aes(color=group),
             size=1.5,
             stroke=0.5,
             alpha=0.5) +
  scale_color_manual(values=c("black","#E89DA0","#88CEE6","#017D6F","#0AAF88","#8AD088","#B686B6","#FDDB76"))+
  new_scale_fill() +
  geom_fruit(
    data=node_info_hgt,
    geom=geom_tile,
    mapping=aes(y=genome1, fill=hgt),
    color=NA,
    width=0.1,
    offset=0.1)+
  scale_fill_manual(
    values=c("#5398CD", "white"),
    guide=guide_legend(keywidth=1, keyheight=1, order=2),
    name="HGT")+
  theme(legend.position = "none")

## 1.2 Figure 1b genome quality

#genome_info <- read.table("dRep999/Widb.csv",sep = ",",header = T)
mean(genome_info$completeness)
mean(genome_info$contamination)

sd(genome_info$completeness)
sd(genome_info$contamination)

plot_comple<-ggplot(genome_info,aes(completeness))+
  geom_histogram(binwidth = 0.5,fill = "#0AAF88",alpha = 0.5)+
  scale_y_reverse()+
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=12))

plot_conta<-ggplot(genome_info,aes(contamination))+geom_histogram(binwidth = 0.2,fill="#0AAF88",alpha = 0.5)+
  coord_flip()+
  scale_y_reverse()+
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=12))

plot_score<-ggplot(genome_info,aes(score))+geom_histogram(binwidth = 0.5)

cluster_group <- gsub("(_).*", "", genome_info$cluster)

genome_group <- data.frame( genome = genome_info$genome, cluster = cluster_group, Value = rep(1, length(cluster_group)))

genome_group1 <- aggregate(Value ~ cluster, data = genome_group, sum)

sorted_data <- genome_group1 %>%
  arrange(desc(Value))

plot_scatter<-ggplot(genome_info, aes(x = completeness, y = contamination)) +
  geom_point(alpha = 0.5,color = "#0AAF88") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=12),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.ticks = element_line(color = "grey", size = 1))  

## 1.3 Figure 1c HGT information ##
# > 500bp
#HGT_info <- read.table("HGT_events/HGT_Clusters_info_filter_conseve_Final.csv",sep=",",header=T,row.names = 1)

HGT_info1 <- HGT_info %>%
  mutate(SNP = Identity)

plot_scatter<-ggplot(HGT_info1, aes(x = SNP, y = Length,color = Identity)) +
  geom_point(alpha = 0.5) + 
  scale_color_gradientn(colors = c("#d8b365", "grey80","#41B6C4"), 
                        values = scales::rescale(c(90, 100)),
                        guide = "none") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=12),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.ticks = element_line(color = "grey", size = 1))  
plot_scatter

plot_Identity<-ggplot(HGT_info1,aes(SNP))+
  geom_histogram(binwidth = 0.01,fill = "#41B6C4", alpha = 0.6)+
  scale_y_reverse()+
  theme_classic() +
  xlab("Identity") +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=12)
  )

plot_Length<-ggplot(HGT_info1,aes(Length))+
  geom_histogram(binwidth = 500,fill="#41B6C4", alpha = 0.6)+
  coord_flip()+
  theme_classic() +
  ylab("SNPs") +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=12))

#### 2 Transfer between species/ Genus / Family & ARG/ BGC tranfer pattern -- Figure2 ####
#data <- read.table("HGT_events/HGT_Clusters_info_filter_conseve_Final.csv", sep=",", header=T,row.names =1)

HGT_cluster <- data.frame(Sequence = data$Sequence_ID,Value=data$Number.of.sequences.in.cluster, Cluster = data$cluster)

HGT_cluster_uniq <- unique(HGT_cluster[, c("Cluster", "Value")])

HGT_cluster_uniq <- HGT_cluster_uniq[order(-HGT_cluster_uniq$Value), ]

HGT_cluster_uniq_filter <- HGT_cluster_uniq[HGT_cluster_uniq$Value >= 20 ,]

table(HGT_cluster_uniq_filter$Value)

#cluster_number_species <- read.table("cd-hit/cluster_number.csv", sep=",", header=T,row.names =NULL)
cluster_number_species_filter <- cluster_number_species[cluster_number_species$species > 6 ,]
#### Species level
ggplot(cluster_number_species_filter, aes(x = reorder(cluster, -species), y = species, fill = (species > 7))) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # 
  geom_text(aes(label=species), vjust=-0.5, size=3.5) +  # 
  scale_fill_manual(values = c("lightgray", "#7DABCF")) +  # 
  xlab("Spanning Species") +
  ylab("Count") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),  # 
    panel.grid.minor = element_line(color = "gray95")   # 
  )
p2
#### Genus level
#cluster_number <- read.table("cd-hit/cluster_number_genus.csv", sep=",", header=T,row.names =NULL)
cluster_number_filter <- cluster_number[cluster_number$species > 4 ,]
p3 <- ggplot(cluster_number_filter, aes(x = reorder(cluster, -species), y = species, fill = (species > 6))) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # 
  geom_text(aes(label=species), vjust=-0.5, size=3.5) +  # 
  scale_fill_manual(values = c("lightgray", "#FF8C00")) +  # 
  xlab("Spanning Genus") +
  ylab("Count") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),  
    panel.grid.minor = element_line(color = "gray95")   
  )
p3

#### Family level
#cluster_number_Family <- read.table("cd-hit/cluster_number_Family.csv", sep=",", header=T,row.names =NULL)
cluster_number_Family_filter <- cluster_number_Family[cluster_number_Family$species > 2 ,]
p4 <- ggplot(cluster_number_Family_filter, aes(x = reorder(cluster, -species), y = species, fill = (species > 3))) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  
  geom_text(aes(label=species), vjust=-0.5, size=3.5) +  
  scale_fill_manual(values = c("lightgray", "#F5D44B")) +  
  xlab("Spanning Family") +
  ylab("Count") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),  
    panel.grid.minor = element_line(color = "gray95")   
  )
p4

#### 2.2 Figure 2a-b HGT transfer pattern Recombinase

#edges_info <- read.csv("HGT_events/HGT_Clusters_info_filter_conseve_Final.csv")

#node_info <- read.csv("HGT_network/HGT_node_info_full.csv",header=T)

node_info %>%
  select(Genus) %>%
  summarise(unique_phylum_count = n_distinct(Genus))

edges_ano <- data.frame(MAG1 = edges_info$MAG.1 , MAG2 = edges_info$MAG.2 ,iden = edges_info$Identity,
                        Cluster = edges_info$cluster, common_tax = edges_info$common_classification, 
                        number_cluster = edges_info$Number.of.sequences.in.cluster)

all_nodes_in_edges <- unique(c(edges_ano$MAG1, edges_ano$MAG2))

all_nodes_in_info <- unique(node_info$genome)
missing_nodes <- setdiff(all_nodes_in_edges, all_nodes_in_info)

# Filter out rows in edges_ano where either MAG1 or MAG2 is in missing_nodes
edges_ano_cleaned <- edges_ano[!edges_ano$MAG1 %in% missing_nodes & !edges_ano$MAG2 %in% missing_nodes,]
g <- graph_from_data_frame(d=edges_ano_cleaned, vertices=node_info, directed=FALSE)

group_graph_tidy <- as_tbl_graph(g)

group_graph_tidy <- group_graph_tidy %>%
  activate(nodes) %>%
  mutate(node_name = name)


# Function to filter and plot a specific cluster
plot_cluster <- function(group_graph_tidy, cluster_name) {
  
  # Filter the edges
  edges_filtered <- group_graph_tidy %>% 
    activate("edges") %>% 
    mutate(from_name = .N()$name[from],
           to_name = .N()$name[to]) %>%
    as_tibble() %>% 
    filter(Cluster == cluster_name)
  
  # Get unique node indices from the filtered edges
  node_indices <- unique(c(edges_filtered$from_name, edges_filtered$to_name))
  
  # Create a subgraph with the filtered nodes and edges
  graph_filtered <- group_graph_tidy %>%
    activate(edges) %>%
    filter(Cluster == cluster_name) %>%
    activate(nodes) %>%
    filter(name %in% node_indices)
  
  
  
  graph_filtered <- as_tbl_graph(graph_filtered)
  
  # Plot the graph
  g <- ggraph(graph_filtered, layout = 'kk') + 
    geom_edge_link(aes(width = iden,label=iden), alpha = 0.6, color="#66c2a5", show.legend = FALSE) +  # Adjust edge width based on iden
    #geom_edges_text(aes(label=iden_label), size=3, vjust = 0.5, hjust = 0.5, color="black", angle=0, fontface="italic") +  # Add iden value labels
    geom_node_point(aes(fill=factor(Family), size = 10), shape=21, color="white", stroke=0.5) +  # Add node border
    #geom_node_text(aes(label=Species), size=3, fontface="bold", repel = TRUE, color="black") +  # Use bolder font for labels
    scale_edge_width_continuous(range=c(0.5,2)) +  # Add some edge width variability
    scale_fill_brewer(palette = "Set2") +  # Change color palette
    labs(title = paste(cluster_name,"network")) +
    theme_graph(base_family="Arial") +  # Choose a clean font family
    theme(legend.background = element_rect(colour = NA),
          panel.background = element_rect(fill = "white",  colour = NA),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),  # Soft background color
          plot.title = element_text(size = 18, face = "bold")) +  # Bolder title
    guides(size=FALSE)
  
  return(g)
}

p_236 <- plot_cluster(group_graph_tidy, "Cluster 236")

## 2.3 Figure 2c ARGs transfer pattern ##
#ARG_table <- read.csv("gene_sum/HGT_ARGs_clusters_all.raw.txt",sep='\t')
Seq2Cluster_map <- setNames(HGT_cluster$Cluster,HGT_cluster$Sequence)
ARG_table$Sequence <- sub("_\\d+$", "", ARG_table$Contig)
ARG_table$cluster <- Seq2Cluster_map[ARG_table$Sequence]

ARG_cluster_uniq <- unique(ARG_table$cluster)

p_902 <- plot_cluster(group_graph_tidy, "Cluster 902")
p_803 <- plot_cluster(group_graph_tidy, "Cluster 803")
p_47 <- plot_cluster(group_graph_tidy, "Cluster 47")
p_1180 <- plot_cluster(group_graph_tidy, "Cluster 1180")
p_117 <- plot_cluster(group_graph_tidy, "Cluster 117")
p_217 <- plot_cluster(group_graph_tidy, "Cluster 217")
p_633 <- plot_cluster(group_graph_tidy, "Cluster 633")
p_113 <- plot_cluster(group_graph_tidy, "Cluster 113")
p_4 <- plot_cluster(group_graph_tidy, "Cluster 4")
p_446 <- plot_cluster(group_graph_tidy, "Cluster 446")


cluster_name <- ARG_cluster_uniq

edges_filtered <- group_graph_tidy %>% 
  activate("edges") %>% 
  mutate(from_name = .N()$name[from],
         to_name = .N()$name[to]) %>%
  as_tibble() %>% 
  filter(Cluster %in% cluster_name)

# Get unique node indices from the filtered edges
node_indices <- unique(c(edges_filtered$from_name, edges_filtered$to_name))

# Create a subgraph with the filtered nodes and edges
graph_filtered <- group_graph_tidy %>%
  activate(edges) %>%
  filter(Cluster %in% cluster_name) %>%
  activate(nodes) %>%
  filter(name %in% node_indices)



graph_filtered <- as_tbl_graph(graph_filtered)

my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
               "#1a55FF", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", 
               "#666666", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", 
               "#66a61e", "#e6ab02")

my_node_colors <- c("Firmicutes_A" = "#1f77b4","Firmicutes_C"= "#7570b3",
                    "Firmicutes" = "#1b9e77","Bacteroidota"= "#ff7f0e")

my_edge_colors <- c("Cluster 902" = "#17becf", "Cluster 803" = "#a6cee3", "Cluster 47" = "#b2df8a",
                    "Cluster 1180" = "#7f7f7f", "Cluster 117" = "#fb9a99", "Cluster 217" = "#e31a1c",
                    "Cluster 4" = "#fdbf6f", "Cluster 633" = "#ff7f00", "Cluster 113" = "#cab2d6",'Cluster 446' = '#66a61e')


# Plot the graph
g <- ggraph(graph_filtered, layout = 'kk') + 
  geom_edge_link(aes(width = iden,color = factor(Cluster),label=iden), alpha = 0.6, show.legend = T) +  # Adjust edge width based on iden
  geom_node_point(aes(fill=factor(Phylum), size = 15), shape=21, color="white", stroke=0.5) +  # Add node border
  geom_node_text(aes(label=Genus), size=4, fontface="bold", repel = TRUE, color="black") +  # Use bolder font for labels
  scale_edge_width_continuous(range=c(0.5,2)) +  # Add some edge width variability
  scale_fill_manual(values = my_node_colors) +
  scale_color_manual(values = my_edge_colors) +
  theme_graph(base_family="Arial") +  # Choose a clean font family
  theme(legend.background = element_rect(colour = NA),
        panel.background = element_rect(fill = "white",  colour = NA),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),  # Soft background color
        plot.title = element_text(size = 18, face = "bold")) +  # Bolder title
  guides(size=FALSE)

#### 3 HGT network & Trait Association -- Figure3 ####

### calculate the network centrality

#head(species_network)

g <- graph_from_data_frame(species_network, directed = FALSE)


## Figure 3a Traitar plot
#trait_metadata <- read.table("traitar/trait_metadata.txt",sep="\t",header = T)
#trait_metadata$Phenotype.a. <- gsub(" ",".", trait_metadata$Phenotype.a.)
#trait_metadata$Phenotype.a. <- gsub("-",".", trait_metadata$Phenotype.a.)
#trait_metadata$Phenotype.a. <- gsub('\\(',".", trait_metadata$Phenotype.a.)
#trait_metadata$Phenotype.a. <- gsub(")",".", trait_metadata$Phenotype.a.)
#trait_metadata$Phenotype.a. <- gsub("°",".", trait_metadata$Phenotype.a.)
trait_metadata <- trait_metadata %>%
  filter(phypat.PGL.b. >= 0.9)
#bac_trait_data <- read.table("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/traitar/predictions_majority-vote_combined_hgt.txt",row.names = 1,sep = "\t",header = T)
bac_trait <- bac_trait_data[rownames(bac_trait_data) %in% centrality_df$Species_group, ]
bac_trait <- bac_trait[,colnames(bac_trait_data) %in% trait_metadata$Phenotype.a. ]
bac_trait[bac_trait == 1 ] <- 0
bac_trait[bac_trait > 1 ] <- 1
trait_cato <- read.table("traitar/traits.tsv",header=T,sep = "\t")

#trait_cato$accession <- gsub(" ",".", trait_cato$accession)
#trait_cato$accession <- gsub("-",".", trait_cato$accession)
#trait_cato$accession <- gsub('\\(',".", trait_cato$accession)
#trait_cato$accession <- gsub(")",".", trait_cato$accession)
#trait_cato$accession <- gsub("°",".", trait_cato$accession)
#trait_cato <- trait_cato[trait_cato$accession %in% trait_metadata$Phenotype.a.,]
#bac_cato <- read.table("traitar/samples.txt",header=T,sep="\t")

set1_colors <- c("#E59B8D","#A5A0CE","#92C2DD","#ADDB88","#369F2D","#FFE34F",
                 "#1663A9","#C4121A","#DE722A","#4B3D7C","#7F7F7F")
set2_colors <- c("#264653","#2a9d8f","#E59B8D","#A5A0CE","#a8dadc")


trait_colors <- colorRampPalette(set1_colors)(length(unique(trait_cato$category)))
bac_colors <- colorRampPalette(set2_colors)(length(unique(bac_cato$category)))

trait_color_mapping <- setNames(trait_colors, unique(trait_cato$category))
bac_color_mapping <- setNames(bac_colors, unique(factor(bac_cato$category)))

bac_trait_matrix <- as.matrix(bac_trait)

trait_to_category_map <- setNames(trait_cato$category, trait_cato$accession)
trait_categories_for_heatmap <- trait_to_category_map[colnames(bac_trait_matrix)]


bac_to_category_map <- setNames(bac_cato$category, bac_cato$sample_name)
bac_categories_for_heatmap <- bac_to_category_map[rownames(bac_trait_matrix)]

bac_ano$sample_name.2 <- paste0(bac_ano$sample_name.1," (",
                                gsub("Group_","",bac_ano$sample_name)
                                ,")")
bac_ano_map <- setNames(bac_ano$sample_name.2, bac_ano$sample_name)
bac_ano_for_heatmap <- bac_ano_map[rownames(bac_trait_matrix)]


heatmap_annotation_trait <- HeatmapAnnotation(
  Trait = trait_categories_for_heatmap,
  col = list(Trait = trait_color_mapping),
  show_annotation_name = TRUE,
  simple_anno_size = unit(5, "mm"),
  show_legend = T,
  annotation_name_gp = gpar(fontsize = 8)
)

heatmap_annotation_bac <- rowAnnotation(
  Phylum = bac_categories_for_heatmap,
  col = list(Phylum = bac_color_mapping),
  show_annotation_name = TRUE,
  simple_anno_size = unit(5, "mm"),
  show_legend = T
)

rownames(bac_trait_matrix) <- bac_ano_map[rownames(bac_trait_matrix)]
Heatmap(bac_trait_matrix,
        border=T,
        name = "Predictions",
        col = c("1" = "#457b9d", #"2" = "#74c69d", "1" = "#a8dadc",
                "0" = "white"),
        row_km = 3,
        column_km = 3,
        #width = unit(25,"cm"),
        #height = unit(10,"cm"),
        column_names_rot = 45,
        top_annotation = heatmap_annotation_trait,
        right_annotation = heatmap_annotation_bac,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        cluster_rows = TRUE,  # Clustering rows
        cluster_columns = TRUE,  # Clustering columns
        show_row_names = TRUE,
        show_column_names = TRUE)

#### Figure 3b HGT network annotated the degree centrality
degree_centrality <- igraph::degree(g) / (vcount(g) - 1)
closeness_centrality <- closeness(g)
betweenness_centrality <- betweenness(g)

centrality_df <- data.frame(
  Species = names(V(g)),
  Degree_Centrality = degree_centrality,
  Closeness_Centrality = closeness_centrality,
  Betweenness_Centrality = betweenness_centrality
)

head(genus_data)

species_to_genus <- setNames(genus_data$Phylum, genus_data$primary_cluster)
genus_list <- unique(genus_data$Phylum)

set_colors <- c( "#a8dadc","#1663A9","#92C2DD","#ADDB88","#92C2DD","#ADDB88",
                 "#E59B8D","black","#DE722A","#4B3D7C","#8C6D31","#264653",
                 "#2a9d8f","#E59B8D","#A5A0CE","#a8dadc","#E59B8D","#9467BD",
                 "#C4121A","#7F7F7F","#7F7F7F","#7F7F7F")
colors <- scales::hue_pal()(length(genus_list))
genus_to_color <- setNames(set_colors, genus_list)

genus_to_color <- setNames(set_colors, genus_list)

node_colors <- sapply(V(g)$name, function(x) genus_to_color[species_to_genus[x]])

V(g)$Degree_Centrality <- degree_centrality
V(g)$Closeness_Centrality <- closeness_centrality
V(g)$Betweenness_Centrality <- betweenness_centrality

V(g)$color <- sapply(V(g)$name, function(x) genus_to_color[species_to_genus[x]])

V(g)$family <- sapply(V(g)$name, function(x) species_to_genus[x])

centrality_df$Species_group <- paste0("Group_",centrality_df$Species)

species_mobile <- bac_trait$Motile == 1
species_gelatin_hydrolysis <- bac_trait$Gelatin.hydrolysis == 1
species_Sucrose <- bac_trait$Sucrose == 1

V(g)$trait_type <- ifelse(species_gelatin_hydrolysis, "gelatin_hydrolysis", "none")

V(g)$Mobile <- ifelse(species_mobile, "Mobile","None") 

ggraph(g, layout = "kk") + 
  geom_edge_link(color = "gray", width = 0.5, alpha = 0.2) +
  geom_node_point(aes(size = Betweenness_Centrality, fill = family, shape = trait_type),show.legend = T) +
  geom_node_point(aes(size = Betweenness_Centrality, color = Mobile, shape = trait_type),stroke = 1) +
  geom_node_text(aes(label = name), vjust = 1.8, size = 3) +
  scale_fill_manual(values = genus_to_color) + 
  scale_color_manual(values = c("red","grey")) +
  #scale_alpha_manual(values = c(0.5, 0.8)) + 
  scale_shape_manual(values = c("both" = 22, "Sucrose" = 24, "gelatin_hydrolysis" = 25, "none" = 21)) +
  theme_graph() 

ggraph(g, layout = "fr") + 
  geom_edge_link(color = "gray", width = 0.5, alpha = 0.2) +
  geom_node_point(aes(size = Betweenness_Centrality, color = Mobile, shape = trait_type),stroke = 1.5,show.legend = T) +
  geom_node_point(aes(size = Betweenness_Centrality, fill = family, shape = trait_type),alpha = 0.6) +
  geom_node_text(aes(label = name), vjust = 1.8, size = 3) +
  scale_fill_manual(values = set_colors) + 
  scale_color_manual(values = c("red","white")) +
  #scale_alpha_manual(values = c(0.5, 0.8)) + 
  scale_shape_manual(values = c( "gelatin_hydrolysis" = 25, "none" = 21)) +
  theme_graph() 

ggraph(g, layout = 'fr') + 
  geom_edge_link(color = "gray", width = 0.5, alpha = 0.2) +
  geom_node_point(aes(size = Degree_Centrality, color = trait_type, shape = Mobile),stroke = 1,show.legend = T) +
  geom_node_point(aes(size = Degree_Centrality, fill = family, shape = Mobile),alpha = 0.6) +
  geom_node_text(aes(label = name), vjust = 1.8, size = 3) +
  scale_fill_manual(values = genus_to_color) + 
  scale_color_manual(values = c("red","white")) +
  #scale_alpha_manual(values = c(0.5, 0.8)) + 
  scale_shape_manual(values = c( "Mobile" = 25, "None" = 21)) +
  theme_graph() 

#### Figure 3c Association between Centrality and Traits

row.names(centrality_df) <- centrality_df$Species_group
centrality_df$Species_group <- NULL
abun_df <- data.frame(species_name=genus_data$Group_cluster,
                      abun_relative=genus_data$abun_relative,
                      abun_cata=genus_data$abun_cata)
row.names(abun_df) <- abun_df$species_name
abun_df$species_name <- NULL
abun_df$abun_relative <- NULL

abun_df$abun_binary[abun_df$abun_cata == "Low"] = 0
abun_df$abun_binary[abun_df$abun_cata == "Medium"] = 1
abun_df$abun_binary[abun_df$abun_cata == "High" ] = 2
abun_df$abun_cata <- NULL

bac_trait1 <- merge(bac_trait,abun_df,by="row.names")
row.names(bac_trait1) <- bac_trait1$Row.names
bac_trait1$Row.names <- NULL
bac_net_tra <- merge(bac_trait1,centrality_df,by="row.names")


#distance_species <- read.csv("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/Tree_HGT_MAGs/Species_distances_bac.csv")
#distance_species$Group1 <- paste0("Group_",distance_species$Group1)
#distance_species$Group2 <- paste0("Group_",distance_species$Group2)
df <- distance_species

species <- unique(c(df$Group1, df$Group2))

matrix <- matrix(NA, length(species), length(species))

rownames(matrix) <- species
colnames(matrix) <- species

for(i in 1:nrow(df)) {
  matrix[rownames(matrix) == df$Group1[i], colnames(matrix) == df$Group2[i]] <- df$Average_Distance[i]
  matrix[rownames(matrix) == df$Group2[i], colnames(matrix) == df$Group1[i]] <- df$Average_Distance[i]
}

diag(matrix) <- 0


distance_matrix <- matrix

#### phylo model
centrality_measures <- c("Degree_Centrality", "Betweenness_Centrality")

results <- data.frame(
  Centrality = character(), 
  Trait = character(), 
  intercept = numeric(), 
  Centrality_number = numeric(), 
  Trait_number = numeric(), 
  coef_trait = numeric(), 
  t_value_trait = numeric(), 
  p_value_trait = numeric(), 
  se_trait = numeric(),
  adj_p_value_trait = numeric(), # New column for adjusted p-values
  stringsAsFactors = FALSE
)
for (p in centrality_measures) {
  temp_results <- data.frame() # Temporary data frame for each phenotype
  for (i in colnames(bac_trait)){
    data_model <- drop_na(bac_net_tra)
    
    model_formula <- as.formula(paste(p, " ~", i,"+ abun_binary + (1 | Row.names)"))
    model_formula0 <- as.formula(paste(p, " ~ abun_binary + (1 | Row.names) "))
    # Using tryCatch to handle errors
    result <- tryCatch({
      model0 <- relmatLmer(model_formula0, data_model, relmat = list(ID = distance_matrix), REML = FALSE)
      model <- relmatLmer(model_formula,data_model, relmat = list(ID = distance_matrix), REML = FALSE)
      
      summa <- summary(model)
      coef_trait <- summa$coefficients[[i,"Estimate"]]
      
      summa_p <- lrtest(model,model0)
      p_value_trait <- summa_p$`Pr(>Chisq)`[2]
      
      # Create a data frame with the results
      temp_results_current <- data.frame(
        Centrality = p, 
        Trait = i, 
        Centrality_number =  sum(data_model[p] != 0 & !is.na(data_model[p])),
        Trait_number = sum(data_model[i] != 0 & !is.na(data_model[i])),
        t_value_trait = t_value_trait,
        coef_trait = coef_trait, 
        p_value_trait = p_value_trait, 
        se_trait = se_trait,
        adj_p_value_trait = NA, # Initialize with NA
        stringsAsFactors = FALSE
      )
      temp_results_current
      
    }, error = function(e) {
      # Return NA values in case of an error
      return(data.frame(
        Centrality = p, 
        Trait = i, 
        Centrality_number =  sum(data_model[p] != 0 & !is.na(data_model[p])),
        Trait_number = sum(data_model[i] != 0 & !is.na(data_model[i])),
        t_value_trait = NA,
        coef_trait = NA, 
        p_value_trait = NA, 
        se_trait = NA,
        adj_p_value_trait = NA,
        stringsAsFactors = FALSE
      ))
    })
    
    temp_results <- rbind(temp_results, result)
  }
  temp_results$adj_p_value_trait <- p.adjust(temp_results$p_value_trait, method = "fdr")
  
  results <- rbind(results, temp_results)
}

## results analyisis
trait_centrality_association <- results
trait_centrality_association <- trait_centrality_association %>%
  filter(Trait_number > 10)
trait_centrality_association$z_value <- -qnorm(trait_centrality_association$p_value_trait/2)*sign(trait_centrality_association$coef_trait)
trait_centrality_association$BH_p_value_trait <- trait_centrality_association$adj_p_value_trait

# Reshape the data to wide format for z-values
heatmap_data <- dcast(trait_centrality_association, Trait ~ Centrality, value.var = "z_value")
# Remove columns where all values are NA
row.names(heatmap_data) <- heatmap_data$Trait
heatmap_data$Trait <- NULL
heatmap_data$Closeness_Centrality <- NULL
heatmap_data <- heatmap_data[, colSums(is.na(heatmap_data)) != nrow(heatmap_data)]
heatmap_data <- heatmap_data[rowSums(is.na(heatmap_data)) != ncol(heatmap_data), ]

# Create a matrix for the significance levels
significance_matrix <- dcast(trait_centrality_association, Trait ~ Centrality, value.var = "BH_p_value_trait")
row.names(significance_matrix) <- significance_matrix$Trait
significance_matrix$Trait <- NULL
significance_matrix <- significance_matrix[rownames(heatmap_data),colnames(heatmap_data)]

significance_matrix_converted <- apply(significance_matrix, c(1, 2), function(x) {
  if (is.na(x)) {
    "" 
  } else if (x < 0.01) {
    "**"
  } else if (x < 0.05) {
    "*"
  } else {
    ""
  }
})

color_mapping <- colorRamp2(c(min(heatmap_data, na.rm = TRUE), 0, max(heatmap_data, na.rm = TRUE)), 
                            c("#8491B4FF","#EAEFF6", "#91D1C2FF"))

color_mapping <- colorRamp2(c(-max(heatmap_data, na.rm = TRUE), 0, max(heatmap_data, na.rm = TRUE)), 
                            c("#8491B4FF","#EAEFF6", "#91D1C2FF"))
heatmap_data[is.na(heatmap_data)] <- 0
p_phy_trait<-Heatmap(heatmap_data,
                     name = "z-value",
                     width = unit(3,"cm"),
                     height = unit(10,"cm"),
                     col = color_mapping,
                     rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                     cluster_rows = T,
                     cluster_columns = F,
                     column_names_rot = 45,
                     na_col = "gray",  # Correct way to specify color for NA values
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       if (!is.na(significance_matrix_converted[i, j])) {
                         grid.text(significance_matrix_converted[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))
                       }
                     },
                     show_row_names = TRUE,
                     show_column_names = F,
                     row_names_gp = gpar(fontsize = 12),
                     column_names_gp = gpar(fontsize = 12),
)

p_phy_trait

trait_centrality_association$size <- -log10(trait_centrality_association$BH_p_value_trait)
trait_centrality_association1 <- trait_centrality_association %>%
  filter(Centrality %in% c("Betweenness_Centrality"))  %>%
  group_by(Trait) %>%
  filter(sum(is.na(z_value)) < n_distinct(Centrality)) %>%
  ungroup()

trait_centrality_association1 <- trait_centrality_association1 %>%
  mutate(shape = case_when(
    BH_p_value_trait <= 0.01 ~ 18,
    BH_p_value_trait <= 0.05 ~ 17,
    TRUE ~ 16
  ))
trait_centrality_association1$color <- ifelse(trait_centrality_association1$BH_p_value_trait <= 0.1, trait_centrality_association1$z_value, NA)
sorted_data <- trait_centrality_association1 %>%
  arrange(desc(z_value))
ggplot(sorted_data, aes(x = Centrality, y = Trait, size = size, color = z_value)) +
  geom_point(aes(shape = factor(shape)), alpha = 0.8) + 
  scale_shape_manual(values = c(16, 17,18,15)) +  
  scale_color_gradient2(low = "#104674", high = "#B45A56", mid = "white", na.value = "grey80") +  
  scale_size(range = c(3, 10)) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right"
  ) +
  labs(x = "Centrality", y = "Trait", size = "P-value (-log10)", color = "Z-value") 

#### 4 HGT rates & Species Correlation & Species phylogenetic distance -- Figure4 ####

#Table_HGT_distance <- read.csv("Tree_HGT_MAGs/Species_HGT_distance_filter_Final.csv",sep=",",header=1,row.names = 2)
#Table_HGT_distance$X <- NULL

## Figure 4a HGT ~ Phylo distance

bin_width <- 0.2
breaks <- seq(0, max(Table_HGT_distance$Average_Distance), by = bin_width)  

Table_HGT_distance$bin <- cut(Table_HGT_distance$Average_Distance, breaks = breaks, labels = FALSE, right = TRUE, include.lowest = TRUE)
Table_HGT_distance <- Table_HGT_distance %>% filter(!is.na(bin))
avg_HGT_by_bin <- aggregate(HGT_Rates ~ bin, data = Table_HGT_distance, FUN = mean)
avg_HGT_by_bin$bin_center <- (as.numeric(as.character(avg_HGT_by_bin$bin)) - 0.5) * bin_width
avg_HGT_by_bin$bin_center <- as.factor(avg_HGT_by_bin$bin_center)

Table_HGT_distance$bin_center <- as.factor(Table_HGT_distance$bin_center)
model_cor <- lm(HGT_Rates ~ Average_Distance,Table_HGT_distance)
model_cor_sum <- summary(model_cor)
Estimate <- model_cor_sum$coefficients[rownames(model_cor_sum$coefficients)[2],1]
Pvalue <- model_cor_sum$coefficients[rownames(model_cor_sum$coefficients)[2],4]

ggplot() +
  geom_bar(data = avg_HGT_by_bin, aes(x = bin_center, y = HGT_Rates), 
           stat = "identity", fill = "#8491B4FF", alpha = 0.5) +
  geom_boxplot(data = Table_HGT_distance, aes(x = bin_center, y = HGT_Rates), 
               fill = "#8491B4FF", alpha = 0.5, width = 0.3) +
  labs(
    x = "Phylogenetics Distance",
    y = "HGT Rates",
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  ) + 
  annotate("text", 
           x = Inf, y = Inf,
           label = paste("p-value : ", 
                         format(Pvalue, digits=4), sep=""),
           hjust = 1.1, vjust = 2, size = 4)

## Figure 4b Correlation ~ phylogeny distance
Table_Cor_distance <- read.csv("Tree_HGT_MAGs/Species_Cor_distance_0.01Final_relative.csv",sep=",",header=1,row.names = 1)
Table_Cor_distance$Correlation_all <- as.factor(Table_Cor_distance$Correlation_all)

bin_width <- 0.2
breaks <- seq(0, max(Table_Cor_distance$Average_Distance), by = bin_width)

Table_Cor_distance$bin <- cut(Table_Cor_distance$Average_Distance, breaks = breaks, labels = FALSE, right = TRUE, include.lowest = TRUE)

avg_cor_by_bin <- aggregate(correlation_value ~ bin, data = Table_Cor_distance, FUN = mean)

avg_cor_by_bin$bin_center <- (as.numeric(as.character(avg_cor_by_bin$bin)) - 1) * bin_width + (bin_width / 2)

avg_cor_by_bin$bin_center <- as.factor(avg_cor_by_bin$bin_center)
Table_Cor_distance$bin_center <- as.factor(Table_Cor_distance$bin_center)
Table_Cor_distance <- Table_Cor_distance %>%
  filter(!is.na(bin_center))
model_cor <- lm(correlation_value ~ Average_Distance,Table_Cor_distance)
model_cor_sum <- summary(model_cor)
Estimate <- model_cor_sum$coefficients[rownames(model_cor_sum$coefficients)[2],1]
Pvalue <- model_cor_sum$coefficients[rownames(model_cor_sum$coefficients)[2],4]

ggplot() +
  geom_bar(data = avg_cor_by_bin, aes(x = bin_center, y = correlation_value), 
           stat = "identity", fill = "#91D1C2FF", alpha = 0.5) +
  geom_boxplot(data = Table_Cor_distance, aes(x = bin_center, y = correlation_value), 
               fill = "#91D1C2FF", alpha = 0.5, width = 0.3) +
  labs(
    x = "Phylogenetics Distance",
    y = "Correlation",
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  ) + annotate("text", 
               x = Inf, y = Inf,
               label = paste("p-value : ", 
                             format(Pvalue, digits=4), sep=""),
               hjust = 1.1, vjust = 2, size = 4)

##  HGT ~ Correlation + Phylo distance
Table_HGT_distance <- read.csv("Tree_HGT_MAGs/Species_HGT_distance_filter_Final.csv",sep=",",header=1)

Merged_Table <- merge(Table_HGT_distance, Table_Cor_distance, by="Species_Pair", all=TRUE)

standardize_pair <- function(gene1, gene2) {
  return(paste0(gene1, "-", gene2))
}
Merged_Table <- mutate(Merged_Table, HGT = ifelse(is.na(common_classification), "nHGT", "HGT"))
Merged_Table <- mutate(Merged_Table, COR = ifelse(is.na(correlation_value), "U", ifelse(correlation_value > 0, "P", "N")))
Merged_Table$HGT_COR <- mapply(standardize_pair, Merged_Table$HGT, Merged_Table$COR)
Merged_Table$Average_Distance <- ifelse(is.na(Merged_Table$Average_Distance.y), 
                                        Merged_Table$Average_Distance.x, 
                                        Merged_Table$Average_Distance.y)

ggplot(data = Merged_Table, aes(x = Average_Distance.y, y = HGT_Rates)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") + 
  ggtitle("Linear Relationship between HGT Rates and Average Distance") +
  xlab("Average_Distance.y") +
  ylab("HGT_Rates")

ggplot(data = Merged_Table, aes(x = Average_Distance.y, y = HGT_Rates, color = Correlation_all)) +
  geom_point() +
  geom_smooth(method = "lm", aes(color = Correlation_all, group = Correlation_all), se=FALSE) +  # 线性回归线
  scale_color_manual(values = c("N-N" = "blue", "P-P" = "green")) + 
  ggtitle("Linear Relationship between HGT Rates and Average Distance") +
  xlab("Average Distance") +
  ylab("HGT Rates") +
  theme(legend.title = element_blank()) 

library(broom)

#### 3 class
Merged_Table1 <- Merged_Table %>%
  filter(!is.na(HGT_Rates)) %>%
  mutate(Correlation_group = case_when(
    Correlation_all == "N-N" ~ "N-N",
    Correlation_all == "P-P" ~ "P-P",
    TRUE ~ "NA"
  ))


HGT_COR_abun_model <- lm(HGT_Rates ~ correlation_value + Average_Distance, data = Merged_Table1)
HGT_COR_abun_model_sum <- summary(HGT_COR_abun_model)

pvalue_dis <- HGT_COR_abun_model_sum$coefficients["Average_Distance",4]
pvalue_cor <- HGT_COR_abun_model_sum$coefficients["correlation_value",4]

plot3 <- ggplot(data = Merged_Table1, aes(x = correlation_value, y = HGT_Rates)) +
  geom_point(aes(color = Correlation_group, alpha = 0.9)) +
  geom_smooth(data = subset(Merged_Table1, Correlation_group != "NA"), 
              method = "lm", 
              aes(color = HGT, group = HGT), se=T) +
  scale_color_manual(values = c("N-N" = "#8491B4FF", "P-P" = "#B45A56"),
                     labels = c("Negative","Positive")) +
  labs(
    y = "HGT Rates",
    x = "Correlation",
    color = "Correlation"
  ) +
  theme_classic() + 
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  ) + 
  guides(alpha = "none")  + annotate("text", 
                                     x = Inf, y = Inf,
                                     label = "HGT_Rates ~ correlation_value + Average_Distance",
                                     hjust = 1.1, vjust = 2, size = 4) +
  guides(alpha = "none")  + annotate("text", 
                                     x = Inf, y = Inf,
                                     label = paste("p-value (Phylogenetics Distance): ", 
                                                   format(pvalue_dis, digits=4), sep=""),
                                     hjust = 1.1, vjust = 6, size = 4) +
  annotate("text", 
           x = Inf, y = Inf,
           label = paste("p-value (Correlation): ", 
                         format(pvalue_cor, digits=4), sep=""),
           hjust = 1.1, vjust = 4, size = 4)

plot3

## HGT ~ Cor
###new version mantel test
Table_HGT_distance1 <- read.csv("Tree_HGT_MAGs/Species_HGT_distance_filter_Final.csv",sep=",",header=1,row.names = 1)
species <- rownames(distance_matrix)
HGT_matrix <- matrix(NA, nrow = length(species), ncol = length(species))
rownames(HGT_matrix) <- colnames(HGT_matrix) <- species
for (i in 1:nrow(Table_HGT_distance1)) {
  species1 <- Table_HGT_distance1$Species1_Group[i]
  species2 <- Table_HGT_distance1$Species2_Group[i]
  HGT_rate <- Table_HGT_distance1$HGT_Rates[i]
  HGT_matrix[species1, species2] <- HGT_rate
  HGT_matrix[species2, species1] <- HGT_rate
}
cor_matrix <- matrix(NA, nrow = length(species), ncol = length(species))
rownames(cor_matrix) <- colnames(cor_matrix) <- species
Table_Cor_distance1 <- Table_Cor_distance %>%
  filter(!Genome1.x %in% c("Group_1","Group_2")) %>%
  filter(!Genome2.x %in% c("Group_1","Group_2"))
for (i in 1:nrow(Table_Cor_distance1)) {
  species1 <- Table_Cor_distance1$Genome1.x[i]
  species2 <- Table_Cor_distance1$Genome2.x[i]
  cor_value <- Table_Cor_distance1$correlation_value[i]
  cor_matrix[species1, species2] <- cor_value
  cor_matrix[species2, species1] <- cor_value
}

HGT_bin_matrix <- HGT_matrix
HGT_bin_matrix[is.na(HGT_bin_matrix)] <- 0

ensure_symmetric <- function(mat) {
  mat <- (mat + t(mat)) / 2
  diag(mat) <- NA
  return(mat)
}
HGT_bin_matrix_clean <- ensure_symmetric(HGT_bin_matrix)
dist_cor <- 1 - cor_matrix
dist_HGT <- 1 - HGT_matrix
dis_HGT_bin <- 1 - HGT_bin_matrix_clean

library(ncf)
result <- partial.mantel.test(M1 = dist_cor, M2 = dis_HGT_bin, M3 = distance_matrix, resamp = 9999)
print(result)
###old version
Table_HGT_distance1 <- read.csv("Tree_HGT_MAGs/Species_HGT_distance_filter_Final.csv",sep=",",header=1,row.names = 1)
Table_HGT_distance1$HGT_Rates[is.na(Table_HGT_distance1$HGT_Rates)] <- 0
Merged_Table <- merge(Table_HGT_distance1, Table_Cor_distance, by="Species_Pair", all=TRUE)
standardize_pair <- function(gene1, gene2) {
  return(paste0(gene1, "-", gene2))
}
Merged_Table <- mutate(Merged_Table, HGT = ifelse(is.na(common_classification), "nHGT", "HGT"))
Merged_Table <- mutate(Merged_Table, COR = ifelse(is.na(correlation_value), "U", ifelse(correlation_value > 0, "P", "N")))
Merged_Table$HGT_COR <- mapply(standardize_pair, Merged_Table$HGT, Merged_Table$COR)
Merged_Table$Average_Distance <- ifelse(is.na(Merged_Table$Average_Distance.y), 
                                        Merged_Table$Average_Distance.x, 
                                        Merged_Table$Average_Distance.y)

Merged_Table$HGT_cato[Merged_Table$HGT_Rates >= 0] <- 1
Merged_Table$HGT_cato[is.na(Merged_Table$HGT_Rates)] <- 0
Merged_Table$correlation_cato[Merged_Table$correlation_value > 0] <- 1
Merged_Table$correlation_cato[Merged_Table$correlation_value < 0] <- 0

model_formula <- as.formula(paste("HGT_cato ~  correlation_value + Average_Distance"))
model <- glm(model_formula, family = binomial(),data = Merged_Table)
model_sum <- summary(model)

model_formula1 <- as.formula(paste("HGT_cato ~  correlation_value"))
model1 <- glm(model_formula1, family = binomial(),data = Merged_Table)
model1_sum <- summary(model1)

model_formula2 <- as.formula(paste("HGT_cato ~ correlation_cato + Average_Distance"))
model2 <- glm(model_formula2, family = binomial(),data = Merged_Table)
model2_sum <- summary(model2)

model_formula3 <- as.formula(paste("HGT_Rates ~ correlation_cato + Average_Distance"))
model3 <- glm(model_formula3, family = binomial(),data = Merged_Table)
model3_sum <- summary(model3)

Merged_Table %>% 
  ggplot(aes(x = correlation_value)) +
  geom_histogram() +
  facet_wrap(~ HGT_cato, ncol = 1)

tidy(model) %>% 
  knitr::kable()

predict_data <- data.frame(
  correlation_value = seq(from = min(Merged_Table$correlation_value, na.rm = TRUE), 
                          to = max(Merged_Table$correlation_value, na.rm = TRUE), length.out = 100)
)

predict_data <- predict_data %>%
  mutate(
    fit = predict(model1, newdata = ., type = "response"),
    ll = predict(model1, newdata = ., type = "response", se.fit = TRUE)$fit - 1.96 * predict(model1, newdata = ., type = "response", se.fit = TRUE)$se.fit,
    ul = predict(model1, newdata = ., type = "response", se.fit = TRUE)$fit + 1.96 * predict(model1, newdata = ., type = "response", se.fit = TRUE)$se.fit
  )
glimpse(predict_data)

predict_data %>%
  ggplot(aes(x = correlation_value, y = fit)) +
  geom_ribbon(aes(ymin = ll, ymax = ul),
              alpha = 1/2) +
  geom_line(aes(y=fit)) + 
  geom_jitter(data = Merged_Table,
              aes(y = HGT_cato),
              size = 1/4, alpha = 1/2, height = 0.05) +
  scale_color_manual(values = c("1" = "#009E73", "0" = "#D55E00"),
                     labels = c("HGT","no-HGT")) +
  scale_y_continuous("probability of HGT", 
                     expand = expansion(mult = 0.01))
library("ggdist")
#Merged_Table$HGT_cato <- as.numeric(Merged_Table$HGT_cato)
predict_data %>%
  ggplot(aes(x = correlation_value, y = fit)) +
  geom_ribbon(aes(ymin = ll, ymax = ul),
              alpha = 1/2) +
  geom_line() + 
  stat_dots(data = Merged_Table,
            aes(y = HGT_cato, side = ifelse(HGT_cato == 0, "top", "bottom")),
            scale = 1/3) +
  scale_y_continuous("probability of HGT", 
                     expand = expansion(mult = 0.01))

library(gghalves)  

p3.5 <- predict_data %>%
  ggplot(aes(y = correlation_value)) +
  geom_half_violin(data = Merged_Table %>% 
                     mutate(binary = factor(binary, levels = c("HGT", "no-HGT"))),
                   aes(y = correlation_value, x = HGT_cato, fill = binary),
                   position = position_nudge(x = .1, y = 0),
                   adjust=1.5, trim=FALSE, colour=NA, side = 'r') +
  geom_boxplot(data = Merged_Table %>% 
                 mutate(binary = factor(binary, levels = c("HGT", "no-HGT"))),
               aes(y = correlation_value, x = HGT_cato, fill = binary), 
               width = 0.05, outlier.shape = NA) + 
  scale_fill_manual("Condition", values = c(alpha("#B45A56", 0.5), alpha("#8491B4FF", 0.5))) +
  scale_x_continuous("HGT Occurance", expand = c(0, 0)) +
  scale_y_continuous("Correlation", expand = c(0, 0)) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  ) +
  annotate("text", x = Inf, y = Inf,
           label = "HGT_cato ~ correlation_value + Average_Distance",
           hjust = 1.1, vjust = 4, size = 4) +
  annotate("text", x = Inf, y = Inf,
           label = paste("p-value (Correlation): ", 
                         format(Pvalue_model_cor, digits=4), sep=""),
           hjust = 1.1, vjust = 6, size = 4) +
  annotate("text", x = Inf, y = Inf,
           label = paste("p-value (Phylogenetics Distance): ", 
                         format(Pvalue_model_dis, digits=4), sep=""),
           hjust = 1.1, vjust = 8, size = 4) +
  coord_flip()

## Chi-square test

count_all <- Merged_Table %>% 
  filter(!is.na(Correlation_all)) %>% 
  group_by(Correlation_all) %>% 
  summarize(count = n())

count_both_not_na <- Merged_Table %>% 
  filter(!is.na(Correlation_all) & !is.na(HGT_Rates)) %>% 
  group_by(Correlation_all) %>% 
  summarize(count = n())

chi_table <- matrix(c(count_all$count[1], count_all$count[2],
                      count_both_not_na$count[1], count_both_not_na$count[2]), 
                    nrow = 2, byrow = TRUE)


colnames(chi_table) <- c("N-N", "P-P")
rownames(chi_table) <- c("ALL", "HGT")

chisq_result <- chisq.test(chi_table)
fisher_result <- fisher.test(chi_table)

chisq_result

## Plot

plot_data <- data.frame(
  Correlation_all = c("N-N", "P-P", "N-N", "P-P"),
  Category = c(rep("ALL", 2), rep("HGT", 2)),
  Count = c(chi_table[1,], chi_table[2,])
)

plot_data$Count[plot_data$Category == "ALL" & plot_data$Correlation_all == "N-N"] <- plot_data$Count[plot_data$Category == "ALL" & plot_data$Correlation_all == "N-N"] - plot_data$Count[plot_data$Category == "HGT" & plot_data$Correlation_all == "N-N"]
plot_data$Count[plot_data$Category == "ALL" & plot_data$Correlation_all == "P-P"] <- plot_data$Count[plot_data$Category == "ALL" & plot_data$Correlation_all == "P-P"] - plot_data$Count[plot_data$Category == "HGT" & plot_data$Correlation_all == "P-P"]

p_value <- chisq_result$p.value

ggplot(data = plot_data, aes(x = Correlation_all, y = Count, fill = interaction(Correlation_all, Category, lex.order = TRUE))) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("N-N.ALL" = scales::alpha("#8491B4FF", 0.5), "N-N.HGT" = "#8491B4FF",
                               "P-P.ALL" = scales::alpha("#B45A56", 0.5), "P-P.HGT" = "#B45A56"),
                    labels = c("Negative No-HGT", "Negative HGT", "Positive no-HGT", "Positive HGT")) +
  labs(
    #title = "Distribution of Correlation_all with HGT_Rates",
    x = "Correlation",
    y = "Count",
    fill = "Category"
  ) +
  theme_classic() + 
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) + geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3.5) +
  scale_x_discrete(labels = c("N-N", "P-P")) +
  annotate("text", 
           x = 1.5, y = max(count_all$count), 
           label = paste("Chi-squared p-value:", format(p_value, digits=4)), 
           hjust = 0.5, vjust = 0, size = 4)


#### 5 Microbial Community Guild -- Figure5 ####

## Figure 5a-b Microbial communities

#data <- read.table("HGT_events/Cor_HGT_network/HGT_COR_0.01_Final.csv",header=T,sep = ",")
head(data)
clean_data <- data %>%
  mutate(edge_weight = ifelse(is.na(correlation_value), 0, correlation_value)) %>%
  select(Genome1, Genome2, edge_weight, HGT_COR)

graph <- graph_from_data_frame(clean_data, directed = FALSE, vertices = NULL)
E(graph)$weight <- clean_data$edge_weight

leiden_cluster <- cluster_leiden(graph, weights = E(graph)$weight,
                                 resolution_parameter = 0.001,
                                 objective_function = c("CPM"),)

membership(leiden_cluster)

V(graph)$grp <- as.character(membership(leiden_cluster))

transformed_weights <- 1 - ((E(graph)$weight + 1) / 2)

layout <- layout_with_kk(graph, weights = transformed_weights)

graph_tbl <- as_tbl_graph(graph)

groups_of_interest <- c("1","2", "3")
selected_nodes <- V(graph)[V(graph)$grp %in% groups_of_interest]

abun_df1<- data.frame(abun = genus_data$abun_cata, name=genus_data$primary_cluster)

abun_network <- setNames(abun_df1$abun, abun_df1$name)

V(graph_tbl)$abun <- abun_network[as.character(V(graph_tbl)$name)]

library(ggforce)
ggraph(graph_tbl, layout = "manual", x = layout[, 1], y = layout[, 2]) +
  geom_edge_link(aes(color = HGT_COR), width = 0.2, alpha = 0.3) +
  geom_node_point(aes(fill = grp,size=abun), shape = 21) +
  geom_mark_hull(
    aes(x = x, y = y, group = grp1, fill = grp1),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0
  ) +
  scale_size_manual(values = c("High" = 4, "Medium" = 3, "Low" = 2)) +
  scale_edge_size_manual(values = c("HGT-N" = 1, "HGT-P" = 1, "HGT-U" = 1, 
                                    "nHGT-N" = 0.0001, "nHGT-P" = 0.0001))+
  scale_edge_colour_manual(values = c(#"HGT-N" = "#4DAF4A", "HGT-P" = "#4DAF4A", "HGT-U" = "#4DAF4A",
    "HGT-N" = "transparent", "HGT-P" = "transparent", "HGT-U" =  "transparent",
    #"nHGT-N" = "transparent", "nHGT-P" = "transparent")) +
    "nHGT-N" = "#377EB8", "nHGT-P" = "#E41A1C")) +
  scale_fill_brewer(palette = "Set1") +
  theme_graph() +
  theme(legend.position = "none") +
  facet_nodes(~grp)

##### Figure5c The genes enriched in each community
cluster_membership <- membership(leiden_cluster)

module_MAGs <- data.frame(node = names(cluster_membership), Guild = cluster_membership)

annotation_table <- read.csv("HGT_network/nodes_group_tax_table.csv",header=T)
annotation_table <- dplyr::rename(annotation_table, node = primary_cluster)

module_MAGs <- merge(annotation_table,module_MAGs,by = "node")
module_MAGs$Group_cluster <- paste0("Group_",module_MAGs$node)
module_MAGs

#### AMP

#data <- read.table("MAGs_Func/result/merged_AMP_files_70.txt", header = T, sep = "\t", row.names = 1)

MAGs_AMP <- data
module_MAGs <- subset(module_MAGs, ! Group_cluster %in% "")
module_MAGs <- subset(module_MAGs, ! Group_cluster %in% c("Group_13","Group_12","Group_11"))
module_MAGs <- subset(module_MAGs, Guild %in% c("1","2","3"))

module_MAGs$name <- module_MAGs$Group_cluster
module_MAGs$Guild <- as.factor(module_MAGs$Guild)
#gene_number_adj <- read.table("MAGs_Func/result/Species_Gene_Number.txt",sep = '\t',header = T,row.names = 1)
annotate_table <- data.frame(row.names = module_MAGs$name, groups=module_MAGs$Guild)

annotation_row_names <- rownames(annotate_table)

MAGs_AMP_colnames <- colnames(MAGs_AMP)

missing_colnames <- setdiff(annotation_row_names, MAGs_AMP_colnames)

for (colname in missing_colnames) {
  MAGs_AMP[, colname] <- 0
}

MAGs_AMP <- MAGs_AMP[, annotation_row_names]

sum_row <- colSums(MAGs_AMP, na.rm = TRUE)

MAGs_AMP1 <- rbind(MAGs_AMP, sum_row)

rownames(MAGs_AMP1)[nrow(MAGs_AMP1)] <- "AMP"


column_ha = HeatmapAnnotation(df = annotate_table)

row_ha = rowAnnotation(df = annotate_table)

#### primary_cluster annotation

primary_cluster_df <- data.frame(primary_cluster = as.factor(module_MAGs$Phylum), row.names = module_MAGs$name)

unique_clusters <- unique(primary_cluster_df$primary_cluster)
library(RColorBrewer)
color_palette <- c("#017D6F","#0AAF88","#FDDB76","#8AD088","#B686B6","#88CEE6","#E89DA0")
color_mapping <- setNames(color_palette, unique_clusters)

MAG_phylum_map <- setNames(module_MAGs$Phylum, as.factor(module_MAGs$name))


library(circlize)
col_fun = colorRamp2(c(0,max(MAGs_AMP)), c("white", "#E41A1C"))
col_fun(seq(-3, 3))

plot_AMP<-Heatmap(t(MAGs_AMP), , show_row_names=F,col=col_fun,
                  row_split = annotate_table,
                  left_annotation = row_ha ,
                  cluster_columns = TRUE
)

plot_AMP

col_fun = colorRamp2(c(0,max(MAGs_AMP1["AMP",])), c("white", "#E41A1C"))
col_fun(seq(-3, 3))
AMP_sum <- t(MAGs_AMP1["AMP",])

plot_AMP1<-Heatmap(t(MAGs_AMP1["AMP",]), show_row_names=F,col=col_fun,
                   row_split = annotate_table,
                   #left_annotation = row_ha ,
                   right_annotation = primary_cluster_annotation,
                   cluster_columns = TRUE
                   #row_split = row.names(MAGs_cog)
)
plot_AMP1

AMP_sum_adj <- merge(AMP_sum,gene_number_adj,by= "row.names")
rownames(AMP_sum_adj) <- AMP_sum_adj$Row.names
AMP_sum_adj$Row.names <- NULL
AMP_sum_adj <- merge(AMP_sum_adj,annotate_table,by= "row.names")
rownames(AMP_sum_adj) <- AMP_sum_adj$Row.names
AMP_sum_adj$Row.names <- NULL

AMP_sum_adj <-  AMP_sum_adj %>%
  mutate(AMP1 = (AMP*1000/number) )
AMP_heat <- AMP_sum_adj[, "AMP1", drop = FALSE]


annotate_table1 <- annotate_table[match(rownames(AMP_heat), rownames(annotate_table)), ]

col_fun = colorRamp2(c(0,max(ARO_heat)/4,max(AMP_heat)), c("white", "white","#E41A1C"))
col_fun(seq(-3, 3))

MAG_phy_forheat <- MAG_phylum_map[rownames(AMP_heat)]

MAG_phy_ano <- rowAnnotation(
  Phylum = MAG_phy_forheat,
  col = list(Phylum = color_mapping),
  width = unit(1, "cm") 
)

plot_AMP1<-Heatmap(AMP_heat, show_row_names=F,col=col_fun,
                   row_names_gp = gpar(fontsize = 6),
                   row_split = annotate_table1,
                   cluster_columns = F
)

plot_AMP1

AMP_data <- t(MAGs_AMP1["AMP",])


##### VFs 

#data <- read.table("MAGs_Func/result/merged_VFs_files_50_70.txt", header = T, sep = "\t", row.names = 1)

#MAGs_VFs <- data

annotation_row_names <- rownames(annotate_table)

MAGs_VFs_colnames <- colnames(MAGs_VFs)

missing_colnames <- setdiff(annotation_row_names, MAGs_VFs_colnames)

for (colname in missing_colnames) {
  MAGs_VFs[, colname] <- 0
}

MAGs_VFs <- MAGs_VFs[, annotation_row_names]

sum_row <- colSums(MAGs_VFs, na.rm = TRUE)

MAGs_VFs1 <- rbind(MAGs_VFs, sum_row)

rownames(MAGs_VFs1)[nrow(MAGs_VFs1)] <- "VFs"

column_ha = HeatmapAnnotation(df = annotate_table)
row_ha = rowAnnotation(df = annotate_table)

library(circlize)
col_fun = colorRamp2(c(0,max(MAGs_VFs)), c("white", "#377EB8"))
col_fun(seq(-3, 3))

plot_VFs<-Heatmap(t(MAGs_VFs), show_row_names=F,col=col_fun,
                  row_split = annotate_table,
                  left_annotation = row_ha ,
                  cluster_columns = TRUE
)
plot_VFs

col_fun = colorRamp2(c(0,max(MAGs_VFs1["VFs",])), c("white", "#377EB8"))
col_fun(seq(-3, 3))

plot_VFs1<-Heatmap(t(MAGs_VFs1["VFs",]), show_row_names=F,col=col_fun,
                   #row_split = annotate_table,
                   #left_annotation = row_ha ,
                   #right_annotation = primary_cluster_annotation,
                   cluster_columns = TRUE
                   #row_split = row.names(MAGs_cog)
)

plot_VFs1

VF_sum <- t(MAGs_VFs1["VFs",])
VF_sum_adj <- merge(VF_sum,gene_number_adj,by= "row.names")
rownames(VF_sum_adj) <- VF_sum_adj$Row.names
VF_sum_adj$Row.names <- NULL
VF_sum_adj <- merge(VF_sum_adj,annotate_table,by= "row.names")
rownames(VF_sum_adj) <- VF_sum_adj$Row.names
VF_sum_adj$Row.names <- NULL

VF_sum_adj <-  VF_sum_adj %>%
  mutate(VFs1 = (VFs*1000/number) )
VF_heat <- VF_sum_adj[, "VFs1", drop = FALSE]
annotate_table1 <- annotate_table[match(rownames(VF_heat), rownames(annotate_table)), ]

col_fun = colorRamp2(c(0,max(VF_heat)/4,max(VF_heat)), c("white","white", "#377EB8"))
col_fun(seq(-3, 3))

plot_VFs1<-Heatmap(VF_heat, show_row_names=F,col=col_fun,
                   row_split = annotate_table1,
                   cluster_columns = TRUE
)

plot_VFs1

plot_AMP1 + plot_VFs1
VF_data <- t(MAGs_VFs1["VFs",])
#### ARO

#data <- read.table("MAGs_Func/result/merged_ARO_files_strict.txt", header = T, sep = "\t", row.names = 1)


#MAGs_AROs <- data

MAGs_AROs_colnames <- colnames(MAGs_AROs)

missing_colnames <- setdiff(annotation_row_names, MAGs_AROs_colnames)

for (colname in missing_colnames) {
  MAGs_AROs[, colname] <- 0
}

MAGs_AROs <- MAGs_AROs[, annotation_row_names]

#ARG_class <- read.table("MAGs_Func/ARO/class_CARD_id.txt",header = T, sep = "\t", row.names = NULL)

sum_row <- colSums(MAGs_AROs, na.rm = TRUE)

MAGs_AROs$category <- ARG_class$Drug.Class[match(row.names(MAGs_AROs), ARG_class$ARO.Accession)]

sum_df_ARG <- aggregate(. ~ category, MAGs_AROs, sum)
rownames(sum_df_ARG) <- sum_df_ARG[,1]
sum_df_ARG<- sum_df_ARG[,-1]

MAGs_AROs1 <- rbind(MAGs_AROs, sum_row)

rownames(MAGs_AROs1)[nrow(MAGs_AROs1)] <- "ARG"


column_ha = HeatmapAnnotation(df = annotate_table)
row_ha = rowAnnotation(df = annotate_table)

library(circlize)
MAGs_AROs <- MAGs_AROs[,colnames(MAGs_AROs) != "category"]
col_fun = colorRamp2(c(0,max(MAGs_AROs)), c("white", "#4DAF4A"))
col_fun(seq(-3, 3))

plot_AROs<-Heatmap(t(MAGs_AROs), show_row_names=F,col=col_fun,
                   row_split = annotate_table,
                   left_annotation = row_ha ,
                   cluster_columns = TRUE
)
plot_AROs

MAGs_AROs1 <- rbind(MAGs_AROs, sum_row)

rownames(MAGs_AROs1)[nrow(MAGs_AROs1)] <- "ARG"


column_ha = HeatmapAnnotation(df = annotate_table)
row_ha = rowAnnotation(df = annotate_table)
col_fun = colorRamp2(c(0,max(t(MAGs_AROs1["ARG",]))), c("white", "#4DAF4A"))
col_fun(seq(-3, 3))

plot_AROs1<-Heatmap(t(MAGs_AROs1["ARG",]), show_row_names=F,col=col_fun,
                    row_split = annotate_table,
                    cluster_columns = TRUE
)

ARO_sum <- t(MAGs_AROs1["ARG",])
ARO_sum_adj <- merge(ARO_sum,gene_number_adj,by= "row.names")
rownames(ARO_sum_adj) <- ARO_sum_adj$Row.names
ARO_sum_adj$Row.names <- NULL
ARO_sum_adj <- merge(ARO_sum_adj,annotate_table,by= "row.names")
rownames(ARO_sum_adj) <- ARO_sum_adj$Row.names
ARO_sum_adj$Row.names <- NULL

ARO_sum_adj <-  ARO_sum_adj %>%
  mutate(ARG1 = (ARG*1000/number) )
ARO_heat <- ARO_sum_adj[, "ARG1", drop = FALSE]
annotate_table1 <- annotate_table[match(rownames(ARO_heat), rownames(annotate_table)), ]
col_fun = colorRamp2(c(0,max(ARO_heat)), c("white", "#4DAF4A"))
col_fun(seq(-3, 3))

MAG_phy_forheat <- MAG_phy_forheat[rownames(ARO_heat)]

MAG_phy_ano <- rowAnnotation(
  Phylum = MAG_phy_forheat,
  col = list(Phylum = color_mapping),
  width = unit(1, "cm") 
)

plot_AROs1<-Heatmap(ARO_heat, show_row_names=F,col=col_fun,
                    row_split = annotate_table1,
                    cluster_columns = TRUE
)
plot_AROs1

ARO_data <- t(MAGs_AROs1["ARG",])
### class ARGs
col_fun = colorRamp2(c(0,max(sum_df_ARG)), c("white", "#4DAF4A"))
col_fun(seq(-3, 3))

plot_AROs2<-Heatmap(t(sum_df_ARG), show_row_names=F,col=col_fun,
                    row_split = annotate_table,
                    cluster_columns = TRUE,
                    column_names_rot = -45
)
plot_AROs2

plot_AMP1 + plot_AROs1 + plot_VFs1

#### CAZy
#data <- read.table("MAGs_Func/result/cazy_all.txt", header = T, sep = "\t", row.names = 1)

#MAGs_CAZys <- data


annotate_table <- data.frame(row.names = module_MAGs$name, groups=module_MAGs$X__leidenCluster)
annotate_table$groups <- factor(annotate_table$groups)

annotation_row_names <- rownames(annotate_table)

MAGs_CAZys_colnames <- colnames(MAGs_CAZys)

missing_colnames <- setdiff(annotation_row_names, MAGs_CAZys_colnames)

for (colname in missing_colnames) {
  MAGs_CAZys[, colname] <- 0
}

# 以annotation_table的行名顺序重新排序MAGs_CAZys的列
MAGs_CAZys <- MAGs_CAZys[, annotation_row_names]

### add classification

category_df <- data.frame(
  gene = c("CE1", "CE2", "CE4", "CE6", "CE7", "GH10", "GH11", "GH115", "GH43", "GH51", "GH67", "GH3", "GH5",
           "GH1", "GH44", "GH48", "GH8", "GH9", "GH3", "GH5", "GH32", "GH91",
           "GH1", "GH2", "GH3", "GH4", "GH18", "GH19", "GH20", "GH29", "GH33", "GH38", "GH58", "GH79", "GH84",
           "GH85", "GH88", "GH89", "GH92", "GH95", "GH98", "GH99", "GH101", "GH105", "GH109", "GH110", "GH113",
           "PL6", "PL8", "PL12", "PL13", "PL21", "CE12", "CE8", "GH28", "PL1", "PL9", "GH13", "GH31", "GH97"), # 我在这里添加了缺失的基因
  category = rep(c("Arabinoxylan-related", "cellulose-related", "inulin-related",
                   "mucin-related", "pectin-related", "starch-related"), 
                 times = c(13, 7, 2, 30, 5, 3))
)

MAGs_CAZys$category <- category_df$category[match(row.names(MAGs_CAZys), category_df$gene)]

sum_df <- aggregate(. ~ category, MAGs_CAZys, sum)
rownames(sum_df) <- sum_df[,1]
sum_df<- sum_df[,-1]

percent_df <- apply(sum_df, MARGIN = 2, function(x) x*100 / sum(x))

percent_df[is.nan(percent_df)] <- 0

column_ha = HeatmapAnnotation(df = annotate_table)
row_ha = rowAnnotation(df = annotate_table)

library(circlize)
col_fun = colorRamp2(c(0,10), c("white", "#FF7F00"))
col_fun(seq(-3, 3))
MAGs_CAZys <- MAGs_CAZys[,colnames(MAGs_CAZys) != "category"]
annotate_table1 <- annotate_table[match(rownames(t(MAGs_CAZys)), rownames(annotate_table)), ]
plot_CAZys<-Heatmap(t(MAGs_CAZys), show_row_names=F,col=col_fun,
                    row_split = annotate_table1,
                    left_annotation = row_ha ,
                    cluster_columns = TRUE
)
plot_CAZys

col_fun = colorRamp2(c(0,max(percent_df)), c("white", "#FF7F00"))
col_fun(seq(-3, 3))

annotate_table1 <- annotate_table[match(rownames(t(percent_df)), rownames(annotate_table)), ]

plot_CAZys1<-Heatmap(t(percent_df), show_row_names=F,col=col_fun,
                     row_split = annotate_table1,
                     cluster_columns = TRUE,
                     column_names_rot = -45,
                     row_names_gp = gpar(fontsize = 4)
)
plot_AMP1  + plot_VFs1 + plot_AROs1 + plot_CAZys1
plot_CAZys1 + plot_AMP1 + plot_AROs1 +  plot_VFs1

plot_AROs1 + plot_AMP1 + plot_VFs1 + plot_CAZys1


data_anno <- data.frame( VFs =  t(MAGs_VFs1["VFs",]), AROs =  t(MAGs_AROs1["ARG",]), AMP = t(MAGs_AMP1["AMP",]) )


COM_Fun_table <- merge(VF_heat, AMP_heat, by = "row.names")
rownames(COM_Fun_table) <- COM_Fun_table$Row.names
COM_Fun_table$Row.names <- NULL
COM_Fun_table <- merge(COM_Fun_table, ARO_heat, by = "row.names")
rownames(COM_Fun_table) <- COM_Fun_table$Row.names
COM_Fun_table$Row.names <- NULL

COM_Fun_table <- merge(COM_Fun_table,abun_df,by="row.names")
row.names(COM_Fun_table) <- COM_Fun_table$Row.names
COM_Fun_table$Row.names <- NULL

COM_Fun_table <- merge(COM_Fun_table,annotate_table,by="row.names")
row.names(COM_Fun_table) <- COM_Fun_table$Row.names
COM_Fun_table$Row.names <- NULL
COM_Fun_table1 <- merge(COM_Fun_table,module_MAGs,by.x="row.names",by.y="Group_cluster")
row.names(COM_Fun_table1) <- COM_Fun_table1$Row.names

COM_Fun_table1 <- COM_Fun_table1 %>%
  mutate(groups = case_when(
    groups == "2" ~ "2",
    groups == "3" ~ "1",
    TRUE ~ as.character(groups)
  ))
## Cazy merge
CAZy_heat <- t(percent_df)
COM_Fun_table1 <- merge(COM_Fun_table1,CAZy_heat,by="row.names")
COM_Fun_table1$Row.names <- NULL


## phylo 
distance_matrix
COM_Fun_table1$Phylum <- as.factor(COM_Fun_table1$Phylum)
COM_Fun_table1$groups <- as.factor(COM_Fun_table1$groups)
Func_cato <- c("VFs1", "AMP1","ARG1")
Func_cato1 <- c("VFs1", "AMP1","ARG1","`Arabinoxylan-related`","`cellulose-related`","`inulin-related`","`mucin-related`","`pectin-related`","`starch-related`")

model_formula <- as.formula(paste("`starch-related`", " ~ groups + abun_binary + Phylum + (1 | Row.names)"))
model_formula0 <- as.formula(paste("`starch-related`", " ~ abun_binary + Phylum  +(1 | Row.names) "))

model0 <- relmatLmer(model_formula0, COM_Fun_table1, relmat = list(ID = distance_matrix), REML = FALSE)
model <- relmatLmer(model_formula,COM_Fun_table1, relmat = list(ID = distance_matrix), REML = FALSE)
lrtest(model,model0)
summa <- summary(model)
coef_trait <- summa$coefficients[[i,"Estimate"]]
t_value_trait <- summa$coefficients[[i,"t value"]]
se_trait <- summa$coefficients[[i,"Std. Error"]]

summa_p <- lrtest(model,model0)
p_value_trait <- summa_p$`Pr(>Chisq)`[2]

results <- data.frame(
  Name = character(), 
  p_value = numeric(), 
  adj_p_value = numeric(), # New column for adjusted p-values
  stringsAsFactors = FALSE
)
for (p in Func_cato1) {
  temp_results <- data.frame() # Temporary data frame for each phenotype
  data_model <- drop_na(bac_net_tra)
  
  model_formula <- as.formula(paste(p, " ~ groups + abun_binary + Phylum + (1 | Row.names)"))
  model_formula0 <- as.formula(paste(p, " ~ abun_binary + Phylum  +(1 | Row.names) "))
  
  # Using tryCatch to handle errors
  result <- tryCatch({
    model0 <- relmatLmer(model_formula0, COM_Fun_table1, relmat = list(ID = distance_matrix), REML = FALSE)
    model <- relmatLmer(model_formula,COM_Fun_table1, relmat = list(ID = distance_matrix), REML = FALSE)
    
    summa_p <- lrtest(model,model0)
    p_value_trait <- summa_p$`Pr(>Chisq)`[2]
    
    # Create a data frame with the results
    temp_results_current <- data.frame(
      Name = p, 
      p_value = p_value_trait, 
      adj_p_value = NA, # Initialize with NA
      stringsAsFactors = FALSE
    )
    temp_results_current
    
  }, error = function(e) {
    # Return NA values in case of an error
    return(data.frame(
      Name = p, 
      p_value = NA, 
      adj_p_value = NA,
      stringsAsFactors = FALSE
    ))
  })
  
  temp_results <- rbind(temp_results, result)
  
  temp_results$adj_p_value <- p.adjust(temp_results$p_value, method = "bonferroni")
  
  results <- rbind(results, temp_results)
}
results$adj_p_value <- p.adjust(results$p_value, method = "fdr")





#### 6 Mobile Elements transmission --Figure6 ####

#hgt_elements <- read.table("ME_connect/ME_connect_strict_Final.csv",sep=",",header=T,row.names = 1)
#hgt_elements[hgt_elements > 0] <- 1
#na_proportion <- rowMeans(is.na(hgt_elements))
#hgt_elements <- hgt_elements[na_proportion <= 0.6, ]
#hgt_elements_info <- read.table("HGT_events/elements_info.csv",sep=",",header=T,row.names = 1)

hgt_elements_info <- hgt_elements_info %>% 
  mutate(row_name = row.names(.)) %>% 
  filter(row_name %in% row.names(hgt_elements)) %>%
  select(-row_name)

colnames(hgt_elements) <- gsub("^X", "", colnames(hgt_elements))
colnames(hgt_elements) <- gsub("LLD_", "", colnames(hgt_elements))

#metadata <- read.table("/Users/mac/Desktop/Fu\ Group/1--Project1/2--DAG3_HGT/metadata.txt",header = T)

cleaned_expression_matrix <- hgt_elements

#phenotype_full <- read.table("/Users/mac/Desktop/Fu\ Group/1--Project1/2--DAG3_HGT/Phenotype/data_pheno_LLD_base_fup_338pairs_62pheno_18med_17disease_log_min_10participants.txt")

LLD_colnames <- sapply(colnames(cleaned_expression_matrix), function(x) {
  metadata$LLD_GoNL_all_id[which(metadata$LLD_id == x)]
})

LLD_FU_colnames <- sapply(colnames(cleaned_expression_matrix), function(x) {
  metadata$LLD_GoNL_all_id[which(metadata$LLD_followup_id == x)]
})

HGT_seq_all <- cleaned_expression_matrix

colnames(HGT_seq_all) <- LLD_FU_colnames

HGT_seq_FU <- HGT_seq_all[,colnames(HGT_seq_all) %in% metadata$LLD_GoNL_all_id]

HGT_seq_all <- cleaned_expression_matrix

colnames(HGT_seq_all) <- LLD_colnames

HGT_seq_BL <- HGT_seq_all[,colnames(HGT_seq_all) %in% metadata$LLD_GoNL_all_id]

colnames(HGT_seq_FU) <- paste0(colnames(HGT_seq_FU),"_F")

HGT_seq_BL$row_name <- rownames(HGT_seq_BL)
HGT_seq_FU$row_name <- rownames(HGT_seq_FU)

HGT_seq_full <- merge(HGT_seq_BL, HGT_seq_FU, by = "row_name", all = TRUE)
rownames(HGT_seq_full) <- HGT_seq_full$row_name
HGT_seq_full<-HGT_seq_full[,-1]

HGT_seq_FU$row_name <- NULL
HGT_seq_BL$row_name <- NULL

## Statistics & Heatmap
HGT_seq_FU1 <- HGT_seq_FU
colnames(HGT_seq_FU1) <- gsub("_F", "", colnames(HGT_seq_FU1))

HGT_seq_FU1 <- HGT_seq_FU1[, colnames(HGT_seq_BL)]

result_df <- data.frame(gene = rownames(HGT_seq_BL), u_statistic = NA, p_value = NA)

for (i in seq_len(nrow(HGT_seq_BL))) {
  gene_BL <- as.numeric(HGT_seq_BL[i, ])
  gene_FU <- as.numeric(HGT_seq_FU1[i, ])
  
  valid_indices <- !is.na(gene_BL) & !is.na(gene_FU)
  gene_BL <- gene_BL[valid_indices]
  gene_FU <- gene_FU[valid_indices]
  
  if(length(gene_BL) > 0 && length(gene_FU) > 0) {
    mw_test_result <- wilcox.test(gene_BL, gene_FU, paired = TRUE)
    result_df[i, "u_statistic"] <- mw_test_result$statistic
    result_df[i, "p_value"] <- mw_test_result$p.value
  } else {
    result_df[i, "u_statistic"] <- NA
    result_df[i, "p_value"] <- NA
  }
}

result_df$p_value_adj <- p.adjust(result_df$p_value, method = "BH")
result_df_BH <- subset(result_df, result_df$p_value_adj <= 0.05)

HGT_seq_BL_Filter <- subset(HGT_seq_BL, rownames(HGT_seq_BL) %in% result_df_BH$gene)
HGT_seq_FU_Filter <- subset(HGT_seq_FU1, rownames(HGT_seq_FU1) %in% result_df_BH$gene)

difference_matrix <- HGT_seq_FU1 - HGT_seq_BL
filtered_matrix <- difference_matrix
filtered_matrix[is.na(filtered_matrix)] <- 0.01

## reads compare

sample_select <- sub("_F","",colnames(filtered_matrix))

calculate_delta_pheno <- function(data_df) {
  samples <- rownames(data_df)
  paired_samples <- samples[grepl("_F$", samples)]
  base_samples <- gsub("_F$", "", paired_samples)
  
  result <- list()
  for (i in seq_along(base_samples)) {
    base_row <- data_df[base_samples[i], , drop = FALSE]
    paired_row <- data_df[paired_samples[i], , drop = FALSE]
    
    delta <- numeric(length(base_row))
    for (j in seq_along(base_row)) {
      if (is.numeric(base_row[1, j]) && is.numeric(paired_row[1, j])) {
        # 0/1 check
        if(all(base_row[1, j] %in% c(0, 1)) && all(paired_row[1, j] %in% c(0, 1))) {
          delta[j] <- paired_row[1, j] + base_row[1, j] 
        } else {
          delta[j] <- paired_row[1, j] - base_row[1, j] 
        }
      } else {
        delta[j] <- NA
      }
    }
    result[[base_samples[i]]] <- delta
  }
  
  delta_df <- do.call(rbind, result)
  delta_df <- as.data.frame(delta_df)
  rownames(delta_df) <- base_samples
  colnames(delta_df) <- colnames(data_df)
  
  delta_df <- delta_df[, !apply(is.na(delta_df), 2, all)]
  
  return(delta_df)
}

delta_phenotype <- calculate_delta_pheno(phenotype_full)

delta_phenotype$SampleType <- ifelse(rownames(delta_phenotype) %in% sample_select, "Selected", "Not Selected")
data_for_plot <- delta_phenotype[, c("clean_reads", "SampleType")]

## Follow filter HGT
difference_matrix <- filtered_matrix

significant_genes <- result_df[result_df$p_value_adj <= 0.01, ]$gene

annotation_matrix <- matrix("", nrow = nrow(difference_matrix), ncol = 1)
rownames(annotation_matrix) <- rownames(difference_matrix)

annotation_matrix[rownames(annotation_matrix) %in% significant_genes, 1] <- "p_adj <= 0.05"


star_anno <- rowAnnotation(df = data.frame(Stars = annotation_matrix),
                           col = list(Stars = c("p_adj <= 0.05" = "black")),
                           show_annotation_name = FALSE,
                           simple_anno_size = unit(2, "mm"))

annotation_matrix <- matrix("", nrow = nrow(difference_matrix), ncol = 2)
rownames(annotation_matrix) <- rownames(difference_matrix)

annotation_matrix <- merge(annotation_matrix, elements_info_unique[, c("HGT_ID", "cluster", "MAG_info")], by.x = "row.names", by.y = "HGT_ID", all.x = TRUE)

set1_colors <- brewer.pal(9, "Set1")
set2_colors <- brewer.pal(8, "Set2")

cluster_colors <- colorRampPalette(set1_colors)(length(unique(annotation_matrix$cluster)))
mag_info_colors <- colorRampPalette(set2_colors)(length(unique(annotation_matrix$MAG_info)))

cluster_color_mapping <- setNames(cluster_colors, unique(annotation_matrix$cluster))
mag_info_color_mapping <- setNames(mag_info_colors, unique(factor(annotation_matrix$MAG_info)))

cluster_category_map <- setNames(annotation_matrix$cluster, annotation_matrix$Row.names)
cluster_categories_for_heatmap <- cluster_category_map[rownames(difference_matrix)]
mag_category_map <- setNames(annotation_matrix$MAG_info, annotation_matrix$Row.names)
mag_categories_for_heatmap <- mag_category_map[rownames(difference_matrix)]

trait_colors <- colorRampPalette(set1_colors)(length(unique(trait_cato$category)))

bac_color_mapping <- setNames(bac_colors, unique(factor(bac_cato$category)))

colnames(difference_matrix) <- gsub("_F","",colnames(difference_matrix))

map_drug <- function(drug) {
  category_map <- setNames(delta_phenotype[[drug]], row.names(delta_phenotype))
  drug_categories_for_heatmap <- category_map[colnames(difference_matrix)]
  return(drug_categories_for_heatmap)
}

PPI <- map_drug("PPI")
angII_receptor_antagonist <- map_drug("angII_receptor_antagonist")
beta_blockers <- map_drug("beta_blockers")
anti_histamine <- map_drug("anti_histamine")
oral_contraceptive <- map_drug("oral_contraceptive")
statin <- map_drug("statin")
ACE_inhibitor <- map_drug("ACE_inhibitor")
#angII_receptor_antagonist <- map_drug("angII_receptor_antagonist")
platelet_aggregation_inhibitor <- map_drug("platelet_aggregation_inhibitor")
clean_reads <- map_drug("clean_reads")
clean_reads_1 <- clean_reads/10000000
category_colors <- c("0" = "white", "1" = "grey", "2" = "black")
reads_colors_category <- colorRamp2(c(min(clean_reads_1, na.rm = TRUE), 0, max(clean_reads_1, na.rm = TRUE)), 
                                    c("#6C96CC", "#EAEFF6", "#C92321"))
dend = as.dendrogram(hclust(dist(t(difference_matrix))))
dend = color_branches(dend, 
                      k = 3, 
                      col = c('#8491B4FF','#F39B7FFF','grey'))
#col = c('#F9C3BF',"#BFE4EE","grey"))
Heatmap(
  difference_matrix,
  name = "Change",
  show_row_names = F,
  cluster_columns = dend,
  column_split = 3,row_names_side = "right",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 4),
  show_column_names = F,
  #"#8491B4FF","#91D1C2FF","#F39B7FFF"
  col = c("-1" = "#8491B4FF", "0" = "white", "1" = "#F39B7FFF","0.01"="grey90"),
  right_annotation = rowAnnotation(
    df = data.frame(
      #Cluster = cluster_categories_for_heatmap,
      Species_info = mag_categories_for_heatmap
    ),
    col = list(
      #Cluster = cluster_color_mapping,
      Species_info = mag_info_color_mapping
    ),
    show_annotation_name = F,
    simple_anno_size = unit(2, "mm"),
    show_legend = F,
    annotation_name_gp = gpar(fontsize = 8)
  ))

Heatmap(
  difference_matrix,
  name = "Change",
  show_row_names = F,
  cluster_columns = dend,
  column_split = 3,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 4),
  show_column_names = F,
  #"#8491B4FF","#91D1C2FF","#F39B7FFF"
  col = c("-1" = "#8491B4FF", "0" = "white", "1" = "#F39B7FFF","0.01"="grey90"),
  top_annotation = HeatmapAnnotation(
    PPI = PPI, 
    angII_receptor_antagonist = angII_receptor_antagonist, 
    beta_blockers = beta_blockers,
    anti_histamine = anti_histamine,
    oral_contraceptive = oral_contraceptive,
    statin = statin,
    ACE_inhibitor = ACE_inhibitor,
    platelet_aggregation_inhibitor = platelet_aggregation_inhibitor,
    #clean_reads = clean_reads_1,
    col = list(
      PPI = category_colors,
      angII_receptor_antagonist = category_colors,
      beta_blockers = category_colors,
      anti_histamine = category_colors,
      oral_contraceptive = category_colors,
      statin = category_colors,
      ACE_inhibitor = category_colors,
      platelet_aggregation_inhibitor = category_colors#,
      #clean_reads = reads_colors_category
    ),
    show_annotation_name = TRUE,
    height = unit(0.1, "cm"),
    simple_anno_size = unit(1, "mm"),
    show_legend = F,
    annotation_name_gp = gpar(fontsize = 6)
  ),
  bottom_annotation = HeatmapAnnotation(
    clean_reads = clean_reads_1,
    col = list(clean_reads = reads_colors_category),
    show_annotation_name = F,
    show_legend = F,
    height = unit(0.1, "cm"),
    simple_anno_size = unit(2, "mm"),
    annotation_name_gp = gpar(fontsize = 6)
  ),
  right_annotation = rowAnnotation(
    df = data.frame(
      #Cluster = cluster_categories_for_heatmap,
      Species_info = mag_categories_for_heatmap
    ),
    col = list(
      #Cluster = cluster_color_mapping,
      Species_info = mag_info_color_mapping
    ),
    show_annotation_name = F,
    simple_anno_size = unit(2, "mm"),
    show_legend = T,
    annotation_name_gp = gpar(fontsize = 8)
  )
)
clusters <- cutree(dend , k=3)

group_all<-as.data.frame(clusters)

## strainphlan + anpan | adjust strain transmission

library(ape)

#tree <- read.tree("/Users/mac/Desktop/Fu\ Group/1--Project1/2--DAG3_HGT/1--Recent_HGT/strainphlan/112/RAxML_bestTree.t__SGB4933_group.StrainPhlAn4.tre")
metadata <- group_all
metadata$Time <- "Baseline"
metadata1 <- metadata
row.names(metadata1) <- paste0(row.names(metadata1), "_F")
metadata1$Time <- "Followup"

metadata_new <- rbind(metadata,metadata1)
metadata_new$ID <- row.names(metadata_new)

df <- metadata_new %>%
  filter(!clusters %in% c("1"))  %>%
  mutate(type = case_when(
    clusters == 2 & Time == "Baseline" ~ "strain1",
    clusters == 2 & Time == "Followup" ~ "strain2",
    clusters == 3 & Time == "Baseline" ~ "strain2",
    clusters == 3 & Time == "Followup" ~ "strain1",
    TRUE ~ NA_character_ # Default case if none of the above conditions are met
  ))
df <- df %>%
  filter(!ID %in% c("LLDeep_1055","LLDeep_1055_F",
                    "LLDeep_0431","LLDeep_0431_F",
                    "LLDeep_1110","LLDeep_1110_F",
                    "LLDeep_0784","LLDeep_0784_F",
                    "LLDeep_0319","LLDeep_0319_F",
                    "LLDeep_0961","LLDeep_0961_F"))
filtered_tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% df$ID])

has_followup <- sapply(filtered_tree$tip.label, function(id) {
  id_base <- gsub("_F$", "", id)
  id_followup <- paste0(id_base, "_F")
  id_base %in% filtered_tree$tip.label & id_followup %in% filtered_tree$tip.label
})

filtered_tree <- drop.tip(filtered_tree, filtered_tree$tip.label[!has_followup])

df_filtered <- df %>% filter(ID %in% filtered_tree$tip.label)

dist_matrix <- cophenetic(filtered_tree)

remove_close_pairs <- function(tree, df, threshold) {
  to_remove <- character(0) 
  for (id in df$ID) {
    if (grepl("_F$", id)) {
      
      id_base <- gsub("_F$", "", id)
      if (id_base %in% tree$tip.label) {
        
        distance <- dist_matrix[id, id_base]
        if (distance < threshold) {
          to_remove <- c(to_remove, id, id_base)
        }
      }
    }
  }
  tree <- drop.tip(tree, to_remove)
  return(tree)
}

threshold <- 0.1

df_filtered$sample_id <- df_filtered$ID
df_filtered$outcome <- as.factor(df_filtered$clusters)
df_filtered$Time <- as.factor(df_filtered$Time)
df_filtered$type <- as.factor(df_filtered$type)

filtered_tree1 <- remove_close_pairs(filtered_tree, df_filtered, threshold)

df_filtered1 <- df_filtered %>% filter(ID %in% filtered_tree1$tip.label)

followup_df1 <- df_filtered1 %>% filter(Time == "Followup")

followup_tree1 <- drop.tip(filtered_tree1, filtered_tree1$tip.label[!filtered_tree1$tip.label %in% followup_df1$ID])

baseline_df1 <- df_filtered1 %>% filter(Time == "Baseline")

baseline_tree1 <- drop.tip(filtered_tree1, filtered_tree1$tip.label[!filtered_tree1$tip.label %in% baseline_df1$ID])

library(ggstar)
library(ggtreeExtra)

ggtree_all <- ggtree(filtered_tree) + geom_fruit(data = df_filtered, geom = geom_star, mapping = aes(y =`ID`,fill = type,size=10, starshape = Time),position = "identity", starstroke = 0.1)+  
  #scale_size_continuous(range = c(1, 3), guide = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(starshape = 15), order = 2)) + 
  scale_fill_manual(values = c("#8dd3c7", "#ffed6f")) +    
  scale_starshape_manual(values = c(1, 1), guide = "none") +
  geom_tiplab(align = F, linetype = 3, linesize = 0.5, size = 4, offset = 0.005,color = "black") +
  theme(legend.position = "top")
ggtree_all

ggtree_filter <- ggtree(filtered_tree1) + geom_fruit(data = df_filtered1, geom = geom_star, mapping = aes(y =`ID`,fill = type,size=10, starshape = Time),position = "identity", starstroke = 0.1)+  
  #scale_size_continuous(range = c(1, 3), guide = guide_legend(keywidth = 0.5, keyheight = 0.5, override.aes = list(starshape = 15), order = 2)) + 
  scale_fill_manual(values = c("#8dd3c7", "#ffed6f")) +     
  scale_starshape_manual(values = c(1, 1), guide = "none") +
  geom_tiplab(align = F, linetype = 3, linesize = 0.5, size = 4, offset = 0.005,color = "black") +
  theme(legend.position = "top")
ggtree_filter

ggtree_all <- as.grob(ggtree_all)
ggtree_filter <- as.grob(ggtree_filter)

layout <- rbind(c(1,2))
grid.arrange(ggtree_all, ggtree_filter,layout_matrix = layout)


## ME transmission

HGT_seq_BL_t <-as.data.frame(HGT_seq_BL)
HGT_seq_BL_t$prevalence <- apply(HGT_seq_BL_t, 1, function(row) {
  non_na_row <- row[!is.na(row)]
  sum(non_na_row != 0) / length(row) * 100
})

histogram(HGT_seq_BL_t$prevalence)
HGT_seq_BL_t_df <- data.frame(name = rownames(HGT_seq_BL_t), prevalence = HGT_seq_BL_t$prevalence, stringsAsFactors = FALSE)

histogram_BL <- ggplot(HGT_seq_BL_t_df, aes(x=prevalence)) +
  geom_histogram(binwidth = 3, color = "black", fill = "#8491B4FF") +
  labs(
    #title = "Histogram of Prevalence 5 Years Later",
    x = "Prevalence",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  ) #+ scale_x_continuous(limits = c(-5, 80)) + scale_y_continuous(limits = c(0,65))
histogram_BL

HGT_seq_FU_t <-as.data.frame(HGT_seq_FU)
HGT_seq_FU_t$prevalence <- apply(HGT_seq_FU_t, 1, function(row) {
  non_na_row <- row[!is.na(row)]
  sum(non_na_row != 0) / length(row) * 100
})

#HGT_seq_FU_t$prevalence <- apply(HGT_seq_FU_t, 1, function(row) sum(row != 0) / length(row) * 100)

# Convert list to data.frame
HGT_seq_FU_t_df <- data.frame(name = rownames(HGT_seq_FU_t), prevalence = HGT_seq_FU_t$prevalence, stringsAsFactors = FALSE)

histogram_FU <- ggplot(HGT_seq_FU_t_df, aes(x=prevalence)) +
  geom_histogram(binwidth = 3, color = "black", fill = "#91D1C2FF") + 
  labs(
    #title = "Histogram of Prevalence 5 Years Later",
    x = "Prevalence",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  )#+ scale_x_continuous(limits = c(-5, 80)) + scale_y_continuous(limits = c(0,65))

histogram_FU

prevalence_plot <- rbind(
  data.frame(prevalence = HGT_seq_BL_t$prevalence, Timepoint = "Baseline"),
  data.frame(prevalence = HGT_seq_FU_t$prevalence, Timepoint = "Followup")
)

combined_histogram <- ggplot(prevalence_plot, aes(x = prevalence, fill = Timepoint)) +
  geom_histogram(data = subset(prevalence_plot, Timepoint == "Baseline"), binwidth = 5, color = "black", alpha = 0.5) +
  geom_histogram(data = subset(prevalence_plot, Timepoint == "Followup"), binwidth = 5, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("Baseline" = "#8491B4FF", "Followup" = "#91D1C2FF")) +
  labs(x = "Prevalence", y = "Number") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  ) + 
  guides(fill="none", color="none")

combined_histogram

## boxplot
data_plot <- rbind(
  data.frame(contig = rownames(HGT_seq_BL_t),Prevalence = HGT_seq_BL_t$prevalence, Time = "Baseline"),
  data.frame(contig = rownames(HGT_seq_FU_t),Prevalence = HGT_seq_FU_t$prevalence, Time = "Followup")
)

Compare_plot <- ggplot(data_plot, aes(x = Time, y = Prevalence, fill = Time)) + 
  #geom_boxplot(outlier.shape = NA) +
  #geom_boxplot() + 
  geom_violin(width=1,adjust = 2,fill = "white", aes(colour = Time),size=1) +
  geom_boxplot(width=0.4, fill = "white", aes(colour = Time),size=1) +
  #scale_colour_manual(values=c("#AF7EC0","#B0A875")) +
  #scale_fill_manual(values=c("#AF7EC0","#B0A875")) +
  scale_colour_manual(values=c("#8491B4FF","#91D1C2FF")) +
  scale_fill_manual(limits=c("Baseline","Followup"), 
                    values=c("#8491B4FF","#91D1C2FF")) +
  theme_classic() +
  labs( x= NULL ,y = "Prevalence") +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  guides(fill="none", color="none") +
  stat_compare_means(method = "wilcox.test",paired = TRUE, 
                     comparisons=list(c("Baseline", "Followup")))

Compare_plot

data_BL_FU <- data.frame(name_BL = rownames(HGT_seq_BL_t),name_FU = rownames(HGT_seq_FU_t)
                         ,Prevalence_BL = HGT_seq_BL_t$prevalence,Prevalence_FU = HGT_seq_FU_t$prevalence)

wilcox.test(Prevalence~ Time,data=data_plot,, paired=TRUE)

data_BL_FU$delta <- data_BL_FU$Prevalence_FU - data_BL_FU$Prevalence_BL

Delta_prevalence <- ggplot(data_BL_FU, aes(x=delta)) +
  geom_histogram(binwidth = 3, color = "black", fill = "#F39B7FFF", alpha = 0.7) + 
  labs(
    y = "Number",
    x = "Delta Frequency"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  )
Delta_prevalence

## Individual distance

HGT_seq_BL1 <- HGT_seq_BL[, !(names(HGT_seq_BL) %in% c("prevalence","row_name") )]

HGT_seq_FU1 <- HGT_seq_FU[, !(names(HGT_seq_FU) %in% c("prevalence","row_name"))]


jaccard_BL <- vegdist(t(HGT_seq_BL1), method = "jaccard",na.rm = TRUE)

jaccard_FU <- vegdist(t(HGT_seq_FU1), method = "jaccard",na.rm = TRUE)

jaccard_BL <- as.matrix(jaccard_BL)
jaccard_FU <- as.matrix(jaccard_BL)

diag(jaccard_BL) <- NA
diag(jaccard_FU) <- NA

common_order <- paste0(colnames(HGT_seq_BL1),"_F")

HGT_seq_FU1 <- HGT_seq_FU1[, common_order]

paired_distance <- numeric(ncol(HGT_seq_BL1))
for(i in 1:ncol(HGT_seq_BL1)) {
  paired_distance[i] <- vegdist(rbind(HGT_seq_BL1[, i], HGT_seq_FU1[, i]), method = "jaccard", na.rm = TRUE)
}

distance_data <- data.frame(
  Baseline = as.vector(jaccard_BL),
  Followup = as.vector(jaccard_FU),
  Paired = paired_distance
)

compare_means <- list(
  c("Baseline", "Paired"),
  c("Followup", "Paired")
)

group1 <- as.numeric(as.vector(as.matrix(jaccard_BL)))
group1 <- na.omit(group1)
group2 <- as.numeric(as.vector(as.matrix(jaccard_FU)))
group2 <- na.omit(group2)
group3 <- as.numeric(paired_distance)

data <- c(group1, group2)
group <- factor(c(rep("Group1", length(group1)), rep("Group2", length(group2))))
library(coin)
test_result12 <- wilcox_test(data ~ group, distribution = approximate(nresample = 9999))

data <- c(group1, group3)
group <- factor(c(rep("Group1", length(group1)), rep("Group2", length(group3))))
test_result13 <- wilcox_test(data ~ group, distribution = approximate(nresample = 9999))

data <- c(group2, group3)
group <- factor(c(rep("Group1", length(group2)), rep("Group2", length(group3))))
test_result23 <- wilcox_test(data ~ group, distribution = approximate(nresample = 9999))

test_result12
test_result13
test_result23

p <- ggplot(melt(distance_data), aes(x = variable, y = value, fill = variable)) + 
  geom_violin(width=1,adjust = 5,fill = "white", aes(colour = variable),size=1) +
  geom_boxplot(width=0.5, fill = "white", aes(colour = variable),size=1) +
  scale_colour_manual(values=c("#8491B4FF","#91D1C2FF","#F39B7FFF")) +
  theme_classic() +
  labs(x = NULL, y = "jaccard distance") +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  guides(fill="none", color="none")

p_individual <- p + stat_compare_means(comparisons = compare_means, method="wilcox", label = "p.format")
p_individual

layout <- rbind(c(1,1),
                c(3,4))
Delta_prevalence
#CairoPDF(file="HGT_events/plot5/Prevalence_change_elements.pdf", width=10, height=15)
grid.arrange(combined_histogram,Compare_plot,p_individual,layout_matrix = layout)
dev.off()
save.image("HGT_pipeline.RData")

## ----Figure S5 Fingerprint-----

#hgt_elements <- read.table("ME_connect/ME_connect_strict_Final.csv",sep=",",header=T,row.names = 1)
hgt_elements[hgt_elements > 0] <- 1
na_proportion <- rowMeans(is.na(hgt_elements))
hgt_elements <- hgt_elements[na_proportion <= 0.8, ]
#hgt_elements_info <- read.table("HGT_events/elements_info.csv",sep=",",header=T,row.names = 1)

hgt_elements_info <- hgt_elements_info %>% 
  mutate(row_name = row.names(.)) %>% 
  filter(row_name %in% row.names(hgt_elements)) %>%
  select(-row_name)

colnames(hgt_elements) <- gsub("^X", "", colnames(hgt_elements))
colnames(hgt_elements) <- gsub("LLD_", "", colnames(hgt_elements))

#metadata <- read.table("/Users/mac/Desktop/Fu\ Group/1--Project1/2--DAG3_HGT/metadata.txt",header = T)

cleaned_expression_matrix <- hgt_elements

#phenotype_full <- read.table("/Users/mac/Desktop/Fu\ Group/1--Project1/2--DAG3_HGT/Phenotype/data_pheno_LLD_base_fup_338pairs_62pheno_18med_17disease_log_min_10participants.txt")

LLD_colnames <- sapply(colnames(cleaned_expression_matrix), function(x) {
  metadata$LLD_GoNL_all_id[which(metadata$LLD_id == x)]
})

LLD_FU_colnames <- sapply(colnames(cleaned_expression_matrix), function(x) {
  metadata$LLD_GoNL_all_id[which(metadata$LLD_followup_id == x)]
})

HGT_seq_all <- cleaned_expression_matrix

colnames(HGT_seq_all) <- LLD_FU_colnames

HGT_seq_FU <- HGT_seq_all[,colnames(HGT_seq_all) %in% metadata$LLD_GoNL_all_id]

HGT_seq_all <- cleaned_expression_matrix

colnames(HGT_seq_all) <- LLD_colnames

HGT_seq_BL <- HGT_seq_all[,colnames(HGT_seq_all) %in% metadata$LLD_GoNL_all_id]

colnames(HGT_seq_FU) <- paste0(colnames(HGT_seq_FU),"_F")

HGT_seq_BL$row_name <- rownames(HGT_seq_BL)
HGT_seq_FU$row_name <- rownames(HGT_seq_FU)

HGT_seq_full <- merge(HGT_seq_BL, HGT_seq_FU, by = "row_name", all = TRUE)
rownames(HGT_seq_full) <- HGT_seq_full$row_name
HGT_seq_full<-HGT_seq_full[,-1]

HGT_seq_FU$row_name <- NULL
HGT_seq_BL$row_name <- NULL

## Statistics & Heatmap
HGT_seq_FU1 <- HGT_seq_FU
colnames(HGT_seq_FU1) <- gsub("_F", "", colnames(HGT_seq_FU1))

HGT_seq_FU1 <- HGT_seq_FU1[, colnames(HGT_seq_BL)]

na_percentage <- colMeans(is.na(HGT_seq_BL)) * 100
na_percentage1 <- colMeans(is.na(HGT_seq_FU1)) * 100

HGT_seq_BL_filtered <- HGT_seq_BL[, na_percentage <= 70]
HGT_seq_FU1_filtered <- HGT_seq_FU1[, na_percentage1 <= 70]

common_columns <- intersect(colnames(HGT_seq_BL_filtered), colnames(HGT_seq_FU1_filtered))

HGT_seq_BL_filtered <- HGT_seq_BL_filtered[, common_columns]
HGT_seq_FU1_filtered <- HGT_seq_FU1_filtered[, common_columns]

# Step 1: Calculate Jaccard distances for each individual in baseline against all individuals in follow-up
paired_distance_all <- matrix(NA, ncol = ncol(HGT_seq_BL_filtered), nrow = ncol(HGT_seq_BL_filtered))

for(i in 1:ncol(HGT_seq_BL_filtered)) {
  for(j in 1:ncol(HGT_seq_FU1_filtered)) {
    # Calculate Jaccard distance between individual i (BL) and individual j (FU)
    paired_distance_all[i, j] <- vegdist(rbind(HGT_seq_BL_filtered[, i], HGT_seq_FU1_filtered[, j]), method = "bray", na.rm = TRUE)
  }
}

# Step 2: Compare each individual's self-distance with all distances and calculate the percentile
percentile_rank <- numeric(ncol(HGT_seq_BL_filtered))

for(i in 1:ncol(HGT_seq_BL_filtered)) {
  # Get the self distance (distance between the same individual in BL and FU)
  self_distance <- paired_distance_all[i, i]
  
  # Calculate the percentile rank of self-distance compared to all distances of the individual
  all_distances <- paired_distance_all[i, ]  # All distances for individual i
  percentile_rank[i] <- sum(all_distances < self_distance, na.rm = TRUE) / sum(!is.na(all_distances)) * 100
}

# Step 3: Create a table with results
results_table <- data.frame(
  Individual = colnames(HGT_seq_BL_filtered),
  Self_Distance = diag(paired_distance_all),
  Percentile = percentile_rank
)

# Display the table
print(results_table)

ggplot(results_table, aes(x = Percentile)) +
  geom_histogram(binwidth = 1, fill = "#8491B4FF", color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

sum(results_table$Percentile < 5)
sum(results_table$Percentile == 0)
p2 + p1

#### 7 ME presence/ absence | delta Association -- Figure7a-c ####

hgt_elements <- read.table("ME_connect/ME_connect_strict_Final.csv",sep=",",header=T,row.names = 1)
hgt_elements[hgt_elements > 0] <- 1
colnames(hgt_elements) <- gsub("^X", "", colnames(hgt_elements))
colnames(hgt_elements) <- gsub("LLD_", "", colnames(hgt_elements))

metadata <- read.table("/Users/mac/Desktop/Fu\ Group/1--Project1/2--DAG3_HGT/metadata.txt",header = T)

cleaned_expression_matrix <- hgt_elements

#phenotype_full <- read.table("/Users/mac/Desktop/Fu Group/1--Project1/2--DAG3_HGT/Phenotype/data_pheno_LLD_base_fup_338pairs_18med_3life_log_min_10participants.txt")

LLD_colnames <- sapply(colnames(cleaned_expression_matrix), function(x) {
  metadata$LLD_GoNL_all_id[which(metadata$LLD_id == x)]
})

LLD_FU_colnames <- sapply(colnames(cleaned_expression_matrix), function(x) {
  metadata$LLD_GoNL_all_id[which(metadata$LLD_followup_id == x)]
})

HGT_seq_all <- cleaned_expression_matrix

colnames(HGT_seq_all) <- LLD_FU_colnames

HGT_seq_FU <- HGT_seq_all[,colnames(HGT_seq_all) %in% metadata$LLD_GoNL_all_id]

HGT_seq_all <- cleaned_expression_matrix

colnames(HGT_seq_all) <- LLD_colnames

HGT_seq_BL <- HGT_seq_all[,colnames(HGT_seq_all) %in% metadata$LLD_GoNL_all_id]

colnames(HGT_seq_FU) <- paste0(colnames(HGT_seq_FU),"_F")

HGT_seq_BL$row_name <- rownames(HGT_seq_BL)
HGT_seq_FU$row_name <- rownames(HGT_seq_FU)

HGT_seq_full <- merge(HGT_seq_BL, HGT_seq_FU, by = "row_name", all = TRUE)
rownames(HGT_seq_full) <- HGT_seq_full$row_name
HGT_seq_full<-HGT_seq_full[,-1]

HGT_seq_association <- HGT_seq_full
phenotype_full <- phenotype_full %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(
    TimePoint = ifelse(str_detect(SampleID, "_F"), "Follow-up", "Baseline"),
    IndividualID = str_replace(SampleID, "_F", "")
  )
rownames(phenotype_full) <- phenotype_full$SampleID 
phenotype_full<-phenotype_full[,-1]

HGT_names <- rownames(HGT_seq_association)
HGT_seq_association <- as.data.frame(t(HGT_seq_association))
phenotype_full <- phenotype_full %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(
    TimePoint = ifelse(str_detect(SampleID, "_F"), "Follow-up", "Baseline"),
    IndividualID = str_replace(SampleID, "_F", "")
  )
rownames(phenotype_full) <- phenotype_full$SampleID 
phenotype_full<-phenotype_full[,-1]
phenotype_names <- colnames(phenotype_full)
phenotype_names <- phenotype_names[!(phenotype_names %in% "SampleID")] # exclude the SampleID column

## delta calculation
calculate_delta <- function(data_df) {
  samples <- rownames(data_df)
  paired_samples <- samples[grepl("_F$", samples)]
  base_samples <- gsub("_F$", "", paired_samples)
  
  result <- list()
  for (i in seq_along(base_samples)) {
    base_row <- data_df[base_samples[i], , drop = FALSE]
    paired_row <- data_df[paired_samples[i], , drop = FALSE]
    
    delta <- numeric(length(base_row))
    for (j in seq_along(base_row)) {
      if (is.na(base_row[1, j]) || is.na(paired_row[1, j])) {
        delta[j] <- NA
      } else {
        #if (base_row[1, j] == 1 && paired_row[1, j] == 1) {
        #if(all(base_row[1, j] %in% c(0, 1)) && all(paired_row[1, j] %in% c(0, 1))) {
        delta[j] <- paired_row[1, j] - base_row[1, j]
        #} else {
        #  delta[j] <- paired_row[1, j] - base_row[1, j]
        #}
      }
    }
    result[[base_samples[i]]] <- delta
  }
  
  delta_df <- do.call(rbind, result)
  delta_df <- as.data.frame(delta_df)
  rownames(delta_df) <- base_samples
  colnames(delta_df) <- colnames(data_df)
  
  delta_df <- delta_df[, !apply(is.na(delta_df), 2, all)]
  
  return(delta_df)
}

delta_HGT <- calculate_delta(HGT_seq_association)
delta_HGT_t <- as.data.frame(t(delta_HGT))

delta_HGT_t1 <- delta_HGT_t[apply(delta_HGT_t, 1, function(row) sum( !is.na(row)) >= 10), ]
delta_HGT <- as.data.frame(t(delta_HGT_t1))

data_df <- phenotype_full
calculate_delta_pheno <- function(data_df) {
  samples <- rownames(data_df)
  paired_samples <- samples[grepl("_F$", samples)]
  base_samples <- gsub("_F$", "", paired_samples)
  
  result <- list()
  for (i in seq_along(base_samples)) {
    base_row <- data_df[base_samples[i], , drop = FALSE]
    paired_row <- data_df[paired_samples[i], , drop = FALSE]
    
    delta <- numeric(length(base_row))
    for (j in seq_along(base_row)) {
      if (is.numeric(base_row[1, j]) && is.numeric(paired_row[1, j])) {
        # 0/1 check
        if(all(base_row[1, j] %in% c(0, 1)) && all(paired_row[1, j] %in% c(0, 1))) {
          if (paired_row[1, j] == 1 && base_row[1, j] == 1) {
            delta[j] <- paired_row[1, j] + base_row[1, j] 
          }
          else {
            delta[j] <- paired_row[1, j] - base_row[1, j] 
          }
        } else {
          delta[j] <- paired_row[1, j] - base_row[1, j] 
        }
      } else {
        delta[j] <- NA
      }
    }
    result[[base_samples[i]]] <- delta
  }
  
  delta_df <- do.call(rbind, result)
  delta_df <- as.data.frame(delta_df)
  rownames(delta_df) <- base_samples
  colnames(delta_df) <- colnames(data_df)
  
  delta_df <- delta_df[, !apply(is.na(delta_df), 2, all)]
  
  return(delta_df)
}

delta_phenotype <- calculate_delta_pheno(phenotype_full)
delta_phenotype_t <- as.data.frame(t(delta_phenotype))
delta_phenotype_t1 <- delta_phenotype_t[apply(delta_phenotype_t, 1, function(row) sum(row != 0 & !is.na(row)) >= 20), ]
delta_phenotype <- as.data.frame(t(delta_phenotype_t1))

phenotype_baseline <- subset(phenotype_full, TimePoint %in% "Baseline")

delta_phenotype$antrop_gender.F1M2 <- phenotype_baseline$antrop_gender.F1M2
delta_phenotype$antrop_age <- phenotype_baseline$antrop_age
delta_phenotype$antrop_height <- phenotype_baseline$antrop_height

HGT_names <- colnames(delta_HGT)
phenotype_names <- colnames(delta_phenotype)
phenotype_names <- phenotype_names[!(phenotype_names %in% "SampleID")] # exclude the SampleID column

results <- data.frame(
  HGT_name = character(), 
  phenotype_name = character(), 
  HGT_number =  numeric(),
  noHGT_noPhe = numeric(),
  noHGT_Phe = numeric(),
  HGT_noPhe = numeric(),
  HGT_Phe = numeric(),
  phenotype_number = numeric(),
  intercept = numeric(), 
  coef_HGT = numeric(), 
  r_squared = numeric(), 
  adj_r_squared = numeric(), 
  p_value_HGT = numeric(), 
  se_HGT = numeric(),
  adj_p_value_HGT = numeric(), # New column for adjusted p-values
  stringsAsFactors = FALSE
)

for (p in phenotype_names) {
  if (p %in% c( "antrop_age", "IndividualID","antrop_height", "TimePoint", "clean_reads", "antrop_gender.F1M2")) {
    next
  }
  temp_results <- data.frame() # Temporary data frame for each phenotype
  for (i in HGT_names){
    data_HGT <- dplyr::select(delta_HGT, i)
    data_HGT$SampleID <- rownames(data_HGT)
    delta_phenotype$SampleID <- rownames(delta_phenotype)
    data_model <- left_join(data_HGT, delta_phenotype, by = "SampleID")
    data_model <- drop_na(data_model)
    colnames(data_model)[1] <- "HGT"
    data_model$HGT <- as.numeric(data_model$HGT)
    data_model$SampleID <- NULL
    
    model_formula <- as.formula(paste("HGT ~", p, "+ clean_reads + antrop_gender.F1M2 +antrop_age"))
    
    # Using tryCatch to handle errors
    result <- tryCatch({
      model <- lm(model_formula, data = data_model)
      summa <- summary(model)
      
      intercept <- coef(model)["(Intercept)"]
      coef_HGT <- coef(model)[p]
      r_squared <- summa$r.squared
      adj_r_squared <- summa$adj.r.squared
      p_value_HGT <- coef(summary(model))[[p,"Pr(>|t|)"]]
      se_HGT <- coef(summary(model))[[p,"Std. Error"]]
      
      # Create a data frame with the results
      temp_results_current <- data.frame(
        HGT_name = i, 
        phenotype_name = p, 
        HGT_number =  sum(!is.na(data_model$HGT) & !is.na(data_model$HGT)),
        phenotype_number = sum(!is.na(data_model[p]) & !is.na(data_model[p])),
        noHGT_noPhe = sum(data_model$HGT == 0 & data_model[[p]] == 0, na.rm = TRUE),
        noHGT_Phe = sum(data_model$HGT == 0 & data_model[[p]] != 0, na.rm = TRUE),
        HGT_noPhe = sum(data_model$HGT != 0 & data_model[[p]] == 0, na.rm = TRUE),
        HGT_Phe = sum(data_model$HGT != 0 & data_model[[p]] != 0, na.rm = TRUE),
        intercept = intercept, 
        coef_HGT = coef_HGT, 
        r_squared = r_squared, 
        adj_r_squared = adj_r_squared, 
        p_value_HGT = p_value_HGT, 
        se_HGT = se_HGT,
        adj_p_value_HGT = NA, # Initialize with NA
        stringsAsFactors = FALSE
      )
      temp_results_current
      
    }, error = function(e) {
      # Return NA values in case of an error
      return(data.frame(
        HGT_name = i, 
        phenotype_name = p, 
        HGT_number =  sum(!is.na(data_model$HGT) & !is.na(data_model$HGT)),
        phenotype_number = sum(!is.na(data_model[p]) & !is.na(data_model[p])),
        noHGT_noPhe = sum(data_model$HGT == 0 & data_model[[p]] == 0, na.rm = TRUE),
        noHGT_Phe = sum(data_model$HGT == 0 & data_model[[p]] != 0, na.rm = TRUE),
        HGT_noPhe = sum(data_model$HGT != 0 & data_model[[p]] == 0, na.rm = TRUE),
        HGT_Phe = sum(data_model$HGT != 0 & data_model[[p]] != 0, na.rm = TRUE),
        intercept = NA, 
        coef_HGT = NA, 
        r_squared = NA, 
        adj_r_squared = NA, 
        p_value_HGT = NA, 
        se_HGT = NA,
        adj_p_value_HGT = NA,
        stringsAsFactors = FALSE
      ))
    })
    
    temp_results <- rbind(temp_results, result)
  }
  # Adjust the p-values for the current phenotype
  temp_results$adj_p_value_HGT <- p.adjust(temp_results$p_value_HGT, method = "BH")
  
  # Append the temp_results to the main results data frame
  results <- rbind(results, temp_results)
}

HGT_elements_result <- results
HGT_elements_result$BH_p_value_HGT <- p.adjust(HGT_elements_result$p_value_HGT, method = "BH")
HGT_elements_result$z_value <- -qnorm(HGT_elements_result$p_value_HGT/2)*sign(HGT_elements_result$coef_HGT)
Delta_HGT_elements <- subset(HGT_elements_result, adj_p_value_HGT < 0.05)

write.csv(Delta_HGT_elements,"Association/Delta_HGT_elements_result_binary_disease_0328.csv")

delta_phenotype$SampleType <- NULL
matrix_phe <- cor(delta_phenotype)
matrix_phe[is.na(matrix_phe)] <- 0
pheatmap(matrix_phe)

## Plot association ##
plot_delta_model <- function(i, phenotype_name) {
  data_HGT <- dplyr::select(delta_HGT, i)
  data_HGT$SampleID <- rownames(data_HGT)
  delta_phenotype$SampleID <- rownames(delta_phenotype)
  
  data_model <- left_join(data_HGT, delta_phenotype, by = "SampleID") %>%
    drop_na() %>%
    na.omit()
  data_model[[phenotype_name]] <- as.factor(data_model[[phenotype_name]])
  
  ggplot(data_model, aes_string(x = phenotype_name, y = i, fill = phenotype_name)) + 
    geom_violin(width=0.5,adjust = 1,fill = "white", aes_string(colour = phenotype_name),size=1) +
    geom_boxplot(alpha =0.5,size=0.5,outlier.shape = NA,width = 0.3)+
    scale_color_manual(limits=c("-1","0","1","2"), 
                       values=c("#1a759f","#99d98c","#1a759f","#1a759f"))+
    scale_fill_manual(limits=c("-1","0","1","2"), 
                      values=c("#1a759f","#99d98c","#1a759f","#1a759f")) +
    geom_jitter(alpha = 0.3,size=3, aes_string(fill=phenotype_name),shape=21)+
    theme_classic() + 
    theme(panel.grid =element_blank(),
          axis.text = element_text(size = 10,colour = "black"),
          legend.position = 'top') +
    xlab(phenotype_name) +
    ylab(i)
}

plot_mixed_model <- function(i, phenotype_name) {
  data_HGT <- dplyr::select(HGT_seq_association, i)
  data_HGT$SampleID <- rownames(data_HGT)
  phenotype_full$SampleID <- rownames(phenotype_full)
  
  data_model <- left_join(data_HGT, phenotype_full, by = "SampleID") %>%
    drop_na() %>%
    na.omit()
  #data_model[[i]] <- as.factor(data_model[[i]])
  data_model[[phenotype_name]] <- as.factor(data_model[[phenotype_name]])
  data_model[[i]] <- log10(data_model[[i]])
  
  ggplot(data_model, aes_string(x = phenotype_name, y = i, fill = phenotype_name)) + 
    geom_violin(width=0.5,adjust = 1,fill = "white", aes_string(colour = phenotype_name),size=1) +
    geom_boxplot(alpha =0.5,size=0.5,outlier.shape = NA,width = 0.3)+
    scale_color_manual(limits=c("0","1"), 
                       values=c("#99d98c","#1a759f"))+
    scale_fill_manual(limits=c("0","1"), 
                      values=c("#99d98c","#1a759f")) +
    geom_jitter(alpha = 0.3,size=3, aes_string(fill=phenotype_name),shape=21)+
    theme_classic() + 
    theme(panel.grid =element_blank(),
          axis.text = element_text(size = 10,colour = "black"),
          legend.position = 'top') +
    xlab(phenotype_name) +
    ylab(i)
}


plot_mixed_model("NODE_3_length_167504_cov_12.261256_21426_22997","AlcoholGlassPerDay")

p_delta_PPI <- plot_delta_model("NODE_53_length_16945_cov_28.913440_140_2344_1","PPI")

p_mix_alcohol <- plot_mixed_model("NODE_113_length_5820_cov_10.897485_5160_5714_1","AlcoholGlassPerDay")
p_mix_alcohol1 <- plot_mixed_model("NODE_49_length_7239_cov_15.926096_5467_3718","AlcoholGlassPerDay")


p_mixed_PPI <- plot_mixed_model("NODE_1_length_198559_cov_32.861927_82794_85322_3","PPI")
p_mixed_PPI.1 <- plot_mixed_model("NODE_1_length_198559_cov_32.861927_82794_85322_4","PPI")
p_mixed_PPI1 <- plot_mixed_model("NODE_53_length_16945_cov_28.913440_140_2344_2","PPI")
p_mixed_PPI2 <- plot_mixed_model("NODE_68_length_15908_cov_9.869867_10865_14777_1","PPI")

p_mixed_PPI

layout <- rbind(c(1),
                c(2))

CairoPDF(file="Association/plot/PPI_association_0713.2.pdf", width=5, height=10)
grid.arrange(p_mixed_PPI,p_mixed_PPI1,layout_matrix = layout)
dev.off()

delta_1023_oral_contraceptive <- plot_delta_model("NODE_141_length_5765_cov_8.752539_208_1252","oral_contraceptive")
delta_1023_1_oral_contraceptive <- plot_delta_model("NODE_124_length_4742_cov_8.943461_1206_2250","oral_contraceptive")
delta_275_1_oral_contraceptive <- plot_delta_model("NODE_2_length_59542_cov_12.998537_9181_13794","oral_contraceptive")
delta_275_oral_contraceptive <- plot_delta_model("NODE_1_length_85092_cov_15.429919_23915_19303","oral_contraceptive")
delta_1306_antagonist <- plot_delta_model("NODE_242_length_2336_cov_16.505480_153_791","angII_receptor_antagonist")
delta_732_antagonist <- plot_delta_model("NODE_20_length_23101_cov_9.725853_20719_22362","angII_receptor_antagonist")
delta_732_1_antagonist <- plot_delta_model("NODE_20_length_23101_cov_9.725853_22335_20745","angII_receptor_antagonist")

delta_574_ACE_inhibitor <- plot_delta_model("NODE_114_length_4728_cov_6.051359_2492_4653","ACE_inhibitor")


mixed_118_PPI <- plot_mixed_model("NODE_19_length_61895_cov_29.057956_230_8368","PPI")
mixed_253_PPI <- plot_mixed_model("NODE_23_length_45707_cov_22.926838_510_5504","PPI")
mixed_480_PPI <- plot_mixed_model("NODE_22_length_38093_cov_85.732136_627_3412","PPI")
mixed_1359_oral_contraceptive <- plot_mixed_model("NODE_28_length_18934_cov_7.162456_15325_14738","oral_contraceptive")
mixed_494_beta_blockers <- plot_mixed_model("NODE_91_length_9374_cov_8.055800_26_2759","beta_blockers")
mixed_1002_oral_ACE_inhibitor <- plot_mixed_model("NODE_96_length_9162_cov_5.391457_1637_548","ACE_inhibitor")


layout <- rbind(c(1,2,3,4),
                c(5,6,7,14),
                c(8,9,10,NA),
                c(11,12,13,NA))

CairoPDF(file="Association/plot/all_ratio_association.pdf", width=15, height=15)
grid.arrange(delta_1023_oral_contraceptive,delta_1023_1_oral_contraceptive,delta_275_1_oral_contraceptive,
             delta_275_oral_contraceptive,delta_1306_antagonist,delta_732_antagonist,delta_732_1_antagonist,
             mixed_118_PPI,mixed_253_PPI,mixed_480_PPI,mixed_1359_oral_contraceptive,mixed_494_beta_blockers,mixed_1002_oral_ACE_inhibitor,p_delta_PPI,layout_matrix = layout)
dev.off()

#### 8 HGT function enrich -- Figure7F ####
HGT_table <- fread("Func_enrich/Species_pfam_emapper_filter.tsv", sep="\t", header=TRUE, fill=TRUE)

HGT_table <- HGT_table %>%
  group_by(query) %>%
  filter(sum_score == max(sum_score)) %>%
  ungroup()

HGT_frequency <- HGT_table %>%
  group_by(PFAMs) %>%
  summarize(count_HGT = n()) %>%
  mutate(freq_HGT = count_HGT / nrow(HGT_table)) %>%
  arrange(-freq_HGT)

AllGenes_table <- fread("Func_enrich/all_pfam_emapper.tsv", sep="\t", header=TRUE, fill=TRUE)

AllGenes_table <- AllGenes_table %>%
  group_by(query) %>%
  filter(sum_score == max(sum_score)) %>%
  ungroup()

AllGenes_frequency <- AllGenes_table %>%
  group_by(PFAMs) %>%
  summarize(count_AllGenes = n()) %>%
  mutate(freq_AllGenes = count_AllGenes / nrow(AllGenes_table)) %>%
  arrange(-freq_AllGenes)

FoldChange_data <- merge(HGT_frequency, AllGenes_frequency, by = "PFAMs") %>%
  mutate(FoldChange = freq_HGT / freq_AllGenes) %>%
  arrange(-FoldChange)
FoldChange_data <- FoldChange_data %>%
  mutate(log2_FoldChange = log2(FoldChange + 1)) # +1 为了避免log(0)的情况

top_20 <- FoldChange_data %>%
  arrange(-log2_FoldChange) %>%
  head(50)

ggplot(top_20, aes(x=reorder(PFAMs, -log2_FoldChange), y=log2_FoldChange, fill=log2_FoldChange)) +
  geom_bar(stat="identity", width=0.7) +
  scale_fill_gradientn(colors=c("#1a759f", "#52b69a", "#99d98c")) +
  labs(#title="Top 20 Log2 Fold Change of Functions in HGT", 
    x="PFAM", 
    y="FoldChange") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "none"
  )

top_20_mat <- top_20 %>%
  arrange(-log2_FoldChange) %>%
  select(PFAMs, log2_FoldChange) %>%
  column_to_rownames(var="PFAMs") %>%
  as.matrix()

FoldChange_data <- FoldChange_data %>%
  rowwise() %>%
  mutate(p_value = fisher.test(matrix(c(count_HGT, (nrow(HGT_table) - count_HGT), 
                                        count_AllGenes, (nrow(AllGenes_table) - count_AllGenes)), 
                                      ncol = 2, byrow = TRUE))$p.value)

FoldChange_data$p_adj <- p.adjust(FoldChange_data$p_value,method = "BH")
FoldChange_data <- FoldChange_data %>%
  mutate(neg_log10_pval = -log10(p_adj))


filtered_data <- FoldChange_data %>%
  filter(neg_log10_pval >= 8, abs(log2_FoldChange) >= 3) %>%
  arrange(-log2_FoldChange) %>%
  head(20)

ggplot(filtered_data, aes(x = log2_FoldChange, y = reorder(PFAMs, log2_FoldChange)
                          , size = neg_log10_pval)) +
  geom_point(shape=21, aes(fill=freq_HGT)) +  # 使用shape变量
  scale_fill_gradientn(colors=c("#1a759f", "#52b69a", "#99d98c")) +
  scale_size(range = c(3, 10)) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1),
    plot.title = element_text(size = 18, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right"
  ) +
  labs(x = "Log2FoldChange", y = "Pfam", size = "FDR (-log10)", fill = "Freq_HGT") 

write.csv(FoldChange_data,"Func_enrich/enrichment_func.csv")