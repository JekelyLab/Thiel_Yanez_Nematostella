#code to generate Figure 6 and Figure 6-figure supplement 1 of the Thiel, Yanez-Guerra et al. paper on Nematostella peptide-GPCRs
#Gaspar Jekely 2023

library(tidyverse)
library(tidygraph)
library(Seurat)
library(SeuratObject)
library(sp)
library(igraph)
library(visNetwork)
library(networkD3)
library(webshot2)
library(RColorBrewer)
library(rgexf)
library(leiden)

# From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
  "#CC79A7", "#000000"
)

Tol_muted <- c(
  '#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933',
  '#CC6677', '#882255', '#AA4499', '#DDDDDD'
)
blues <- brewer.pal(9, 'Blues')
bluepurple <- brewer.pal(9, 'BuPu')
oranges <- brewer.pal(9, 'YlOrRd')

GPCR <- readRDS("data/GPC.data_new_names.RDS")

#some methods
GPCR@meta.data
GetAssayData(GPCR)
#FindAllMarkers(GPCR)
#rownames(GPCR)
#GPCR@assays
# Print unique identities in the 'ident' slot
#unique(Idents(GPCR))
# Print unique values in the relevant metadata column
#print(unique(GPCR@meta.data$lifehistory))

# Split Seurat object based on the 'lifehistory' metadata column
seurat_list <- SplitObject(GPCR, split.by = "lifehistory")

# Access the Seurat object for 'DevSubset'
GPCR_DevSubset <- seurat_list$DevSubset

# Access the Seurat object for 'AdultSubset'
GPCR_AdultSubset <- seurat_list$AdultSubset

genes <- rownames(GPCR_DevSubset)

DotPlot(
  object = GPCR_DevSubset,
  features = genes,
  col.min = -2,
  col.max = 2
)


#GetAssayData(GPCR)
#GPCR@assays$RNA@counts
#GPCR@assays$RNA@scale.data
#GPCR@assays$RNA@data

# Calculate average expression of gene1 across all cells
avg_expression_genes <- AverageExpression(object = GPCR, features = genes)

# Convert to data frame and reshape to long format tibble
# Convert to tibble and reshape to long format
avg_expression_df <- as.data.frame(avg_expression_genes)


avg_expression_tibble <- avg_expression_df %>%
  rownames_to_column(var = "genes") %>%
  pivot_longer(cols = -genes, names_to = "celltype", values_to = "avg_expression") %>%
  mutate(celltype = str_replace_all(celltype, "RNA\\.", ""))

celltypes <- unique(avg_expression_tibble %>%
  select(celltype) %>%
  pull())

pNPs <- c(
  "GLWL-a : NV2.23205",
  "HIRa : NV2.8166",
  "HIRa : NV2.8166",
  "LRWa-1 : NV2.10311",
  "LRWa-1 : NV2.10311",
  "LRWa-1 : NV2.10311",
  "LRWa-1 : NV2.10311",
  "PRGa : NV2.16299", "PRGa : NV2.16299",
  "PRGa : NV2.16299", "PRGa : NV2.16299",
  "PRGa : NV2.16299", "PRGa : NV2.16299",
  "PRGa : NV2.16299", "PRGa : NV2.16299",
  "PRGa : NV2.16299", "PRGa : NV2.16299",
  "PRGa : NV2.16299", "PRGa : NV2.16299",
  "PRGa : NV2.16299",
  "QGRFa : NV2.8437", "QGRFa : NV2.8437",
  "QITRFa : NV2.23834",
  "QWa : NV2.4017",
  "FLRNa : NV2.1448",
  "FLRNa : NV2.1448",
  "FLRNa : NV2.1448",
  "FLRNa : NV2.1448",
  "FLRNa : NV2.1448",
  "VRHa : NV2.15165"
)

GPCRs <- c(
  "GLWL-R.18a : NV2.13124",
  "HIRa-R.21 : NV2.23510",
  "HIRa-R.29 : NV2.24982",
  "LRWa-R.193 : NV2.10074",
  "LRWa1-R.19 : NV2.20941",
  "LRWa2-R.213 : NV2.23402",
  "LRWa3-R.204 : NV2.10076",
  "PRGa-R.198 : NV2.22132", "PRGa-R.199 : NV2.21302",
  "PRGa-R.200 : NV2.11164", "PRGa-R.202 : NV2.21659",
  "PRGa-R.210 : NV2.3314", "PRGa-R.211 : NV2.1942", 
  "PRGa-R.219 : NV2.511","PRGa-R.220 : NV2.512",  
  "PRGa-R.221 : NV2.1941", "PRGa-R.222 : NV2.1945", 
  "PRGa-R.223 : NV2.2307", "PRGa-R.28 : NV2.1943",  
  "PRGa-R.32 : NV2.22136",
  "QGRFa-R.234 : NV2.1352", "QGRFa-R.70 : NV2.9922",
  "QITRFa-R.196 : NV2.9618",
  "QWa-R.69 : NV2.9996",
  "FLRNa-R.197-i2 : NV2.9939",
  "FLRNa-R.230 : NV2.16286",
  "FLRNa-R.26 : NV2.22838",
  "FLRNa-R.197-i1 : NV2.14466",
  "FLRNa-R.187 : NV2.23421",
  "VRHa-R.186 : NV2.3564"
)

#EC50 values for the peptide-GPCR pairs
EC50 <- c(
  9e-07, 
  2.1e-9, 
  1.9e-7, 
  6.6e-11, 
  9.1e-8, 
  5.4e-8, 
  7.4e-8, 
  6.9e-10, 5.2e-10, 
  3.1e-9, 4.6e-9, 
  1.6e-10, 2.4e-10, 
  1.5e-8, 5.1e-9, 
  4.5e-8, 1.4e-9, 
  2e-9, 1e-8, 
  3.4e-8, 
  4.7e-9, 2.1e-6, 
  7.7e-11, 
  8e-7, 
  1.1e-8, 
  4.7e-9, 
  7.4e-8, 
  1.1e-8, 
  7e-8, 
  2.7e-10
)

length(pNPs)
length(GPCRs)
length(EC50)


#pNP GPCR pairs
pNP_GPCR <- tibble(peptide = pNPs,
                   receptor = GPCRs,
                   EC50 =EC50)
hex_colors <- c(Okabe_Ito, oranges, Tol_muted, bluepurple[3:9])
length(hex_colors)

#define function to generate graph
create_pNP_GPCR_graph <- function(pNP, GPCR, expressions, colors, cell_types, N_pairs) {
  # create graph with cell types as nodes
  node_IDs <- data.frame(name = cell_types)
  graph.funct <- tbl_graph(nodes = node_IDs)
  
  for (i in 1:N_pairs) {
    # pNP and GPCR name
    pNP_name <- pNP[i]
    GPCR_name <- GPCR[i]
    absEC50_log <- abs(log10(EC50[i]))

    # define sources, targets and expression for pNP and GPCR
    sources <- expressions %>%
      filter(genes == pNP_name) %>%
      select(celltype) %>%
      pull()

    pNP_level <- expressions %>%
      filter(genes == pNP_name) %>%
      select(avg_expression) %>%
      pull()

    targets <- expressions %>%
      filter(genes == GPCR_name) %>%
      select(celltype) %>%
      pull()

    GPCR_level <- expressions %>%
      filter(genes == GPCR_name) %>%
      select(avg_expression) %>%
      pull()
    
    # calculate strength for each edge, defined as the geometric mean of pNP and GPCR expression
    x <- rep(pNP_level, each = length(targets))
    y <- rep(GPCR_level, length(sources))
    geometric_mean <- sqrt(x * y)*absEC50_log

    edge_color <- colors[i]

    # edges are connecting pNP-expressing source nodes with GPCR-expressing target nodes
    graph_edges <- data.frame(
      from = rep(sources, each = length(targets)),
      to = rep(targets, length(sources)),
      value = geometric_mean,
      color = edge_color,
      peptide = pNP_name,
      receptor = GPCR_name
    )
    
    # add edges to graph
    graph.funct <- graph.funct %>%
      bind_edges(graph_edges)
    
   
  }
  return(graph.funct)
}



# networks for ID.separate for dev subset ------------------------------------------

GPCR.separate <- SetIdent(GPCR_DevSubset, value = "ID.separate")

unique(Idents(GPCR.separate))

# set the new IDs in the Seurat object
#Idents(GPCR.separate) <- cleaned_ids

# Calculate average expression of gene1 across all cells
avg_expression_genes.separate <- AverageExpression(object = GPCR.separate, features = genes)

# Convert to tibble and reshape to long format
avg_expression_df.separate <- as.data.frame(avg_expression_genes.separate)

# Replace column names using regex
{
names(avg_expression_df.separate) <- sub("RNA\\.", "", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("neurogland\\.all\\.", "", names(avg_expression_df.separate))
#names(avg_expression_df.separate) <- sub("\\.\\.", ".", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("endomesoderm", "edMes", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("ectomesoderm", "ecMes", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("ectoderm", "ect", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("gastrodermis", "gderm", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("embryonic", "emb", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("spermatagonia", "sperm", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("cnidocyte", "cnido", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("planula", "pla", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("muscle", "mus", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("epithelia", "epith", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("pharyngeal", "pha", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("gland", "gld", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("early\\.state", "early", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("retractor", "retr", names(avg_expression_df.separate))
names(avg_expression_df.separate) <- sub("gland\\.mucous", "gld\\.muc", names(avg_expression_df.separate))

names(avg_expression_df.separate)
}

#rename one GPCR rowname
rownames(avg_expression_df.separate) <- sub("RYAR-like-1 : NV2.3564", "VRHa-R.186 : NV2.3564", rownames(avg_expression_df.separate))
rownames(avg_expression_df.separate) <- sub("QRFPR-like-3 : NV2.23421", "FLRNa-R.187 : NV2.23421", rownames(avg_expression_df.separate))
rownames(avg_expression_df.separate) <- sub("NPY2R-like-1 : NV2.3638", "FRPa-RPa-R.248 : NV2.3638", rownames(avg_expression_df.separate))

avg_expression_tibble.separate <- avg_expression_df.separate %>%
  rownames_to_column(var = "genes") %>%
  pivot_longer(cols = -genes, names_to = "celltype", values_to = "avg_expression")

avg_expression_tibble.separate %>%
  filter(genes == unlist(pNPs) | genes == unlist(GPCRs)) %>%
  ggplot(aes(x=celltype, y=genes, color=sqrt(avg_expression), size=sqrt(avg_expression))) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),
        axis.text = element_text(size=6)) +
  scale_color_gradientn(colors = c("white", Okabe_Ito[c(6,5,5,5,8,8,8)]))

celltypes.separate <- unique(avg_expression_tibble.separate %>%
                      select(celltype) %>%
                      pull())

avg_expression_tibble.separate %>%
  select(avg_expression) %>%
  ggplot(aes(x=avg_expression)) +
  geom_histogram() +
  scale_x_log10()
  
# do graph with function ------------------------------------

graph.separate <- create_pNP_GPCR_graph(
  pNPs, GPCRs, avg_expression_tibble.separate, 
  hex_colors[1:length(pNPs)], celltypes.separate, length(pNPs)
)

# remove zero weight edges
graph.separate <- graph.separate %>%
  activate(edges) %>%
  filter(value > 0) 

# Save `tbl_graph` as a serialized binary file
saveRDS(graph.separate, "source_data/Figure6_source_data_1.rds")

# filter by edge weight
graph.separate <- graph.separate %>%
  activate(edges) %>%
  filter(value > 30) 

#number of connections after filtering
graph.separate %>%
  activate(edges) %>%
  select(value) %>%
  pull() %>%
  length()

# Get the subgraph of connected nodes
connected_nodes <- components(graph.separate)$membership == which.max(components(graph.separate)$csize)
graph.separate.connected <- induced_subgraph(graph.separate, which(connected_nodes))

#add weighted degree to nodes
graph.separate.connected <- graph.separate.connected %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(value = centrality_degree(
    weights = NULL,
    mode = "all",
    loops = TRUE,
    normalized = FALSE
  )
)

# Calculate the Leiden clustering and modularity
lc <- leiden(graph.separate.connected)
max(lc)

# visNetwork plot with peptides colored --------------------------------

## convert to VisNetwork
Conn_graph.visn <- toVisNetworkData(graph.separate.connected)

#assing module id to nodes
Conn_graph.visn$nodes$group <- lc

#change node color to grey
Conn_graph.visn$nodes$color <- "grey"

Conn_graph.visn$nodes$fontSize <- Conn_graph.visn$nodes$value

visNet <- visNetwork(Conn_graph.visn$nodes, Conn_graph.visn$edges) %>%
  visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 23, type = "full") %>%
  visPhysics(
    solver = "hierarchicalRepulsion", 
    hierarchicalRepulsion = list(
      nodeDistance = 200, 
      centralGravity = 0.1, 
      springConstant = 0.00005
    )
  ) %>%
  visEdges(
    smooth = list(style = "curve-style: haystack;"),
    scaling = list(min = 5, max = 35),
    color = list(inherit = TRUE, opacity = 0.7),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 1.5, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(border = "black"),
    opacity = 1,
    scaling = list(min = 65, max = 120),
    shape = "dot",
    font = list(color = "black", strokeWidth = 2, background = "white", size = 75),
  ) %>%
  visOptions(highlightNearest = TRUE, width = 2000, height = 1600)


saveNetwork(visNet, "pictures/peptidergic_networks_full_peptides_dev.html")
user_agent <- "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:109.0) Gecko/20100101 Firefox/112.0"

webshot2::webshot(
  url = "pictures/peptidergic_networks_full_peptides_dev.html",
  file = "pictures/peptidergic_networks_full_peptides_dev.png", useragent = user_agent,
  vwidth = 2000, vheight = 1600, # define the size of the browser window
  cliprect = c(100, 120, 1800, 1480), zoom = 2, delay = 100
)

library(rgexf)
# Serialise the graph to Gephi format
gexf_data <- rgexf::igraph.to.gexf(as.igraph(as_tbl_graph(Conn_graph.visn)))
# Write the network into a gexf (Gephi) file
write.gexf(gexf_data, output = "source_data/Figure6_source_data2.gexf")


# visNetwork plot with modules highlighted --------------------------------

## convert to VisNetwork
Conn_graph.visn <- toVisNetworkData(graph.separate.connected)

#assing module id to nodes
Conn_graph.visn$nodes$group <- lc

Conn_graph.visn$edges$color <- c()

# Define colors for each group
color_lookup <- c("1" = Okabe_Ito[1], "2" = Okabe_Ito[2], "3" = Okabe_Ito[3], "4" = Okabe_Ito[5], "5" = Okabe_Ito[7])

# Assign colors based on group
Conn_graph.visn$nodes$color <- sapply(Conn_graph.visn$nodes$group, function(group) {
  return(color_lookup[group])
})

unique(Conn_graph.visn$edges$peptide)
Conn_graph.visn$nodes$fontSize <- Conn_graph.visn$nodes$value

visNet <- visNetwork(Conn_graph.visn$nodes, Conn_graph.visn$edges) %>%
  visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 23, type = "full") %>%
  visPhysics(
    solver = "hierarchicalRepulsion", 
    hierarchicalRepulsion = list(
      nodeDistance = 200, 
      centralGravity = 0.1, 
      springConstant = 0.00005
      )
    ) %>%
  visEdges(
    smooth = list(style = "curve-style: haystack;"),
    scaling = list(min = 5, max = 35),
    color = list(inherit = TRUE, opacity = 0.7),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 1.5, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(border = "black"),
    opacity = 1,
    scaling = list(min = 65, max = 120),
    shape = "dot",
    font = list(color = "black", strokeWidth = 2, background = "white", size = 75),
  ) %>%
  visOptions(highlightNearest = TRUE, width = 2000, height = 1600)


saveNetwork(visNet, "pictures/peptidergic_networks_full_modules_dev.html")
user_agent <- "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:109.0) Gecko/20100101 Firefox/112.0"

webshot2::webshot(
  url = "pictures/peptidergic_networks_full_modules_dev.html",
  file = "pictures/peptidergic_networks_full_modules_dev.png", useragent = user_agent,
  vwidth = 2000, vheight = 1600, # define the size of the browser window
  cliprect = c(100, 120, 1800, 1480), zoom = 2, delay = 100
)

# export graphs for each peptide ------------------------------------------

## convert to VisNetwork
Conn_graph.visn <- toVisNetworkData(graph.separate.connected)

getIDColor <- function(id) {
  if(grepl("N_", id)) {
    return("#56B4E9") 
  } else if(grepl("N_gland", id)) {
    return("#E69F00")
  } else if(grepl("mucous", id)) {
    return("#CC79A7") 
  } else if(grepl("PGCs", id)) {
    return("#F0E442")
  } else if(grepl("endomesoderm", id)) {
    return("#0072B2")
  } else if(grepl("muscle", id)) {
    return("#111111")
  } else {
    return("#CCCCCC") # (default color)
  }
}

colors <- sapply(Conn_graph.visn$nodes$id, getIDColor) # apply getIDColor to each ID
Conn_graph.visn$nodes$color <- colors


unique_pNPs <- unique(Conn_graph.visn$edges$peptide)
unique_pNPs

for (i in 1:length(unique_pNPs)) {

#filter connections
connections_to_show <- !sapply(
  Conn_graph.visn$edges$peptide, 
  grepl, 
  pattern = unique_pNPs[i]
)

#add FALSE or TRUE to edges$hidden variable (will toggle edges on/off)
Conn_graph.visn$edges$hidden <- connections_to_show
Conn_graph.visn$edges$hidden
visNet <- visNetwork(Conn_graph.visn$nodes, Conn_graph.visn$edges) %>%
  visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 23, type = "full") %>%
  visPhysics(
    solver = "hierarchicalRepulsion", 
    hierarchicalRepulsion = list(
      nodeDistance = 200, 
      centralGravity = 0.1, 
      springConstant = 0.00005
    )
  ) %>%
  visEdges(
    smooth = list(style = "curve-style: haystack;"),
    scaling = list(min = 5, max = 35),
    color = list(inherit = TRUE, opacity = 0.7),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 1.5, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(border = "black"),
    opacity = 1,
    scaling = list(min = 65, max = 120),
    shape = "dot",
    font = list(color = "black", strokeWidth = 2, background = "white", size = 75),
  ) %>%
  visOptions(highlightNearest = TRUE, width = 2000, height = 1600)

filename <- paste("pictures/peptidergic_networks_full_peptide_dev_", unique_pNPs[i], ".html", sep = "")
filename2 <- paste("pictures/peptidergic_networks_full_peptide_dev_", unique_pNPs[i], ".png", sep = "")
saveNetwork(visNet, filename)

webshot2::webshot(
  url = filename,
  file = filename2, useragent = user_agent,
  vwidth = 2000, vheight = 1600, # define the size of the browser window
  cliprect = c(160, 120, 1750, 1480), zoom = 2, delay = 100
)

}

# networks for ID.separate for adult subset ------------------------------------------

GPCR.separate <- SetIdent(GPCR_AdultSubset, value = "ID.separate")

DotPlot(
  object = GPCR.separate,
  features = genes,
  col.min = -2,
  col.max = 2
) +
  theme(axis.text = element_text(size = 4))

unique(Idents(GPCR.separate))

# Calculate average expression of gene1 across all cells
avg_expression_genes.separate <- AverageExpression(object = GPCR.separate, features = genes)

# Convert to tibble and reshape to long format
avg_expression_df.separate <- as.data.frame(avg_expression_genes.separate)

# Replace column names using regex
{
  names(avg_expression_df.separate) <- sub("RNA\\.", "", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("neurogland\\.all\\.", "", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("endomesoderm", "edMes", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("ectomesoderm", "ecMes", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("ectoderm", "ect", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("gastrodermis", "gderm", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("embryonic", "emb", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("spermatagonia", "sperm", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("cnidocyte", "cnido", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("planula", "pla", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("muscle", "mus", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("epithelia", "epith", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("pharyngeal", "pha", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("gland", "gld", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("early\\.state", "early", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("retractor", "retr", names(avg_expression_df.separate))
  names(avg_expression_df.separate) <- sub("gland\\.mucous", "gld\\.muc", names(avg_expression_df.separate))
  
  names(avg_expression_df.separate)
}

rownames(avg_expression_df.separate) <- sub("RYAR-like-1 : NV2.3564", "VRHa-R.186 : NV2.3564", rownames(avg_expression_df.separate))
rownames(avg_expression_df.separate) <- sub("QRFPR-like-3 : NV2.23421", "FLRNa-R.187 : NV2.23421", rownames(avg_expression_df.separate))
rownames(avg_expression_df.separate) <- sub("NPY2R-like-1 : NV2.3638", "FRPa-RPa-R.248 : NV2.3638", rownames(avg_expression_df.separate))

avg_expression_tibble.separate <- avg_expression_df.separate %>%
  rownames_to_column(var = "genes") %>%
  pivot_longer(cols = -genes, names_to = "celltype", values_to = "avg_expression")

avg_expression_tibble.separate %>%
  select(avg_expression) %>%
  max()

avg_expression_tibble.separate %>%
  filter(genes == unlist(pNPs) | genes == unlist(GPCRs)) %>%
  ggplot(aes(x=celltype, y=genes, color=sqrt(avg_expression), size=sqrt(avg_expression))) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),
        axis.text = element_text(size=6)) +
  scale_color_gradientn(colors = c("white", Okabe_Ito[c(6,5,5,5,8,8,8)]))

celltypes.separate <- unique(avg_expression_tibble.separate %>%
                               select(celltype) %>%
                               pull())


# create and grow graph with function ------------------------------------

graph.separate <- create_pNP_GPCR_graph(
  pNPs, GPCRs, avg_expression_tibble.separate, 
  hex_colors[1:length(pNPs)], celltypes.separate, length(pNPs)
)

# remove zero weight edges
graph.separate <- graph.separate %>%
  activate(edges) %>%
  filter(value > 0) 

# Save `tbl_graph` as a serialized binary file
saveRDS(graph.separate, "source_data/Figure6_source_data_3.rds")

#plot edge weight distribution
graph.separate %>%
  activate(edges) %>%
  as.data.frame() %>%
  ggplot(aes(x=(value))) +
  geom_histogram() +
  scale_x_log10() +
  theme_minimal()


# filter by edge weight
graph.separate <- graph.separate %>%
  activate(edges) %>%
  filter(value > 15) 

#number of connections after filtering
graph.separate %>%
  activate(edges) %>%
  select(value) %>%
  pull() %>%
  length()

# Get the subgraph of connected nodes
connected_nodes <- components(graph.separate)$membership == which.max(components(graph.separate)$csize)
graph.separate.connected <- induced_subgraph(graph.separate, which(connected_nodes))

#add weighted degree to nodes
graph.separate.connected <- graph.separate.connected %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(value = centrality_degree(
    weights = NULL,
    mode = "all",
    loops = TRUE,
    normalized = FALSE
  )
  )

#plot edge weight distribution
graph.separate.connected %>%
  activate(edges) %>%
  as.data.frame() %>%
  ggplot(aes(x=(value))) +
  geom_histogram() +
  scale_x_log10() +
  theme_minimal()

# Calculate the Leiden clustering and modularity
lc <- leiden(graph.separate.connected)
max(lc)

# visNetwork plot with peptides colored --------------------------------

## convert to VisNetwork
Conn_graph.visn <- toVisNetworkData(graph.separate.connected)

#assing module id to nodes
Conn_graph.visn$nodes$group <- lc

#change node color to grey
Conn_graph.visn$nodes$color <- "grey"

Conn_graph.visn$nodes$fontSize <- Conn_graph.visn$nodes$value

visNet <- visNetwork(Conn_graph.visn$nodes, Conn_graph.visn$edges) %>%
  visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 23, type = "full") %>%
  visPhysics(
    solver = "hierarchicalRepulsion", 
    hierarchicalRepulsion = list(
      nodeDistance = 200, 
      centralGravity = 0.1, 
      springConstant = 0.00005
    )
  ) %>%
  visEdges(
    smooth = list(style = "curve-style: haystack;"),
    scaling = list(min = 5, max = 35),
    color = list(inherit = TRUE, opacity = 0.7),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 1.5, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(border = "black"),
    opacity = 1,
    scaling = list(min = 65, max = 120),
    shape = "dot",
    font = list(color = "black", strokeWidth = 2, background = "white", size = 75),
  ) %>%
  visOptions(highlightNearest = TRUE, width = 2000, height = 1600)


saveNetwork(visNet, "pictures/peptidergic_networks_full_peptides.html")
user_agent <- "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:109.0) Gecko/20100101 Firefox/112.0"

webshot2::webshot(
  url = "pictures/peptidergic_networks_full_peptides.html",
  file = "pictures/peptidergic_networks_full_peptides.png", useragent = user_agent,
  vwidth = 2000, vheight = 1600, # define the size of the browser window
  cliprect = c(100, 120, 1800, 1480), zoom = 2, delay = 100
)

library(rgexf)
# Serialise the graph to Gephi format
gexf_data <- rgexf::igraph.to.gexf(as.igraph(as_tbl_graph(Conn_graph.visn)))
# Write the network into a gexf (Gephi) file
write.gexf(gexf_data, output = "source_data/Figure6_source_data4.gexf")


# visNetwork plot with modules highlighted --------------------------------

## convert to VisNetwork
Conn_graph.visn <- toVisNetworkData(graph.separate.connected)

#assing module id to nodes
Conn_graph.visn$nodes$group <- lc

Conn_graph.visn$edges$color <- c()

# Define colors for each group
color_lookup <- c("1" = Okabe_Ito[1], "2" = Okabe_Ito[2], "3" = Okabe_Ito[3], "4" = Okabe_Ito[5], "5" = Okabe_Ito[7])

# Assign colors based on group
Conn_graph.visn$nodes$color <- sapply(Conn_graph.visn$nodes$group, function(group) {
  return(color_lookup[group])
})

unique(Conn_graph.visn$edges$peptide)
Conn_graph.visn$nodes$fontSize <- Conn_graph.visn$nodes$value

visNet <- visNetwork(Conn_graph.visn$nodes, Conn_graph.visn$edges) %>%
  visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 23, type = "full") %>%
  visPhysics(
    solver = "hierarchicalRepulsion", 
    hierarchicalRepulsion = list(
      nodeDistance = 200, 
      centralGravity = 0.1, 
      springConstant = 0.00005
    )
  ) %>%
  visEdges(
    smooth = list(style = "curve-style: haystack;"),
    scaling = list(min = 5, max = 35),
    color = list(inherit = TRUE, opacity = 0.7),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 1.5, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(border = "black"),
    opacity = 1,
    scaling = list(min = 65, max = 120),
    shape = "dot",
    font = list(color = "black", strokeWidth = 2, background = "white", size = 75),
  ) %>%
  visOptions(highlightNearest = TRUE, width = 2000, height = 1600)


saveNetwork(visNet, "pictures/peptidergic_networks_full_modules.html")
user_agent <- "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:109.0) Gecko/20100101 Firefox/112.0"

webshot2::webshot(
  url = "pictures/peptidergic_networks_full_modules.html",
  file = "pictures/peptidergic_networks_full_modules.png", useragent = user_agent,
  vwidth = 2000, vheight = 1600, # define the size of the browser window
  cliprect = c(100, 120, 1800, 1480), zoom = 2, delay = 100
)


# export graph for each peptide ------------------------------------------

## convert to VisNetwork
Conn_graph.visn <- toVisNetworkData(graph.separate.connected)

getIDColor <- function(id) {
  if(grepl("N_", id)) {
    return("#56B4E9") 
  } else if(grepl("N_gland", id)) {
    return("#E69F00")
  } else if(grepl("mucous", id)) {
    return("#CC79A7") 
  } else if(grepl("PGCs", id)) {
    return("#F0E442")
  } else if(grepl("endomesoderm", id)) {
    return("#0072B2")
  } else if(grepl("muscle", id)) {
    return("#111111")
  } else {
    return("#CCCCCC") # (default color)
  }
}

colors <- sapply(Conn_graph.visn$nodes$id, getIDColor) # apply getIDColor to each ID
Conn_graph.visn$nodes$color <- colors


unique_pNPs <- unique(Conn_graph.visn$edges$peptide)
unique_pNPs

for (i in 1:length(unique_pNPs)) {
  
  #filter connections
  connections_to_show <- !sapply(
    Conn_graph.visn$edges$peptide, 
    grepl, 
    pattern = unique_pNPs[i]
  )
  
  #add FALSE or TRUE to edges$hidden variable (will toggle edges on/off)
  Conn_graph.visn$edges$hidden <- connections_to_show
  Conn_graph.visn$edges$hidden
  visNet <- visNetwork(Conn_graph.visn$nodes, Conn_graph.visn$edges) %>%
    visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 23, type = "full") %>%
    visPhysics(
      solver = "hierarchicalRepulsion", 
      hierarchicalRepulsion = list(
        nodeDistance = 200, 
        centralGravity = 0.1, 
        springConstant = 0.00005
      )
    ) %>%
    visEdges(
      smooth = list(style = "curve-style: haystack;"),
      scaling = list(min = 5, max = 35),
      color = list(inherit = TRUE, opacity = 0.7),
      arrows = list(to = list(
        enabled = TRUE,
        scaleFactor = 1.5, type = "arrow"
      ))
    ) %>%
    visNodes(
      borderWidth = 0.3,
      color = list(border = "black"),
      opacity = 1,
      scaling = list(min = 65, max = 120),
      shape = "dot",
      font = list(color = "black", strokeWidth = 2, background = "white", size = 75),
    ) %>%
    visOptions(highlightNearest = TRUE, width = 2000, height = 1600)
  
  filename <- paste("pictures/peptidergic_networks_full_peptide_", unique_pNPs[i], ".html", sep = "")
  filename2 <- paste("pictures/peptidergic_networks_full_peptide_", unique_pNPs[i], ".png", sep = "")
  saveNetwork(visNet, filename)
  
  webshot2::webshot(
    url = filename,
    file = filename2, useragent = user_agent,
    vwidth = 2000, vheight = 1600, # define the size of the browser window
    cliprect = c(160, 120, 1750, 1480), zoom = 2, delay = 100
  )
  
}

# plot color scale of edges -----------------------------------------------

graph.separate.connected %>%
  activate(edges) %>%
  as_tibble() %>%
  select(color,peptide) %>%
  unique() %>%
  ggplot(aes(y = peptide, fill = color)) +
  geom_bar(position = "stack") +
  guides(fill = "none") +
  theme_minimal() +
  scale_x_discrete() +
  theme(axis.title = element_blank())


# assemble figure ---------------------------------------------------------

library(cowplot)
library(patchwork)
library(png)
library(ggplot2)

# Figure 6  ------------------

layout_Fig6 <- "
AB"

panel_pep_dev <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_peptides_dev.png"
  )
) +
  draw_label("developmental subset", x = 0.1, y = 0.98, size = 14, hjust = 0)

panel_pep_ad <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_peptides.png"
  )
) +
  draw_label("adult subset", x = 0.12, y = 0.98, size = 14, hjust = 0)


Fig6 <- panel_pep_dev + panel_pep_ad +
  plot_layout(design = layout_Fig6, guides = "collect", widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "plain"))

ggsave("figures/Figure6.png",
       limitsize = FALSE,
       units = c("px"), Fig6, width = 3500, height = 1580, bg = "white"
)

ggsave("figures/Figure6.pdf",
       limitsize = FALSE,
       units = c("px"), Fig6, width = 3500, height = 1580
)


# Fig6 fig suppl 1 ----------------------
panel_all_dev <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_modules_dev.png"
  )
) +
  draw_label("Leiden modules\ndev. subset", x = 0.12, y = 0.98, size = 14, hjust = 0)


panel_LRWa_dev <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_peptide_dev_LRWa-1 : NV2.10311.png"
  )
) +
  draw_label("LRWa\ndev. subset", x = 0.1, y = 0.98, size = 14, hjust = 0)

panel_PRGa_dev <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_peptide_dev_PRGa : NV2.16299.png"
  )
) +
  draw_label("PRGa\ndev. subset", x = 0.1, y = 0.98, size = 14, hjust = 0)

panel_FLRNa_dev <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_peptide_dev_FLRNa : NV2.1448.png"
  )
) +
  draw_label("FLRNa\ndev. subset", x = 0.1, y = 0.98, size = 14, hjust = 0)



#adult panels

panel_all <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_modules.png"
    )
  ) +
  draw_label("Leiden modules\nadult subset", x = 0.12, y = 0.98, size = 14, hjust = 0)

panel_HIRa <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_peptide_HIRa : NV2.8166.png"
    )
  ) +
  draw_label("HIRa\nadult subset", x = 0.1, y = 0.98, size = 14, hjust = 0)

panel_LRWa <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_peptide_LRWa-1 : NV2.10311.png"
    )
  ) +
  draw_label("LRWa\nadult subset", x = 0.1, y = 0.98, size = 14, hjust = 0)

panel_PRGa <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_peptide_PRGa : NV2.16299.png"
    )
  ) +
  draw_label("PRGa\nadult subset", x = 0.1, y = 0.98, size = 14, hjust = 0)

panel_QGRFa <- ggdraw() + draw_image(
  readPNG(
    "pictures//peptidergic_networks_full_peptide_QGRFa : NV2.8437.png"
    )
  ) +
  draw_label("QGRFa\nadult subset", x = 0.1, y = 0.98, size = 14, hjust = 0)

panel_QWa <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_peptide_QWa : NV2.4017.png"
    )
  ) +
  draw_label("QWa\nadult subset", x = 0.1, y = 0.98, size = 14, hjust = 0)

panel_FLRNa <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_peptide_FLRNa : NV2.1448.png"
    )
  ) +
  draw_label("FLRNa\nadult subset", x = 0.1, y = 0.98, size = 14, hjust = 0)

panel_VRHa <- ggdraw() + draw_image(
  readPNG(
    "pictures/peptidergic_networks_full_peptide_VRHa : NV2.15165.png"
    )
  ) +
  draw_label("VRHa\nadult subset", x = 0.1, y = 0.98, size = 14, hjust = 0)

layout <- "
ABCD
####
EFGH
####
IJKL
"

Fig6_fig_suppl1 <- panel_all_dev +  panel_LRWa_dev + panel_FLRNa_dev + panel_PRGa_dev +
  panel_all +  panel_LRWa + panel_FLRNa + panel_PRGa + 
  panel_QGRFa + panel_HIRa + panel_VRHa + panel_QWa +
  plot_layout(design = layout, guides = "collect", heights = c(1, 0.05, 1, 0.05, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "plain"))

ggsave("figures/Fig6_fig_suppl1.png",
       limitsize = FALSE,
       units = c("px"), Fig6_fig_suppl1, width = 1750*4, height = 1480*3+150, bg = "white"
)

ggsave("figures/Fig6_fig_suppl1.pdf",
       limitsize = FALSE,
       units = c("px"), Fig6_fig_suppl1, width = 1750*4, height = 1480*3+150
)

