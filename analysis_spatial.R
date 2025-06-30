library(Rcpp)
library(ggplot2)
sourceCpp("~/project/LUNG/analysis/moran.cpp")
analyze_spatial <- function(coords, size, k=15, type=c("global", "local")) {
  # Input validation
  if (nrow(coords) != length(size)) 
    stop("Number of rows in coordinate matrix must match expression vector length")
  
  # Convert input types
  coords <- as.matrix(coords)
  size <- as.numeric(size)
  
  # Call C++ core function
  result <- cpp_moran(coords, size, k=k, global=(type[1]=="global"))
  
  # Process results
  if (type[1] == "global") {
    cat(sprintf("Global Moran's I = %.4f\n", result$global))
    return(result$global)
  } else {
    df <- data.frame(
      x = coords[,1],
      y = coords[,2],
      size = size,
      local_I = result$local
    )
    
    # Calculate significance
    df$p_val <- 2 * pnorm(-abs(df$local_I))
    df$sig <- ifelse(df$p_val < 0.05, 
                     ifelse(df$size > mean(df$size), 
                            ifelse(df$local_I > 0, "High-High", "High-Low"),
                            ifelse(df$local_I < 0, "Low-Low", "Low-High")),
                     "Non-sig")
    return(df)
  }
}


# Generate simulated data for 600,000 cells
set.seed(123)
n_cells <- 60000
coords <- cbind(runif(n_cells, 0, 1000), runif(n_cells, 0, 1000))
size <- 10 * exp(-((coords[,1]-500)^2 + (coords[,2]-500)^2)/(200^2)) + rnorm(n_cells)
# Global analysis
system.time({
  global_I <- analyze_spatial(coords, size, k=15, type="global")
})
# Visualization function
plot_moran <- function(result_df, alpha=0.5) {
  ggplot(result_df, aes(x, y)) +
    geom_point(aes(color=sig, size=abs(local_I)), alpha=alpha) +
    scale_color_manual(values=c("High-High"="#FFFF00", "Low-Low"="#00FFFF",
                                "High-Low"="#FF00FF", "Low-High"="#00FF00",
                                "Non-sig"="#cccccc")) +
    scale_size_continuous(range=c(0.1, 2)) +
    labs(title="Local Spatial Autocorrelation",
         x="X Coordinate", y="Y Coordinate",
         color="Cluster Type", size="|Local Moran's I|") +
    theme_minimal() +
    coord_fixed()
}
library(tidyverse)
cell_size <- read_csv('~/project/LUNG/analysis/sup9/cells.csv.gz')
cell_segment <- read_csv('~/project/LUNG/analysis/sup9/cell_boundaries.csv.gz')
cells_meta <- left_join(cell_segment, cell_size)
cell_cluster <- read_csv('~/project/LUNG/analysis/sup9/gene_expression_kmeans_8_clusters/clusters.csv')
cells_meta <- inner_join(cells_meta, cell_cluster, by=c('cell_id'='Barcode'))
cells_meta$vertex_x <- -cells_meta$vertex_x

# Calculate global Moran's I for all clusters
global_I_list <- pbmcapply::pbmclapply(1:8, function(x){
  cells_meta_C <- subset(cells_meta, Cluster %in% x)
  cell_data <- cells_meta_C
  coords <- cbind(cell_data$vertex_x, cell_data$vertex_y)
  size <- cell_data$cell_area
  global_I <- analyze_spatial(coords, size, k=20, type="global")
  return(global_I)
}, mc.cores = 10)

# Calculate local Moran's I for all clusters
local_I_list <- pbmcapply::pbmclapply(1:8, function(x){
  cells_meta_C <- subset(cells_meta, Cluster %in% x)
  cell_data <- cells_meta_C
  coords <- cbind(cell_data$vertex_x, cell_data$vertex_y)
  size <- cell_data$cell_area
  sample_idx <- sample(1:length(size), 50000)
  local_I <- analyze_spatial(coords, size, k=25, type="local")
  local_I$cluster <- x
  return(local_I)
}, mc.cores = 10)

# Combine and filter results
local_I_dat <- Reduce(rbind, local_I_list)
local_I_dat <- filter(local_I_dat, !sig %in% 'Non-sig')

# Save results visualization
pdf('~/project/LUNG/analysis/figures/fig8_Moran_bar.pdf', width=10, height=10)
ggplot(local_I_dat, aes(x = factor(cluster), fill=sig)) +
  geom_bar()
dev.off()

# Generate Moran visualization plots
pdf('~/project/LUNG/analysis/figures/fig8_Moran.pdf', width=10, height=10)
for (i in 1:8) {
  p <- plot_moran(local_I_list[[i]]) +
    ggtitle(paste("Cluster", i, "Spatial Pattern of Cell Size")) +
    theme(
      panel.background = element_rect(fill="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill="white")
    )
  print(p)
}
dev.off()

# Save all data
save(list=ls(), file='~/project/LUNG/analysis/fig8_Moran_alldata.rda')

# Domain analysis
Selection_1 <- read_csv('~/project/LUNG/analysis/sup9/Selection_1_cells_stats.csv') %>% data.frame()
Selection_2 <- read_csv('~/project/LUNG/analysis/sup9/Selection_2_cells_stats.csv') %>% data.frame()
Selection_1$Domain <- 'Domain1'
Selection_2$Domain <- 'Domain2'
Domain_Cluster3 <- rbind(Selection_1, Selection_2)

# Statistical test (Mann-Whitney U test)
wilcox.test(Area ~ Domain, data=Domain_Cluster3, exact=FALSE)  # p < 0.0001

# Boxplot visualization
library(ggpattern)
library(ggpubr)
pdf('~/project/LUNG/analysis/figures/sup9/Domain_Cluster3_Area_boxplot.pdf', width=3.2, height=5)
ggboxplot(Domain_Cluster3, x="Domain", y="Area",
          color="Domain",
          palette=c("#E60012", "#2EA7E0"),
          add="jitter") +
  stat_compare_means(comparisons=list(c('Domain1','Domain2'))) + 
  stat_compare_means(label.y=50)   
dev.off()

# Generate final Moran plots
pdf('~/project/LUNG/analysis/figures/sup9/fig9_Moran.pdf', width=10, height=10)
for (i in 1:8) {
  p <- plot_moran(local_I_list[[i]]) +
    ggtitle(paste("Cluster", i, "Spatial Pattern of Cell Size")) +
    theme(
      panel.background = element_rect(fill="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill="white")
    )
  print(p)
}
dev.off()