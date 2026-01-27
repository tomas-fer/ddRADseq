#!/usr/bin/env Rscript
#Plots piecharts on the map
#Runs within LEAmakePieMaps.sh
#Tomas Fer, 2025-2026
#tomas.fer@natur.cuni.cz

library(maps)
library(plotrix)
library(RColorBrewer)
args <- commandArgs()
clusters <- as.numeric(args[6])
tbl <- read.table(paste0("K",clusters,".txt"), header=T)
my_data <- as.data.frame(tbl)
pie_cols <- names(tbl[,grep("c", colnames(tbl))])
K <- length(pie_cols)
pie_matrix <- my_data[, pie_cols]
totals <- rowSums(pie_matrix)
max_total <- max(totals)
scaling_factor <- 0.15
my_radii <- (totals / max_total) * scaling_factor
pie_colors <- RColorBrewer::brewer.pal(K, "Set2")
pdf(paste0("K",clusters,".pdf"))
#plot(my_data$x, my_data$y,type = "n",xlim = c(min(my_data$x) - 1, max(my_data$x) + 1),ylim = c(48,53),xlab="",ylab="",main=paste0("K ",K))
map(database= "world", ylim = c(min(my_data$y) - 1, max(my_data$y) + 1), xlim = c(min(my_data$x) - 1, max(my_data$x) + 1), col="grey80", fill=F, lwd=2)

for (i in 1:nrow(my_data)) {
  # Get the data for the current pie
  current_pie_data <- as.numeric(pie_matrix[i, ])
  # Check if data is all zero to avoid errors
  if(sum(current_pie_data) > 0) {
    floating.pie(
      xpos = my_data$x[i],   # X position
      ypos = my_data$y[i],   # Y position
      x = current_pie_data,  # The data values for the slices
      radius = my_radii[i],  # Use the scaled radius
      col = pie_colors       # Apply the colors
    )
  }
}
legend("topright",legend = pie_cols,fill = pie_colors,title = "Cluster",cex = 0.8)
title(paste0("K",clusters))
dev.off()
