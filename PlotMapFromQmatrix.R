library(maps)
args <- commandArgs()
K <- as.numeric(args[5])
fullname <- readLines("fullname")
print(paste("Sample full name is: ", fullname))
print(paste("K: ", K))

source("POPSutilitiesCorrectedConstraines.r")
coord = read.table("coordinates.txt")
asc.raster = "Europe.asc"
grid = createGridFromAsciiRaster(asc.raster)
constraints = getConstraintsFromAsciiRaster(asc.raster,cell_value_min=-10)

myColorGradients = list( 
	c("gray95",brewer.pal(9,"Greens")),
	c("gray95",brewer.pal(9,"Reds")),
	c("gray95",brewer.pal(9,"Blues")),
	c("gray95",brewer.pal(9,"YlOrBr")),
	c("gray95",brewer.pal(9,"RdPu")),
	c("gray95",brewer.pal(9,"Greys")),
	c("gray95",brewer.pal(9,"Purples"))
)

#Loop over different K (PDF width and height need to be adjusted based on xlim and ylim to keep "correct projection", i.e. appearance of the map)
for (x in 2:K) {
  if (x == 2) {
    Qmatrix = read.table(paste0("k", x, ".Q"))
  } else {
    Qmatrix = read.table(paste0("k", x, "new.Q"))
  }
  pdf(file = paste0("K", x, "_regions.pdf"), width = 12, height = 9)
  maps(matrix = Qmatrix, coord, grid, constraints, method = "max", main = paste0(fullname, ", K=", x), xlab = "Longitude", ylab = "Latitude", xlim = c(12, 19), ylim = c(48.5, 51.1), colorGradientsList = myColorGradients, cex.main = 2, cex.axis = 1)
  map(add = T, interior = T, col = "black", lwd=3)
  dev.off()
  #Reorder Q matrix for x+1
  if (x == 2) {
    qA = read.table(paste0("k", x, ".Q"))
  } else {
    qA = read.table(paste0("k", x, "new.Q"))
  }
  if (x == K) { break () }
  qB = read.table(paste0("k", x+1, ".Q"))
  cor_matrix <- cor(qA, qB)
  best_matches <- apply(cor_matrix, 1, which.max)
  all_cols <- 1:ncol(qB)
  unassigned_cols <- setdiff(all_cols, best_matches)
  new_column_order <- c(best_matches, unassigned_cols)
  reordered_qB <- qB[, new_column_order]
  write.table(reordered_qB, file=paste0("k", x+1, "new.Q"))
}

print("Finished...")
