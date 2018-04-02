library(raster)

# Example Data
set.seed(1)
r <- raster(matrix(sample(1:10, 100, replace=T), 10, 10))

# Calculate a weights matrix, and reset elements to 0s and 1s 
# rather than true weights
fw <- focalWeight(r, 0.2, 'circle')
fw <- ifelse(fw == 0, NA, 1)

# Neighbourhood richness
richness <- function(x, ...) {
  length(unique(na.omit(x)))
}
richOut <- focal(r, fw, fun=richness, pad=T)

# Neighbourhood Shannon Diversity Index
shannon <- function(x, ...) {
  cnts <- table(x)
  cnts <- cnts / sum(cnts)
  -sum(cnts * log(cnts))
}
shanOut <- focal(r, fw, fun=shannon, pad=T)