library(fastverse)
library(rayshader)
library(elevatr)
library(terra)
library(rnaturalearth)

# d = fread("data/adj_2026-07-08T22-44-23-438Z.csv")

# https://www.paulamoraga.com/book-spatial/r-packages-to-download-open-spatial-data.html#elevation

# https://libguides.mit.edu/gis/elevation

# Create a bounding box with these corners:
d = data.frame(long = c(-70.928, -71.181),
               lat =   c(42.316,  42.481))
