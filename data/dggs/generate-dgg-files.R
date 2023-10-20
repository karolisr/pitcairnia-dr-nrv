library(fs)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(tmaptools)
library(tmap)
library(dggridR)

setwd("/home/karolis/scratch/pop-gen-01-scripts/data/dggs")

# for (i in 1:10) {
#
#   grid.resolution <- i
#   out.file.name <- paste0(
#     'dgg_world_triangle_res_', sprintf('%02d', grid.resolution), '.shp')
#
#   if (file_exists(out.file.name)) {
#     next
#   }
#
#   message('Generating:', out.file.name)
#
# dgearthgrid(
#   dgconstruct(
#     res=grid.resolution,
#     projection='ISEA',
#     aperture=4,
#     topology='TRIANGLE'),
#   savegrid=out.file.name)
#
# }

# -----------------------------------------------------------------------------

ht <- ne_countries(
  scale = "large", country = "Haiti", returnclass = "sf"
)

dr <- ne_countries(
  scale = "large", country = "Dominican Republic", returnclass = "sf"
)

hisp <- rbind(ht, dr)
st_write(hisp, "hispanola.shp", append = FALSE)

grid_resolution <- 09

dggs <- dgconstruct(
  res = grid_resolution,
  precision = 30,
  projection = "ISEA",
  aperture = 4,
  topology = "TRIANGLE"
)

# hisp.grid <- dgshptogrid(dggs, 'hispanola.shp', cellsize=0.01)
# plot(hisp.grid)

out_file_name <- paste0(
  "dgg_hisp_triangle_res_", sprintf("%02d", grid_resolution), ".shp"
)
hisp_grid_file_name_as_written <- dgshptogrid(
  dggs, "hispanola.shp",
  cellsize = 0.01, savegrid = out_file_name
)

outline_1 <- read.csv("../dominican_republic_outline_coords_1", header = FALSE)
outline_2 <- read.csv("../dominican_republic_outline_coords_2", header = FALSE)

# world.grid <- read_sf('dgg_world_triangle_res_03.shp')
hisp_grid <- read_sf("dgg_hisp_triangle_res_11.shp")

# qtm(c('Haiti', 'Dominican Republic'))
bb_ht_dr <- bb_poly(bb(c(dr$geometry, ht$geometry)))
qtm(c(bb_ht_dr, hisp_grid$geometry, ht$geometry, dr$geometry))
