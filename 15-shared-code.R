rm(list=ls())

if ("RStudioGD" %in% names(dev.list())) {
  dev.off(dev.list()["RStudioGD"])
}

suppressPackageStartupMessages(library(fs, quietly=T))
suppressPackageStartupMessages(library(readxl, quietly=T))
suppressPackageStartupMessages(library(tidyverse, quietly=T, verbose=F))
suppressPackageStartupMessages(library(sp, quietly=T))

# Paths
path.base <- path_real("~/SyncThing/nrv-dr-pitcairnia/KR")
path.scpt <- path(path.base, "pop-gen-01-scripts")
path.rslt <- path(path.base, "pop-gen-02-results")
path.resamp.snps <- path(path.rslt, '29-resampled-snps')
path.data <- path(path.scpt, "data")
path.data.rdata <- path(path.data, "RData")
path.data.clim <- path(path.data, "climate")

dir_create(c(path.base, path.scpt, path.rslt, path.data, path.data.clim))

path.sample.data <- path(path.data, "nrv-sample-metadata.xlsx")
sample.dat <- read_xlsx(path.sample.data, sheet='geo')
sample.dat <- sample.dat %>% arrange(id)

sample.dat.raster.overlap.fix <- sample.dat
sample.dat.raster.overlap.fix[16,]$lon <- sample.dat.raster.overlap.fix[16,]$lon - 0.005
sample.dat.raster.overlap.fix[16,]$lat <- sample.dat.raster.overlap.fix[16,]$lat - 0.005
sample.dat.raster.overlap.fix[17,]$lon <- sample.dat.raster.overlap.fix[17,]$lon - 0.003
sample.dat.raster.overlap.fix[17,]$lat <- sample.dat.raster.overlap.fix[17,]$lat - 0.000
sample.dat.sp <- SpatialPointsDataFrame(
  select(sample.dat.raster.overlap.fix, lon, lat), sample.dat)
