#!/usr/bin/env Rscript

source("15-shared-code.R")

library(readxl)
library(tidyverse)
library(geosphere)
library(fs)
library(parallel)
library(foreach)
library(doParallel)
library(conStruct)

setwd(path.scpt)

# -----------------------------------------------------------------------------
theta <- 0.001
min.samp <- 15
dataset <- 1
n.iter <- 20e3
k.min <- 1
k.max <- 6
train.prop <- 0.90
n.nodes <- 4
# -----------------------------------------------------------------------------
pop.col.name <- 'pop1'
# -----------------------------------------------------------------------------

# ---------------------------------------------------------------------------
dataset.dir <- file.path(path.resamp.snps,
                         paste0('theta-', sprintf('%0.3f', theta)),
                         paste0('min-samp-', sprintf('%02d', min.samp)))
# ---------------------------------------------------------------------------
snps.f.path <- file.path(
  dataset.dir,
  paste0(sprintf('%03d', dataset), '_theta-', sprintf('%0.3f', theta),
         '_min-samp-', sprintf('%02d', min.samp), '.RData'))
# -----------------------------------------------------------------------------
load(snps.f.path)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
geo.coords.df <- select(sample.dat, id, !!pop.col.name, lon, lat) %>%
  filter(!is.na(!!as.symbol(pop.col.name)))

# geo.coords.df %>% filter(id!='08') -> geo.coords.df

geo.coords.df <- geo.coords.df %>% group_by((!!as.symbol(pop.col.name))) %>%
  summarise_at(vars(-id), funs(mean(., na.rm=F)))

geo.coords <- as.matrix(select(geo.coords.df, lon, lat))
rownames(geo.coords) <- geo.coords.df[[pop.col.name]]

pw.geo.dist <- distm(geo.coords, fun=distGeo)
colnames(pw.geo.dist) <- geo.coords.df[[pop.col.name]]
rownames(pw.geo.dist) <- geo.coords.df[[pop.col.name]]

chosen.snps.pop.df <- left_join(
  select(sample.dat, id, !!pop.col.name),
  as_tibble(chosen.snps, rownames='id'), by='id') %>%
  filter(!is.na(!!as.symbol(pop.col.name)))
# chosen.snps.pop.df <- chosen.snps.pop.df %>%
#   select(-sample, -location, -region, -lat, -lon, -lab_adj_x, -lab_adj_y)
# chosen.snps.pop.df %>% filter(id!='08') -> chosen.snps.pop.df
chosen.snps.pop.df <- chosen.snps.pop.df %>% group_by((!!as.symbol(pop.col.name))) %>%
  summarise_at(vars(-id), funs(mean(., na.rm=F)))

chosen.snps <- as.matrix(chosen.snps.pop.df %>% select(-(!!as.symbol(pop.col.name))))
rownames(chosen.snps) <- chosen.snps.pop.df[[pop.col.name]]
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
chosen.snps.tmp <- t(chosen.snps)
tmp <- apply(chosen.snps.tmp, 1, unique)
tmp <- lapply(tmp, na.omit)
n.variants <- sapply(tmp, length)
rm(tmp)
print(table(n.variants))
chosen.snps.invars <- chosen.snps.tmp[n.variants <= 1,]
chosen.snps.vars <- chosen.snps.tmp[n.variants > 1,]
chosen.snps <- t(chosen.snps.vars)
sanity.check <- (nrow(chosen.snps.tmp) - (nrow(chosen.snps.invars) +
                                            nrow(chosen.snps.vars))) == 0
message('           snps: ', nrow(chosen.snps.tmp))
message('invariable snps: ', nrow(chosen.snps.invars))
message('  variable snps: ', nrow(chosen.snps.vars))
message('   sanity check: ', sanity.check)
if (sanity.check != TRUE) {
  stop('Abort: Sanity check failed!')
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
conStruct.dir <- paste0('41-conStruct-cv')
cross.val.run.id <- paste0(sprintf('%02d', min.samp), 'k', k.min, 'k', k.max,
                           'snps')
cross.val.sub.dir <- paste0('cross-val-s', cross.val.run.id)
cross.val.dir.path <- file.path(path.rslt, conStruct.dir, pop.col.name,
                                cross.val.sub.dir)
cross.val.pdf.path <- file.path(path.rslt, conStruct.dir, pop.col.name,
                                paste0(cross.val.sub.dir, '.pdf'))
cross.val.rda.path <- file.path(path.rslt, conStruct.dir, pop.col.name,
                                paste0(cross.val.sub.dir, '.RData'))
cross.val.gt.path <- file.path(path.rslt, conStruct.dir, pop.col.name,
                               paste0(cross.val.sub.dir, '.snps.RData'))
dir_create(cross.val.dir.path)
save(chosen.snps, file=cross.val.gt.path)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
setwd(cross.val.dir.path)
# -----------------------------------------------------------------------------
clust <- makeCluster(n.nodes, type='FORK')
registerDoParallel(clust)
cross.val <- x.validation(
  train.prop=train.prop,
  n.reps=n.nodes,
  K=k.min:k.max,
  freqs=chosen.snps,
  geoDist=pw.geo.dist,
  coords=geo.coords,
  prefix=cross.val.run.id,
  n.iter=n.iter,
  make.figs=T,
  save.files=T,
  parallel=T,
  n.nodes=n.nodes,
  algorithm='NUTS',
  control=list(
    adapt_engaged=FALSE,
    metric='diag_e',
    stepsize=0.0175,
    max_treedepth=10
  ),
)
stopCluster(clust)
save(cross.val, file=cross.val.rda.path)
# -----------------------------------------------------------------------------
# read in results from text files
sp.results <- as.matrix(
  read.table(paste0(cross.val.run.id, '_sp_xval_results.txt'),
             header=T,
             stringsAsFactors=F))
nsp.results <- as.matrix(
  read.table(paste0(cross.val.run.id, '_nsp_xval_results.txt'),
             header=T,
             stringsAsFactors=F))
# -----------------------------------------------------------------------------
# t.test(sp.results[4,],sp.results[3,],paired=TRUE,alternative='greater')
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
setwd(path.scpt)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
pdf(cross.val.pdf.path,
    width=16, height=8,
    pointsize=12,
    bg='transparent',
    colormodel='srgb')
# -----------------------------------------------------------------------------
# first, get the 95% confidence intervals for the spatial and nonspatial
#   models over values of K (mean +/- 1.96 the standard error)
sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

par(mfrow=c(1,3))

plot(rowMeans(sp.results),
     pch=19,col='blue',
     ylab='predictive accuracy',xlab='values of K',
     ylim=range(sp.results,nsp.results),
     main='cross-validation results')
points(rowMeans(nsp.results),col='green',pch=19)

plot(rowMeans(sp.results),
     pch=19,col='blue',
     ylab='predictive accuracy',xlab='values of K',
     ylim=range(sp.CIs),
     main='spatial cross-validation results')
segments(x0=1:nrow(sp.results),
         y0=sp.CIs[1,],
         x1=1:nrow(sp.results),
         y1=sp.CIs[2,],
         col='blue',lwd=2)

plot(rowMeans(nsp.results),
     pch=19,col='green',
     ylab='predictive accuracy',xlab='values of K',
     ylim=range(nsp.CIs),
     main='non-spatial cross-validation results')
segments(x0=1:nrow(nsp.results),
         y0=nsp.CIs[1,],
         x1=1:nrow(nsp.results),
         y1=nsp.CIs[2,],
         col='green',lwd=2)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

par(mfrow=c(n.nodes, 1), mar=c(0, 0, 0, 0))

for (i in 1:n.nodes) {
  # Loop through output files generated by conStruct
  #   runs with k.min through k.max and calculate the
  #   layer contributions for each layer in each run
  layer.contributions <- matrix(NA, nrow=k.max, ncol=k.max)

  # load the conStruct.results.Robj and data.block.Robj
  #   files saved at the end of a conStruct run

  f1 <- file.path(cross.val.dir.path,
                  paste0(cross.val.run.id,
                         '_sp_rep', i, 'K', k.min, '_conStruct.results.Robj'))
  f2 <- file.path(cross.val.dir.path,
                  paste0(cross.val.run.id,
                         '_sp_rep', i, 'K', k.min, '_data.block.Robj'))

  load(f1)
  load(f2)

  # calculate layer contributions
  layer.contributions[,k.min] <- c(
    calculate.layer.contribution(conStruct.results[[1]],data.block),
    rep(0,k.max-k.min))
  tmp <- conStruct.results[[1]]$MAP$admix.proportions

  for(k in 2:k.max){

    f1 <- file.path(cross.val.dir.path,
                    paste0(cross.val.run.id,
                           '_sp_rep', i, 'K', k, '_conStruct.results.Robj'))
    f2 <- file.path(cross.val.dir.path,
                    paste0(cross.val.run.id,
                           '_sp_rep', i, 'K', k, '_data.block.Robj'))

    # load the conStruct.results.Robj and data.block.Robj
    #   files saved at the end of a conStruct run
    load(f1)
    load(f2)

    # match layers up across runs to keep plotting colors consistent
    #   for the same layers in different runs
    tmp.order <- match.layers.x.runs(tmp,conStruct.results[[1]]$MAP$admix.proportions)

    # calculate layer contributions
    layer.contributions[,k] <- c(
      calculate.layer.contribution(conStruct.results=conStruct.results[[1]],
                                   data.block=data.block,
                                   layer.order=tmp.order), rep(0,k.max-k))
    tmp <- conStruct.results[[1]]$MAP$admix.proportions[,tmp.order]
  }

  barplot(layer.contributions,
          col=c('blue', 'red', 'goldenrod1', 'forestgreen', 'darkorchid1', 'pink'),
          xlab=NULL,
          ylab=NULL,
          axes=F,
          names.arg=paste0('K=', 1:k.max)
  )
}
# -----------------------------------------------------------------------------
dev.off()
# -----------------------------------------------------------------------------
