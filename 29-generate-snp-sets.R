#!/usr/bin/env Rscript

source("15-shared-code.R")

library(fs)
library(readr)
library(dplyr)

# -----------------------------------------------------------------------------
# args <- commandArgs(trailingOnly=TRUE)
# if (length(args) != 1) {
#   stop('min.samp?')
# }
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
theta <- c('0.001')
# min.samp <- c(args[1])
# min.samp <- c('01')
min.samp <- c('01', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19')
# -----------------------------------------------------------------------------

set.seed(10241983)

for (ms in min.samp) {
  # ---------------------------------------------------------------------------
  message(paste0('--- ', ms, ' ----------------------------------------------'))
  gt.dir <- paste0('17-freebayes-theta-', theta, '-R')
  gt.file <- paste0('nrv-pitcarn.vcf.filt.', ms, '.vcf.02.gt.RData')
  gt.file.path <- file.path(path.rslt, gt.dir, gt.file)
  load(gt.file.path)
  # ---------------------------------------------------------------------------

  # ---------------------------------------------------------------------------
  dataset.dir <- file.path(path.resamp.snps,
                           paste0('theta-', theta),
                           paste0('min-samp-', ms))
  dir_create(dataset.dir)
  # ---------------------------------------------------------------------------

  # ---------------------------------------------------------------------------
  cn <- gsub(x=rownames(gt), pattern='(\\:\\d+)(.*$)', replacement='', perl=T)
  cn.uniq <- unique(cn)
  rm(cn)
  message('                   loci: ', length(cn.uniq))
  # ---------------------------------------------------------------------------

  # ---------------------------------------------------------------------------
  gt[gt=='.'] <- NaN
  gt[gt=='0/1'] <- 0.5
  gt[gt=='1/1'] <- 1
  gt[gt=='0/0'] <- 0
  gt <- matrix(
    as.numeric(gt), ncol=ncol(gt), dimnames=list(rownames(gt), colnames(gt)))
  message('                  sites: ', nrow(gt))
  tmp <- apply(gt, 1, unique)
  tmp <- lapply(tmp, na.omit)
  n.variants <- sapply(tmp, length)
  rm(tmp)
  gt <- gt[n.variants != 1,]
  message('          variable snps: ', nrow(gt))
  # ---------------------------------------------------------------------------

  # ---------------------------------------------------------------------------
  cn <- gsub(x=rownames(gt), pattern='(\\:\\d+)(.*$)', replacement='', perl=T)
  cn.uniq <- unique(cn)
  rm(cn)
  message('loci with variable snps: ', length(cn.uniq))
  # ---------------------------------------------------------------------------

  # ---------------------------------------------------------------------------
  for (i in 1:5) {
    message(sprintf('------- %03d -------', i))
    chosen.snps <- NULL
    for (contig in cn.uniq) {
      contig.snps <- gt[grep(pattern=contig, x=rownames(gt)),, drop=F]
      chosen.idx <- sample.int(nrow(contig.snps), size=1)
      chosen.snp <- contig.snps[chosen.idx,, drop=F]
      chosen.snps <- rbind(chosen.snps, chosen.snp)
    }
    rm(contig.snps)
    rm(chosen.snp)
    message('            chosen snps: ', nrow(chosen.snps))
    chosen.snps <- t(chosen.snps)
    print(sort(unique(as.vector(chosen.snps)), na.last=T))
    snps.f.path <- file.path(
      dataset.dir,
      paste0(sprintf('%03d', i), '_theta-', theta, '_min-samp-', ms, '.RData'))
    save(chosen.snps, file=snps.f.path, ascii=F, compress=T)

    #####    #####    #####    #####    #####    #####    #####    #####    #####
    # load(snps.f.path, verbose = TRUE)
    expected.sample.order <- c(
      "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12",
      "13", "14", "15", "16", "17", "18", "19")
    names.ok <- assertthat::are_equal(row.names(chosen.snps),
                                      expected.sample.order)
    if (!names.ok) {
      message('Sample names are not as expected and/or are not in an expected order.')
    } else {
      message('Sample names are OK and in order.')
      chosen.snps.feems <- as_tibble(chosen.snps * 2)
      snps.f.path.feems <- file.path(
        dataset.dir,
        paste0(sprintf('%03d', i), '_theta-', theta, '_min-samp-', ms,
               '_feems.csv'))
      write_csv(chosen.snps.feems, file = snps.f.path.feems, col_names = F)
    }
    #####    #####    #####    #####    #####    #####    #####    #####    #####
  }
  # ---------------------------------------------------------------------------
  rm(cn.uniq)
}
