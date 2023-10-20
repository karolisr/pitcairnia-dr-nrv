#!/usr/bin/env Rscript

# Andy Raduski (KR Edited)
# R code to calculate pairwise seq. divergence between samples given a vcf file.

calc.nuc.div <- function(vcf.file.path, f.prefix) {

  f.gt.sn <- paste0(f.prefix, '.02.gt.RData')

  if (!exists('gt') & file.exists(f.gt.sn)) {
    load(f.gt.sn)
  } else {

    # import vcf file
    message('Reading VCF file.')
    vcf <- readVcf(vcf.file.path)
    message('Reading VCF file done.')

    save(vcf, file=paste0(f.prefix, '.01.vcf.RData'), ascii=F, compress=T)

    # get genotypes from vcf file
    gt <- geno(vcf)$GT
    rm(vcf)

    # get names of sites in vcf file
    sn <- gsub(x=rownames(gt), pattern='(\\:\\d+)(.*$)', replacement='\\1', perl=T)
    sn.dupes <- sn[duplicated(sn)]
    sn.dupes.count <- length(sn.dupes)
    if (sn.dupes.count > 0) {
      stop(paste0('Count of duplicate site names: ', sn.dupes.count))
    }
    message('Number of sites: ', length(sn))

    # get names of samples
    samples <- colnames(gt)
    samples <- gsub('NAV', '0', samples, perl=T)
    samples <- gsub('NA', '', samples, perl=T)
    colnames(gt) <- samples
    rownames(gt) <- sn
    samples.sorted <- sort(samples)
    gt <- gt[, samples.sorted]

    save(gt, sn, file=f.gt.sn, ascii=F, compress=T)
  }

  # get names of samples
  samples <- colnames(gt)

  # get names of contigs in vcf file
  cn <- unique(gsub(x=sn, pattern='\\:.*$', replacement='', perl=T))
  # get number of contigs in vcf file
  n.cn <- length(cn)
  message('Number of contigs: ', n.cn)

  # get numbers of samples
  n.samples <- length(samples)
  # get all possible combinations of samples in vcf file
  combos <- combn(x=samples, m=2)
  # get number of combinations
  n.combos <- ncol(combos)

  f.pw.seq.d <- paste0(f.prefix, '.03.pw.seq.d.RData')
  if (!exists('pw.seq.d') & file.exists(f.pw.seq.d)) {
    load(f.pw.seq.d)
  } else {
    # create list of length n.cn to store dataframes for each cn/gene
    pw.seq.d <- vector(mode='list', length=n.cn)
    # create empty dataframes in each element of pw.seq.d list
    for (i in 1:n.cn) {
      pw.seq.d[[i]] <- as.data.frame(matrix(data=0, nrow=n.samples, ncol=n.samples))
      rownames(pw.seq.d[[i]]) <- samples
      colnames(pw.seq.d[[i]]) <- samples
    }
    names(pw.seq.d) <- cn

    # create a vector of full contig length (largest possible denominator for each contig)
    # cn.len <- numeric(length=n.cn)
    # names(cn.len) <- cn

    # for each contig calculate seq divergence
    for (i in 1:n.cn) {
      message('Starting contig ', i, ' of ', n.cn, ': ', cn[i])
      # grab only the genotypes for that contig
      gt.tmp <- gt[grep(pattern=cn[i], x=rownames(gt)),, drop=FALSE]
      # cn.len[i] <- nrow(gt.tmp)
      # for all combinations of samples at that contig
      for (j in 1:n.combos) {
        # grab genotypes for each pairwise comparison
        geno_to_compare <- gt.tmp[, combos[, j], drop=FALSE]
        # turn missing genotypes to NA
        geno_to_compare[geno_to_compare == '.'] <- NA
        # only take sites that are genotyped for both samples
        complete_geno_to_compare <-
          geno_to_compare[complete.cases(geno_to_compare), ]
        # Throws an error if <=1 sites are left
        if (length(complete_geno_to_compare) == 2) {
          complete_geno_to_compare <- t(complete_geno_to_compare)
        }
        if (length(complete_geno_to_compare) >= 2) {
          # get number of sites that are found in both samples
          ns <- nrow(complete_geno_to_compare)
          # make vector to store site based divergence values
          site_div <- numeric(length=ns)
          # calculate pi across sites of contig
          for (k in 1:ns) {
            # Extract site info
            t <- strsplit(x=complete_geno_to_compare[k, ], '/')
            #######################################################################
            a <- t[[1]][1]
            b <- t[[1]][2]
            c <- t[[2]][1]
            d <- t[[2]][2]
            site_div[k] <- 1 - (sum(c(a==c, a==d, b==c, b==d)) / 4)
            #######################################################################
            # a <- t[[1]]
            # b <- t[[2]]
            # la <- length(unique(a))
            # lb <- length(unique(b))
            #
            # #1 -- samples are homozygous for same nucleotide
            # if (identical(a, b) & la == 1) {
            #   # site_div[k] <- 0
            #   # print('1')
            #   next
            # }
            #
            # #2 -- samples are homozygous (or heterozygous) for different nucleotides
            # if (sum(a %in% b, b %in% a) == 0) {
            #   site_div[k] <- 1
            #   # print('2')
            #   next
            # }
            #
            # #3 -- both samples are heterozygous for the same pair of nucleotides
            # if (sum(a %in% b) == 2 & la == 2 & lb == 2) {
            #   site_div[k] <- 0.5
            #   # print('3')
            #   next
            # }
            #
            # #4 -- one sample is het, one is hom and they share one allele
            # if ( ( (la == 1 & lb == 2) | (lb == 1 & la == 2) ) &
            #      (  sum(a %in% b) == 1 | sum(b %in% a) == 1 ) ) {
            #   site_div[k] <- 0.5
            #   # print('4')
            #   next
            # }
            #
            # #5 -- both samples are het and only share one allele
            # if (((la == 2 & lb == 2)) & (sum(a %in% b) == 1 & sum(b %in% a) == 1)) {
            #   site_div[k] <- 0.75
            #   # print('5')
            #   next
            # }
            #######################################################################
          }
          # sum up divergence across the contig
          sum_div <- sum(site_div)
          # plug divergence sum into appropriate element in big list; below the diagonal
          pw.seq.d[[i]][combos[2, j], combos[1, j]] <- sum_div
          # plug length of contig that I have comapred the two sequences over above the diagonal
          pw.seq.d[[i]][combos[1, j], combos[2, j]] <- ns
        }
      }
    }
    save(pw.seq.d, file=f.pw.seq.d, ascii=F, compress=T)
  }

  m.zero <-as.data.frame(matrix(data=0, nrow=n.samples, ncol=n.samples), row.names=samples)
  m.ones <-as.data.frame(matrix(data=1, nrow=n.samples, ncol=n.samples), row.names=samples)
  names(m.zero) <- samples
  names(m.ones) <- samples

  lt <- lower.tri(m.ones)
  ut <- upper.tri(m.ones)
  dg <- as.data.frame(diag(NA, n.samples, n.samples), row.names=samples)
  names(dg) <- samples

  ntdiv.counts <- m.zero
  ntdiv.sums <- m.zero
  for (i in 1:length(pw.seq.d)) {
    m <- pw.seq.d[[i]]
    c <- m
    c[c > 0] <- 1
    ntdiv.counts <- ntdiv.counts + c

    m.t <- t(m)
    m[ut] <- 0
    m.t[ut] <- 1
    m.t <- m.t + dg

    m.div <- m / m.t
    m.div[is.na(m.div)] <- 0

    ntdiv.sums <- ntdiv.sums + m.div
  }

  ntdiv.counts[lt] <- 1
  ntdiv.counts <- t(ntdiv.counts) + dg

  ntdiv.all <- ntdiv.sums / ntdiv.counts
  ntdiv.all[ut] <- NA

  colnames(ntdiv.all) <- samples
  rownames(ntdiv.all) <- samples

  x <- ntdiv.all
  y <- t(x)

  x[is.na(x)] <- 0
  y[is.na(y)] <- 0

  ntdiv.all <- x + y

  save(ntdiv.all, file=paste0(f.prefix, '.04.pw.seq.diff.RData'),
       ascii=F, compress=T)
  write.csv(ntdiv.all, paste0(f.prefix, '.05.pw.seq.diff.csv'))

}

# -----------------------------------------------------------------------------

library(fs)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  stop('Please provide a path to a VCF file.')
}

p <- args[1]

p <- path_expand(p)
p <- path_norm(p)

p.dir <- path_dir(p)
p.file <- path_file(p)

o.dir <- file.path(
  path_dir(p.dir),
  paste0(sub('16-', '17-', path_file(p.dir)), '-R'))

dir_create(o.dir)

# -----------------------------------------------------------------------------

if (file_exists(p) && is_file(p) && dir_exists(o.dir)) {
  message(p)
  setwd(o.dir)

  library(VariantAnnotation)

  vcf.file.path <- p
  f.prefix <- p.file

  rm(p, p.file, p.dir, o.dir)

  calc.nuc.div(vcf.file.path, f.prefix)
}
