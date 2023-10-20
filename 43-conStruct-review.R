#!/usr/bin/env Rscript

source("15-shared-code.R")

suppressPackageStartupMessages(library(conStruct, quietly=T))
suppressPackageStartupMessages(library(fs, quietly=T))

###############################################################################

conStruct.dir.path <- file.path(path.rslt, '42-conStruct', 'pop1')

# -----------------------------------------------------------------------------
setwd(conStruct.dir.path)
# -----------------------------------------------------------------------------
prefix1 <- 'conStruct-run_ms_15-ds_001-k_3-sp-nc_04-ni_0020000'
prefix2 <- 'conStruct-run_ms_15-ds_001-k_4-sp-nc_04-ni_0020000'

load(file.path(prefix1, 'conStruct.results.RData'))
load(file.path(prefix1, 'data.block.RData'))
r1 <- conStruct.results
d1 <- data.block

load(file.path(prefix2, 'conStruct.results.RData'))
load(file.path(prefix2, 'data.block.RData'))
r2 <- conStruct.results
d2 <- data.block

rm(conStruct.results, data.block)
# -----------------------------------------------------------------------------
setwd(path.scpt)
# -----------------------------------------------------------------------------
conStruct::compare.two.runs(r1, d1, r2, d2, 'ZZZ')
# admix.prps.r1.c1 <- r1$chain_1$MAP$admix.proportions
# admix.prps.r1.c2 <- r1$chain_2$MAP$admix.proportions
# admix.prps.r1.c3 <- r1$chain_3$MAP$admix.proportions
# admix.prps.r1.c4 <- r1$chain_4$MAP$admix.proportions
# conStruct::match.layers.x.runs(admix.prps.r1.c1, admix.prps.r1.c2)
# -----------------------------------------------------------------------------

###############################################################################
