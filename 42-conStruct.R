#!/usr/bin/env Rscript

source("15-shared-code.R")

# -----------------------------------------------------------------------------
# suppressPackageStartupMessages(library(argparse))
#
# parser <- ArgumentParser()
# parser$add_argument('-t', '--theta', type='double', default=0.001,
#                     help='theta', metavar='double')
# parser$add_argument('-ms', '--min.samp', type='integer', default=1,
#                     help='min.samp', metavar='integer')
# parser$add_argument('-ds', '--dataset', type='integer', default=1,
#                     help='dataset', metavar='integer')
# parser$add_argument('-k', '--k', type='integer', default=1,
#                     help='k', metavar='integer')
# parser$add_argument('-nc', '--n.chains', type='integer', default=4,
#                     help='n.chains', metavar='integer')
# parser$add_argument('-ni', '--n.iter', type='integer', default=10000,
#                     help='n.iter', metavar='integer')
# args <- parser$parse_args()
#
# theta <- args$theta
# min.samp <- args$min.samp
# dataset <- args$dataset
# k <- args$k
# n.chains <- args$n.chains
# n.iter <- args$n.iter
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
theta <- 0.001
min.samp <- 15
dataset <- 1
# ks <- c(1, 2, 3, 4)
ks <- c(3)
n.chains <- 4
n.iter <- 20e3
# sp.nsp <- c(TRUE, FALSE)
sp.nsp <- c(TRUE)
# -----------------------------------------------------------------------------
pop.col.name <- 'pop1'
# -----------------------------------------------------------------------------

for (spatial in sp.nsp) {
  for (k in ks) {

    message()
    message(
      # 'theta=', theta,
      ' min.samp=', min.samp,
      ' dataset=', dataset,
      ' k=', k,
      ' spatial=', spatial,
      ' n.chains=', n.chains,
      ' n.iter=', n.iter)
    message()

    suppressPackageStartupMessages(library(fs, quietly=T))

    # -----------------------------------------------------------------------------
    dataset.dir <- file.path(path.resamp.snps,
                             paste0('theta-', sprintf('%0.3f', theta)),
                             paste0('min-samp-', sprintf('%02d', min.samp)))
    snps.f.path <- file.path(
      dataset.dir,
      paste0(sprintf('%03d', dataset),
             '_theta-', sprintf('%0.3f', theta),
             '_min-samp-', sprintf('%02d', min.samp), '.RData'))
    # -----------------------------------------------------------------------------

    message('Loading: ', snps.f.path)
    message()

    if (!file_exists(snps.f.path)) {
      stop('SNPs file not found.')
    }

    load(snps.f.path)

    suppressPackageStartupMessages(library(readxl, quietly=T))
    suppressPackageStartupMessages(library(tidyverse, quietly=T, verbose=F))
    suppressPackageStartupMessages(library(geosphere, quietly=T))
    suppressPackageStartupMessages(library(conStruct, quietly=T))

    # -----------------------------------------------------------------------------
    setwd(path.scpt)
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    geo.coords.df <- select(sample.dat, id, !!pop.col.name, lon, lat) %>%
      filter(!is.na(!!as.symbol(pop.col.name)))

    # geo.coords.df %>% filter(id!='08') -> geo.coords.df

    geo.coords.df <- geo.coords.df %>% group_by((!!as.symbol(pop.col.name))) %>%
      summarise_at(vars(-id), list(~ mean(., na.rm=F)))

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
      summarise_at(vars(-id), list(~ mean(., na.rm=F)))

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
    conStruct.dir <- paste0('42-conStruct')
    conStruct.run.id <- sprintf('ms_%02d-ds_%03d-k_%01d-%s-nc_%02d-ni_%07d',
                                min.samp, dataset, k,
                                ifelse(spatial==T, 'sp', 'nsp'),
                                n.chains, n.iter)
    conStruct.dir.path <- file.path(path.rslt, conStruct.dir, pop.col.name,
                                    paste0('conStruct-run_', conStruct.run.id))
    dir_create(conStruct.dir.path)
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    message(conStruct.run.id)
    setwd(conStruct.dir.path)
    # -----------------------------------------------------------------------------

    # --- OVERLOADING THE conStruct FUNCTION --------------------------------------
    # conStruct <- function(spatial=T, K, freqs, geoDist=NULL, coords,
    #                       n.chains=1, n.iter=1000,
    #                       make.figs=T, save.files=T, ...)
    # {
    #   call.check <- conStruct:::check.call(args=as.list(environment()))
    #   freq.data <- conStruct:::process.freq.data(freqs)
    #   data.block <- conStruct:::make.data.block(K, freq.data, coords, spatial,
    #                                             geoDist)
    #   if (save.files) {
    #     save(data.block,
    #          file=paste0('data.block.RData'), ascii=F, compress=T)
    #   }
    #
    #   stan.model <- conStruct:::pick.stan.model(spatial, K)
    #
    #   model.fit <- rstan::sampling(
    #     object=conStruct:::stanmodels[[stan.model]],
    #     refresh=min(n.iter/10, 500),
    #     data=data.block,
    #     iter=n.iter,
    #     chains=n.chains,
    #     # thin=ifelse(n.iter/500 > 1, n.iter/500, 1),
    #     sample_file=paste0('sample_file'),
    #     diagnostic_file=paste0('diagnostic_file'),
    #     ...
    #   )
    #
    #   conStruct.results <- conStruct:::get.conStruct.results(data.block, model.fit,
    #                                                          n.chains)
    #   data.block <- conStruct:::unstandardize.distances(data.block)
    #   if (save.files) {
    #     save(data.block,
    #          file=paste0('data.block.RData'), ascii=F, compress=T)
    #     save(model.fit,
    #          file=paste0('model.fit.RData'), ascii=F, compress=T)
    #     save(conStruct.results,
    #          file=paste0('conStruct.results.RData'), ascii=F, compress=T)
    #   }
    #   if (make.figs) {
    #     conStruct:::make.all.the.plots(conStruct.results, data.block, '',
    #                                    layer.colors=NULL)
    #   }
    #   return(conStruct.results)
    # }
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    # warmup.fraction <- 1/4
    # # thin <- max(10 ** (floor(log10(n.iter * (1-warmup.fraction))) - 3) * 2, 1)
    # # thin <- floor((1 - warmup.fraction) * n.iter / 500)
    # thin <- 5

    # conStruct.results <- conStruct(
    #   spatial=spatial,
    #   K=as.numeric(k),
    #   freqs=chosen.snps,
    #   geoDist=pw.geo.dist,
    #   coords=geo.coords,
    #   n.chains=as.numeric(n.chains),
    #   n.iter=as.numeric(n.iter),
    #   make.figs=TRUE,
    #   save.files=TRUE,
    #   cores=as.numeric(n.chains),
    #   # cores=1,
    #   verbose=TRUE,
    #
    #   init='random',
    #   seed=1,
    #
    #   # algorithm='HMC',
    #   # control=list(int_time=1)
    #
    #   algorithm='NUTS',
    #   control=list(
    #     adapt_engaged=FALSE,
    #     metric='diag_e',
    #     stepsize=0.0175,
    #     # stepsize_jitter=1,
    #     # adapt_gamma=0.05,
    #     # adapt_delta=0.8,
    #     # adapt_kappa=0.75,
    #     # adapt_t0=10,
    #     max_treedepth=10
    #   ),
    #
    #   save_warmup=F,
    #   warmup=floor(n.iter * warmup.fraction),
    #   thin=thin
    # )

    #
    # control
    #
    # A named list of parameters to control the sampler's behavior. It defaults to
    # NULL so all the default values are used. First, the following are adaptation
    # parameters for sampling algorithms. These are parameters used in Stan with
    # similar names here.
    #
    # adapt_engaged (logical)
    # adapt_gamma (double, positive, defaults to 0.05)
    # adapt_delta (double, between 0 and 1, defaults to 0.8)
    # adapt_kappa (double, positive, defaults to 0.75)
    # adapt_t0 (double, positive, defaults to 10)
    # adapt_init_buffer (integer, positive, defaults to 75)
    # adapt_term_buffer (integer, positive, defaults to 50)
    # adapt_window (integer, positive, defaults to 25)
    #
    # In addition, algorithm HMC (called 'static HMC' in Stan) and NUTS
    # share the following parameters:
    #
    # stepsize (double, positive, defaults to 1) Note: this controls the initial stepsize only, unless adapt_engaged=FALSE.
    # stepsize_jitter (double, [0,1], defaults to 0)
    # metric (string, one of "unit_e", "diag_e", "dense_e", defaults to "diag_e")
    #
    # For algorithm NUTS, we can also set:
    # max_treedepth (integer, positive, defaults to 10)
    #
    # For algorithm HMC, we can also set:
    # int_time (double, positive)
    #
    # For test_grad mode, the following parameters can be set:
    # epsilon (double, defaults to 1e-6)
    # error (double, defaults to 1e-6)
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    # file_move(
    #   paste0(conStruct.run.id, '_conStruct.results.Robj'),
    #   paste0(conStruct.run.id, '_conStruct.results.RData'))
    # file_move(
    #   paste0(conStruct.run.id, '_data.block.Robj'),
    #   paste0(conStruct.run.id, '_data.block.RData'))
    # file_move(
    #   paste0(conStruct.run.id, '_model.fit.Robj'),
    #   paste0(conStruct.run.id, '_model.fit.RData'))
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    # load(paste0(conStruct.run.id, '_conStruct.results.RData'))
    # load(paste0(conStruct.run.id, '_data.block.RData'))
    load('conStruct.results.RData')
    load('data.block.RData')
    r1 <- conStruct.results
    d1 <- data.block

    r1.c1 <- r1$chain_1
    r1.c2 <- r1$chain_2
    r1.c3 <- r1$chain_3
    r1.c4 <- r1$chain_4

    admix.prps.r1.c1 <- r1.c1$MAP$admix.proportions
    admix.prps.r1.c2 <- r1.c2$MAP$admix.proportions
    admix.prps.r1.c3 <- r1.c3$MAP$admix.proportions
    admix.prps.r1.c4 <- r1.c4$MAP$admix.proportions
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    r.matched <- r1

    m <- conStruct::match.layers.x.runs(admix.prps.r1.c1, admix.prps.r1.c2)
    r.matched$chain_2$posterior$layer.params[[1]] <- r1.c2$posterior$layer.params[[m[1]]]
    r.matched$chain_2$posterior$layer.params[[2]] <- r1.c2$posterior$layer.params[[m[2]]]
    r.matched$chain_2$posterior$layer.params[[3]] <- r1.c2$posterior$layer.params[[m[3]]]
    r.matched$chain_2$MAP$layer.params[[1]] <- r1.c2$MAP$layer.params[[m[1]]]
    r.matched$chain_2$MAP$layer.params[[2]] <- r1.c2$MAP$layer.params[[m[2]]]
    r.matched$chain_2$MAP$layer.params[[3]] <- r1.c2$MAP$layer.params[[m[3]]]
    r.matched$chain_2$MAP$admix.proportions <- r1.c2$MAP$admix.proportions[,m]

    m <- conStruct::match.layers.x.runs(admix.prps.r1.c1, admix.prps.r1.c3)
    r.matched$chain_3$posterior$layer.params[[1]] <- r1.c3$posterior$layer.params[[m[1]]]
    r.matched$chain_3$posterior$layer.params[[2]] <- r1.c3$posterior$layer.params[[m[2]]]
    r.matched$chain_3$posterior$layer.params[[3]] <- r1.c3$posterior$layer.params[[m[3]]]
    r.matched$chain_3$MAP$layer.params[[1]] <- r1.c3$MAP$layer.params[[m[1]]]
    r.matched$chain_3$MAP$layer.params[[2]] <- r1.c3$MAP$layer.params[[m[2]]]
    r.matched$chain_3$MAP$layer.params[[3]] <- r1.c3$MAP$layer.params[[m[3]]]
    r.matched$chain_3$MAP$admix.proportions <- r1.c3$MAP$admix.proportions[,m]

    m <- conStruct::match.layers.x.runs(admix.prps.r1.c1, admix.prps.r1.c4)
    r.matched$chain_4$posterior$layer.params[[1]] <- r1.c4$posterior$layer.params[[m[1]]]
    r.matched$chain_4$posterior$layer.params[[2]] <- r1.c4$posterior$layer.params[[m[2]]]
    r.matched$chain_4$posterior$layer.params[[3]] <- r1.c4$posterior$layer.params[[m[3]]]
    r.matched$chain_4$MAP$layer.params[[1]] <- r1.c4$MAP$layer.params[[m[1]]]
    r.matched$chain_4$MAP$layer.params[[2]] <- r1.c4$MAP$layer.params[[m[2]]]
    r.matched$chain_4$MAP$layer.params[[3]] <- r1.c4$MAP$layer.params[[m[3]]]
    r.matched$chain_4$MAP$admix.proportions <- r1.c4$MAP$admix.proportions[,m]
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    admix.props.1 <- r.matched$chain_1$MAP$admix.proportions
    admix.props.2 <- r.matched$chain_2$MAP$admix.proportions
    admix.props.3 <- r.matched$chain_3$MAP$admix.proportions
    admix.props.4 <- r.matched$chain_4$MAP$admix.proportions

    admix.props <- (admix.props.1 + admix.props.2 + admix.props.3 + admix.props.4) / 4
    # -----------------------------------------------------------------------------

    conStruct:::make.all.the.plots(r.matched, d1, '', layer.colors=NULL)

    # -----------------------------------------------------------------------------
    pdf(paste0(conStruct.run.id, '.pdf'),
        width=11, height=8,
        pointsize=12,
        bg='transparent',
        colormodel='srgb')
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    # make a STRUCTURE plot using the
    #   maximum a posteriori (MAP) estimates
    make.structure.plot(admix.proportions=admix.props.1,
                        sample.names=row.names(d1$coords),
                        mar=c(2, 4, 0, 0))

    make.structure.plot(admix.proportions=admix.props.2,
                        sample.names=row.names(d1$coords),
                        mar=c(2, 4, 0, 0))

    make.structure.plot(admix.proportions=admix.props.3,
                        sample.names=row.names(d1$coords),
                        mar=c(2, 4, 0, 0))

    make.structure.plot(admix.proportions=admix.props.4,
                        sample.names=row.names(d1$coords),
                        mar=c(2, 4, 0, 0))

    make.structure.plot(admix.proportions=admix.props,
                        sample.names=row.names(d1$coords),
                        mar=c(2, 4, 0, 0))

    # add pie plot to an existing map
    # make the desired map
    # -----------------------------------------------------------------------------
    maps::map(xlim=range(d1$coords[,1]) + c(-0.15, 0.15),
              ylim=range(d1$coords[,2]) + c(-0.15, 0.15),
              col='gray')

    # add the admixture pie plot
    make.admix.pie.plot(admix.proportions=admix.props.1,
                        coords=d1$coords,
                        add=T)
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    maps::map(xlim=range(d1$coords[,1]) + c(-0.15, 0.15),
              ylim=range(d1$coords[,2]) + c(-0.15, 0.15),
              col='gray')

    # add the admixture pie plot
    make.admix.pie.plot(admix.proportions=admix.props.2,
                        coords=d1$coords,
                        add=T)
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    maps::map(xlim=range(d1$coords[,1]) + c(-0.15, 0.15),
              ylim=range(d1$coords[,2]) + c(-0.15, 0.15),
              col='gray')

    # add the admixture pie plot
    make.admix.pie.plot(admix.proportions=admix.props.3,
                        coords=d1$coords,
                        add=T)
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    maps::map(xlim=range(d1$coords[,1]) + c(-0.15, 0.15),
              ylim=range(d1$coords[,2]) + c(-0.15, 0.15),
              col='gray')

    # add the admixture pie plot
    make.admix.pie.plot(admix.proportions=admix.props.4,
                        coords=d1$coords,
                        add=T)
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    maps::map(xlim=range(d1$coords[,1]) + c(-0.15, 0.15),
              ylim=range(d1$coords[,2]) + c(-0.15, 0.15),
              col='gray')

    # add the admixture pie plot
    make.admix.pie.plot(admix.proportions=admix.props,
                        coords=d1$coords,
                        add=T)
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    par(mfrow=c(1, length(r.matched)), mar=c(2, 2, 1, 0))
    breaks <- 50

    for(c in 1:length(r.matched)) {
      hist(r.matched[[c]]$posterior$gamma,
           main=paste0('gamma, chain: ', c), breaks=breaks)
    }

    for(l in 1:length(r.matched[[c]]$posterior$layer.params)) {

      for(c in 1:length(r.matched)) {
        hist(r.matched[[c]]$posterior$layer.params[[l]]$alpha0,
             main=paste0('alpha0, chain: ', c, ', layer: ', l), breaks=breaks)
      }

      for(c in 1:length(r.matched)) {
        hist(r.matched[[c]]$posterior$layer.params[[l]]$alphaD,
             main=paste0('alphaD, chain: ', c, ', layer: ', l), breaks=breaks)
      }

      for(c in 1:length(r.matched)) {
        hist(r.matched[[c]]$posterior$layer.params[[l]]$alpha2,
             main=paste0('alpha2, chain: ', c, ', layer: ', l), breaks=breaks)
      }

      for(c in 1:length(r.matched)) {
        hist(r.matched[[c]]$posterior$layer.params[[l]]$phi,
             main=paste0('phi, chain: ', c, ', layer: ', l), breaks=breaks)
      }
    }

    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    tmp <- dev.off()
    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    setwd(path.scpt)
    # -----------------------------------------------------------------------------

  }
}
