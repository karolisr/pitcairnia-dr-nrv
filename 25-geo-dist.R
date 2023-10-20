#!/usr/bin/env Rscript

source('15-shared-code.R')

library(readxl)
library(rgbif)
library(tidyverse)
library(ggmap)
library(geosphere)
library(ggplot2)
library(sf)
library(conStruct)
library(fs)

setwd(path.data.rdata)

if (!exists('raw.gbif') & file.exists('25-geo-dist_raw.gbif_.RData')) {
  load('25-geo-dist_raw.gbif_.RData')
}

if (!exists('pitcairnia.gbif') & file.exists('25-geo-dist_pitcairnia.gbif_.RData')) {
  load('25-geo-dist_pitcairnia.gbif_.RData')
}

if (!exists('sample.dat') & file.exists('25-geo-dist_sample.dat_.RData')) {
  load('25-geo-dist_sample.dat_.RData')
}

setwd(path.data)

if (!exists('raw.gbif')) {
  raw.gbif <- occ_data(scientificName='Pitcairnia', country=c('DO', 'HT'))
}

pitcairnia.gbif <- rows_append(raw.gbif$DO$data, raw.gbif$HT$data)
pitcairnia.gbif %>%
  rename(lat=decimalLatitude, lon=decimalLongitude) -> pitcairnia.gbif

lim.left <- -72
lim.right <- -69
lim.top <- 20
lim.bottom <- 17.5

pitcairnia.gbif %>%
  filter(lon > lim.left) %>%
  filter(lon < lim.right) %>%
  filter(lat < lim.top) %>%
  filter(lat > lim.bottom) -> pitcairnia.gbif

pitcairnia.gbif %>% arrange(scientificName, lat, lon) -> pitcairnia.gbif

###############################################################################
# FEEMS -----------------------------------------------------------------------
feems_nodes <- read_csv('feems_nodes.csv', col_names='node_id', col_types='i')
feems_node_pos <- read_csv('feems_node_pos.csv', col_names=c('lon', 'lat'), col_types='dd')
feems_edges <- read_csv('feems_edges.csv', col_names=c('n1', 'n2'), col_types='ii')
feems_w <- read_csv('feems_w.csv', col_names='w', col_types='d')

nodes <- add_column(feems_nodes, feems_node_pos)
edges <- add_column(feems_edges, feems_w)

edges_n1 <- left_join(edges, nodes, by=c('n1' = 'node_id'))
edges_n2 <- left_join(edges_n1, nodes, by=c('n2' = 'node_id'))
feems.edges <- select(edges_n2, lon1=lon.x, lat1=lat.x, lon2=lon.y, lat2=lat.y, w)
feems.mid.w = ((max(feems.edges$w) - min(feems.edges$w)) / 2) + min(feems.edges$w)
rm(feems_nodes, feems_node_pos, feems_edges, feems_w, nodes, edges, edges_n1, edges_n2)
# END FEEMS -------------------------------------------------------------------
###############################################################################

###############################################################################
# conStruct -------------------------------------------------------------------
conStruct.dir.path <- file.path(path.rslt, '42-conStruct', 'pop1')
# -----------------------------------------------------------------------------
setwd(conStruct.dir.path)
# -----------------------------------------------------------------------------
prefix1 <- 'conStruct-run_ms_15-ds_001-k_3-sp-nc_04-ni_0020000'
load(file.path(prefix1, 'conStruct.results.RData'))
load(file.path(prefix1, 'data.block.RData'))
r1 <- conStruct.results
d1 <- data.block
rm(conStruct.results, data.block)
# -----------------------------------------------------------------------------
setwd(path.scpt)
# -----------------------------------------------------------------------------
r1.c1 <- r1$chain_1
r1.c2 <- r1$chain_2
r1.c3 <- r1$chain_3
r1.c4 <- r1$chain_4

admix.props.r1.c1 <- r1.c1$MAP$admix.proportions
admix.props.r1.c2 <- r1.c2$MAP$admix.proportions
admix.props.r1.c3 <- r1.c3$MAP$admix.proportions
admix.props.r1.c4 <- r1.c4$MAP$admix.proportions
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
r.matched <- r1
m <- conStruct::match.layers.x.runs(admix.props.r1.c1, admix.props.r1.c2)
r.matched$chain_2$MAP$admix.proportions <- r1.c2$MAP$admix.proportions[,m]
m <- conStruct::match.layers.x.runs(admix.props.r1.c1, admix.props.r1.c3)
r.matched$chain_3$MAP$admix.proportions <- r1.c3$MAP$admix.proportions[,m]
m <- conStruct::match.layers.x.runs(admix.props.r1.c1, admix.props.r1.c4)
r.matched$chain_4$MAP$admix.proportions <- r1.c4$MAP$admix.proportions[,m]
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
admix.props.1 <- r.matched$chain_1$MAP$admix.proportions
admix.props.2 <- r.matched$chain_2$MAP$admix.proportions
admix.props.3 <- r.matched$chain_3$MAP$admix.proportions
admix.props.4 <- r.matched$chain_4$MAP$admix.proportions

admix.props <- (admix.props.1 + admix.props.2 + admix.props.3 + admix.props.4) / 4
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
colnames(admix.props) <- c('L1', 'L2', 'L3')
admix.props <- cbind(admix.props, d1$coords)
admix.props <- as_tibble(admix.props)
sample.dat <- left_join(sample.dat, admix.props, by=c('lon', 'lat'), copy=TRUE)
# END conStruct ---------------------------------------------------------------
###############################################################################

###############################################################################
###############################################################################
# make.conStruct.admix.pie.list.kr.ggplot <- function(
    #     admix.proportions,
#     coords) {
#   if (!inherits(admix.proportions, 'matrix')) {
#     stop('\nyou must specify a matrix of admixture proportions\n')
#   }
#   if (is.null(coords)) {
#     stop('\nyou must specify sampling coordinates\n')
#   }
#   else {
#     if (nrow(coords) != nrow(admix.proportions)) {
#       stop('\nyou must specify one set of coordinates for each row of the admixture proportion matrix\n')
#     }
#   }
#   K <- ncol(admix.proportions)
#   N <- nrow(admix.proportions)
#   # layer.colors <- c('blue', 'red', 'goldenrod1', 'forestgreen', 'darkorchid1', 'deepskyblue', 'darkorange1', 'seagreen2', 'yellow1', 'black')
#   layer.names <- paste0('layer_', 1:K)
#   sample.names <- paste0('sample_', 1:N)
#   # color.tab <- caroline::nv(c(layer.colors[1:K]), layer.names)
#   pie.list <- lapply(1:N, function(i) {
#     caroline::nv(admix.proportions[i, ], layer.names)
#   })
#   names(pie.list) <- sample.names
#   # if (add) {
#   #   graphics::par(new = TRUE)
#   # }
#   # else {
#   #   graphics::par(mar = mar)
#   # }
#   # if (is.null(x.lim)) {
#   #   x.lim <- c(min(coords[, 1]) - 1, max(coords[, 1]) + 1)
#   # }
#   # if (is.null(y.lim)) {
#   #   y.lim <- c(min(coords[, 2]) - 1, max(coords[, 2]) + 1)
#   # }
#   # suppressWarnings(caroline::pies(pie.list, x0 = coords[, 1],
#   #                                 y0 = coords[, 2], color.table = color.tab, border = 'black',
#   #                                 radii = radii, xlab = '', ylab = '', main = '', lty = 1,
#   #                                 density = NULL, xlim = x.lim, ylim = y.lim))
#   # return(invisible(0))
#   return(pie.list)
# }
# make.admix.pie.plot.kr.ggplot(admix.proportions=admix.props, coords=d1$coords)
###############################################################################
###############################################################################v



# library(sp)
# library(terra)
# library(geodata)
# library(dismo)
# library(scatterpie)
# library(ggnewscale)
# library(raster)



# if (!file.exists('25-geo-dist-map-samples.pdf')) {

# setwd(path.data.rdata)
# if (file.exists('25-geo-dist_plt.map.samples_.RData')) {
#   load('25-geo-dist_plt.map.samples_.RData')
# } else {

feems.colors <- c(
  '#994000', '#CC5800', '#FF8F33', '#FFAD66', '#FFCA99', '#FFE6CC', '#FBFBFB',
           '#CCFDFF', '#99F8FF', '#66F0FF', '#33E4FF', '#00AACC', '#007A99')

# construct.colors <- c('blue', 'red', 'goldenrod1', 'forestgreen', 'darkorchid1', 'deepskyblue', 'darkorange1', 'seagreen2', 'yellow1', 'black')
# construct.colors.desaturated <- shades::saturation(construct.colors, scalefac(0.5))
# fair.cols <- c('#BF1B0B', '#FFC465', '#66ADE5', '#252A52', '#38170B')
# fair.cols.desaturated <- shades::saturation(fair.cols, scalefac(0.5))
construct.colors <- c('#66c2a5', '#fc8d62', '#8da0cb')

bbox.pitcairnia <- c(left=lim.left, right=lim.right,
                     bottom=lim.bottom, top=lim.top)

# -----------------------------------------------------------------------------
hti.elev <- geodata::worldclim_country(country='HTI', var='elev', path=path.data.clim)
dom.elev <- geodata::worldclim_country(country='DOM', var='elev', path=path.data.clim)
# -----------------------------------------------------------------------------
hti.elev.df <- raster::as.data.frame(hti.elev, xy=TRUE)
dom.elev.df <- raster::as.data.frame(dom.elev, xy=TRUE)
# -----------------------------------------------------------------------------

hisp <- st_read('/Users/karolis/SyncThing/nrv-dr-pitcairnia/KR/pop-gen-01-scripts/data/dggs/hispanola.shp')

plt.map.samples <- ggplot(hisp) +
  # xlim(bbox.pitcairnia[1], bbox.pitcairnia[2]) +
  # ylim(bbox.pitcairnia[3], bbox.pitcairnia[4]) +
  geom_tile(data=dom.elev.df, aes(x=x, y=y, fill=DOM_wc2.1_30s_elev)) +
  geom_tile(data=hti.elev.df, aes(x=x, y=y, fill=HTI_wc2.1_30s_elev)) +
  scale_fill_gradientn(colors=c('white', 'black'), guide='none') +
  geom_sf(color='black', fill=NA, linewidth=0.5) +

  # geom_segment(
  #   mapping=aes(x=lon1, xend=lon2, y=lat1, yend=lat2, color=w),
  #   data=feems.edges,
  #   linewidth=0.4) +
  #
  # scale_color_gradientn(colors=feems.colors) +

  # geom_point(pitcairnia.gbif,
  #            color='red',
  #            size=2,
  #            shape=3,
  #            mapping=aes(x=lon, y=lat)) +

  # new_scale_fill() +
  # geom_scatterpie(data=sample.dat,
  #                 cols=c('L1', 'L2', 'L3'),
  #                 color='white', alpha=1.0, linewidth=0.3,
  #                 legend_name='Layer',
  #                 mapping=aes(x=lon, y=lat, group=id, r=0.06)) +
  #
  # scale_fill_manual(values=construct.colors, guide='none') +
  # scale_fill_brewer(palette='Set1') +

  # geom_label(sample.dat,
  #            color='black',
  #            fill='NA',
  #            alpha=0.75,
  #            nudge_x=sample.dat$lab_adj_x * 0.056,
  #            nudge_y=sample.dat$lab_adj_y * 0.046,
  #            label.size=0,
  #            size=3,
  #            mapping=aes(x=lon, y=lat, label=id)) +

  # geom_point(sample.dat,
  #            color='#EEEEEE',
  #            fill='black',
  #            size=2,
  #            shape=21,
  #            mapping=aes(x=lon, y=lat)) +

  coord_sf() +
  labs(x=NULL, y=NULL) +
  # theme_classic()
theme_nothing()
plt.map.samples

# stamen.maptypes <- c(
#   'terrain',
#   'terrain-background',
#   'terrain-labels',
#   'terrain-lines',
#   'toner',
#   'toner-2010',
#   'toner-2011',
#   'toner-background',
#   'toner-hybrid',
#   'toner-labels',
#   'toner-lines',
#   'toner-lite',
#   'watercolor')

# domr <- map_data('world', c('.'), resolution=0) %>%
#   select(lon=long, lat, group, id=region) %>% filter(group==2)

# ggplot(data=domr, aes(lon, lat)) +
# xlim(bbox.pitcairnia[1]-3, bbox.pitcairnia[2]+1) +
# ylim(bbox.pitcairnia[4], bbox.pitcairnia[3])

# ggmap(get_stamenmap(bbox.pitcairnia, zoom=9, color='bw', force=TRUE,
#                     maptype=stamen.maptypes[2]))

# geom_polygon(aes(group=group)) +
# geom_segment(aes(x=lon1, xend=lon2, y=lat1, yend=lat2, color=w), data=feems.edges) +
# scale_color_gradientn(colours=colorspace::heat_hcl(7)) +
# labs(x=NULL, y=NULL) +
# theme_nothing()

# plt.map.samples <- ggmap(
# plt.map.samples <- ggplot(


# get_map(bbox.pitcairnia,
#         maptype='terrain-background')

# get_stamenmap(bbox.pitcairnia,
#               zoom=9,
#               color='bw',
#               force=TRUE,
#               maptype=stamen.maptypes[2])

#   ) +
#
# geom_segment(aes(x=lon1, xend=lon2, y=lat1, yend=lat2, color=w), data=feems.edges) +
# scale_color_gradientn(colours=colorspace::heat_hcl(7)) +

# scale_color_gradient2(midpoint=feems.mid.w,
#                       low='red',
#                       mid='white',
#                       high='blue') +

# labs(x=NULL, y=NULL) +
# theme_classic()
# theme_nothing()

# if (!file.exists('25-geo-dist_plt.map.samples_.RData')) {
#   save(plt.map.samples, file='25-geo-dist_plt.map.samples_.RData')
# }
# }

# setwd(path.data)

# plot(plt.map.samples)

ggsave('25-geo-dist-map-samples.pdf',
       plt.map.samples,
       width=9, height=6,
       scale=1.25)

# }

# setwd(path.data.rdata)
#
# save(raw.gbif, file='25-geo-dist_raw.gbif_.RData')
# save(pitcairnia.gbif, file='25-geo-dist_pitcairnia.gbif_.RData')
# save(sample.dat, file='25-geo-dist_sample.dat_.RData')

# -----------------------------------------------------------------------------
# setwd(path.data)
# pw.geo.dist <- distm(select(sample.dat, lon, lat), fun=distGeo)
# colnames(pw.geo.dist) <- sample.dat$id
# rownames(pw.geo.dist) <- sample.dat$id
# write.csv(pw.geo.dist, 'pw-geo-dist.csv')

# pdf('25-geo-dist-ntdiff-vs-geodist.pdf',
#   width=8, height=8,
#   pointsize=12,
#   bg='transparent',
#   colormodel='srgb')

# thetas <- c('0.001', '0.00001')
# theta <- thetas[1]
# min.samp.v <- c(1, 10:19)

# for (i in 1:length(min.samp.v)) {

#   min.samp <- formatC(min.samp.v[i], mode='integer', format='d', width=2, flag='0')
#   r.dir <- paste0('17-freebayes-theta-', theta, '-R')
#   r.file <- paste0('nrv-pitcarn.vcf.filt.', min.samp, '.vcf.05.pw.seq.diff.csv')
#   pw.nuc.div.f <- file.path('..', '..', 'pop-gen-02-results', r.dir, r.file)

#   pw.nuc.div <- as.matrix(read.csv(pw.nuc.div.f, row.names=1))
#   colnames(pw.nuc.div) <- sample.dat$id
#   rownames(pw.nuc.div) <- sample.dat$id

#   lt <- lower.tri(pw.nuc.div)

#   fit <- lm(pw.nuc.div[lt] ~ pw.geo.dist[lt])
#   print(summary(fit))

#   rsqr <- paste('R^{2}=', round(summary(fit)$r.squared, digits=2), sep='')
#   message(i, ' :: ', rsqr, ' :: ', cor(pw.geo.dist[lt], pw.nuc.div[lt]))

#   aspect.ratio <- (
#     (max(pw.geo.dist[lt]) - min(pw.geo.dist[lt])) /
#       (max(pw.nuc.div[lt]) - min(pw.nuc.div[lt]))
#   )

#   fg <- rgb(0, 0, 0, 0.15)

#   par(pty='s')
#   plot(pw.geo.dist[lt], pw.nuc.div[lt],
#        # ylim=c(min(pw.nuc.div[lt]), max(pw.nuc.div[lt])),
#        # xlim=c(min(pw.geo.dist[lt]), max(pw.geo.dist[lt])),
#        asp=aspect.ratio,
#        pch=21,
#        col='black',
#        bg=fg,
#        xlab='Distance (m)',
#        ylab='Pairwise nuc. diff. per base pair',
#        main=paste0('Minimum samples included: ', min.samp))

#   abline(fit)

#   # The average numbers of pairwise nucleotide differences per base pair at
#   # 4-fold degenerate sites (𝜋) between samples within regions were calculated
#   # to estimate genetic diversity.
# }
# dev.off()
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
setwd(path.scpt)
# -----------------------------------------------------------------------------
