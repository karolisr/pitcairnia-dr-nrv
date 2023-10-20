source("15-shared-code.R")

library(SNPRelate)
library(gdsfmt)
library(tidyverse)
library(plotly)
library(PCAviz)
library(ggpubr)
library(vegan)
library(readxl)
library(geosphere)
library(elevatr)
library(fs)
library(sf)
library(rgdal)

# -----------------------------------------------------------------------------
setwd(path.data)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
geo.coords <- as.matrix(select(sample.dat, lon, lat))
rownames(geo.coords) <- sample.dat$id
# write.table(geo.coords, "sample_coords", quote=F, sep=",", row.names=F, col.names=F)
geo.coords.df <- as_tibble(geo.coords, rownames='sample')
elev <- get_elev_point(as.data.frame(geo.coords), prj="EPSG:4326", src='aws')
geo.coords.df %>% add_column(elev=elev$elevation) -> geo.coords.df
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# library(SNPRelate)
# snpgdsVCF2GDS('nrv-pitcarn.vcf.filt.01.vcf', 'nrv-pitcarn.vcf.filt.01.gds')
genofile <- snpgdsOpen(filename='/Users/karolis/Library/CloudStorage/GoogleDrive-kraman2@uic.edu/Shared drives/nrv-dr-pitcairnia/KR/pop-gen-02-results/16-freebayes-theta-0.001/nrv-pitcarn.vcf.filt.01.gds', allow.duplicate=T)
# -----------------------------------------------------------------------------

set.seed(seed=39)

snpset <- snpgdsLDpruning(
  genofile,
  sample.id=NULL,
  ld.threshold=0.0,
  missing.rate=1/19 + 1e-10,
  autosome.only=F,
  verbose=F)

snpset.id <- unlist(snpset)

# snpgdsPCA {SNPRelate}
pca.results <- snpgdsPCA(
  genofile,
  snp.id=snpset.id,
  autosome.only=F,
  maf=0.05)

##### !!! #####
# snpgdsAdmixProp(eigobj, groups, bound=FALSE)
##### !!! #####

# -----------------------------------------------------------------------------
snpgdsClose(genofile)
# -----------------------------------------------------------------------------

pca.results.df <- tibble(
  sample=pca.results$sample.id,
  pc1=pca.results$eigenvect[,1],
  pc2=pca.results$eigenvect[,2],
  pc3=pca.results$eigenvect[,3],
  pc4=pca.results$eigenvect[,4],
  pc5=pca.results$eigenvect[,5],
  pc6=pca.results$eigenvect[,6],
  pc7=pca.results$eigenvect[,7],
  pc8=pca.results$eigenvect[,8],
  pc9=pca.results$eigenvect[,9]
)

# get names of samples
samples <- pca.results.df$sample
samples <- gsub('NAV', '0', samples, perl=T)
samples <- gsub('NA', '', samples, perl=T)
pca.results.df$sample <- samples
pca.results.df %>% arrange(sample) -> pca.results.df
pca.results.df %>% left_join(geo.coords.df, by='sample') %>%
  dplyr::select(sample, lon, lat, elev,
                pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9) -> pca.results.df

# pca.results.df <- add_column(pca.results.df, bc@presence)

pca.plt.1 <- ggplot(data=pca.results.df,
       aes(x=pc1,
           y=pc2,
           label=sample,
           color='red',
           lon=lon,
           lat=lat,
           xmin=0.5,
           xmax=-0.5,
           ymin=0.5,
           ymax=-0.5,
           )) +
  # geom_point(size=5, alpha=0.95) +
  geom_text(size=5, alpha=0.95, color='black') +
  theme(plot.subtitle=element_text(vjust=1),
        plot.caption=element_text(vjust=1)) +
  labs(x=paste('PC1 (', round(pca.results$varprop[1] * 100, 2), '%)', sep=''),
  y=paste('PC2 (', round(pca.results$varprop[2] * 100, 2), '%)', sep='')) +
  labs(colour=NULL) +
  theme_bw() + coord_fixed(ratio=1)
# pca.plt.1
ggsave('26-pca-a.pdf', plot=pca.plt.1, width=8, height=8, units='in')

pca.results.pcaviz <- pcaviz(dat=pca.results.df)
# plot(pca.results.pcaviz, show.legend=F)

###############################################################################
# i <- pca.results.df$pc1

# # v <- c(colnames(bc@presence), 'lon', 'lat')

# for (x in v) {
#   # message('#####################################################################')
#   # print(x)

#   # j <- pca.results.df[, x][[1]]
#   j <- scale(pca.results.df[, x][[1]])[, 1] + scale(pca.results.df$lon)[, 1]

#   fit <- lm(j ~ i)
#   # print(summary(fit))

#   rsqr <- paste('R^{2}=', round(summary(fit)$r.squared, digits=2), sep='')
#   message(x, ' :: ', rsqr, ' :: ', cor(i, j))

#   plot(
#     i,
#     j,
#     xlab='i',
#     ylab='j',
#     # xlim=c(0, 14),
#     # ylim=c(0, 14),
#     main=bquote(.(x) : R^2 == .(round(summary(fit)$r.squared, digits=2))),
#     pch=19)

#   abline(fit)
# }
###############################################################################
# i <- pca.results.df$pc1
# j <- scale(pca.results.df$bio17)[,1] + scale(pca.results.df$lat)[,1]
# j <- scale(pca.results.df$lon, center=T)[,1] +
#   scale(pca.results.df$lat, center=T)[,1] +
#   scale(pca.results.df$elev, center=min(pca.results.df$elev))[,1]
# cor(i, j)
# fit <- lm(j ~ i)
# print(summary(fit))
###############################################################################

# plot(pca.results.pcaviz, draw.points=F, label='sample', color='elev')

plot1 <- plot(pca.results.pcaviz, draw.points=F, label='sample', color='pc1')
plot2 <- plot(pca.results.pcaviz, draw.points=F, label='sample', color='pc3')
plot3 <- plot(pca.results.pcaviz, draw.points=F, label='sample', color='lat')
plot4 <- plot(pca.results.pcaviz, draw.points=F, label='sample', color='elev')
plot5 <- plot(pca.results.pcaviz, draw.points=F, label='sample', color='pc4')
plot6 <- plot(pca.results.pcaviz, draw.points=F, label='sample', color='pc5')
plts <- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol=2, nrow=3,
                  labels='AUTO')

ggsave('26-pca-b.pdf', plot=plts, width=8, height=11.5, units='in')

# pairs(as.data.frame(pca.results.df[c(2:10)]), pch=19, cex=0.5, col='red')

# # rotate
# pca.results.pcaviz %>%
#   pcaviz_rotate(-30) %>%
#   plot(show.legend=F)
#
# plot(pca.results.df$lat, pca.results.df$pc1, pch=19)
# plot(pca.results.df$elev, pca.results.df$pc8, pch=19)
# plot(pca.results.df$lat, pca.results.df$pc2, pch=19)
# plot(pca.results.df$lon, pca.results.df$pc2, pch=19)
#
# ggplotly(plot(pca.results.pcaviz %>% pcaviz_rotate(-90), draw.points=T, color='lat'))

# -----------------------------------------------------------------------------
setwd(path.scpt)
# -----------------------------------------------------------------------------
