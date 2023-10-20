#!/usr/bin/env python3

from os.path import join as opj
from os.path import expanduser as ope

from pprint import pprint

# base
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
# from pandas_plink import read_plink

# viz
# import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz
from feems.cross_validation import run_cv

# matplotlib.use('pdf')

# change matplotlib fonts
# plt.rcParams['font.family'] = 'DejaVu Sans'
# plt.rcParams['font.sans-serif'] = 'DejaVu Sans'

data_dir = ope('~/SyncThing/nrv-dr-pitcairnia/KR/pop-gen-01-scripts/data')

# plink_file_path = opj(data_dir, 'nrv-pitcarn.vcf.filt.01')

gt_file_path = ope(opj(
    '~',
    'SyncThing',
    'nrv-dr-pitcairnia',
    'KR',
    'pop-gen-02-results',
    '29-resampled-snps',
    'theta-0.001',
    'min-samp-15',
    '001_theta-0.001_min-samp-15_feems.csv'))

outline_coords_path = opj(data_dir, 'dominican_republic_outline_coords_2')
sample_coords_path = opj(data_dir, 'sample_coords')

dgg_res = 10

pprint(dgg_res)

# path to discrete global grid
grid_path = opj(data_dir, 'dggs', f'dgg_hisp_triangle_res_{dgg_res:02d}.shp')

# NAV9 NAV6 NAV4 NA10 NA15 NAV5 NAV8 NA17 NAV2 NAV1 NA19 NAV7 NA11 NAV3 NA12 NA18 NA14 NA16 NA13
# 1/1  1/1  1/1  0/1  0/1  0/0  1/1  0/1  1/1  1/1  1/1  1/1  0/1  0/1  0/1  0/1  0/1  0/0  1/1
#  0    0    0    1    1    2    0    1    0    0    0    0    1    1    1    1    1    2    0

# read the genotype data
# (bim, fam, bed) = read_plink(plink_file_path, verbose=False)
# genotypes = np.array(bed)
# genotypes = [tuple(row) for row in genotypes if len(np.unique(row)) > 1]
# genotypes = np.array(genotypes).T
genotypes = np.genfromtxt(gt_file_path, delimiter=',')
# genotypes = np.delete(genotypes, 14, 0)
print(f'n_samples={genotypes.shape[0]}, n_snps={genotypes.shape[1]}')

imp = SimpleImputer(missing_values=np.nan, strategy='mean')
imp = imp.fit(genotypes)
genotypes = imp.transform(genotypes)
genotypes = genotypes.T
genotypes = [tuple(row) for row in genotypes if len(np.unique(row)) > 1]
genotypes = np.array(genotypes).T
print(f'n_samples={genotypes.shape[0]}, n_snps={genotypes.shape[1]}')

sample_coords = np.loadtxt(sample_coords_path, delimiter=',')
# sample_coords = np.delete(sample_coords, 14, 0)

outline_coords = np.loadtxt(outline_coords_path, delimiter=',')

outline_coords, edges, grid, _ = prepare_graph_inputs(
    coord=sample_coords,
    ggrid=grid_path,
    translated=True,
    buffer=0,
    outer=outline_coords)

sp_graph = SpatialGraph(genotypes, sample_coords, grid, edges, scale_snps=True)


def plot_feems_results(sp_graph_local, projection_local, pdf_file_name_local):
    fig_local = plt.figure(dpi=1200)
    ax_local = fig_local.add_subplot(1, 1, 1, projection=projection_local)
    coastline_m = '10m'
    v = Viz(ax_local, sp_graph_local,
            projection=projection_local,
            coastline_m=coastline_m,

            edge_zorder=2,
            obs_node_zorder=3,
            sample_pt_zorder=4,

            edge_color='green',
            edge_alpha=1,
            edge_width=0.5,

            obs_node_color='blue',
            obs_node_alpha=0.25,
            obs_node_linewidth=0.25,
            obs_node_size=11,

            sample_pt_color='red',
            sample_pt_alpha=1,
            sample_pt_linewidth=0.25,
            sample_pt_size=8,

            cbar_font_size=4,
            cbar_ticklabelsize=4,
            cbar_width='20%',
            cbar_height='5%',
            cbar_loc='upper left')

    def draw_map(ax_local_local):
        ax_local_local.add_feature(cfeature.LAND, facecolor="#ffffff", zorder=0)
        ax_local_local.coastlines(
            coastline_m,
            color="#000000",
            linewidth=0.5,
            zorder=0
        )

    draw_map(ax_local)
    v.draw_samples()
    v.draw_edges(use_weights=True)
    v.draw_obs_nodes(use_ids=False)
    v.draw_edge_colorbar()

    # fig.show()
    fig_local.savefig(f'{pdf_file_name_local}.pdf', format='pdf')


# ----------------------------------------------------------------------------
projection = ccrs.EquidistantConic(central_longitude=-75, central_latitude=17)
# ----------------------------------------------------------------------------

lambda_vals = (00.01, 00.05, 00.10, 00.25, 00.50, 00.75, 01.00, 10.00)

for lambda_val in lambda_vals:
    sp_graph.fit(lamb=lambda_val, verbose=True)
    plot_feems_results(
        sp_graph, projection,
        f'feems_dgg_hisp_triangle_res_{dgg_res:02d}_lambda_{lambda_val:05.2f}')

# ----------------------------------------------------------------------------

# define grids
# reverse the order of lambdas and alphas for warmstart
lamb_grid = np.geomspace(1e-6, 1e2, 30)[::-1]
pprint(lamb_grid)

# run cross-validation
cv_err = run_cv(sp_graph, lamb_grid, n_folds=sp_graph.n_observed_nodes, factr=1e10)
# pprint(cv_err)

# average over folds
mean_cv_err = np.mean(cv_err, axis=0)
pprint(mean_cv_err)

# argmin of cv error
lamb_cv = float(lamb_grid[np.argmin(mean_cv_err)])
pprint(lamb_cv)

# ----------------------------------------------------------------------------

fig, ax = plt.subplots(dpi=600)
ax.plot(np.log10(lamb_grid), mean_cv_err, ".")
ax.set_xlabel("log10(lambda)")
ax.set_ylabel("L2 CV Error")
ax.axvline(np.log10(lamb_cv), color="orange")

pdf_file_name = f'feems_dgg_hisp_triangle_res_{dgg_res:02d}_cv'
fig.savefig(f'{pdf_file_name}.pdf', format='pdf')

# ----------------------------------------------------------------------------

lambda_val = lamb_cv
sp_graph.fit(lamb=lambda_val, verbose=True)
plot_feems_results(
    sp_graph, projection,
    f'feems_dgg_hisp_triangle_res_{dgg_res:02d}_lambda_{lambda_val:05.2f}'
)

# ----------------------------------------------------------------------------

feems_nodes = sp_graph.nodes
pd.DataFrame(feems_nodes).to_csv('feems_nodes.csv', header=False, index=False)

feems_node_pos = sp_graph.node_pos
pd.DataFrame(feems_node_pos).to_csv('feems_node_pos.csv', header=False, index=False)

feems_edges = sp_graph.edges
pd.DataFrame(feems_edges).to_csv('feems_edges.csv', header=False, index=False)

feems_w = sp_graph.w
pd.DataFrame(feems_w).to_csv('feems_w.csv', header=False, index=False)
