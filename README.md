# Pacific trench segmentation

Code and data for the paper:

> P. Lemenkova, A. L. Piskarev. *Along-Strike Segmentation of Pacific Trenches from Axial Bathymetry.*

Detects along-strike segment boundaries on 20 Pacific trench traces from open
bathymetry and tests their coincidence with subducting fracture zones, plateaus,
ridges and seamounts against two null models.

## Layout

```
./          analysis and figure scripts
data/       trench axes, per-node series, fabric contacts
results/    boundaries, trench summary, coincidence and sensitivity tests
figures/    figure output (workflow.tex kept)
gmt/        GMT gravity map and trench traces
```

### Scripts

| file | stage | purpose |
|------|-------|---------|
| `extract_series.py` | 1a | cross-axis profiles, Viterbi thalweg pick, per-node depth, width and wall rises → `series_XX.csv` |
| `fabric_contacts.py` | 1b | fracture zones, LIP plateaus/ridges and seamounts projected onto each axis → `contacts.csv` |
| `analysis.py` | 2 | trench-presence trimming, PELT segmentation, 54-run sensitivity grid, coincidence tests → `results/` |
| `figures.py` | 3 | all data figures → `figures/fig_*.pdf/.png` |
| `gmt/gravity_fabric.sh` | 3 | free-air gravity map with the trench traces |
| `figures/workflow.tex` | – | TikZ workflow diagram |

### data/

| file | description |
|------|-------------|
| `trench_axes_full.gmt` | 20 PB2002 trench traces |
| `series/series_XX.csv` | per-node series: `s_km, u_pick_km, d_m, rise_land_m, rise_sea_m, W_km, relief_outer_m` |
| `fabric_contacts.csv` | fabric contacts per trench |
| `fabric_geometry.gmt` | fabric geometry for maps |

### results/

| file | description |
|------|-------------|
| `trench_summary.csv` | analysed window, depth range, segments per trench |
| `boundary_catalogue.csv` | 479 boundaries: position, stability, step sizes |
| `contacts_in_windows.csv` | contacts inside the analysed windows |
| `coincidence_by_class.csv` | distance and tolerance tests, uniform and circular-shift nulls |
| `coincidence_by_trench.csv` | per-trench tests |
| `coincidence_loto.csv` | leave-one-trench-out tests |
| `null_distributions.json` | null distributions of the mean distance |
| `sensitivity_grid.csv` | boundary counts over the parameter grid |

## Requirements

Python 3.11 with `numpy`, `ruptures`, `shapely`, `matplotlib`, `basemap`,
`basemap-data`, `pillow`; Nimbus Sans font for the figures. GMT 6 for stage 1a
and the gravity map; GDAL (`ogr2ogr`) for stage 1b; `pdflatex` + TikZ for
`workflow.tex`.

```
pip install -r requirements.txt
```

## Running

Stages 2–3 run from `data/` alone:

```
python analysis.py
python figures.py
cd gmt && bash gravity_fabric.sh
```

Stage 1 rebuilds `data/` from the source datasets (GMT remote grid
`@earth_gebco_30s`; in the working folder `SHP/GSFML_SF_FZ_KM.shp`,
`IgneousProvinces/Whittaker_2015/SHP/Whittaker_etal_2015_LIPs.shp` and `kw.txt`
from the Kim & Wessel census):

```
python extract_series.py data/trench_axes_full.gmt 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
python fabric_contacts.py data/trench_axes_full.gmt
```

Move `series_XX.csv` to `data/series/` and `contacts.csv` to
`data/fabric_contacts.csv`.

## Data sources

GEBCO; PB2002; GSFML; Whittaker et al. (2015) LIPs; Kim & Wessel (2011) seamount
census; ETOPO1 and GSHHG (basemap-data); IGPP free-air gravity (GMT data server).
Fabric databases not redistributed here.

## License

MIT (see `LICENSE`).
