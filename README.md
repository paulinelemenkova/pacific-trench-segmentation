# Pacific trench segmentation

Code for the paper:

> P. Lemenkova. *Along-Strike Segmentation of Pacific Trenches from Axial Bathymetry.*

Detects along-strike segment boundaries on 20 Pacific trenches from open
bathymetry and tests their coincidence with subducting fabric.

## Layout

```
scripts/   analysis and figure code
data/      trench axes and derived outputs
```

### scripts/ — analysis

| file | purpose |
|------|---------|
| `results.py` | series and boundaries → `trench_summary.csv`, `boundary_catalogue.csv`, `coincidence_by_class.csv` |
| `make_fabric_intersections.py` | fabric–axis overlay → `fabric_intersections.csv` |

### scripts/ — figures

One `.py` per figure; a matching `.sh` is its GMT/PyGMT driver.

| file(s) | figure |
|---------|--------|
| `studyarea.py` / `.sh` | study-area map |
| `fabric_overview.py` / `.sh` | fabric overview |
| `axis_extraction.py` / `.sh` | axis extraction |
| `changepoint_model.py` / `.sh` | changepoint model |
| `alongstrike_panel.py` / `.sh` | along-strike profiles |
| `boundary_map.py` / `.sh` | boundaries by class |
| `sensitivity.py` / `.sh` | boundary sensitivity |
| `segment_contrast.py` | cross-boundary contrast |
| `coincidence.py` | coincidence, class fractions |
| `workflow.tex` | TikZ workflow diagram |

### data/

| file | description |
|------|-------------|
| `trench_axes_full.gmt` | trench axes (GMT, 5 km) |
| `trench_summary.csv` | per-trench summary |
| `boundary_catalogue.csv` | detected boundaries |

## Requirements

`numpy`, `pandas`, `matplotlib`, PyGMT + GMT 6; `ruptures` optional; `pdflatex`
+ TikZ for `workflow.tex`.

```
pip install -r requirements.txt
```

## Running

Axes read from `axes_full/<n>_<Name>.txt` (`lon lat`, 0–360), split from
`data/trench_axes_full.gmt`.

```
python scripts/results.py
python scripts/make_fabric_intersections.py
python scripts/results.py
bash scripts/studyarea.sh
python scripts/coincidence.py
python scripts/segment_contrast.py
```

## Data sources

GEBCO 2023; GSFML, EarthByte, seamount census, PB2002, MORVEL. Fabric databases
not redistributed here.

## License

MIT (see `LICENSE`).
