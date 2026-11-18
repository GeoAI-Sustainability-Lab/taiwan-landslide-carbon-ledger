# Living-biomass carbon lost to and re-accumulated after landslides in Taiwan, 2004–2025

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24350620.svg)](https://doi.org/10.5281/zenodo.24350620)

Code, derived tables and figures for a national annual ledger of forest living-biomass carbon removed by catalogued landslides
in Taiwan and re-accumulated on the scars, 2004–2025, and for its comparison with the disturbance-loss term of the national
greenhouse gas inventory. The ledger follows the gain–loss method of the IPCC guidelines on a 30 m grid: the event-based
landslide catalogue gives the activity data; the stand stock of the fourth national forest inventory, scaled by the measured
pre-failure canopy stature of the failing cells, gives the stock removed; and canopy-height recovery measured on dated
canopy-height imagery gives the re-accumulation. Three independent data sets check the three key steps.

## Main results

All values in Mt CO<sub>2</sub>; intervals are 95% block-bootstrap intervals (300 replicates). Source: `results/ledger_summary.json`, `results/parameters.json`.

| Quantity | Value |
|---|---|
| Living-biomass loss 2004–2025 | 14.72 (12.88–16.80) |
| Loss in 2009 (Typhoon Morakot year) | 8.27 (7.22–9.39), 56.2% of the total |
| Share of the four largest loss years (2009, 2007, 2024, 2016) | 80.2% |
| Re-accumulation on scars to 2025 | 3.97 (2.68–5.43), 27.0% of the loss |
| Unrecovered stock at the end of 2025 | 10.75 (8.74–12.60) |
| Pre-failure relative biomass of closed-canopy failing cells | 0.55–0.81 of the local reference forest |
| Recovery rate of relative canopy height (% yr<sup>−1</sup>): broadleaf, conifer, mixed, bamboo | 0.78, 0.92, 1.15, 4.06 |
| Disturbance loss in the national inventory 2004–2024 | 3.04; ledger loss over the same years 14.35 |
| Twelve alternative settings | loss 11.7–22.8; re-accumulation 0.7–5.2 |

## Repository layout

| Path | Content |
|---|---|
| `src/` | Analysis pipeline, numbered in run order (table below); `src/config.py` holds the input and working locations |
| `src/prep/` | Preparation of inputs: rasterisation of the fourth-inventory forest-type and stand-volume maps, GEDI L4A download through SlideRule, download and compilation of the third-inventory photo points |
| `src/cases/` | Image chips (Planetary Computer STAC), canopy-height windows and the NDVI series of the two case sites |
| `results/` | Derived tables and statistics written by the pipeline (CSV, JSON, bootstrap replicates `bootstrap.npz`) |
| `figures/` | Figures 1–8 and S1–S4 (PNG, 300 dpi, and PDF); Figure 2 is a schematic, all other panels are computed from `results/` |

## Running the pipeline

```bash
pip install -r requirements.txt
export LRC_DATASET=/path/to/input    # input data store, layout below
export LRC_WORK=/path/to/work        # working directory; s02 onwards run with it as the current directory
cd "$LRC_WORK"
python src/q01_catalogue_cells.py "$LRC_WORK/q01"
python src/s01_chm_date_raster.py "$LRC_WORK/s01/chmv2_acqdate_30m.npz"
python src/s02_build_cells.py v1
python src/s05_extend_cells.py
for y in 2017 2018 2019 2020 2021 2022 2023 2024; do python src/s07_annual_state_raster.py $y <annual layer of year y (.shp or .zip)> annual_state/annual_state_$y.npz; done
python src/s07b_join_annual_state.py
python src/s06_gedi_join.py
NBOOT=300 python src/s08_ledger.py
python src/s09_tables.py
python src/s14_nfi3_check.py
python src/cases/fetch_chips.py
python src/cases/ts.py 2005 2012; python src/cases/ts.py 2013 2019; python src/cases/ts.py 2020 2025
python src/s22_case_sites.py
python src/s21_figures.py
```

(Write `src/` as the path to this repository's `src` folder.) The random seeds are fixed, so with the inputs listed below the
scripts reproduce the files in `results/`. The cell-level intermediate tables (several hundred MB) are not included; the
pipeline regenerates them in `v1/` and `v2/` of the working directory.

| Script | Step |
|---|---|
| `q01_catalogue_cells.py` | Catalogue polygons on a 6 m sub-grid, aggregated to the fraction of each 30 m cell per catalogue record |
| `s01_chm_date_raster.py` | Image date of the canopy-height map on the 30 m grid |
| `s02_build_cells.py` | Cell table: landslide fraction, failure dates, forest type, stand volume, tree cover 2000, terrain, canopy height and image date, L-band backscatter; 500 000 never-catalogued reference forest cells |
| `s05_extend_cells.py` | Stand origin, Global Forest Change loss year, polygon size, interior index |
| `s07_annual_state_raster.py`, `s07b_join_annual_state.py` | Official annual landslide maps 2017–2024 on the 30 m grid and in the cell table |
| `s06_gedi_join.py` | GEDI L4A footprints with cell attributes, canopy height and landslide history |
| `s08_ledger.py` | Relative canopy height, pre-failure relative biomass, recovery rates, annual ledger 2004–2025, block bootstrap with conversion-factor uncertainty, twelve alternative settings, projection after 2025 |
| `s09_tables.py` | Supporting tables and statistics (recovery bins, scar-size gradient, annual-map persistence, GEDI checks) |
| `s14_nfi3_check.py` | Independent check of pre-failure stocking with the third-inventory photo points |
| `s22_case_sites.py`, `s21_figures.py` | Figures and case-site statistics |

## Input data

The pipeline reads openly available data. Paths are relative to `LRC_DATASET`.

| Data | Expected location | Source and licence |
|---|---|---|
| Event-based landslide catalogue 2004–2025 (polygons, event, pre- and post-event image dates) | `02_CURATED/Hazard_Events/Landslide_Taiwan/Phase19_Event_Catalog_2004_2025_20260825/geodata/phase19_event_catalog_candidates.gpkg` | Agency of Rural Development and Soil and Water Conservation (ARDSWC), https://tech.ardswc.gov.tw/Results/EventBasedLandslideInventory |
| Forest-type map and stand-volume (carbon storage) map of the fourth national forest inventory, rasterised with `src/prep/` | `02_CURATED/Forestry_Reference/Taiwan_forest_type_map_30m_TWD97/forest_type_30m.npz`, `02_CURATED/Forestry_Reference/Taiwan_4th_NFI_volume_30m_TWD97/nfi4_volume_30m.npz` | Forestry and Nature Conservation Agency; https://data.gov.tw/dataset/57873 and https://data.moa.gov.tw/open_detail.aspx?id=J18; Open Government Data License, Taiwan, v1.0 |
| CHMv2 10 m mean canopy height (tiles 1321230–1321233) and image-date polygons | `01_SOURCE/Terrain_Canopy/CHMv2_10m_mean*/mean/*.tif`, `01_SOURCE/Terrain_Canopy/CHMv2_Meta_WRI_2026/metadata/tile_observation_dates/` | Brandt et al. 2026, https://doi.org/10.48550/arXiv.2603.06382; AWS Open Data `dataforgood-fb-data`; CC BY 4.0 |
| Global Forest Change v1.13 (tree cover 2000, loss year) and SRTM elevation and slope on the 30 m grid | `01_SOURCE/Earth_Observation/Lab_platform_assets/hansen_gfc_2025_30m.tif`, `dem_terrain_30m.tif` | Hansen et al. 2013, https://doi.org/10.1126/science.1244693 (CC BY 4.0); SRTM, NASA/USGS (public domain) |
| GEDI L4A footprint aboveground biomass density v2.1 | `01_SOURCE/Earth_Observation/GEDI_L4A_V2_1_SlideRule_Taiwan_20261111/derived/gedi_l4a_taiwan_shots.npz` | Dubayah et al. 2022, https://doi.org/10.3334/ORNLDAAC/2056, through SlideRule Earth (`src/prep/acquire_gedi_l4a_sliderule.py`) |
| Taiwan landslide inventory, annual layers 2017–2024 | input of `s07_annual_state_raster.py` | ARDSWC, https://tech.ardswc.gov.tw/Results/TaiwanLandslideInventory |
| Photo points of the third national forest inventory, 1990–1994 | `01_SOURCE/Forestry_Reference/Taiwan_NFI3_SamplePoints_FANC_MapServer_20261111/compiled/nfi3_photo_points_layer1_core.csv` | Forestry and Nature Conservation Agency map service MA_Map0210, layer 1 (`src/prep/download_nfi3_photo_points.py`, `src/prep/compile_nfi3_photo_points.py`) |
| ALOS-2 PALSAR-2 2020 global mosaic, HV | `01_SOURCE/Earth_Observation/ALOS2_PALSAR2_Mosaic_CentralTaiwan_2018_2020_20261107/2020/*_HV_*.tif` | JAXA, Shimada et al. 2014, https://doi.org/10.1016/j.rse.2014.04.014 |
| Natural Earth admin-0 countries, 1:10 m | `01_SOURCE/Reference_Boundaries/NaturalEarth_admin_10m/ne_10m_admin_0_countries.shp` | Natural Earth, public domain |
| Landsat Collection 2 Level-2 and Sentinel-2 Level-2A (case sites) | read on demand | Microsoft Planetary Computer, https://planetarycomputer.microsoft.com |
| National greenhouse gas inventory 2026, chapter 6 (forest-land net removals, disturbance losses, conversion factors and their uncertainty) | values written in `s08_ledger.py` | Ministry of Environment, Taiwan, https://www.cca.gov.tw/information-service/publications/national-ghg-inventory-report/34874.html |

Grid: TWD97 / TM2 zone 121 (EPSG:3826), 30 m, upper-left corner (151050, 2798910), 6680 columns × 12558 rows; a cell index is
`row × 6680 + column`.

## Citation

Hsu, Y.-C. (2026). taiwan-landslide-carbon-ledger: national annual ledger of forest living-biomass carbon lost to and
re-accumulated after landslides in Taiwan, 2004–2025 (code, derived tables and figures), version 1.0.0. Zenodo.
https://doi.org/10.5281/zenodo.24350620 (all versions: https://doi.org/10.5281/zenodo.24350619). See also `CITATION.cff`.

## Licence

Code: MIT (`LICENSE`). Derived tables and figures: see `DATA_LICENSE.md` for the attributions that the source licences require.
