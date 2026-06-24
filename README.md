# Habitat Connectivity Analysis for Pine Marten Conservation

Assessment of landscape connectivity between Cairngorms National Park and Loch Lomond National Park, Scotland, to support conservation planning for the pine marten (*Martes martes*).

Completed as part of the GEOG71922 Spatial Ecology module, MSc Geographical Information Science, University of Manchester.

---

## Overview

Pine marten populations in Scotland declined sharply during the 19th and 20th centuries due to deforestation and persecution, leaving isolated populations in the Scottish Highlands. This project evaluates the potential for functional connectivity between two major national parks using cost-surface modelling and graph-theoretic network analysis, identifying both the least-cost movement corridor and the relative importance of 330 conifer patches as stepping stones.

---

## Study Area

- **Cairngorms National Park** (4,528 km²): largest national park in the UK; Caledonian pine forests provide core pine marten habitat
- **Loch Lomond and The Trossachs National Park** (1,865 km²): southernmost extent of the Scottish Highlands; mixed deciduous and coniferous woodland

---

## Methods

### 1. Cost-Surface Resistance Modelling
- Corine Landcover 2018 dataset reclassified into a resistance layer representing movement difficulty for pine martens across 34 land cover types
- Resistance values sourced from Stringer et al. (2018), ranging from forests (low resistance) to urban areas and large water bodies (high resistance)

### 2. Least-Cost Corridor Delineation
- Cumulative cost surfaces generated from the centroids of each national park using the `accCost` function, which applies Dijkstra's algorithm via a transition matrix with 8-neighbour mean conductance values
- The two cumulative cost layers were overlaid; the lowest 7% cost threshold was used to delineate the corridor

### 3. Stepping-Stone Patch Analysis
- 330 conifer patches within the study extent evaluated for their relative importance as stepping stones
- **Betweenness Centrality**: calculated from an adjacency matrix based on a negative exponential function using the pine marten dispersal distance of 8 km (Stringer et al., 2018); identifies patches critical for shortest-path movement across the network
- **Probability of Connectivity (dPC Index)**: measures the reduction in network-wide connectivity when each patch is individually removed; highlights patches whose loss would most severely fragment the landscape

---

## Key Findings

- The area between the two parks offers moderate to low resistance to pine marten movement, largely due to the prevalence of forests, grasslands, and scrubland
- Betweenness Centrality identified patches closely aligned with the least-cost corridor as the most critical nodes for movement
- dPC analysis revealed a smaller subset of patches whose removal would disproportionately reduce overall connectivity, providing prioritisation targets for conservation effort

---

## Files

| File | Description |
|------|-------------|
| `spatial_ecology.R` | Full R script: resistance layer parameterisation, least-cost corridor generation, Betweenness Centrality and dPC analysis |
| `corridor_analysis_report.pdf` | Full project report with methodology, results, maps, and references |

---

## Dependencies

```r
library(raster)
library(gdistance)
library(igraph)
library(sf)
library(tmap)
```

---

## Data Sources

- Corine Landcover 2018 — European Environment Agency (publicly available at [land.copernicus.eu](https://land.copernicus.eu))
- Resistance values — Stringer et al. (2018), *The Feasibility of Reintroducing Pine Martens to the Forest of Dean and Lower Wye Valley*, Gloucestershire Wildlife Trust (values documented in the project report)
- National park boundaries — Scottish Natural Heritage / NatureScot (publicly available at [spatialdata.gov.scot](https://www.spatialdata.gov.scot))

> **Note:** Input data files are not included in this repository. The Corine Landcover and national park boundary datasets are publicly accessible at the links above. The full methodology and resistance parameterisation are documented in `corridor_analysis_report.pdf`.

---

## Limitations

- Resistance values are literature-derived rather than empirically calibrated; observed pine marten movement data were not available for this analysis
- The model does not account for ecological interactions such as predation or competition
- Future work should incorporate GPS-tracked movement data to validate and refine corridor delineation
