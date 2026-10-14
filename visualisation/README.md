# 3D visualisation

An interactive 3D map of the paper's data: Canadian CMAs rendered as extruded
columns, column height mapped to Average AQHI and column color mapped to
Gini coefficient, with AQHI monitoring station locations shown as small
spikes.

Built with [deck.gl](https://deck.gl) (`GeoJsonLayer`, `ColumnLayer`) and
[MapLibre GL JS](https://maplibre.org) for the base map. No API key needed,
no build step, just static files.

## Running it

The data files are already generated in `data/`, so you can just serve this
folder:

```bash
cd visualisation
python3 -m http.server 8000
```

Then open `http://localhost:8000` in a browser. It must be served over HTTP
(not opened as a `file://` URL) because the browser blocks local `fetch()`
of the GeoJSON files otherwise.

Drag to rotate/tilt, scroll to zoom, hover a column for its stats. The
elevation slider rescales column height; the checkbox toggles the station
markers.

## Regenerating the data

`data/cma.geojson` and `data/stations.geojson` are derived from the same
pipeline as the paper (`Code/DataCleaner2.R`): AQHI readings joined to CMA
boundaries, income data joined in, Gini coefficient estimated via the
monotone-spline Lorenz curve method from Appendix A. To rebuild them from
the raw data in `Data/`:

```bash
Rscript visualisation/export_data.R
```

Requires R with the `sf`, `data.table`, `jsonlite`, `dplyr`, `tidyr`, and
`stringr` packages installed. Run from the repo root.
