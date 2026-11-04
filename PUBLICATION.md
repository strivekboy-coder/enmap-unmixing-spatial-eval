# Public snapshot

This repository contains the completed Stage A, Stage A additions and Stage B
code, configurations, measured result tables, figures and small provenance records.
It is exported from local research commit `2352b459e8d6be4402f3a4cd4304d747be2ea3d3` on 2026-11-04.

Public history starts from a clean snapshot: older local commits contain raster
outputs and are not pushed. Original tutorial inputs, GeoTIFFs, caches, downloaded
assets and machine settings are excluded. No scientific result has been changed
for publication. Run the pipeline to regenerate complete raster outputs locally.

Five historical configuration snapshots substitute `data/tutorial_materials.zip`
for the original machine-specific archive path. Numerical settings and seeds are
unchanged. Input checksums remain in `results/input_manifest.json`.

Validation before publication: 12 scientific invariant tests passed; the downloader
passed pinned-archive extraction, exact input selection, repeated-use preservation
and corrupt-file rejection checks. See [results notes](RESULTS_NOTES.md).

Code is [MIT](LICENSE); external data and derived results retain the attribution
and licenses in [DATA_SOURCES.md](DATA_SOURCES.md).
