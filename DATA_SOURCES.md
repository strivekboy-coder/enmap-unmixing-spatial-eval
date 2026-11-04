# Data attribution and licenses

Project code is licensed under [MIT](LICENSE). That permission does not replace
the licenses of external data. Original tutorial data, downloaded assets, GeoTIFFs
and processing caches are not redistributed in this repository.

## EnMAP-Box tutorial and Berlin Urban Gradient

The inputs come from the [EnMAP-Box regression-based urban unmixing tutorial](https://enmap-box.readthedocs.io/en/latest/usr_section/application_tutorials/urban_unmixing/tutorial.html).
Run `python download_data.py` to obtain the [20240122 tutorial archive](https://box.hu-berlin.de/f/e4fa78c198bc4d868d30/?dl=1).
The script pins the archive and three required files by SHA256 and refuses changed
content. The archive checksum is
`22c4f40d89cc329d0fe6fa500b8ced50239c7e318658eee52a08bce6dc4b6c55`.

Please cite: Okujeni, A., van der Linden, S., and Hostert, P. (2016),
*Berlin-Urban-Gradient dataset 2009 — An EnMAP Preparatory Flight Campaign*,
[GFZ Data Services, doi:10.5880/enmap.2016.002](https://doi.org/10.5880/enmap.2016.002),
and the [technical report, doi:10.2312/enmap.2016.002](https://doi.org/10.2312/enmap.2016.002).
Data license: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).
Published result tables and figures derived from these inputs are shared under
CC BY-SA 4.0 with this attribution. They are new analytical outputs, not copies
of the original tutorial rasters or spectral library.

## Sentinel-2 and WorldCover

Stage B uses Copernicus Sentinel-2 L2A observations (2017/2018 selection audit,
2018 analysis) through the public Earth Search STAC catalogue. Contains modified
Copernicus Sentinel data (2017–2018). The spectral response workbook comes from
[ESA SentiWiki](https://sentiwiki.copernicus.eu/web/s2-documents).

Independent land-cover labels use ESA WorldCover 2021 v200, produced by the ESA
WorldCover consortium. WorldCover is distributed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); see the
[WorldCover data access and attribution page](https://esa-worldcover.org/en/data-access).
STAC assets and source URLs, calibration overrides, input checksums and product
years are recorded in the small provenance JSON/CSV files under `results/`.

The 2009 tutorial scene, 2018 Sentinel-2 observations and 2021 WorldCover map
represent different years. Derived Stage B comparisons demonstrate methods;
they do not establish same-date sensor superiority or verified temporal stability.


## Real-data first look (2025)

EnMAP L2A scene `ENMAP01-____L2A-DT0000129412_20250511T104201Z_009_V010502_20250515T040235Z`
was acquired on 11 May 2025 and retrieved from the [DLR Geoservice STAC](https://geoservice.dlr.de/eoc/ogc/stac/v1/collections/ENMAP_HSI_L2A).
The original EnMAP product remains local; DLR/EnMAP data terms apply independently
of the MIT code license. Actual Sentinel-2C L2A `S2C_33UUU_20250512_0_L2A`
comes from Earth Search / Copernicus, with original ESA XML/JP2 calibration checks.
Contains modified Copernicus Sentinel data (2025). Convolution uses the official
ESA Sentinel-2C spectral response functions, version 4.0. The first-look figures
also derive from the attributed 2009 tutorial inputs above and retain CC BY-SA 4.0.

To reproduce, first run Stage A, then the ordered stages in
[REAL_DATA_FIRST_LOOK.md](REAL_DATA_FIRST_LOOK.md). Install the optional isolated
AROSICS environment with `micromamba create -p ../.ar_env -f environment_arosics.yaml`
([official installation guide](https://danschef.git-pages.gfz-potsdam.de/arosics/doc/installation.html));
configure `micromamba_executable` and `arosics_environment` in
`real_enmap_config.yaml`. The public S2 date/extent search is fetched on demand.
Only EnMAP download requires local `ENMAP_USERNAME` and `ENMAP_PASSWORD` user
environment variables (IPS main identity); no credentials belong in this repository.
Detailed input checksums, quality counts, calibration audits and band metrics are
regenerated under ignored `data/real_enmap/20250511/`, not committed.
