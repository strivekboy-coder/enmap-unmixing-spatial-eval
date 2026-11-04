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
