"""Download the two small, public primary-source assets used by the study."""
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).parent
URLS = {
    'data/sentinel2_srf.xlsx': 'https://sentiwiki.copernicus.eu/__attachments/a_ece6183b6698587e7ecd9804974bece3d82e39e151aa5aeca9d8fd95e1f1558c/COPE-GSEG-EOPG-TN-15-0007%20-%20Sentinel-2%20Spectral%20Response%20Functions%202024%20-%204.0.xlsx?cb=05a866b9653a4d2bc442df394ea6cc7a',
    'diagnostics/tutorial_accuracy.png': 'https://enmap-box.readthedocs.io/en/latest/_images/20_accuracies.png',
}
for filename, url in URLS.items():
    path = ROOT / filename
    if path.exists():
        continue
    request = Request(url, headers={'User-Agent': 'EnMAP educational reproducibility study'})
    with urlopen(request, timeout=90) as response:
        content = response.read()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    print(filename, len(content), flush=True)
