"""Download and verify the EnMAP-Box tutorial inputs; standard library only.

Source: EnMAP-Box regression-based urban unmixing tutorial.
Data: Okujeni, van der Linden & Hostert (2009 Berlin Urban Gradient),
https://doi.org/10.5880/enmap.2016.002, CC BY-SA 4.0.
The code's MIT license does not replace the data license.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import tempfile
import urllib.request
import zipfile
from pathlib import Path


SOURCE_URL = 'https://box.hu-berlin.de/f/e4fa78c198bc4d868d30/?dl=1'
ARCHIVE_SHA256 = '22c4f40d89cc329d0fe6fa500b8ced50239c7e318658eee52a08bce6dc4b6c55'
INPUT_SHA256 = {
    'enmap_berlin.tif': 'c4dfc8c0a4749b44923b93534be17f1bf59cbb989e2525452be229171282d113',
    'library_berlin.gpkg': 'f9a0eb1543e3570157021253cb88a8dd78effb49eadd4716c57f98407d497eb6',
    'landcover_berlin.gpkg': '07cf337e149c3e85da41aabbccb33c22d4295ab3972183f6bcf8ba297f68cad0',
}


def sha256(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def check(path: Path, expected: str) -> None:
    if sha256(path) != expected:
        raise ValueError(f'SHA256 mismatch: {path}. Expected the tutorial version 20240122; no inputs were replaced.')


def download(archive: Path, url: str = SOURCE_URL) -> None:
    archive.parent.mkdir(parents=True, exist_ok=True)
    # A unique temporary file avoids overwriting a user's existing files.
    with tempfile.NamedTemporaryFile(dir=archive.parent, suffix='.part', delete=False) as stream:
        temporary = Path(stream.name)
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'enmap-unmixing-reproducibility/1.0'})
            with urllib.request.urlopen(request, timeout=120) as response:
                shutil.copyfileobj(response, stream, length=1024 * 1024)
        except BaseException:
            stream.close()
            temporary.unlink(missing_ok=True)
            raise
    try:
        check(temporary, ARCHIVE_SHA256)
        temporary.replace(archive)
    finally:
        temporary.unlink(missing_ok=True)


def prepare(archive: Path, destination: Path) -> None:
    """Verify the pinned archive and extract only the three required inputs."""
    check(archive, ARCHIVE_SHA256)
    destination.mkdir(parents=True, exist_ok=True)
    # Validate every member before replacing any destination. Never extract paths
    # supplied by a ZIP: only fixed allowlisted basenames are written.
    with tempfile.TemporaryDirectory(dir=destination, prefix='.tutorial-') as staging:
        staging = Path(staging)
        with zipfile.ZipFile(archive) as source:
            for name, expected in INPUT_SHA256.items():
                members = [item for item in source.infolist()
                           if not item.is_dir() and Path(item.filename).name == name]
                if len(members) != 1:
                    raise ValueError(f'Expected exactly one archive member for {name}')
                with source.open(members[0]) as member, (staging / name).open('wb') as target:
                    shutil.copyfileobj(member, target)
                check(staging / name, expected)
        for name, expected in INPUT_SHA256.items():
            target = destination / name
            if target.exists():
                check(target, expected)
        for name in INPUT_SHA256:
            target = destination / name
            if not target.exists():
                (staging / name).replace(target)
            print(f'Verified {name}')


def main() -> None:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=root / 'data/tutorial_materials.zip',
                        help='Existing tutorial ZIP, or download destination if absent')
    parser.add_argument('--destination', type=Path, default=root / 'data',
                        help='Directory for the three verified tutorial inputs')
    args = parser.parse_args()
    if not args.archive.exists():
        print(f'Downloading ~274 MB from the EnMAP-Box tutorial: {SOURCE_URL}')
        download(args.archive)
    prepare(args.archive, args.destination)
    print('Ready. Data attribution/license: https://doi.org/10.5880/enmap.2016.002 (CC BY-SA 4.0).')


if __name__ == '__main__':
    main()
