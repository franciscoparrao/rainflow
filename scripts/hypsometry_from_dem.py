#!/usr/bin/env python3
"""Derive a catchment hypsometric curve from a DEM clipped to its polygon.

Reproduces the `data/camels-cl/<id>_hypsometry.csv` curves consumed by
`rainflow split-sample --hypsometry-file`. Requires GDAL (ogr2ogr, gdalwarp,
gdalbuildvrt) and the Python GDAL bindings.

Inputs:
  - a catchment polygon (GeoJSON/shapefile) in EPSG:4326
  - one or more Copernicus DEM GLO-30 tiles covering it (1x1 deg GeoTIFFs;
    public, no auth: https://copernicus-dem-30m.s3.amazonaws.com/<TILE>/<TILE>.tif)

Output: CSV with 101 knots, "area_fraction,elevation_m".

Example:
  python hypsometry_from_dem.py --polygon 4703002.geojson \\
      --dem S32_W071.tif S33_W071.tif --out 4703002_hypsometry.csv
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from osgeo import gdal

gdal.UseExceptions()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--polygon", required=True, help="catchment polygon (EPSG:4326)")
    ap.add_argument("--dem", required=True, nargs="+", help="DEM tile(s) covering the catchment")
    ap.add_argument("--out", required=True, help="output hypsometry CSV")
    ap.add_argument("--knots", type=int, default=101, help="number of curve knots (default 101)")
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        vrt, clip = tmp / "mosaic.vrt", tmp / "clip.tif"
        # Repair the polygon (CAMELS-CL boundaries have self-intersections).
        poly = tmp / "valid.geojson"
        subprocess.run(["ogr2ogr", "-makevalid", str(poly), args.polygon], check=True)
        subprocess.run(["gdalbuildvrt", "-q", str(vrt), *args.dem], check=True)
        subprocess.run(
            ["gdalwarp", "-q", "-overwrite", "-cutline", str(poly),
             "-crop_to_cutline", "-dstnodata", "-9999", str(vrt), str(clip)],
            check=True,
        )

        band = gdal.Open(str(clip)).GetRasterBand(1)
        a = band.ReadAsArray().astype("float64")
        nd = band.GetNoDataValue()
        z = a[(a != nd) & np.isfinite(a) & (a > -1000)]
        z.sort()

        frac = np.linspace(0.0, 1.0, args.knots)
        elev = np.quantile(z, frac)
        with open(args.out, "w") as f:
            f.write("area_fraction,elevation_m\n")
            for fa, e in zip(frac, elev):
                f.write(f"{fa:.4f},{e:.1f}\n")

    print(f"{args.out}: {len(z)} pixels, {args.knots} knots "
          f"(min={z.min():.0f} median={np.median(z):.0f} max={z.max():.0f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
