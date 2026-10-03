"""Sentinel-2 pilot, step 1: pick clear summer scenes and cut small chips around each station.

Source: Element84 Earth Search STAC catalog (public, no login); Sentinel-2 L2A (Sen2Cor) cloud-optimized
GeoTIFFs on public AWS storage. Only a small window around each station is read, so the download is a few MB.

For each pilot year the clearest July-September dates are chosen (all ten stations covered, lowest tile cloud
cover). For every scene and station a 41 x 41 pixel (410 m) chip is saved for the red, green, near-infrared and
shortwave-infrared bands and the scene classification, as surface reflectance (scale and offset taken from the
STAC metadata, which handles the processing-baseline change in January 2022).

Run with the Python environment that has rasterio, numpy, pandas and requests. Output: data/s2/pilot/
(gitignored via data/s2/).
"""
import json, os, sys, urllib.request
import numpy as np
import rasterio
from rasterio.windows import Window
from rasterio.warp import transform

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "s2", "pilot")
os.makedirs(OUT, exist_ok=True)

STAC = "https://earth-search.aws.element84.com/v1"
YEARS = [2018, 2019, 2020, 2025, 2026]
PER_YEAR = 2
HALF = 20                      # chip is (2*HALF+1) pixels square at 10 m
STATIONS = {                    # CDEC station coordinates (lat, lon)
    "SJR": (37.679290, -121.265110), "MSD": (37.786000, -121.306000), "RRI": (37.963000, -121.365000),
    "PPT": (38.066000, -121.562000), "TWI": (38.096900, -121.669100), "BLP": (38.032452, -121.719040),
    "ANH": (38.017830, -121.802963), "TRN": (37.992800, -121.454200), "OSJ": (38.071100, -121.578900),
    "MOK": (38.106100, -121.571100),
}


def stac_search(body):
    req = urllib.request.Request(STAC + "/search", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=90))


def in_polygon(lon, lat, ring):
    inside = False
    n = len(ring)
    for i in range(n):
        x1, y1 = ring[i]; x2, y2 = ring[(i + 1) % n]
        if (y1 > lat) != (y2 > lat) and lon < (x2 - x1) * (lat - y1) / (y2 - y1) + x1:
            inside = not inside
    return inside


def covers(item, lon, lat):
    g = item["geometry"]
    polys = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]
    return any(in_polygon(lon, lat, p[0]) for p in polys)


def scenes_for_year(year):
    r = stac_search({"collections": ["sentinel-2-l2a"], "bbox": [-121.9, 37.6, -121.2, 38.2],
                     "datetime": f"{year}-07-01T00:00:00Z/{year}-09-30T23:59:59Z", "limit": 300})
    by_date = {}
    for it in r["features"]:
        by_date.setdefault(it["properties"]["datetime"][:10], []).append(it)
    rows = []
    for d, items in by_date.items():
        cover = {}
        for s, (lat, lon) in STATIONS.items():
            c = [it for it in items if covers(it, lon, lat)]
            if c:
                cover[s] = min(c, key=lambda it: it["properties"].get("eo:cloud_cover", 100))
        if len(cover) < 8:
            continue
        cloud = max(it["properties"].get("eo:cloud_cover", 100) for it in cover.values())
        rows.append((cloud, len(cover), d, cover))
    rows.sort(key=lambda x: (x[0] - 2 * x[1] / 10.0, x[2]))
    return rows


def band_scale(item, asset):
    rb = item["assets"][asset].get("raster:bands")
    if rb and "scale" in rb[0]:
        return rb[0]["scale"], rb[0].get("offset", 0.0)
    base = float(item["properties"].get("s2:processing_baseline", "02.00"))
    return 0.0001, (-0.1 if base >= 4.0 else 0.0)


def read_chip(item, asset, lat, lon, half, resolution_factor=1):
    href = item["assets"][asset]["href"]
    with rasterio.open(href) as ds:
        xs, ys = transform("EPSG:4326", ds.crs, [lon], [lat])
        row, col = ds.index(xs[0], ys[0])
        h = max(1, half // resolution_factor)
        win = Window(col - h, row - h, 2 * h + 1, 2 * h + 1)
        data = ds.read(1, window=win, boundless=True, fill_value=0)
    return data


def main():
    manifest = []
    for year in YEARS:
        rows = scenes_for_year(year)
        print(f"{year}: {len(rows)} candidate dates covering >= 8 stations; clearest:",
              [(r[2], round(r[0], 1), r[1]) for r in rows[:4]])
        for cloud, n, d, cover in rows[:PER_YEAR]:
            for s, it in cover.items():
                lat, lon = STATIONS[s]
                chip = {}
                for asset, key in [("red", "red"), ("green", "green"), ("nir", "nir"), ("swir16", "swir"), ("scl", "scl")]:
                    f = 2 if asset in ("swir16", "scl") else 1
                    raw = read_chip(it, asset, lat, lon, HALF, f)
                    if asset == "scl":
                        scl = np.kron(raw, np.ones((2, 2), dtype=raw.dtype))[: 2 * HALF + 1, : 2 * HALF + 1]
                        chip[key] = scl.astype("uint8")
                    else:
                        sc, off = band_scale(it, asset)
                        arr = raw.astype("float32")
                        if f == 2:
                            arr = np.kron(arr, np.ones((2, 2), dtype="float32"))[: 2 * HALF + 1, : 2 * HALF + 1]
                        refl = np.where(arr > 0, arr * sc + off, np.nan)
                        chip[key] = refl.astype("float32")
                name = f"{d}_{s}_{it['properties'].get('grid:code', 'tile')}"
                np.savez_compressed(os.path.join(OUT, name + ".npz"), **chip)
                manifest.append({"date": d, "station": s, "item": it["id"], "datetime_utc": it["properties"]["datetime"],
                                 "tile_cloud": it["properties"].get("eo:cloud_cover"), "file": name + ".npz",
                                 "baseline": it["properties"].get("s2:processing_baseline")})
            print("  saved", d, f"({len(cover)} stations, max tile cloud {cloud:.1f}%)")
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    print("chips:", len(manifest), "; total MB:", round(sum(os.path.getsize(os.path.join(OUT, m['file'])) for m in manifest) / 1e6, 2))


if __name__ == "__main__":
    main()
