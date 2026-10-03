"""Download one Sentinel-2 Level-1C product (SAFE folder) from the public Google Cloud bucket
gcp-public-data-sentinel-2 (anonymous access, no login), into data/s2/l1c/<product>.SAFE/.

ACOLITE needs the full SAFE structure (band JP2 files plus the metadata XML files), so the whole product is
fetched (about 850 MB per tile and date). Already-downloaded files of the right size are skipped.

Usage: python sentinel2/download_l1c.py 10SFH 2019-07-24
   or: python sentinel2/download_l1c.py 10SFH S2B_MSIL1C_20190724T184929_N0500_R113_T10SFH_20230508T233807.SAFE
A date is resolved to the product name through the Earth Search STAC catalog (s2:product_uri of the Level-1C item).
"""
import os, sys, json, urllib.request, urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "s2", "l1c")
BUCKET = "gcp-public-data-sentinel-2"


def resolve_product(tile, date):
    point = {"10SFH": [-121.33, 37.95], "10SFG": [-121.265, 37.679]}.get(tile, [-121.33, 37.95])   # a point inside the tile
    body = {"collections": ["sentinel-2-l1c"], "intersects": {"type": "Point", "coordinates": point},
            "datetime": f"{date}T00:00:00Z/{date}T23:59:59Z", "limit": 20}
    req = urllib.request.Request("https://earth-search.aws.element84.com/v1/search", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    items = [f for f in json.load(urllib.request.urlopen(req, timeout=60))["features"] if f"_{tile}_" in f["id"]]
    if not items:
        raise SystemExit(f"no Level-1C item for tile {tile} on {date}")
    return items[0]["properties"]["s2:product_uri"]


def main(tile, product):
    if len(product) == 10 and product[4] == "-":
        product = resolve_product(tile, product)
    utm, band, sq = tile[:2].lstrip("0"), tile[2], tile[3:]
    prefix = f"tiles/{utm}/{band}/{sq}/{product}/"
    url = ("https://storage.googleapis.com/storage/v1/b/%s/o?" % BUCKET) + urllib.parse.urlencode({"prefix": prefix, "maxResults": 1000})
    objs = json.load(urllib.request.urlopen(url, timeout=90)).get("items", [])
    # Zero-byte objects are the folder placeholders in this bucket
    # (older products also list zero-byte "$folder$" markers, which would confuse ACOLITE's product detection)
    objs = [o for o in objs if not o["name"].endswith("/") and "$folder$" not in o["name"]]
    total = sum(int(o["size"]) for o in objs)
    print(f"{product}: {len(objs)} files, {total / 1e6:.0f} MB", flush=True)
    done = 0
    for o in objs:
        rel = o["name"][len(prefix):]
        dest = os.path.join(OUT, product, rel)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        if os.path.exists(dest) and os.path.getsize(dest) == int(o["size"]):
            done += int(o["size"]); continue
        link = f"https://storage.googleapis.com/{BUCKET}/" + urllib.parse.quote(o["name"])
        with urllib.request.urlopen(link, timeout=300) as r, open(dest + ".part", "wb") as f:
            while True:
                chunk = r.read(1 << 20)
                if not chunk: break
                f.write(chunk)
        os.replace(dest + ".part", dest)
        done += int(o["size"])
        print(f"  {done / 1e6:6.0f} / {total / 1e6:.0f} MB  {rel}", flush=True)
    print("done:", os.path.join(OUT, product))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
