"""Matchups between ACOLITE-derived turbidity and the continuous stations.

Reads the ACOLITE L2R files (surface reflectance, rhos_*) in an output directory, not the L2W files: ACOLITE's own
water mask (rhot at 1610 nm below 0.0215) removed almost every Delta pixel in the test scene, so water is
identified here instead.

Water pixel: NDWI > 0 (green vs near-infrared), near-infrared reflectance below 0.05, red reflectance above zero.
Turbidity: Nechad et al. (2009) red band, T = A * rho / (1 - rho / C) with A = 282.95 FNU and C = 0.1728 (665 nm).
Station value: median over the water pixels within `radius` metres of the station coordinates (the CDEC coordinates
sit on piers and banks, where the 20 m pixel is mixed with land), compared with the median of the station's
readings within one hour of the overpass. Run with ACOLITE's dsf_residual_glint_correction on (sentinel2/run_acolite.py
... glint): without it the near-infrared and shortwave-infrared reflectance of the Delta's water is a flat 0.03-0.07
and almost no pixel qualifies as water.

Usage: python sentinel2/acolite_matchups.py <output_dir_with_L2R_nc_files>
"""
import glob, os, sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pilot_matchups import station_series, FNU_STATIONS, nechad   # same station readers, unit groups and algorithm

STATIONS = {"SJR": (37.679290, -121.265110), "MSD": (37.786000, -121.306000), "RRI": (37.963000, -121.365000),
            "PPT": (38.066000, -121.562000), "TWI": (38.096900, -121.669100), "BLP": (38.032452, -121.719040),
            "ANH": (38.017830, -121.802963), "TRN": (37.992800, -121.454200), "OSJ": (38.071100, -121.578900),
            "MOK": (38.106100, -121.571100)}
RADII = (100, 200, 400)
MIN_PIXELS = 5


def main(outdir):
    rows = []
    for f in sorted(glob.glob(os.path.join(outdir, "*_L2R.nc"))):
        ds = Dataset(f)
        lat, lon = ds["lat"][:].filled(np.nan), ds["lon"][:].filled(np.nan)
        rhos = {k: ds[k][:].filled(np.nan) for k in ["rhos_559", "rhos_665", "rhos_833"]}
        ndwi = (rhos["rhos_559"] - rhos["rhos_833"]) / (rhos["rhos_559"] + rhos["rhos_833"])
        water = (ndwi > 0) & (rhos["rhos_833"] < 0.05) & (rhos["rhos_665"] > 0)
        tur = nechad(rhos["rhos_665"])
        stamp = pd.Timestamp(ds.getncattr("isodate")).tz_localize(None)
        date = stamp.strftime("%Y-%m-%d")
        t_pst = stamp - pd.Timedelta(hours=8)
        print(f"{os.path.basename(f)}: water-like pixels {water.mean():.2%} of the window")
        for st, (la, lo) in STATIONS.items():
            dm = np.sqrt((lat - la) ** 2 + ((lon - lo) * np.cos(np.radians(la))) ** 2) * 111000
            if not np.nanmin(dm) < 40:
                continue                                       # station is outside this tile
            s = station_series(st, int(date[:4]))
            w = s[(s.index >= t_pst - pd.Timedelta(hours=1)) & (s.index <= t_pst + pd.Timedelta(hours=1))]
            row = {"date": date, "station": st, "unit": "FNU" if st in FNU_STATIONS else "NTU",
                   "sonde": float(w.median()) if len(w) >= 4 else np.nan}
            for r in RADII:
                m = (dm <= r) & water & np.isfinite(tur)
                row[f"n{r}"] = int(m.sum())
                row[f"s2_{r}m"] = float(np.nanmedian(tur[m])) if m.sum() >= MIN_PIXELS else np.nan
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, "acolite_matchups.csv"), index=False)
    print("\n" + df.round(2).to_string(index=False))
    for r in RADII:
        ok = df.dropna(subset=[f"s2_{r}m", "sonde"])
        if len(ok) >= 4:
            x, y = ok[f"s2_{r}m"], ok["sonde"]
            print(f"\nradius {r} m: n={len(ok)}  Spearman={ok[[f's2_{r}m', 'sonde']].corr(method='spearman').iloc[0, 1]:.2f}  "
                  f"R2={np.corrcoef(x, y)[0, 1] ** 2:.2f}  median ratio S2/sonde={np.median(x / y):.2f}  "
                  f"median bias={np.median(x - y):.1f}  median abs % error={100 * np.median(np.abs(x / y - 1)):.0f}%")


if __name__ == "__main__":
    main(sys.argv[1])
