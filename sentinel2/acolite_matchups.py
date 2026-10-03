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


def rhos_band(ds, wave):
    """Surface reflectance band closest to `wave` nm. ACOLITE names bands by the sensor's own centre wavelength
    (for example 559 nm on Sentinel-2B but 560 nm on 2A and 2C), so the name differs between satellites."""
    names = [k for k in ds.variables if k.startswith("rhos_")]
    best = min(names, key=lambda k: abs(int(k.split("_")[1]) - wave))
    if abs(int(best.split("_")[1]) - wave) > 15:
        raise KeyError(f"no rhos band within 15 nm of {wave} in {list(names)}")
    return ds[best][:].filled(np.nan)


def main(outdir):
    rows = []
    for f in sorted(glob.glob(os.path.join(outdir, "*_L2R.nc"))):
        ds = Dataset(f)
        lat, lon = ds["lat"][:].filled(np.nan), ds["lon"][:].filled(np.nan)
        green, red, nir = (rhos_band(ds, w) for w in (560, 665, 833))
        ndwi = (green - nir) / (green + nir)
        water = (ndwi > 0) & (nir < 0.05) & (red > 0)
        tur = nechad(red)
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
    # A station inside the overlap of two tiles appears twice on a date: keep the row with more water pixels around it
    df = df.sort_values("n200", ascending=False).drop_duplicates(["date", "station"]).sort_values(["date", "station"]).reset_index(drop=True)
    df.to_csv(os.path.join(outdir, "acolite_matchups.csv"), index=False)
    print("\n" + df.round(2).to_string(index=False))
    for r in RADII:
        ok = df.dropna(subset=[f"s2_{r}m", "sonde"])
        if len(ok) >= 4:
            x, y = ok[f"s2_{r}m"], ok["sonde"]
            print(f"\nradius {r} m: n={len(ok)}  Spearman={ok[[f's2_{r}m', 'sonde']].corr(method='spearman').iloc[0, 1]:.2f}  "
                  f"R2={np.corrcoef(x, y)[0, 1] ** 2:.2f}  median ratio S2/sonde={np.median(x / y):.2f}  "
                  f"median bias={np.median(x - y):.1f}  median abs % error={100 * np.median(np.abs(x / y - 1)):.0f}%")
    # Does the satellite track the sondes across dates (not just across stations)? Radius 200 m.
    col = "s2_200m"
    ok = df.dropna(subset=[col, "sonde"]).assign(year=lambda x: x["date"].str.slice(0, 4).astype(int),
                                                 ls2=lambda x: np.log10(x[col]), lso=lambda x: np.log10(x["sonde"]))
    if ok["date"].nunique() >= 3:
        print(f"\n--- cross-date checks (radius 200 m; {ok['date'].nunique()} dates) ---")
        print("median ratio S2/sonde by year (is the bias stable?):")
        print(ok.groupby("year").apply(lambda x: pd.Series({"pairs": len(x), "ratio": 10 ** np.median(x["ls2"] - x["lso"])}),
                                       include_groups=False).round(2).to_string())
        within = ok.assign(a2=ok["ls2"] - ok.groupby("station")["ls2"].transform("mean"),
                           ao=ok["lso"] - ok.groupby("station")["lso"].transform("mean"))
        print(f"\ndate-to-date tracking within stations (anomalies from each station's own mean, log10): "
              f"Pearson={within['a2'].corr(within['ao']):.2f}, Spearman={within['a2'].corr(within['ao'], method='spearman'):.2f}, n={len(within)}")
        base_year = ok["year"].min()
        late = ok[ok["year"] > base_year].groupby("station")[["ls2", "lso"]].mean()
        early = ok[ok["year"] == base_year].groupby("station")[["ls2", "lso"]].mean()
        chg = (late - early).dropna()
        if len(chg) >= 4:
            print(f"\nchange from {base_year} to the later years, per station, in log10 (+ve = more turbid later); "
                  f"satellite vs sonde: Pearson={chg['ls2'].corr(chg['lso']):.2f}, n={len(chg)}")
            print((100 * (10 ** chg - 1)).round(0).rename(columns={"ls2": "S2 % change", "lso": "sonde % change"}).to_string())


if __name__ == "__main__":
    main(sys.argv[1])
