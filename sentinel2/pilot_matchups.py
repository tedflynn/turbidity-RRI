"""Sentinel-2 pilot, step 2: turbidity from the chips and matchups with the continuous stations.

Turbidity from the L2A (Sen2Cor) red-band surface reflectance with the Nechad et al. (2009) algorithm,
T = A * rho / (1 - rho / C), using the coefficients ACOLITE uses at 665 nm (A = 282.95 FNU, C = 0.1728), as in
Lee et al. (2021) except for the atmospheric correction (Lee used ACOLITE dark spectrum fitting on Level-1C).

Matchups follow Lee et al.: median of the water pixels in a 3 x 3 window around the station, against the median
of the station's readings within one hour of the overpass. Also reports a 5 x 5 window and a local calibration.

Station data: data/<station>_data.csv (to 2025), data/CDEC_<station>_turb_2026.csv (2026), RRI from the 1 m sonde
file (clock values are UTC, shifted 8 h to PST as in turbidity-depth-RRI.R). Stations reporting NTU are
SJR, MSD, RRI, PPT, TWI, ANH; FNU: BLP, TRN, OSJ, MOK.
"""
import json, os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PILOT = os.path.join(DATA, "s2", "pilot")
A_T, C_T = 282.95, 0.1728
FNU_STATIONS = {"BLP", "TRN", "OSJ", "MOK"}


def nechad(rho):
    rho = np.asarray(rho, dtype="float64")
    out = A_T * rho / (1 - rho / C_T)
    out[(rho <= 0) | (rho >= C_T)] = np.nan
    return out


def station_series(st, year):
    """Return a Series of turbidity indexed by PST naive timestamps."""
    if year >= 2026:
        d = pd.read_csv(os.path.join(DATA, f"CDEC_{st}_turb_2026.csv"), dtype=str)
        t = pd.to_datetime(d["DATE TIME"], format="%Y%m%d %H%M", errors="coerce")
        v = pd.to_numeric(d["VALUE"], errors="coerce")
    elif st == "RRI":
        d = pd.read_csv(os.path.join(DATA, "RRI_data_1m.csv"), dtype=str)
        clock = pd.to_datetime(d["date_time"].str.slice(0, 19), format="%Y-%m-%dT%H:%M:%S", errors="coerce")
        t = clock - pd.Timedelta(hours=8)          # clock values are UTC mislabelled -08:00
        v = pd.to_numeric(d["turb_fnu"], errors="coerce")
    else:
        d = pd.read_csv(os.path.join(DATA, f"{st}_data.csv"), dtype=str)
        if "Date_Time" in d.columns:
            t = pd.to_datetime(d["Date_Time"], format="%m/%d/%Y %H:%M:%S", errors="coerce")
            v = pd.to_numeric(d["Turbidity"], errors="coerce")
        else:
            t = pd.to_datetime(d["date"] + " " + d["time"], format="%Y-%m-%d %H:%M:%S", errors="coerce")
            v = pd.to_numeric(d["turbidity"], errors="coerce")
    s = pd.Series(v.values, index=t).dropna()
    s = s[(s > 0) & (s < 1000)]
    return s[~s.index.duplicated()].sort_index()


def chip_turbidity(npz, half_window):
    z = np.load(npz)
    red, green, nir, scl = z["red"], z["green"], z["nir"], z["scl"]
    c = red.shape[0] // 2
    ndwi = (green - nir) / (green + nir)
    water = (ndwi > 0) & ~np.isin(scl, [1, 3, 8, 9, 10, 11]) & np.isfinite(red) & (red > 0) & (nir < 0.2)
    sl = (slice(c - half_window, c + half_window + 1), slice(c - half_window, c + half_window + 1))
    rho = np.where(water, red, np.nan)[sl]
    rho_valid = rho[np.isfinite(rho)]
    if rho_valid.size == 0:
        return np.nan, 0, np.nan
    t = nechad(rho_valid)
    return float(np.nanmedian(t)), int(rho_valid.size), float(np.nanmedian(rho_valid))


def main():
    man = json.load(open(os.path.join(PILOT, "manifest.json")))
    rows = []
    cache = {}
    for m in man:
        st, d = m["station"], m["date"]
        year = int(d[:4])
        key = (st, year >= 2026)
        if key not in cache:
            cache[key] = station_series(st, year)
        s = cache[key]
        t_pst = pd.Timestamp(m["datetime_utc"]).tz_localize(None) - pd.Timedelta(hours=8)
        w = s[(s.index >= t_pst - pd.Timedelta(hours=1)) & (s.index <= t_pst + pd.Timedelta(hours=1))]
        t3, n3, r3 = chip_turbidity(os.path.join(PILOT, m["file"]), 1)
        t5, n5, _ = chip_turbidity(os.path.join(PILOT, m["file"]), 2)
        rows.append({"date": d, "year": year, "station": st, "unit": "FNU" if st in FNU_STATIONS else "NTU",
                     "s2_3x3": t3, "n3": n3, "rho_red": r3, "s2_5x5": t5, "n5": n5,
                     "sonde": float(w.median()) if len(w) >= 4 else np.nan, "n_sonde": len(w)})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(PILOT, "matchups.csv"), index=False)

    ok = df.dropna(subset=["s2_3x3", "sonde"]).query("n3 >= 5")
    print(f"scene-station pairs: {len(df)}; with a sonde value: {df['sonde'].notna().sum()}; "
          f"with >=5 water pixels in 3x3: {(df['n3'] >= 5).sum()}; usable matchups: {len(ok)}")
    print("\nPairs lost, by station (no sonde value / fewer than 5 water pixels):")
    print(df.assign(no_sonde=df["sonde"].isna(), few_pixels=df["n3"] < 5).groupby("station")[["no_sonde", "few_pixels"]].sum().T.to_string())

    def stats(x, label):
        if len(x) < 4:
            print(f"{label}: n={len(x)} (too few)"); return
        lr = np.log10(x["s2_3x3"] / x["sonde"])
        r = np.corrcoef(x["s2_3x3"], x["sonde"])[0, 1]
        rl = np.corrcoef(np.log10(x["s2_3x3"]), np.log10(x["sonde"]))[0, 1]
        rs = x[["s2_3x3", "sonde"]].corr(method="spearman").iloc[0, 1]
        print(f"{label}: n={len(x)}  R2={r*r:.2f}  R2(log)={rl*rl:.2f}  Spearman={rs:.2f}  "
              f"median ratio S2/sonde={10**lr.median():.2f}  RMSE={np.sqrt(((x['s2_3x3']-x['sonde'])**2).mean()):.1f}  "
              f"median abs % error={100*np.median(np.abs(x['s2_3x3']/x['sonde']-1)):.0f}%")

    print("\nDefault Nechad (665 nm) against sonde, 3x3 window:")
    stats(ok, "all")
    for u in ["NTU", "FNU"]:
        stats(ok[ok["unit"] == u], u + " stations")
    stats(ok[ok["sonde"] < 5], "sonde < 5")
    stats(ok[(ok["sonde"] >= 5) & (ok["sonde"] < 10)], "sonde 5-10")
    stats(ok[ok["sonde"] >= 10], "sonde >= 10")

    print("\nBy station (median ratio S2/sonde, Spearman):")
    g = ok.groupby("station").apply(lambda x: pd.Series({
        "n": len(x), "median_sonde": x["sonde"].median(), "median_s2": x["s2_3x3"].median(),
        "ratio": 10 ** np.log10(x["s2_3x3"] / x["sonde"]).median(),
        "spearman": x[["s2_3x3", "sonde"]].corr(method="spearman").iloc[0, 1] if len(x) > 3 else np.nan}), include_groups=False)
    print(g.round(2).to_string())

    print("\nBy year (does the retrieval bias drift?): median ratio S2/sonde")
    print(ok.groupby("year").apply(lambda x: pd.Series({"n": len(x), "ratio": 10 ** np.log10(x["s2_3x3"] / x["sonde"]).median()}),
                                   include_groups=False).round(2).to_string())

    # Local calibration of A with C fixed, in log space, leave-one-date-out
    def fit_A(x):
        rho = x["rho_red"].values
        base = rho / (1 - rho / C_T)
        return 10 ** np.median(np.log10(x["sonde"].values / base))
    ok = ok.copy()
    preds = []
    for d in ok["date"].unique():
        tr, te = ok[ok["date"] != d], ok[ok["date"] == d]
        A = fit_A(tr)
        preds.append(te.assign(s2_cal=A * te["rho_red"] / (1 - te["rho_red"] / C_T)))
    cv = pd.concat(preds)
    lr = np.log10(cv["s2_cal"] / cv["sonde"])
    print(f"\nLocal calibration (A refitted, leave-one-date-out): fitted A over all data = {fit_A(ok):.0f} FNU (default {A_T});")
    print(f"  cross-validated: R2={np.corrcoef(cv['s2_cal'], cv['sonde'])[0,1]**2:.2f}  "
          f"median abs % error={100*np.median(np.abs(cv['s2_cal']/cv['sonde']-1)):.0f}%  "
          f"Spearman={cv[['s2_cal','sonde']].corr(method='spearman').iloc[0,1]:.2f}")

    print("\nS2 window sensitivity (3x3 vs 5x5): Spearman with sonde:",
          round(ok[["s2_3x3", "sonde"]].corr(method="spearman").iloc[0, 1], 2), "vs",
          round(ok.dropna(subset=["s2_5x5"])[["s2_5x5", "sonde"]].corr(method="spearman").iloc[0, 1], 2))


if __name__ == "__main__":
    main()
