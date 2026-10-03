"""Run ACOLITE (Dark Spectrum Fitting) on a downloaded Sentinel-2 Level-1C SAFE product, as in Lee et al. (2021).

Needs the ACOLITE source (git clone --depth 1 https://github.com/acolite/acolite) and a Python environment with
GDAL (conda-forge: gdal, libgdal-jp2openjpeg, libgdal-netcdf, netcdf4, scipy, pyproj, scikit-image, h5py, requests).
Lookup tables are fetched automatically from github.com/acolite/acolite_luts the first time they are needed.

Settings follow Lee et al.: 20 m target resolution, turbidity from the Nechad et al. (2009) red-band algorithm
(665 nm). No ancillary meteorological data are used (that download needs a NASA Earthdata login), so ACOLITE's
default pressure and ozone are assumed; Dark Spectrum Fitting needs no external inputs.

Usage: python sentinel2/run_acolite.py <acolite_dir> <SAFE_dir> <output_dir>
"""
import sys, os

acolite_dir, safe, outdir = sys.argv[1:4]
sys.path.insert(0, acolite_dir)
import acolite as ac

os.makedirs(outdir, exist_ok=True)
settings = {
    "inputfile": safe,
    "output": outdir,
    "limit": [37.60, -122.00, 38.20, -121.15],      # south, west, north, east: the central Delta and the San Joaquin corridor
    "s2_target_res": 20,
    "l2w_parameters": ["rhow_665", "tur_nechad2009_665"],
    "ancillary_data": False,
    "rgb_rhot": False,
    "rgb_rhos": False,
    "map_l2w": False,
    "l1r_delete_netcdf": True,                       # the 200 MB L1R intermediate is not needed after the run
}
# Optional fourth argument "radcor" switches the atmospheric correction to RAdCor (Castagna and Vanhellemont 2025),
# which includes a correction for adjacency effects from land around narrow channels. The default (as in
# Lee et al. 2021) is Dark Spectrum Fitting without an adjacency correction.
if len(sys.argv) > 4 and sys.argv[4] == "radcor":
    settings.update({"atmospheric_correction_method": "radcor"})
# Fourth argument "glint" turns on the residual sun glint correction in Dark Spectrum Fitting (a spectrally flat
# offset estimated from the shortwave-infrared bands and removed from all bands).
if len(sys.argv) > 4 and sys.argv[4] == "glint":
    settings.update({"dsf_residual_glint_correction": True, "dsf_residual_glint_wave_range": [1500, 2400]})
ac.acolite.acolite_run(settings=settings)
