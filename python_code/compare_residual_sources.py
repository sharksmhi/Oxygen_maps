import xarray as xr
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

year = 1970
season = "Autumn"
results_dir = Path(f"/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results_lena_temp/Baltic_Proper/20260320_1659/")

############### netcdf data ###################
pattern = f"Oxygen_{year}_{season}_0*.nc"
files = list((results_dir / "DIVArun").glob(pattern))
if not files:
    raise FileNotFoundError(f"No file found for {pattern}")
if len(files) > 1:
    print("Warning: multiple files found, using first one:", files[0])
# take the first match
observation_nc = files[0]

residuals_nc = results_dir / f"DIVArun/Oxygen_{year}_{season}_residuals.nc"
ds_obs = xr.open_dataset(observation_nc, engine='h5netcdf')
ds_res = xr.open_dataset(residuals_nc, engine='h5netcdf')
obs = ds_obs["Oxygen_data"].values
res = ds_res["Oxygen_residual"].values
lon = ds_obs["obslon"].values
lat = ds_obs["obslat"].values
depth = ds_obs["obsdepth"].values
time = ds_obs["obstime"].values
sel = ~np.isnan(res)

obs_sel = obs[sel]
res_sel = res[sel]

diva_result = obs_sel - res_sel

################ textfile data #############

filepath=results_dir / f"DIVArun/Oxygen_{year}_{season}_residual.txt"
df = pd.read_csv(filepath, sep="\t")

# Extrahera variabler
obs = df["obsval"]
interp = df["diva"]
residual = df["residual"]
depth = df["obsdepth"]


############## plot ##################
plt.scatter(res_sel, residual, alpha=0.5)
plt.xlabel("residuals from netcdf")
plt.ylabel("residuals from textfile")
plt.show()