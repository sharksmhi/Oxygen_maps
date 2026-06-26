import xarray as xr
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import matplotlib.colors as mcolors
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as cm 
import seaborn as sns

year = 1968
season = "Autumn"
results_dir_kt = Path(f"/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results/Kattegat/20260622_1633/full/")
results_dir_bp = Path(f"/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results/Baltic_Proper/20260622_1637/full/")

############### netcdf data ###################
observation_nc_kt = results_dir_kt / f"DIVArun/Oxygen_{year}_{season}_0.2_25000.0_50000.0_0.025_5.0_2.0_Kattegat_varcorrlenz"
observation_nc_bp = results_dir_bp / f"DIVArun/Oxygen_{year}_{season}_0.2_25000.0_50000.0_0.025_5.0_2.0_Baltic_Proper_varcorrlenz"
ds_obs_kt = xr.open_dataset(observation_nc_kt, engine='h5netcdf')
ds_obs_bp = xr.open_dataset(observation_nc_bp, engine='h5netcdf')
obs_bp = ds_obs_bp["Oxygen_data"].values
obs_kt = ds_obs_kt["Oxygen_data"].values

diff = ds_obs_bp["Oxygen_data"]- ds_obs_kt["Oxygen_data"]




diff = obs_bp-obs_kt