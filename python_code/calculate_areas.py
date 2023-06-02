import xarray as xr
import numpy as np

### open netcdf file ###
netcdf_filename = "Oxygen_1960-2020_Autumn_1_gebco_30sec_4"
location = "." # or other location, like havgem path
ds = xr.open_dataset(f"{location}/resultat/nc/O2/{netcdf_filename}.nc")
print(ds)

# for att in ds.attrs:
#     print(att)
# exit()
# Read the bathymetry file using Xarray
bath_file = xr.open_dataset("./bathymetry/gebco_30sec_4.nc")
print(bath_file)
# Extract the required variables
b = bath_file["bat"]

### extract values that are within our limits, save to a new variable and nc-file. ####
# in ml/l
hypox = 2
anox = 0.1
unit = 'ml/l'
# in umol/l
# 1 ml/l of O2 is approximately 43.570 µmol/kg 
# (assumes a molar volume of O2 of 22.392 l/mole and 
# a constant seawater potential density of 1025 kg/m3).
# https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html
hypox = 90
anox = 0.1
unit = 'umol/l'

var_name = "Oxygen"
ds["ANOX"]=xr.where((ds[var_name]<=0.5),ds[var_name]/ds[var_name]*1,ds[var_name]*np.nan,keep_attrs=True)
ds["HYPOX"]=xr.where((ds[var_name]<=hypox),ds[var_name]/ds[var_name]*1,ds[var_name]*np.nan,keep_attrs=True)

threshold = hypox  # Define the threshold value
ds["hypox_depth"]=xr.where((ds[var_name]<=hypox),ds['depth'],ds[var_name]*np.nan,keep_attrs=True)
# Find the minimum threshold depth at each (time, lon, lat) coordinate
ds['min_hypox_depth']=ds['hypox_depth'].min(dim='depth', skipna=True)
threshold = anox  # Define the threshold value
ds["anox_depth"]=xr.where((ds[var_name]<=anox),ds['depth'],ds[var_name]*np.nan,keep_attrs=True)
# Find the minimum threshold depth at each (time, lon, lat) coordinate
ds['min_anox_depth']=ds['anox_depth'].min(dim='depth', skipna=True)

# Print the updated dataset
print(ds)
ds.to_netcdf(f'resultat/nc/modified_{netcdf_filename}.nc') # rewrite to netcdf
