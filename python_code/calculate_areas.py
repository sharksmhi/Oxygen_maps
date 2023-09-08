import xarray as xr
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt 

### Program to find anox and hypox gridcells
### Also to get the depth and min depth for anox and hypox cells
### Saves the results back to a nc-file
### Calculate a corrected area grid to be used when calculationg the total areas of hypox and anox

def haversine(lat1, lon1, lat2, lon2, radius=6371):
    """
    Calculate the great-circle distance between two points
    on the Earth's surface using the Haversine formula.
    """
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    distance = radius * c
    return distance

def calculate_grid_areas(latitudes, longitudes):
    """
    Calculate the approximate area of a grid on a sphere
    using latitude and longitude coordinates.
    """
    radius = 6371  # Radius of the Earth in kilometers

    # Initialize an array to store the grid cell areas
    areas = np.zeros_like(latitudes)

    # Iterate over each grid cell
    for i in range(latitudes.shape[0]-1):
        for j in range(latitudes.shape[1]-1):
            lat1 = latitudes[i, j]
            lon1 = longitudes[i, j]
            lat2 = latitudes[i+1, j+1]
            lon2 = longitudes[i+1, j+1]

            # Calculate the distance between the latitude and longitude points
            lat_distance = haversine(lat1, lon1, lat2, lon1, radius)
            lon_distance = haversine(lat1, lon1, lat1, lon2, radius)

            # Convert the distances from kilometers to meters
            #lat_distance_meters = lat_distance * 1000
            #lon_distance_meters = lon_distance * 1000

            # Calculate the area of the rectangular grid cell
            #area = lat_distance_meters * lon_distance_meters
            area = lat_distance * lon_distance
            # Store the area in the corresponding grid cell
            areas[i, j] = area

    return areas

location = "//winfs-proj/proj/havgem/DIVA/syrekartor/" # or other location, like havgem path
df = pd.DataFrame()
fig, axs = plt.subplots(1, 1, figsize=(10, 8))
fig2, axs2 = plt.subplots(1, 1, figsize=(10, 8))
for season in ['Winter', 'Spring', 'Summer', 'Autumn']:
    ### open netcdf file ###
    netcdf_filename = f"Oxygen_1960-2018_{season}_1_39000.0_gebco_30sec_4"
    ds = xr.open_dataset(f"{location}/resultat/nc/O2/{netcdf_filename}.nc")

    ### extract values that are within our limits, save to a new variable and nc-file. ####
    # 1 ml/l of O2 is approximately 43.570 µmol/kg
    # (assumes a molar volume of O2 of 22.392 l/mole and 
    # a constant seawater potential density of 1025 kg/m3).
    # https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html
    hypox = 90
    anox = 4.357      # Detta är kanske för lågt? borde vara kanske ~4 µmol/l = 0.1 ml/l?
    unit = 'umol/l'

    var_name = "Oxygen"
    ### Find anox gridcells and set them to 1 all others to NaN
    ds["ANOX"]=xr.where((ds[var_name]<=anox),ds[var_name]/ds[var_name]*1,ds[var_name]*np.nan,keep_attrs=True)
    ds["HYPOX"]=xr.where((ds[var_name]<=hypox),ds[var_name]/ds[var_name]*1,ds[var_name]*np.nan,keep_attrs=True)

    threshold = hypox  # Define the threshold value
    ds["hypox_depth"]=xr.where((ds[var_name]<=hypox),ds['depth'],ds[var_name]*np.nan,keep_attrs=True)
    # Find the minimum threshold depth at each (time, lon, lat) coordinate
    ds['min_hypox_depth']=ds['hypox_depth'].min(dim='depth', skipna=True)
    ds['min_hypox_depth_relerr']=ds['hypox_depth'].min(dim='depth', skipna=True)
    threshold = anox  # Define the threshold value
    ds["anox_depth"]=xr.where((ds[var_name]<=anox),ds['depth'],ds[var_name]*np.nan,keep_attrs=True)
    # Find the minimum threshold depth at each (time, lon, lat) coordinate
    ds['min_anox_depth']=ds['anox_depth'].min(dim='depth', skipna=True)

    ### Calculate area of all grid cells
    # Get the latitude and longitude coordinates
    lon, lat = np.meshgrid(ds['lon'], ds['lat'])
    grid_areas = calculate_grid_areas(latitudes=lat, longitudes=lon)
    # assign the calculated areas to the dataset and return the updated dataset
    ds = ds.assign(grid_area = (('lat', 'lon'), grid_areas))

    ### Sum anox and hypox areas per season
    ds["HYPOX_area"]=xr.where((ds['min_hypox_depth']>=0),ds['grid_area'],ds['min_hypox_depth']*np.nan,keep_attrs=True).sum(dim=['lat', 'lon'], skipna=True)
    ds["HYPOX_area_relerr"]=xr.where(((ds['depth']==ds['min_hypox_depth'])),ds['Oxygen_relerr'],ds['min_hypox_depth']*np.nan,keep_attrs=True).sum(dim=['depth'], skipna=True)
    #ds["HYPOX_area_relerr"]=xr.where(((ds['depth']<=ds['min_hypox_depth'])),ds['Oxygen_relerr'],ds['min_hypox_depth']*np.nan,keep_attrs=True).sum(dim=['lat', 'lon'], skipna=True)
    ds["HYPOX_area_relerr_1"] = xr.where(((ds[var_name] <= hypox)), ds['Oxygen_relerr'],ds['min_hypox_depth'] * np.nan, keep_attrs=True)
    ds["HYPOX_area_relerr_1"] = xr.where(((ds[var_name] <= hypox)), ds['Oxygen_relerr'],ds['min_hypox_depth'] * np.nan, keep_attrs=True)

    ### plot the resulting timeseries
    # ds["HYPOX_area"].plot(ax=axs, label = season)
    ds["HYPOX_area"].plot(ax=axs, label = season)
    ds["HYPOX_area_relerr"].plot(ax=axs2, label = season)
    df[season] = ds["HYPOX_area"]
    # df[f"{season}_relerr"] = ds["HYPOX_area_relerr"]
    # save the updated dataset
    ds.to_netcdf(f'{location}/resultat/nc/processed/{netcdf_filename}.nc') # rewrite to netcdf
    print(f'{location}/resultat/nc/processed/{netcdf_filename}.nc has been modified')

df.set_index(ds.time.values, inplace=True)

axs.set_xlabel('Year')
axs.set_ylabel('Area [km3]')
axs.legend()
plt.savefig(f'{location}resultat/figures/hypox_area.png')

fig, axs = plt.subplots(1, 1, figsize=(10, 8))
df.plot.bar(ax = axs,)
plt.savefig(f'{location}resultat/figures/hypox_area_barchart.png')