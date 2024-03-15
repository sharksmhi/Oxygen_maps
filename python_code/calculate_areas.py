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


def area_at_threshold(threshold, ds):
    # Step 1: Create a 3D mask where Oxygen values are below a given threshold
    threshold = threshold  # Example threshold value

    # Step 2: Create a 3D mask where Oxygen values are below the threshold and False everywhere else
    mask_below_threshold = ds['Oxygen'] <= threshold
    ds[f"{threshold}_mask"] = mask_below_threshold
    # Step 3: Get the first layer in the depth dimension where Oxygen is below the threshold
    # First, find the indices where the condition is True along the depth dimension
    first_layer_indices = mask_below_threshold.argmax(dim='depth')

    # Then, use these indices to extract the corresponding depth values
    ds[f"Min_depth_{threshold}"] = ds['depth'].where(mask_below_threshold).isel(depth=first_layer_indices)

    # Step 5: Extract corresponding values for the variable Oxygen_relerr using the mask from step 4
    ds[f'Relerr_per_grid_at_min_{threshold}_depth'] = ds['Oxygen_relerr'].where(mask_below_threshold).isel(
        depth=first_layer_indices)


def calculate_areas(results_dir, file_list, threshold_list):

    area_results = []
    fig, axs = plt.subplots(1, 1, figsize=(10, 8))
    fig2, axs2 = plt.subplots(1, 1, figsize=(10, 8))
    for netcdf_filename in file_list:

        print(netcdf_filename)
        ds = xr.open_dataset(f"{results_dir}nc/O2/{netcdf_filename}")
        ### Calculate area of all grid cells
        # Get the latitude and longitude coordinates
        lon, lat = np.meshgrid(ds['lon'], ds['lat'])
        grid_areas = calculate_grid_areas(latitudes=lat, longitudes=lon)
        # assign the calculated areas to the dataset and return the updated dataset
        ds = ds.assign(grid_area=(('lat', 'lon'), grid_areas))

        var_name = "Oxygen"
        ### Find anox gridcells and set them to 1 all others to NaN
        #To do: Ändra så att vi kan skicka in en lista på gränsvärden.
        for threshold in threshold_list:
            area_at_threshold(threshold,ds)
        #----------------------------
        #Här nedan skall vi rensa. Men det skrivs areor till nc-filen som vi använder för att skapa bar-plotstidsserien. Det får vi fixa nästa gång.
        #
        hypox = 90
        anox = 9      # Detta är kanske för lågt? borde vara kanske 9 µmol/l = ~0.2 ml/l? 18 = ~0.4 ml/l
        relerr_lim = 0.5
        unit = 'umol/l'


        ds["Depth_hypoxia"]=xr.where((ds[var_name]<=hypox),ds['depth'],ds[var_name]*np.nan,keep_attrs=True)
        # Find the minimum threshold depth at each (time, lon, lat) coordinate
        ds['Min_depth_hypoxia']=ds['Depth_hypoxia'].min(dim='depth', skipna=True)


        ds["Depth_anoxia"]=xr.where((ds[var_name]<=anox),ds['depth'],ds[var_name]*np.nan,keep_attrs=True)
        # Find the minimum threshold depth at each (time, lon, lat) coordinate
        ds['Min_depth_anoxia']=ds['Depth_anoxia'].min(dim='depth', skipna=True)

        ### Sum hypox areas per season
        ds["Hypoxic_area"]=xr.where((ds['Min_depth_hypoxia']>=0),ds['grid_area'],ds['Min_depth_hypoxia']*np.nan,keep_attrs=True).sum(dim=['lat', 'lon'], skipna=True)
        #HRelativt fel på alla grundaste djup som har hypoxi.
        ds["Hypoxic_relerr_per_depth"] = xr.where(((ds['depth'] == ds['Min_depth_hypoxia'])), ds['Oxygen_relerr'],
                                           ds['Min_depth_hypoxia'] * 1000., keep_attrs=True)
        #Summerar relativt felet över alla djup till en platt matris, dvs relativt fel, per gridcell och år.
        #ds["Hypoxic_relerr_per_grid"]=xr.where(((ds['depth']==ds['Min_depth_hypoxia'])),ds['Oxygen_relerr'],ds['Min_depth_hypoxia']*np.nan,keep_attrs=True).sum(dim=['depth'], skipna=True)
        ds["Hypoxic_relerr_per_grid"] = ds["Hypoxic_relerr_per_depth"].sum(dim=['depth'], skipna=True)
        ds["Hypoxic_relerr_area"] = xr.where((ds['Hypoxic_relerr_per_grid']>=relerr_lim),ds['grid_area'],ds['Min_depth_hypoxia']*np.nan,keep_attrs=True).sum(dim=['lat', 'lon'], skipna=True)

        ### Sum anoxic areas per season
        ds["Anoxic_area"] = xr.where((ds['Min_depth_anoxia'] >= 0), ds['grid_area'], ds['Min_depth_anoxia'] * np.nan,
                                      keep_attrs=True).sum(dim=['lat', 'lon'], skipna=True)
        # HRelativt fel på alla grundaste djup som har anoxi.
        ds["Anoxic_relerr_per_depth"] = xr.where(((ds['depth'] == ds['Min_depth_anoxia'])), ds['Oxygen_relerr'],
                                                  ds['Min_depth_anoxia'] * np.nan, keep_attrs=True)

        # Summerar relativt felet över alla djup till en platt matris, dvs relativt fel, per gridcell och år.
        ds["Anoxic_relerr_per_grid"] = ds["Anoxic_relerr_per_depth"].sum(dim=['depth'], skipna=True)
        ds["Anoxic_relerr_area"] = xr.where((ds['Anoxic_relerr_per_grid'] >= relerr_lim), ds['grid_area'],
                                             ds['Min_depth_anoxia'] * np.nan, keep_attrs=True).sum(dim=['lat', 'lon'],
                                                                                                    skipna=True)

        # save the updated dataset
        print(f"saving {netcdf_filename}...")
        ds.to_netcdf(f'{results_dir}nc/processed/{netcdf_filename}') # rewrite to netcdf
        print(f'{netcdf_filename} has been modified')


