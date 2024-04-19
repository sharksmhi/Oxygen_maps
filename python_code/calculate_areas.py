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


def area_at_threshold(threshold, ds, df):
    relerr_lim = 0.5
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

    ds[f"Area_{threshold}"] = xr.where((ds[f'Min_depth_{threshold}'] >= 0), ds['grid_area'], ds[f'Min_depth_{threshold}'] * np.nan,
                                      keep_attrs=True).sum(dim=['lat', 'lon'], skipna=True)
    # HRelativt fel på alla grundaste djup som har anoxi.
    ds[f"{threshold}_relerr_per_depth"] = xr.where(((ds['depth'] == ds[f'Min_depth_{threshold}'])), ds['Oxygen_relerr'],
                                                  ds[f'Min_depth_{threshold}'] * np.nan, keep_attrs=True)

    # Summerar relativt felet över alla djup till en platt matris, dvs relativt fel, per gridcell och år.
    ds[f"{threshold}_relerr_per_grid"] = ds[f"{threshold}_relerr_per_depth"].sum(dim=['depth'], skipna=True)
    ds[f"Relerr_area_{threshold}"] = xr.where((ds[f'{threshold}_relerr_per_grid'] >= relerr_lim), ds['grid_area'],
                                             ds[f'Min_depth_{threshold}'] * np.nan, keep_attrs=True).sum(dim=['lat', 'lon'],
                                                                                                    skipna=True)
    df[f"Area_{threshold}_km2"] = ds[f"Area_{threshold}"].round(3)
    df[f"Relerr_area_{threshold}_km2"] = ds[f"Relerr_area_{threshold}"].round(3)


def calculate_areas(results_dir, file_list, threshold_list, save_area_data=False):
    #List for results to txt
    area_results = []
    fig, axs = plt.subplots(1, 1, figsize=(10, 8))
    fig2, axs2 = plt.subplots(1, 1, figsize=(10, 8))
    for netcdf_filename in file_list:
        ds = xr.open_dataset(f"{results_dir}nc/O2/{netcdf_filename}")
        season = ds.attrs['season']
        start_year = ds.attrs['start year']
        end_year = ds.attrs['end year']
        ### Calculate area of all grid cells
        # Get the latitude and longitude coordinates
        lon, lat = np.meshgrid(ds['lon'], ds['lat'])
        grid_areas = calculate_grid_areas(latitudes=lat, longitudes=lon)
        # assign the calculated areas to the dataset and return the updated dataset
        ds = ds.assign(grid_area=(('lat', 'lon'), grid_areas))

        var_name = "Oxygen"
        relerr_lim = 0.5
        ### Find anox gridcells and set them to 1 all others to NaN
        #To do: Ändra så att vi kan skicka in en lista på gränsvärden.
        df = pd.DataFrame()
        for threshold in threshold_list:
            area_at_threshold(threshold,ds,df)
        df["year"] = ds.time.values.astype('datetime64[Y]').astype(int) + 1970
        df["season"] = season
        area_results.append(df)

        # save the updated dataset
        print(f"writing to {results_dir}nc/processed/{netcdf_filename}...")
        ds.to_netcdf(f'{results_dir}nc/processed/{netcdf_filename}') # rewrite to netcdf

    # combing area results from all seasons and saving to a textfile
    if save_area_data:
        pd.concat(area_results).to_csv(f'{results_dir}area_data_{start_year}_{end_year}.txt', sep='\t', index=False)

def area_data(results_dir,file_list):
    area_results=[]
    for netcdf_filename in file_list:
        ds = xr.open_dataset(f"{results_dir}nc/processed/{netcdf_filename}", engine='h5netcdf')
        season = ds.attrs['season']
        start_year = ds.attrs['start year']
        end_year = ds.attrs['end year']

        df = pd.DataFrame()
        df["Hypoxic_area_km2"] = ds["Hypoxic_area"].round(3)
        df["Hypoxic_relerr_area_km2"] = ds["Hypoxic_relerr_area"].round(3)
        df["Anoxic_area_km2"] = ds["Anoxic_area"].round(3)
        df["Anoxic_relerr_area_km2"] = ds["Anoxic_relerr_area"].round(3)
        df["year"] = ds.time.values.astype('datetime64[Y]').astype(int) + 1970
        df["season"] = season
        area_results.append(df)

    # combing area results from all seasons and saving to a textfile
    pd.concat(area_results).to_csv(f'{results_dir}area_data_{start_year}_{end_year}.txt', sep='\t', index=False)

if __name__ == "__main__":
    # Result directory
    results_dir = "//winfs-proj/proj/havgem/DIVA/syrekartor/resultat/"
    # Open the JSON file
    file_list = ["Oxygen_1960-2021_Autumn_0.2_80000_0.05_5.0_2.0_gebco_30sec_4_varcorrlenz.nc","Oxygen_1960-2021_Winter_0.2_80000_0.05_5.0_2.0_gebco_30sec_4_varcorrlenz.nc",
                 "Oxygen_1960-2021_Summer_0.2_80000_0.05_5.0_2.0_gebco_30sec_4_varcorrlenz.nc", "Oxygen_1960-2021_Spring_0.2_80000_0.05_5.0_2.0_gebco_30sec_4_varcorrlenz.nc"]
    # Thresholds to analyse in µmol/l oxygen (0, 2, 4 ml/l)
    threshold_list = [0, 90, 180]
    calculate_areas(results_dir, file_list, threshold_list)
    area_data(results_dir, file_list)