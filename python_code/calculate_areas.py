from pathlib import Path
import shapely
import xarray as xr
import numpy as np
import math
import pandas as pd
import geopandas as gpd

BASINS = gpd.read_file("./data/HELCOM_subbasins_2022_level2/HELCOM_subbasins_2022_level2.shp")
BASINS = BASINS.to_crs("EPSG:4326")
BASINS["basin_name"] = BASINS["level_2"].astype(str)
BASINS["basin_id"] = BASINS["HELCOM_ID"].astype(str)

BASINIDS = ["SEA-014", "SEA-006", "SEA-005", 
             "SEA-007", "SEA-009", "SEA-012", "SEA-010", 
             "SEA-008", "SEA-013", "SEA-011", "SEA-001", "SEA-004", "SEA-003"]
BASINIDS_SUBSET = ["SEA-007", "SEA-009", "SEA-012", "SEA-010"]

### Program to find anox and hypox gridcells
### Also to get the depth and min depth for anox and hypox cells
### Saves the results back to a nc-file
### Calculate a corrected area grid to be used when calculationg the total areas of hypox and anox

def haversine_single_grid(lat1, lon1, lat2, lon2, radius=6371):
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

def haversine_vectorized(lat1, lon1, lat2, lon2, radius=6371):
    """
    Vectorized Haversine distance between two sets of points.
    Inputs are arrays of the same shape.
    Returns distance in km.
    """
    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)

    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad

    a = np.sin(dlat / 2) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distance = radius * c
    return distance

def calculate_grid_areas_single_grids(latitudes, longitudes):
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
            lat_distance = haversine_single_grid(lat1, lon1, lat2, lon1, radius)
            lon_distance = haversine_single_grid(lat1, lon1, lat1, lon2, radius)

            # Calculate the area of the rectangular grid cell
            area = lat_distance * lon_distance
            # Store the area in the corresponding grid cell
            areas[i, j] = area

    return areas

def calculate_grid_areas(latitudes, longitudes, radius=6371):
    """
    Vectorized calculation of grid cell areas in km^2 using the Haversine formula.
    latitudes, longitudes: 2D arrays of shape (n_lat, n_lon)
    Returns: 2D array of shape (n_lat-1, n_lon-1)
    """
    # dy: distance between lat edges
    dy = haversine_vectorized(
        latitudes[:-1, :-1], longitudes[:-1, :-1],
        latitudes[1:, :-1], longitudes[:-1, :-1],
        radius
    )

    # dx: distance between lon edges
    dx = haversine_vectorized(
        latitudes[:-1, :-1], longitudes[:-1, :-1],
        latitudes[:-1, 1:], longitudes[:-1, 1:],
        radius
    )

    # area = dx * dy
    areas = dx * dy
    return areas

def area_at_threshold(threshold, ds, relerr_lim = 0.5):

    # Detta kör vi sen, när vi kollat lite mer.
    if threshold != 0:
        mask_below_threshold = ds['Oxygen'] <= threshold
    else:
        mask_below_threshold = ds['Oxygen'] <= threshold + 9 #4.5

    ds[f"{threshold}_mask"] = mask_below_threshold
    # Get the first layer in the depth dimension where Oxygen is below the threshold
    # First, find the indices where the condition is True along the depth dimension
    first_layer_indices = mask_below_threshold.argmax(dim="depth")
    # Then, use these indices to extract the corresponding depth values
    ds[f"Min_depth_{threshold}"] = ds['depth'].where(mask_below_threshold).isel(depth=first_layer_indices)

    # Then, use these indices to extract the corresponding depth values - not used?
    ds[f"{threshold}_mask_firstlayer"] = ds[f"{threshold}_mask"].where(mask_below_threshold).isel(depth=first_layer_indices)

    # Step 5: Extract corresponding values for the variable Oxygen_relerr using the mask from step 4
    # change Relerr_per_grid_at_min_{threshold}_depth to Relerr_per_grid_at_Min_depth_{threshold} to follow the same pattern
    # 2D array (only lat,lon, time dimensions left)
    ds[f"Relerr_per_grid_at_min_{threshold}_depth"] = (
        ds["Oxygen_relerr"].where(mask_below_threshold).isel(depth=first_layer_indices)
    )
    # 3D array with NaN at all cells except Min_depth_<threhsold>
    ds[f"{threshold}_relerr_per_depth"] = xr.where(
        (ds["depth"] == ds[f"Min_depth_{threshold}"]),
        ds["Oxygen_relerr"],
        ds[f"Min_depth_{threshold}"] * np.nan,
        keep_attrs=True,
    )
    # sums all cells skipping nan, so the same as "Relerr_per_grid_at_min_{threshold}_depth 
    # but can be used if more than one depths are used, like for total volume
    ds[f"{threshold}_relerr_per_grid"] = ds[f"{threshold}_relerr_per_depth"].sum(dim=['depth'], skipna=True)
    
    # Calculate total area by summing
    ds[f"Area_{threshold}"] = xr.where(
        (ds[f"Min_depth_{threshold}"] >= 0),
        ds["grid_area"],
        ds[f"Min_depth_{threshold}"] * np.nan,
        keep_attrs=True,
    ).sum(dim=["lat", "lon"], skipna=True)

    # Calculate area with relative error above the limit 
    ds[f"Relerr_area_{threshold}"] = xr.where(
        (ds[f"Relerr_per_grid_at_min_{threshold}_depth"] >= relerr_lim),
        ds["grid_area"],
        ds[f"Min_depth_{threshold}"] * np.nan,
        keep_attrs=True,
    ).sum(dim=["lat", "lon"], skipna=True)

    # add the calculated areas to the df
    # df[f"Area_{threshold}_km2"] = ds[f"Area_{threshold}"].round(3)
    # df[f"Relerr_area_{threshold}_km2"] = ds[f"Relerr_area_{threshold}"].round(3)

def build_basin_masks(ds: xr.Dataset, basins: gpd.GeoDataFrame, basin_ids: list, basin_id_col: str="basin_id"):
    """
    Create boolean basin masks on a regular lat/lon grid using vectorized shapely operations.
    Returns a DataArray of shape (basin, lat, lon) with pd.Index for basin dimension.
    """

    lon2d, lat2d = np.meshgrid(ds["lon"].values, ds["lat"].values)

    masks = []
    basin_names = []

    for basin_id in basin_ids:
        subset = basins.loc[basins[basin_id_col] == basin_id, "geometry"]
        if subset.empty:
            raise ValueError(f"Basin ID {basin_id} not found in basins GeoDataFrame")
        geom = subset.iloc[0]

        # Vectorized contains check
        mask = shapely.contains_xy(geom, lon2d, lat2d)
        masks.append(mask)

        basin_names.append(
            basins.loc[basins[basin_id_col] == basin_id, "basin_name"].iloc[0]
        )

    # Stack into xarray DataArray
    mask_da = xr.DataArray(
        np.stack(masks),
        dims=("basin", "lat", "lon"),
        coords={
            "basin": pd.Index(basin_ids, name="basin"),
            "basin_name": ("basin", basin_names),
            "lat": ds["lat"],
            "lon": ds["lon"],
        },
        name="basin_mask"
    )

    return mask_da

def basin_area_at_threshold(ds, basin_mask, threshold, relerr_lim=0.5):
    area = xr.where(
        basin_mask & (ds[f"Min_depth_{threshold}"] >= 0),
        ds["grid_area"],
        0.0
    ).sum(dim=["lat", "lon"], skipna=True)

    error_area = xr.where(
        basin_mask & (ds[f"Relerr_per_grid_at_min_{threshold}_depth"] >= relerr_lim),
        ds["grid_area"],
        0.0
    ).sum(dim=["lat", "lon"], skipna=True)

    return area, error_area

def calculate_areas(results_dir, threshold_list, save_area_data=False):
    print("calculating areas...")
    #List for results to txt
    raw_files = Path(results_dir) / "processed"
    df_list = []
    df_BG_list = []
    for netcdf_filepath in list(raw_files.glob('*.nc')):
        netcdf_filename = netcdf_filepath.name
        ds = xr.open_dataset(netcdf_filepath)
        ds = ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=False)
        ds = ds.rio.write_crs("EPSG:4326", inplace=False)
        print(netcdf_filepath)
        print(ds.attrs)
        season = ds.attrs['season']
        start_year = ds.attrs['start year']
        end_year = ds.attrs['end year']
        ### Calculate area of all grid cells
        # Get the latitude and longitude coordinates
        lon, lat = np.meshgrid(ds['lon'], ds['lat'])

        grid_areas = calculate_grid_areas(latitudes=lat, longitudes=lon)
        # assign the calculated areas to the dataset and return the updated dataset
        assert grid_areas.shape == (len(ds.lat)-1, len(ds.lon)-1)
        ds = ds.assign(grid_area=(('lat', 'lon'), np.pad(grid_areas, ((0,1),(0,1)), "edge")))
        ### Find anox gridcells and set them to 1 all others to NaN
        #To do: Ändra så att vi kan skicka in en lista på gränsvärden.
        ds["basin_mask"] = build_basin_masks(
                                            ds,
                                            BASINS,
                                            BASINIDS
                                        )
        basin_names = [BASINS.set_index("basin_id").loc[b, "basin_name"] for b in BASINIDS]
        n_basins = len(BASINIDS)
        n_times = ds.dims.get("time", 1)  # default to 1 if no time dim

        # create array of shape (n_basins, n_times)
        
        for threshold in threshold_list:
            area_at_threshold(threshold, ds)
            area_vals_full = np.zeros((n_basins, n_times))
            error_area_vals_full = np.zeros((n_basins, n_times))

            for i, basin_id in enumerate(BASINIDS):
                basin_mask = ds["basin_mask"].sel(basin=basin_id)
                for t_idx in range(n_times):
                    # slice ds along time if needed
                    ds_time_slice = ds.isel(time=t_idx) if "time" in ds.dims else ds
                    area, error_area = basin_area_at_threshold(ds_time_slice, basin_mask, threshold)
                    area_vals_full[i, t_idx] = float(area.values)
                    error_area_vals_full[i, t_idx] = float(error_area.values)
                        
            # Now create the DataArrays with guaranteed full basin dimension
            ds[f"Area_{threshold}_basin"] = xr.DataArray(
                area_vals_full,                       # shape = (n_basins,)
                dims=("basin", "time"),                   # dimension along the array
                coords={
                    "basin": np.arange(len(BASINIDS)), 
                    "basin_id": ("basin", BASINIDS),       # strings like "SEA-007"
                    "basin_name": ("basin", basin_names),  # human-readable
                    "time": ds.time
                },
                name=f"Area_{threshold}_basin"
            ).set_index(basin="basin_id")

            ds[f"Relerr_area_{threshold}_basin"] = xr.DataArray(
                error_area_vals_full,
                dims=("basin", "time"),
                coords={
                    "basin": np.arange(len(BASINIDS)), 
                    "basin_id": ("basin", BASINIDS),
                    "basin_name": ("basin", basin_names),
                    "time": ds.time
                },
                name=f"Relerr_area_{threshold}_basin"
            ).set_index(basin="basin_id")
        print(ds["Relerr_area_180_basin"].coords)
        ds[f"Area_{180}_basin"].sel(basin="SEA-007")
        # save the updated dataset
        print(f"writing to {results_dir}/processed_test/{netcdf_filename}...")
        ds.to_netcdf(f'{results_dir}/processed_test/{netcdf_filename}') # rewrite to netcdf
        if "Background" in netcdf_filename:
            df_BG = extract_area_data(ds, threshold_list, basin_ids=BASINIDS_SUBSET)
            df_BG_list.append(df_BG)
        else:            
            df = extract_area_data(ds, threshold_list, basin_ids=BASINIDS_SUBSET)
            df_list.append(df)
    pd.concat(df_list).to_csv(f"{results_dir}/test_area_data_{start_year}_{end_year}.txt", sep="\t", index=False)    
    pd.concat(df_BG_list).to_csv(f"{results_dir}/test_BG_area_data_{start_year}_{end_year}.txt", sep="\t", index=False)    


def extract_area_data(ds, threshold_list, basin_ids=None):
    rows = []
    season = ds.attrs['season']

    for t_idx, t_val in enumerate(ds.time.values):
        # Extract year from timestamp
        year = pd.to_datetime(t_val).year
        row = {"year": year}
        row["season": season]

        # Total areas
        for threshold in threshold_list:
            row[f"Area_{threshold}_km2"] = float(ds[f"Area_{threshold}"].isel(time=t_idx).values).round(3)
            row[f"Relerr_area_{threshold}_km2"] = float(ds[f"Relerr_area_{threshold}"].isel(time=t_idx).values).round(3)

        # Basin areas
        if basin_ids:
            for basin_id_val in basin_ids:
                for threshold in threshold_list:
                    row[f"Area_{threshold}_{basin_id_val}_km2"] = float(
                        ds[f"Area_{threshold}_basin"].sel(basin=basin_id_val).isel(time=t_idx).values
                    ).round(3)
                    row[f"Relerr_area_{threshold}_{basin_id_val}_km2"] = float(
                        ds[f"Relerr_area_{threshold}_basin"].sel(basin=basin_id_val).isel(time=t_idx).values
                    ).round(3)

        rows.append(row)

    df = pd.DataFrame(rows)
    return df

if __name__ == "__main__":
    # Result directory
    results_dir = "./results/Baltic_Proper/20260114_1653/"
   # Thresholds to analyse in µmol/l oxygen (0, 2, 4 ml/l)
    threshold_list = [0, 90, 180]
    calculate_areas(results_dir, threshold_list)