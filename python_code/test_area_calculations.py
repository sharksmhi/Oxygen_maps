import xarray as xr
import numpy as np
import pandas as pd
import math

import math

lat1 = 54.
lon1 = 10.
lat2 = 54.05
lon2 = 10.05

# controll measure in Qgis of the above coordinates show the area is ~18 km2. 0.05 degree resoiution

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

def calculate_grid_area(lat1, lon1, lat2, lon2, radius=6371):
    """
    Calculate the approximate area of a grid on a sphere
    using latitude and longitude coordinates.
    """
    # Calculate the distance between the latitude and longitude points
    lat_distance = haversine(lat1, lon1, lat2, lon1, radius)
    lon_distance = haversine(lat1, lon1, lat1, lon2, radius)

    # Convert the distances from kilometers to meters
    lat_distance_meters = lat_distance * 1000
    lon_distance_meters = lon_distance * 1000

    # Calculate the area of the rectangular grid cell
    area = lat_distance_meters * lon_distance_meters

    return area

print(calculate_grid_area(lat1, lon1, lat2, lon2,))

def calculate_grid_areas(latitudes, longitudes):
    """
    Calculate the approximate area of a grid on a sphere
    using latitude and longitude coordinates.
    """
    radius = 6371  # Radius of the Earth in kilometers

    # Initialize an array to store the grid cell areas
    areas = np.zeros_like(latitudes)

    # TODO: add row and column to get areas of last grid cells in the grid
    # print(latitudes[-1])
    # # Get the values of the last row
    # last_row = latitudes[-1, :]
    # print(len(last_row))

    # # Calculate the difference between the last two rows
    # new_row = latitudes[-1] + latitudes[-1] - latitudes[-2]
    # print(new_row, len(new_row))

    # # Create a new row andusing the calculated differences
    # # new_row = last_row + diff_row

    # # Add the new row to the original array
    # new_array = np.concatenate((latitudes, [new_row]), axis=0)
    # print(repr(new_array))
    # # Get the values of the last column
    # last_column = new_array[:, -1]
    # print(len(last_column))

    # # Calculate the difference between the last two columns
    # diff_column = last_column[-1] - last_column[-2]

    # # Calculate the difference between the last two rows
    # new_column = new_array[:,-1] + new_array[: ,-1] - new_array[:, -2]
    # print(new_column, len(new_column))

    # # Add the new row to the original array
    # final_array = np.concatenate((new_array, [new_column]), axis=1)
    # print(final_array[-1])

    # exit()

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
            lat_distance_meters = lat_distance * 1000
            lon_distance_meters = lon_distance * 1000

            # Calculate the area of the rectangular grid cell
            area = lat_distance_meters * lon_distance_meters
            # Store the area in the corresponding grid cell
            areas[i, j] = area

    return areas


from geopy.distance import geodesic

def calculate_grid_area_geopy(lat1, lon1, lat2, lon2):
    # Calculate the distance between two latitude and longitude coordinates
    distance_x1 = geodesic((lat1, lon1), (lat1, lon2)).kilometers
    distance_y1 = geodesic((lat1, lon1), (lat2, lon1)).kilometers
    distance_x2 = geodesic((lat2, lon1), (lat2, lon2)).kilometers
    distance_y2 = geodesic((lat1, lon2), (lat2, lon2)).kilometers

    # lat2, lon1    lat2, lon2
    # lat1, lon1    lat1, lon2
    print(f"{distance_x1=}\n{distance_x2=}\n{distance_y1=}\n{distance_y2=}")

    # Assuming the grid is square, calculate the area
    area = distance_y1*distance_x1

    return area

print(calculate_grid_area_geopy(lat1, lon1, lat2, lon2))

# Read the bathymetry file using Xarray
bath_file = xr.open_dataset("../bathymetry/gebco_30sec_8.nc")
print(bath_file)
# Extract the required variables
b = bath_file["bat"]

## or....
### open netcdf file ###
ds = xr.open_dataset('../resultat/nc/modified_Oxygen_1960-2018_Autumn_1_50000.0_gebco_30sec_4.nc')
print(ds)

# Get the latitude and longitude coordinates
# varför med meshgrid?
lon, lat = np.meshgrid(ds['lon'], ds['lat'])
print(lon[0:5])
print(lat[0:5])

grid_areas = calculate_grid_areas(latitudes=lat, longitudes=lon)
print(f"{grid_areas[0:5,0:5]=}")

# da = xr.DataArray(
#     data=grid_areas,
#     dims=["x", "y"],
#     coords=dict(
#         lon=(["lon"], lon),
#         lat=(["lat"], lat),
#     ),
#     attrs=dict(
#         description="Grid area",
#         units="km2",
#     )
#     )
# print(da)


print(grid_areas.shape)

# assign the calculated areas to the dataset and return the updated dataset
ds2 = ds.assign(grid_area = (('lat', 'lon'), grid_areas))

# Print the updated dataset
print(ds)
print(ds2)

# save the updated dataset
location = "." # or other location like havgem
ds2.to_netcdf(f'{location}./resultat/nc/modified_modified_Oxygen_1960-2018_Autumn_1_50000.0_gebco_30sec_4.nc') # rewrite to netcdf