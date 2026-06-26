import xarray as xr
import numpy as np
import glob
import os

tile_dir = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/data/bathymetry/tiles/"
out_file = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/data/bathymetry/emodnet_bathymetry_merged.nc"

files = sorted(glob.glob(os.path.join(tile_dir, "*.nc")))

print("Antal tiles:", len(files))
print(files)

varname = "elevation"  # byt om variabeln heter något annat

merged = None

for f in files:
    print("Läser:", f)
    ds = xr.open_dataset(f)
    da = ds[varname]

    rename_dict = {}

    if "longitude" in da.dims:
        rename_dict["longitude"] = "lon"
    if "latitude" in da.dims:
        rename_dict["latitude"] = "lat"
    if "longitude" in da.coords:
        rename_dict["longitude"] = "lon"
    if "latitude" in da.coords:
        rename_dict["latitude"] = "lat"

    if rename_dict:
        da = da.rename(rename_dict)

    # Ta bara 2D lon/lat om extra dimensioner finns
    extra_dims = [d for d in da.dims if d not in ["lat", "lon"]]
    for d in extra_dims:
        da = da.isel({d: 0})

    da = da.load()

    # Sortera
    da = da.sortby("lat")
    da = da.sortby("lon")

    # Ta bort dubletter inom varje tile
    _, lat_idx = np.unique(da["lat"].values, return_index=True)
    _, lon_idx = np.unique(da["lon"].values, return_index=True)

    da = da.isel(
        lat=np.sort(lat_idx),
        lon=np.sort(lon_idx)
    )

    print("  dims:", da.dims)
    print("  shape:", da.shape)
    print("  lon:", float(da.lon.min()), float(da.lon.max()))
    print("  lat:", float(da.lat.min()), float(da.lat.max()))

    if merged is None:
        merged = da
    else:
        merged = merged.combine_first(da)

# Slutlig sortering
merged = merged.sortby("lat").sortby("lon")

# Begränsa till lon 9-31 och lat 53-66.5
merged = merged.sel(
    lon=slice(9, 31),
    lat=slice(53, 66.5)
)

# Maska bort övre vänstra hörnet:
# lon 9-16.5 och lat 61.5-66.5
corner_mask = (
    (merged["lon"] >= 9)
    & (merged["lon"] <= 16.5)
    & (merged["lat"] >= 61.5)
    & (merged["lat"] <= 66.5)
)

merged = merged.where(~corner_mask, other=np.nan)

# Maska bort land
merged = merged.where(merged <= 0)

# Spara
ds_out = merged.to_dataset(name="bat")
ds_out.to_netcdf(out_file)

print("Sparade:", out_file)
print(ds_out)