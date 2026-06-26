import xarray as xr
from pathlib import Path
import numpy as np

# Input files
filename_kat = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results/Kattegat/20260616_0837/full/processed/Oxygen_2020_Autumn_0.2_25000.0_0.025_5.0_2.0_Kattegat_varcorrlenz.nc"
filename_bp = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results/Baltic_Proper/20260601_2334/full/processed/Oxygen_2020_Autumn_0.2_50000.0_0.025_5.0_2.0_Baltic_Proper_varcorrlenz.nc"

regional_file = Path(filename_kat)       # Kattegat-resultatet
main_file = Path(filename_bp)      # Huvudfilen som ska patchas
outfile = Path("baltic_proper_patched_with_kattegat.nc")

ds_reg = xr.open_dataset(regional_file)
ds_main = xr.open_dataset(main_file)

ds_out = ds_main.copy(deep=True)

vars_to_patch = [
    "Oxygen",
    "Oxygen_L1",
    "Oxygen_L2",
    "Oxygen_relerr",
    "0_mask",
    "90_mask",
    "180_mask",
    "Min_depth_0",
    "Min_depth_90",
    "Min_depth_180",
    "0_mask_firstlayer",
    "90_mask_firstlayer",
    "180_mask_firstlayer",
    "Relerr_per_grid_at_min_0_depth",
    "Relerr_per_grid_at_min_90_depth",
    "Relerr_per_grid_at_min_180_depth",
    "0_relerr_per_grid",
    "90_relerr_per_grid",
    "180_relerr_per_grid",
]

# Mask defining where the Kattegat solution exists
regional_2d_mask = (
    ds_reg["Oxygen"]
    .isel(time=0, depth=0)
    .interp(
        lon=ds_main.lon,
        lat=ds_main.lat,
        method="nearest",
    )
    .notnull()
)

for varname in vars_to_patch:

    if varname not in ds_reg or varname not in ds_main:
        print(f"Skipping {varname}: missing in one dataset")
        continue

    print(f"Patching {varname}")

    reg_var = ds_reg[varname]
    main_var = ds_main[varname]

    interp_coords = {
        "lon": ds_main.lon,
        "lat": ds_main.lat,
    }

    # ------------------------------------------------------------------
    # Variables containing depth dimension
    # ------------------------------------------------------------------
    if "depth" in reg_var.dims and "depth" in main_var.dims:

        common_depths = ds_main.depth.where(
            ds_main.depth.isin(ds_reg.depth),
            drop=True,
        )

        deep_depths = ds_main.depth.where(
            ds_main.depth > ds_reg.depth.max(),
            drop=True,
        )

        # --------------------------------------------------------------
        # Replace common depth levels
        # --------------------------------------------------------------
        if len(common_depths) > 0:

            reg_selected = reg_var.sel(depth=common_depths)

            if reg_selected.dtype == bool:

                reg_interp = (
                    reg_selected.astype(float)
                    .interp(
                        interp_coords,
                        method="nearest",
                    )
                    .fillna(0)
                    .astype(bool)
                )

            else:

                reg_interp = reg_selected.interp(
                    interp_coords,
                    method="linear",
                )

            mask = reg_interp.notnull()

            ds_out[varname].loc[dict(depth=common_depths)] = xr.where(
                mask,
                reg_interp,
                ds_out[varname].sel(depth=common_depths),
            )

        # --------------------------------------------------------------
        # Set deeper layers to NaN inside Kattegat patch area
        # --------------------------------------------------------------
        if len(deep_depths) > 0:

            deep_mask = regional_2d_mask.broadcast_like(
                ds_out[varname].sel(depth=deep_depths)
            )

            if np.issubdtype(ds_out[varname].dtype, np.bool_):

                ds_out[varname].loc[dict(depth=deep_depths)] = xr.where(
                    deep_mask,
                    False,
                    ds_out[varname].sel(depth=deep_depths),
                )

            else:

                ds_out[varname].loc[dict(depth=deep_depths)] = xr.where(
                    deep_mask,
                    np.nan,
                    ds_out[varname].sel(depth=deep_depths),
                )

    # ------------------------------------------------------------------
    # Variables with lat/lon but no depth
    # ------------------------------------------------------------------
    elif "lat" in reg_var.dims and "lon" in reg_var.dims:

        if reg_var.dtype == bool:

            reg_interp = (
                reg_var.astype(float)
                .interp(
                    interp_coords,
                    method="nearest",
                )
                .fillna(0)
                .astype(bool)
            )

        else:

            reg_interp = reg_var.interp(
                interp_coords,
                method="linear",
            )

        mask = reg_interp.notnull()

        ds_out[varname] = xr.where(
            mask,
            reg_interp,
            ds_out[varname],
        )

    else:

        print(f"Skipping {varname}: no lat/lon dimensions")

# ----------------------------------------------------------------------
# Metadata
# ----------------------------------------------------------------------
ds_out.attrs["comment_patch"] = (
    "Kattegat region replaced with regional Kattegat result. "
    "Only regional depth levels were patched; deeper levels within the "
    "regional patch area were removed."
)

# ----------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------
ds_out.to_netcdf(outfile)

print(f"Saved patched file: {outfile}")