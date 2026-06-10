#= using JSON
using CSV
using DataFrames

cv_settings = JSON.parse(read("cv_settings.json", String))
@show cv_settings

stations = CSV.read(
    "standard_stations.txt",
    DataFrame;
    delim = '\t'
)
@show stations
tol = cv_settings["standard_stations"]["match_tolerance_deg"]
@show tol
for station in eachrow(stations)
    station_mask =
        (abs.(11 .- station.lon) .< tol) .&
        (abs.(57 .- station.lat) .< tol)
end =#

using NCDatasets
year = "2015"
season = "Autumn"
results_dir = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results/Baltic_Proper/20260601_2334/full/DIVArun"
results_dir = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results_lena_temp/Baltic_Proper/20260408_1932_high_res_1960_2025/DIVArun"

nc_filename = "Oxygen_$(year)_$(season)_0.2_50000.0_0.025_5.0_2.0_Baltic_Proper_varcorrlenz.nc"
nc_filename = "Oxygen_2015_Autumn_0.2_50000.0_0.025_5.0_2.0_bat_elevation_Baltic_Sea_masked_varcorrlenz.nc"
nc_filename_res = "Oxygen_$(year)_$(season)_residuals.nc"
nc_filepath = joinpath("$(results_dir)", nc_filename)
nc_filepath_res = joinpath("$(results_dir)", nc_filename_res)

ds = NCDataset(nc_filepath)

println(keys(ds))
println(ds.dim)

close(ds)