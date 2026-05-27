using JSON
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
end