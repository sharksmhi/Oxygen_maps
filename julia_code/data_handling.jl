#=
syrekartor_data_proc:
- Julia version: 
- Author: a001109
- Date: 2023-08-30

Program to handle O2 data before DIVAnd.
Emodnet, SHARKweb, ICES and SYKE data
-Load data
-Removes duplicates, with the help of DIVA-funktions.
-Handles both BTL and CTD-data. BTL-data is prioritiz
-Merge to one dataset
# ### Lena Viktorsson & Martin Hansson
=#
#using Pkg
#Pkg.add("PyCall")

# #### Add necessary packages
using DIVAnd
using PyPlot
using NCDatasets
using Dates
using Statistics
using DelimitedFiles
using DataStructures
using Printf
using Missings
using DataFrames
using CSV
using PyCall

function remove_matching_rows(data::DataFrame, bad_data::DataFrame, columns::Vector{Symbol})
    """
    Removes rows from `data` DataFrame where the specified `columns` match in the `to_remove` DataFrame.

    Args:
        data (DataFrame): Original dataset
        to_remove (DataFrame): Dataset describing rows to be removed
        columns (Vector{Symbol}): List of column names to use for comparison

    Returns:
        DataFrame: Filtered DataFrame with rows removed
    """
    # Perform a left join with indicator column
    merged = leftjoin(data, bad_data, on=columns, makeunique=true, indicator=:merge_status)

    # Filter out rows where `merge_status` indicates a match
    filtered_data = merged[ismissing.(merged.merge_status), :]

    # Drop the indicator column
    select!(filtered_data, Not(:merge_status))
    return filtered_data
end


# ## Configuration
# * Define variabel and set horizontal, vertical and temporal resolutions.
### Oxygen ###
varname = "Oxygen"
savevar = "O2"
sdnp35 = "SDN:P35::EPC00002"
sdnp02 = "SDN:P02::DOXY"
doi = "10.6092/B5BB9EA1-4BDA-48A1-92B6-03645AC12FAE"
unit = "umol/l";

# ## Where to save the result. Create path if not there.
# File name based on the variable (but all spaces are replaced by _) _varlenz
# data-files
#location = "//winfs-proj/proj/havgem/DIVA/syrekartor/"
location = "C:/Work/DIVAnd/Oxygen_maps/"
outputdir = joinpath(location,"data/");
if !isdir(outputdir)
    mkpath(outputdir)
end
# Figures
figdir = joinpath(location,"resultat/figures/");

if !isdir(figdir)
    mkpath(figdir)
    #mkpath(joinpath(figdir, "test"))
end

# ## Load data big files
# emodnet BTL + CTD
@show("Loading SHARK BTL/CTD...")
fname_shark = joinpath(location, "data/all_baltic/sharkweb_btlctd_02_241107.txt")
@time obsval_shark,obslon_shark,obslat_shark,obsdepth_shark,obstime_shark,obsid_shark = loadbigfile(fname_shark);

@show("Loading EMODNET BTL...")
datafile_emod_btl = joinpath(location, "data/EMODNET_2024/BTL_data_from_EMODnet_Eutrophication_European_2024_unrestricted.txt")
@time obsval_emod_btl,obslon_emod_btl,obslat_emod_btl,obsdepth_emod_btl,obstime_emod_btl,obsid_emod_btl = ODVspreadsheet.load(Float64,[datafile_emod_btl],
                           ["Water body dissolved oxygen concentration"]; nametype = :localname );
# obsid_emod_btl = LOCAL_CDI_ID + EDMO_CODEz
obsid_emod_btl = string.("emod_btl-",obsid_emod_btl)

@show("Loading EMODNET CTD...")
datafile_emod_ctd = joinpath(location, "data/EMODNET_2024/CTD_data_from_EMODnet_Eutrophication_European_2024_unrestricted.txt")
@time obsval_emod_ctd,obslon_emod_ctd,obslat_emod_ctd,obsdepth_emod_ctd,obstime_emod_ctd,obsid_emod_ctd = ODVspreadsheet.load(Float64,[datafile_emod_ctd],
                           ["Water body dissolved oxygen concentration"]; nametype = :localname );
obsid_emod_ctd = string.("emod_ctd-",obsid_emod_ctd)
@show(obsid_emod_ctd[1])

@show("Loading ICES...")
datafile_ices_btlctd = joinpath(location, "data/all_baltic/ICES_btl_lowres_ctd_02_241107.txt")
@time obsval_ices_btlctd,obslon_ices_btlctd,obslat_ices_btlctd,obsdepth_ices_btlctd,obstime_ices_btlctd,obsid_ices_btlctd = loadbigfile(datafile_ices_btlctd);
obsid_ices_btlctd = string.("ICES-", obsid_ices_btlctd)
@show(obsid_ices_btlctd[1])

@show("Loading SYKE data...")
datafile_syke_btlctd = joinpath(location, "data/all_baltic/syke_data_no_header_241107.txt")
@time obsval_syke_btlctd,obslon_syke_btlctd,obslat_syke_btlctd,obsdepth_syke_btlctd,obstime_syke_btlctd,obsid_syke_btlctd = loadbigfile(datafile_syke_btlctd);
@show(obsid_syke_btlctd[1])

@show("Loading BAD-data file...data to be removed")
datafile_bad_data = joinpath(location, "data/EMODNET_2024/bad_data.txt")
df_bad_data = CSV.read(datafile_bad_data, DataFrame)

@show(length(obsval_emod_btl));
@show(length(obsval_emod_ctd));
@show(length(obsval_shark));
@show(length(obsval_ices_btlctd));
@show(length(obsval_syke_btlctd));

# ## EMODNET: Remove low res CTD when BTL is available.
# ## Criteria (can be adapted according to the application):
# Horizontal distance: 0.01 degree (about 1km)
xy_dist = 0.01
# Vertical separation: 0.01 m depth
depth_dist= 0.01
#Time separation: 60 minute.
time_sep = 60
#obsval difference: 90 µmol/l correspond to ~2 ml/l
obsval_diff = 90

@time dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon_emod_btl,obslat_emod_btl,obsdepth_emod_btl,obstime_emod_btl), obsval_emod_btl,
    (obslon_emod_ctd,obslat_emod_ctd,obsdepth_emod_ctd,obstime_emod_ctd),obsval_emod_ctd,
    (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);

# ## Find the indices of the possible duplicates:
index = findall(.!isempty.(dupl));
ndupl = length(index);
pcdupl = round(ndupl / length(obslon_emod_ctd) * 100; digits=2);
@info("Number of possible duplicates emodnet BTL/CTD: $ndupl")
@info("Percentage of duplicates: $pcdupl%")
# ## If you decide to combine the 2 (or more) datasets:
newpoints = isempty.(dupl);
@info("Number of new points: $(sum(newpoints)))")

obslon_emod = [obslon_emod_btl; obslon_emod_ctd[newpoints]];
obslat_emod = [obslat_emod_btl; obslat_emod_ctd[newpoints]];
obsdepth_emod = [obsdepth_emod_btl; obsdepth_emod_ctd[newpoints]];
obstime_emod = [obstime_emod_btl; obstime_emod_ctd[newpoints]];
obsval_emod = [obsval_emod_btl; obsval_emod_ctd[newpoints]];
obsid_emod = [obsid_emod_btl; obsid_emod_ctd[newpoints]];

# ## Remove SHARK_data when EMODnet_data is available.
# Remove true duplicates, hence when exactly the same data is found in both datasets.
# ## Criteria (can be adapted according to the application):
# Horizontal distance: 0.01 degree (about 1km)
xy_dist = 0.01
# Vertical separation: 0.01 m depth
depth_dist= 0.01
#Time separation: 10 minute.
time_sep = 10
#obsval difference: 1 µmol/l correspond to ~pyttelite ml/l
obsval_diff = 1

@time dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon_emod,obslat_emod,obsdepth_emod,obstime_emod), obsval_emod,
    (obslon_shark,obslat_shark,obsdepth_shark,obstime_shark),obsval_shark,
    (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);

# ## Find the indices of the possible duplicates:
index = findall(.!isempty.(dupl));
ndupl = length(index);
pcdupl = round(ndupl / length(obslon_shark) * 100; digits=2);
@info("Number of possible duplicates in emodnet/SHARKweb: $ndupl")
@info("Percentage of duplicates: $pcdupl%")
# ## If you decide to combine the 2 (or more) datasets:
newpoints_shark = isempty.(dupl);
@info("Number of new points: $(sum(newpoints)))")

obslon_emodshark = [obslon_emod; obslon_shark[newpoints_shark]];
obslat_emodshark = [obslat_emod; obslat_shark[newpoints_shark]];
obsdepth_emodshark = [obsdepth_emod; obsdepth_shark[newpoints_shark]];
obstime_emodshark = [obstime_emod; obstime_shark[newpoints_shark]];
obsval_emodshark = [obsval_emod; obsval_shark[newpoints_shark]];
obsid_emodshark = [obsid_emod; obsid_shark[newpoints_shark]];

# ## Remove ICES data when EMODnet_SHARK is available.
# Remove true duplicates, hence when exactly the same data is found in both datasets.
# ## Criteria (can be adapted according to the application):
# Horizontal distance: 0.01 degree (about 1km)
xy_dist = 0.05
# Vertical separation: 0.01 m depth
depth_dist= 1
#Time separation: 10 minute.
time_sep = 60
#obsval difference: 1 µmol/l correspond to ~pyttelite ml/l
obsval_diff = 1

@time dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon_emodshark,obslat_emodshark,obsdepth_emodshark,obstime_emodshark), obsval_emodshark,
    (obslon_ices_btlctd,obslat_ices_btlctd,obsdepth_ices_btlctd,obstime_ices_btlctd),obsval_ices_btlctd,
    (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);

# ## Find the indices of the possible duplicates:
index = findall(.!isempty.(dupl));
ndupl = length(index);
pcdupl = round(ndupl / length(obslon_ices_btlctd) * 100; digits=2);
@info("Number of possible duplicates in emodnetSHARKweb/ICES: $ndupl")
@info("Percentage of duplicates: $pcdupl%")
# ## If you decide to combine the 2 (or more) datasets:
newpoints_ICES = isempty.(dupl);
@info("Number of new points: $(sum(newpoints)))")

obslon_emodsharkices = [obslon_emodshark; obslon_ices_btlctd[newpoints_ICES]];
obslat_emodsharkices = [obslat_emodshark; obslat_ices_btlctd[newpoints_ICES]];
obsdepth_emodsharkices = [obsdepth_emodshark; obsdepth_ices_btlctd[newpoints_ICES]];
obstime_emodsharkices = [obstime_emodshark; obstime_ices_btlctd[newpoints_ICES]];
obsval_emodsharkices = [obsval_emodshark; obsval_ices_btlctd[newpoints_ICES]];
obsid_emodsharkices = [obsid_emodshark; obsid_ices_btlctd[newpoints_ICES]];

# ## Remove SYKE data when EMODnet_SHARK_ICES data is available.
# Remove true duplicates, hence when exactly the same data is found in both datasets.
# ## Criteria (can be adapted according to the application):
# Horizontal distance: 0.01 degree (about 1km)
xy_dist = 0.05
# Vertical separation: 0.01 m depth
depth_dist= 1
#Time separation: 10 minute.
time_sep = 60
#obsval difference: 1 µmol/l correspond to ~pyttelite ml/l
obsval_diff = 1

@time dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon_emodsharkices,obslat_emodsharkices,obsdepth_emodsharkices,obstime_emodsharkices), obsval_emodsharkices,
    (obslon_syke_btlctd,obslat_syke_btlctd,obsdepth_syke_btlctd,obstime_syke_btlctd),obsval_syke_btlctd,
    (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);

# ## Find the indices of the possible duplicates:
index = findall(.!isempty.(dupl));
ndupl = length(index);
pcdupl = round(ndupl / length(obslon_ices_btlctd) * 100; digits=2);
@info("Number of possible duplicates in emodnetSHARKwebICES/SYKE: $ndupl")
@info("Percentage of duplicates: $pcdupl%")
# ## If you decide to combine the 2 (or more) datasets:
newpoints_SYKE = isempty.(dupl);
@info("Number of new points: $(sum(newpoints)))")

obslon = [obslon_emodsharkices; obslon_syke_btlctd[newpoints_SYKE]];
obslat = [obslat_emodsharkices; obslat_syke_btlctd[newpoints_SYKE]];
obsdepth = [obsdepth_emodsharkices; obsdepth_syke_btlctd[newpoints_SYKE]];
obstime = [obstime_emodsharkices; obstime_syke_btlctd[newpoints_SYKE]];
obsval = [obsval_emodsharkices; obsval_syke_btlctd[newpoints_SYKE]];
obsid = [obsid_emodsharkices; obsid_syke_btlctd[newpoints_SYKE]];


# ## Do a range check and remove outliers µomol/l
min_value = -500.0
max_value = 600.0

# Perform the range check
in_range = (obsval .>= min_value) .& (obsval .<= max_value)

# Extract values within the range (valid data)
valid_obsval = obsval[in_range]

# Extract values outside the range (invalid data)
out_of_range = .!in_range
invalid_obsval = obsval[out_of_range]

# Store the invalid values in a separate file
invalid_values_file = joinpath(location, "data/all_baltic/invalid_obsval.txt")
open(invalid_values_file, "w") do file
    for value in invalid_obsval
        write(file, "$value\n")
    end
end

# Save the valid data back to obsval_shark
obsval = valid_obsval
obslon = obslon[in_range];
obslat = obslat[in_range];
obsdepth = obsdepth[in_range];
obstime = obstime[in_range];
obsid = obsid[in_range];

@show("sets all neg o2 to 0.44661 µmol/l")
obsval[obsval .<= 0] .= 0.44661

# Output some statistics
@show sum(in_range)
@show sum(out_of_range)

# Now obsval_shark contains only the valid values


# ## Hantera negativa värden i obsdepth_shark och gör dem positiva
obsdepth_shark = abs.(obsdepth_shark)


# ## Create a plot showing the additional data points:
figure("Additional-Data")
aspectratio = 1/cos(mean(54:0.05:61) * pi/180)
ax = subplot(1,1,1)
ax.tick_params("both",labelsize=6)
ylim(53, 66);
xlim(8, 31);
#contourf(bx, by, permutedims(Float64.(mask_edit[:,:,1]),[2,1]),
#    levels=[-1e5,0],cmap="binary");
plot(obslon_emod, obslat_emod, "bo", markersize=.2, label="Emodnet")
plot(obslon_shark[newpoints_shark], obslat_shark[newpoints_shark], "go",
   markersize=.2, label="Additional data\nfrom SHARKweb")
plot(obslon_ices_btlctd[newpoints_ICES], obslat_ices_btlctd[newpoints_ICES], "ro", mfc="none",
   markersize=.2, label="Additional data\nfrom ICES BTL lowres CTD")
plot(obslon_syke_btlctd[newpoints_SYKE], obslat_syke_btlctd[newpoints_SYKE], "yo", mfc="none",
   markersize=.2, label="Additional data\nfrom SYKE")


legend(loc=4, fontsize=4)
gca().set_aspect(aspectratio)
figname = "additional_data.png"
@show joinpath(figdir,"$(figname)")
PyPlot.savefig(joinpath(figdir,"$(figname)"), dpi=300);
PyPlot.close_figs()

df  = DataFrame(obslon=obslon,obslat=obslat,obsval=obsval,obsdepth=obsdepth,obsdepth1=obsdepth,obsdepth2=obsdepth,obsdepth3=obsdepth,obsdepth4=obsdepth,obsdepth5=obsdepth,obstime=obstime,obsid=obsid)

@show("Removing bad data...from det bad-data file..")
columns_to_match = [:obslon, :obslat, :obsdepth, :obsval, :obstime, :obsid]
df = remove_matching_rows(df,df_bad_data,columns_to_match)

filename = "EMODNET_SHARK_ICES_SYKE_241204"
CSV.write(joinpath(outputdir, "$(filename).txt"), df, delim="\t", writeheader=false)
#DIVAnd.saveobs(joinpath(outputdir, "$(filename).nc"),varname, obsval, (obslon,obslat,obsdepth,obstime),obsid)

# Formaterar det till önskat format för Johannas/fredriks QC....
#formatted_times = [Dates.format(time, "yyyymmddHHMM") for time in obstime]
#obsval = obsval.* (22.414./1000) #Tillbaka till ml/l för QC?
#QC = ones(Int, length(formatted_times))
#NN = fill(NaN, length(formatted_times))
##                   tid,                lat,    lat_qc, lon,          lon_qc,    djup,       djup_qc,      tryck,     tryck_qc, temp, temp_qc, salt, salt_qc,chl,  chl_qc, turb, turb_qc, syre,       syre_qc, PAr, PAR_qc, ljudhast, ljudhast_qc, ID (Skall eg inte va med)
#               40014	8002	88002	8003	88003	8178	88178	8168	88168	8179	88179	8181	88181	8063	88063	8174	88174	8191	88191	8175	88175	8187	88187
#df_QC  = DataFrame(formatted_times=formatted_times, obslat=obslat, QC1=QC, obslon=obslon, QC2=QC, obsdepth=obsdepth, QC3=QC, obsdepth_tryck=obsdepth, QC4=QC, NN5=NN, QC5=QC, NN6=NN, QC6=QC, NN7=NN, QC7=QC, NN8=NN, QC8=QC, obsval=obsval, QC9=QC, NN10=NN, QC10=QC, NN11=NN, QC11=QC, obsid=obsid)
#filename = "EMODNET_SHARK_ICES_SYKE_241107_QC"
#CSV.write(joinpath(outputdir, "$(filename).txt"), df_QC, delim="\t", writeheader=false)