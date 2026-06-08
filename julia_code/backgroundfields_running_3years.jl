# ### Oxygen maps
# ### 2023
# ### Lena Viktorsson & Martin Hansson
# ### Based on copy from Karin Wesslander and DIVA workshop notebooks

# calculates background fields for ten year periods, winter, spring, summer, autumn. Saves to netcdf
# Run seperate for each basin...

# #### Add necessary packages
#using Pkg
#Pkg.add("JLD")

using DIVAnd
using PyPlot
using NCDatasets
using Dates
using Statistics
using DelimitedFiles
using DataStructures
using Printf
using Missings
using JLD
using JSON
using DataFrames
using CSV

function rdiag_stats(
    year,
    obsid,
    obslon,
    obslat,
    obsdepth,
    obstime,
    rdiag,
    epsilon_weighted
)

    return DataFrame(
        year = [year],

        nobs = [length(obsid)],
        nstations = [length(unique(obsid))],

        nlon_unique = [length(unique(obslon))],
        nlat_unique = [length(unique(obslat))],
        ndepth_unique = [length(unique(obsdepth))],

        depth_min = [minimum(obsdepth)],
        depth_max = [maximum(obsdepth)],

        ndates = [length(unique(Date.(obstime)))],

        rdiag_min = [minimum(rdiag)],
        rdiag_mean = [mean(rdiag)],
        rdiag_median = [median(rdiag)],
        rdiag_p95 = [quantile(rdiag, 0.95)],
        rdiag_max = [maximum(rdiag)],

        eps_min = [minimum(epsilon_weighted)],
        eps_mean = [mean(epsilon_weighted)],
        eps_max = [maximum(epsilon_weighted)]
    )
end

function get_weights(year)
    if year < 1980
        return (0.1, 0.1, 50.0, 7.0)
    elseif year < 2000
        return (0.1, 0.1, 20.0, 5.0)
    else
        return (0.1, 0.1, 5.0, 2.0)
    end
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
# NC-files
location = "C:/LenaV/code/DIVAnd/Oxygen_maps/"
location = "C:/Work/DIVAnd/Oxygen_maps/"
location = ""

freja_location = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/"
freja_location = "freja/inprut_dur"
freja_location = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/"

on_freja = false
if isdir(freja_location)
    location = freja_location
    on_freja = true
end

outputdir = joinpath(location, "data/");
if !isdir(outputdir)
    @info("directory path for data does not exist $(outputdir)")
end

# ## Load data big files created by program "data_handling"
data_fname = "SHARK_SYKE_IOW_EMODNET_ICES_260325"
#data_fname = "SHARK_SYKE_IOW_EMODNET_ICES_250619"
#data_fname = "mat_file_1960_2024_reordered"
@time obsval,obslon,obslat,obsdepth,obstime,obsid = loadbigfile(joinpath(location, "data/$data_fname.txt"));

weighting_dir = joinpath(outputdir, "weighting", "background_fields", data_fname)
mkpath(weighting_dir)

background_dir = joinpath(outputdir, "background_fields", data_fname)
mkpath(background_dir)

#Bottniska viken
#basin = "Gulf_of_Bothnia"

#Eg Östersjön o Kattegatt
basin ="Baltic_Proper"

#Kattegatt
#basin = "Kattegat"

# Läs in filens innehåll som en sträng
json_content = read("./settings.json", String)

# Parsar strängen som JSON
settings = JSON.parse(json_content)

# Exempel på att hämta specifika värden
lonr_range = settings[basin]["lonr"]
latr_range = settings[basin]["latr"]
dx = Float64.(settings[basin]["dx"])
dy = Float64.(settings[basin]["dy"])
# Skapa intervall (range) i Julia med angivet dx
lonr_min = lonr_range[1]
lonr_max = lonr_range[2]
latr_min = latr_range[1]
latr_max = latr_range[2]
lonr = lonr_min:dx:lonr_max  # 14:0.05:31
latr = latr_min:dy:latr_max

layers = settings[basin]["layers"]
depthr = Float64[layer[1] for layer in layers]
lenz_ = Float64[layer[2] for layer in layers]

#depthr = Float64.(settings[basin]["depthr"])
#lenz_ = Float64.(settings[basin]["lenz_"])
lenf = Float64.(settings[basin]["lenf_background"])
#yearlist_json = settings[basin]["yearlist_background"]
# Konvertera varje par i yearlist till ett intervall (range) i Julia
#year_list = [year[1]:year[2] for year in yearlist_json]
#years = settings[basin]["years"]
threshold_list = settings[basin]["threshold_list"]

# Ange vilket intervall du vill ha på bakgrundfältet
start_year = 1959
end_year   = 2025

# Skapa listan med rullande treårsintervall # [[1959,1960,1961],[], osv...]
year_list = [[y-1, y, y+1] for y in start_year:end_year]  

@show(year_list)

#BACKGROUND
#month_list = [ [11,12,1,2], [3,4,5], [6,7,8], [8,9,10]];  # Seasonal climatology
month_list = [[1,2,3,4,5,6,7,8,9,10,11,12]];  # 3 whole years climatology

#seasons=["Winter","Spring","Summer","Autumn"]
#months=["(Nov-Feb)","(Mar-May)","(June-Aug)","(Aug-Oct)"];
seasons=["All"]
months=["(Jan-Dec)"];

# Time origin for the NetCDF file
timeorigin = DateTime(1900,1,1,0,0,0);
aspect_ratio = 1/cos(mean(latr) * pi/180);

"BATHYMETRY"
bathname = joinpath(location, "data", "bat_elevation_Baltic_Sea_masked.nc")
bath_file_name = split(bathname,"/")[end]
bathisglobal = true;
@time bx,by,b = DIVAnd.extract_bath(bathname,bathisglobal,lonr,latr);

# ## MASK
@time xmask,ymask,mmask = load_mask(bathname,true,lonr,latr,depthr);

label = DIVAnd.floodfill(mmask);
new_mask = (label .== 1);

# #### 4) mask for Skagerrak
grid_bx = [i for i in xmask, j in ymask];
grid_by = [j for i in xmask, j in ymask];
mask_edit = copy(new_mask);
#mask_edit = copy(mmask);                                                           #Maska bort:
sel_mask1 = (grid_by .>= 57.75) .& (grid_bx .<= 12.2);                              #Skagerrak
sel_mask2 = (grid_by .>= 57.4) .& (grid_by .< 57.75) .& (grid_bx .<= 10.4);         #Skagerrak
sel_mask3 = (grid_by .>= 57.) .& (grid_by .< 57.4) .& (grid_bx .<= 10.);            #Skagerrak

## If basin is....exclude other basins.
if basin == "Gulf_of_Bothnia"
    sel_mask4 = (grid_by .<= latr_min)        # keep all except north of min_lat (Åland Sea)
elseif basin == "Baltic_Proper"
    sel_mask4 = (grid_by .>= latr_max)        # keep all except north of max_lat (Åland sea)
elseif basin == "Kattegat"
    sel_mask4 = (grid_bx .>= lonr_max);        # keep all except east of lonr_max (Arkona/Bornholm)
end
new_mask = mask_edit .* .!sel_mask1 .* .!sel_mask2 .* .!sel_mask3 .* .!sel_mask4;

# ## Analysis parameters
sz = (length(lonr), length(latr), length(depthr));
lenx = fill(lenf,sz)   # 200 km
leny = fill(lenf,sz)   # 200 km
lenz =  [lenz_[k] for i = 1:sz[1], j = 1:sz[2], k = 1:sz[3]];
len = (lenx, leny, lenz);

lx = lenf
ly = lenf

epsilon = Float64.(settings["Global"]["epsilon_background"])

w_depth = 5.
w_days = 2.
w_lon = 0.1
w_lat = 0.1

# Settings for DIVAnd-------------------------------------------------------------------------------
error_thresholds = [("L1", 0.3), ("L2", 0.5)];

# One metadata set up per season
metadata=Array{DataStructures.OrderedDict{String,Any}}(undef,length(seasons));
stats_file = joinpath(weighting_dir, "rdiag_statistics_background.txt")

# Loop over the 3 years and make a BG-field for those years to be used for the analysis of the mid year
# BG field for [1960,1961,1962] [All month] to be used for year 1961 in analysis. 
# BG field for [1961,1962,1963] [All month] to be used for year 1962 in analysis. 
for year_list_index in 1:length(year_list)
    window_years = year_list[year_list_index]
    year = year_list[year_list_index][2] # Year to the analysis +/- 1 år
    @show(year, window_years)    # Lägg in subsetting av data här så vi bara använder de tre åren i fönstret
    year_mask = in.(Dates.year.(obstime), Ref(window_years))
    obsval_sub   = obsval[year_mask]
    obslon_sub   = obslon[year_mask]
    obslat_sub   = obslat[year_mask]
    obsdepth_sub = obsdepth[year_mask]
    obstime_sub  = obstime[year_mask]
    obsid_sub    = obsid[year_mask]
    @info "nobs = $(length(obslon_sub))"
    @info "unique lons = $(length(unique(obslon_sub)))"
    @info "unique lats = $(length(unique(obslat_sub)))"
    @info "unique depth = $(length(unique(obsdepth_sub)))"
    # räkna rdiag för treårsperioden och vikta med endast de data
    # w_lon, w_lat, w_depth, w_days = get_weights(year)
    rdiag_jldfile = joinpath(weighting_dir, "$(data_fname)_$(year)_weighted_$(w_depth)_$(w_days).jld")
    if isfile(rdiag_jldfile)
        @load rdiag_jldfile rdiag
        @info "Loading saved rdiag file for $(year) with w_depth = $(w_depth) and w_days = $(w_days)"
    else
        @info "Calculating rdiag for $(year) with w_depth = $(w_depth) and w_days = $(w_days)!"
        @time rdiag = 1 ./ DIVAnd.weight_RtimesOne((obslon_sub, obslat_sub, obsdepth_sub, float.(Dates.dayofyear.(obstime_sub))), (w_lon, w_lat, w_depth, w_days));
        @save rdiag_jldfile rdiag
    end
    # spara rdiag filen, behövs en per år till bg istället för en för hela perioden
    @show maximum(rdiag),mean(rdiag)
    epsilon_weighted = epsilon * rdiag;
    @show minimum(epsilon_weighted), maximum(epsilon_weighted),mean(epsilon_weighted)

    stats_df = rdiag_stats(
        year,
        obsid_sub,
        obslon_sub,
        obslat_sub,
        obsdepth_sub,
        obstime_sub,
        rdiag,
        epsilon_weighted
    )
    if isfile(stats_file)
        CSV.write(stats_file, stats_df; append=true)
    else
        CSV.write(stats_file, stats_df)
    end

    for monthlist_index in 1:length(month_list)
        season = seasons[monthlist_index]
        # File name based on the variable (but all spaces are replaced by _)
        nc_filename = "Background_$(replace(varname,' '=>'_'))_$(year)_$(season)_$(epsilon)_$(lx)_$(dx)_$(w_depth)_$(w_days)_$(basin).nc"
        nc_filepath = joinpath(background_dir, nc_filename)

        if isfile(nc_filepath)
            @info "Backgroundfield already exists. Skip to next year"
            continue
            # rm(nc_filepath) # delete the previous analysis
        end
        
        @info("Creating metadata dicitonary for the season $(season)")
        metadata_season = OrderedDict(
            # set attributes for DVIA run from our settings
            "threshold_list" => "$threshold_list",
            "season" => season,
            "epsilon" => "$epsilon",
            "horizontal correlation length m" => "$lx",
            "start year" => string(year_list[year_list_index][1]), 
            "end year" => string(year_list[year_list_index][3]),
            # Name of the project (SeaDataCloud, SeaDataNet, EMODNET-chemistry, ...)
            "project" => "EMODNET-chemistry",
            # URN code for the institution EDMO registry,
            # e.g. SDN:EDMO::1579
            "institution_urn" => "SDN:EDMO::545",
            # Production group
            "production" => "SMHI",
            # Name and emails from authorsseasons
            "Author_e-mail" => ["Martin Hansson <martin.hansson@smhi.se>"],
            # Source of the observation
            "source" => "observational data from SeaDataNet/EMODnet Chemistry Data Network",
            # Additional comment
            "comment" => "Every year of the time dimension corresponds to a 1-year centred average for one season.",

            # SeaDataNet Vocabulary P35 URN
            # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p35
            # example: SDN:P35::WATERTEMP
            "parameter_keyword_urn" => "$sdnp35", # Water body phosphate

            # List of SeaDataNet Parameter Discovery Vocabulary P02 URNs
            # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p02
            # example: ["SDN:P02::TEMP"]
            "search_keywords_urn" => ["$sdnp02"], # Water body phosphate

            # List of SeaDataNet Vocabulary C19 area URNs
            # SeaVoX salt and fresh water body gazetteer (C19)
            # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=C19
            # example: ["SDN:C19::3_1"]
            "area_keywords_urn" => ["SDN:C19::2"],
            "product_version" => "2.0",
            "product_code" => "SMHI-Baltic Sea-$(replace(varname,' '=>'_'))-v2023-ANA",

            # bathymetry source acknowledgement
            # see, e.g.
            # * EMODnet Bathymetry Consortium (2016): EMODnet Digital Bathymetry (DTM).
            # https://doi.org/10.12770/c7b53704-999d-4721-b1a3-04ec60c87238
            #
            # taken from
            # http://www.emodnet-bathymetry.eu/data-products/acknowledgement-in-publications
            #
            # * The GEBCO Digital Atlas published by the British Oceanographic Data Centre on behalf of IOC and IHO, 2003
            #
            # taken from
            # https://www.bodc.ac.uk/projects/data_management/international/gebco/gebco_digital_atlas/copyright_and_attribution/

            "bathymetry_source" => "The GEBCO 30sec Digital Atlas published by the British Oceanographic Data Centre on behalf of IOC and IHO, 2003",

            # NetCDF CF standard name
            # http://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
            # example "standard_name" = "sea_water_temperature",
            "netcdf_standard_name" => "$(replace(varname,' '=>'_'))",
            "netcdf_long_name" => "$varname",
            "netcdf_units" => "$unit",

            # Abstract for the product
            #"abstract" => "...",

            # This option provides a place to acknowledge various types of support for the
            # project that produced the data
            "acknowledgement" => "Aggregated data products are generated by SMHI with the support of SWAM",
            "documentation" => "https://doi.org/10.6092/A8CFB472-10DB-4225-9737-5A60DA9AF523",

            # Digital Object Identifier of the data product
            "doi" => "$doi",

            "DIVAnd_source" => "https://github.com/gher-ulg/DIVAnd.jl",
            "DIVAnd_version" => "2.7.5",
            "DIVA_code_doi" => "10.5281/zenodo.4715361",
            "DIVA_references" => "Barth, A.; Beckers, J.-M.; Troupin, C.; Alvera-Azcárate, A. & Vandenbulcke, L.
            divand-1.0: n-dimensional variational data analysis for ocean observations
            Geoscientific Model Development, 2014, 7, 225-241. DOI: :10.5194/gmd-7-225-2014");

        metadata[monthlist_index]=metadata_season
        ncglobalattrib = metadata_season
        ncvarattrib = OrderedDict(
            "standard_name" => "$(replace(varname,' '=>'_'))",
            "long_name" => "$varname",
            "units" => "$unit")

        @info("starting DIVAnd computations for $(seasons[monthlist_index])")
        @info(Dates.now())

        # Time selection for the analyse. This was already defined together with yearlist, month_list, seasons
        TS = DIVAnd.TimeSelectorYearListMonthList([year_list[year_list_index]],month_list[monthlist_index:monthlist_index])
        @show(TS)
        @info("$(month_list[monthlist_index:monthlist_index])")

        @info("Calling diva3d. Save background field in $nc_filepath")
        # create attributes for the netcdf file (need an internet connexion) and does sometimes not work. We fixed this by doing our own attributes.
        #ncglobalattrib,ncvarattrib = SDNMetadata(metadata_season,nc_filepath,varname,lonr,latr)

        dbinfo = @time diva3d((lonr,latr,depthr,TS),
            (obslon_sub,obslat_sub,obsdepth_sub,obstime_sub),
            obsval_sub,
            len,
            epsilon_weighted,
            nc_filepath, varname,
            bathname=bathname,
            bathisglobal = bathisglobal,
            ncvarattrib = ncvarattrib,
            ncglobalattrib = ncglobalattrib,
            timeorigin = timeorigin,
            mask = new_mask,
            solver = :direct,
            #fitcorrlen = false, # inte med i vanliga
            niter_e = 1,
            error_thresholds = error_thresholds,
            surfextend = true,
            alphabc = 0,
            minfield = 0.44662,  
            maxfield = 600,
            #stat_per_timeslice = true,
            MEMTOFIT = 250 # tidigare 120 dvs ej samma som analysen
        );

        # Save the observation metadata in the NetCDF file
        DIVAnd.saveobs(nc_filepath, "Oxygen_data", obsval_sub, (obslon_sub, obslat_sub, obsdepth_sub, obstime_sub), obsid_sub, used = dbinfo[:used])
    end    
end
