# ### Oxygen maps
# ### 2023-2024
# ### Lena Viktorsson & Martin Hansson
# ### Based on copy from Karin Wesslander and DIVA workshop notebooks

# read dataset
# set up metadata info
# run DIVAnd and save to netcdf.

# #### Add necessary packages
# using Pkg
# Pkg.add("CSV")
# Pkg.add("DataFrames")

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
using PyCall
using CSV
using DataFrames


function rdiag_stats(
    year,
    season,
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
        season = [season],

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

function subset_year_month_mask(years, months, year_subset, month_subset)

    return (years .== year_subset) .& in.(months, Ref(month_subset))
end

function apply_cv_scenario(
    scenario,
    obslon,
    obslat,
    obsdepth,
    obsid
)
    nobs = length(obsid)
    keep_mask = trues(nobs)
    scenario_type = scenario["type"]

    if scenario_type == "remove_area"
        remove_mask =
            (obslon .>= scenario["lonmin"]) .&
            (obslon .<= scenario["lonmax"]) .&
            (obslat .>= scenario["latmin"]) .&
            (obslat .<= scenario["latmax"])
        keep_mask[remove_mask] .= false

    elseif scenario_type == "keep_selected_stations"

        keep_mask .= false
        stations = CSV.read(
            "standard_stations.txt",
            DataFrame;
            delim = '\t'
        )
        tol = scenario["match_tolerance_deg"]

        for station in eachrow(stations)
            station_mask =
                (abs.(obslon .- station.lon) .< tol) .&
                (abs.(obslat .- station.lat) .< tol)
            keep_mask[station_mask] .= true
        end
    else
        error("Unknown scenario type: $scenario_type")
    end
    @info "Keeping $(sum(keep_mask)) / $(length(keep_mask)) observations"
    return keep_mask
end

args= ARGS
input_dir = args[1]
results_dir = args[2]
data_fname = args[3]
year_list = JSON.parse(args[4])
month_list = JSON.parse(args[5])
seasons = JSON.parse(args[6])
lenf = JSON.parse(args[7])
epsilon = JSON.parse(args[8])
dx = JSON.parse(args[9])
bath_file_name = args[10]
w_depth = JSON.parse(args[11])
w_days = JSON.parse(args[12])
depthr = JSON.parse(args[13])
lenz_ = JSON.parse(args[14])
lonr = args[15]
latr = args[16]
basin = args[17]
threshold_list = JSON.parse(args[18])
epsilon_background = JSON.parse(args[19])
lenf_background = JSON.parse(args[20])
cv_mode = JSON.parse(args[21])
scenario_name = JSON.parse(args[22])

@show "read all arguments from python in julia"

weighting_dir = joinpath(input_dir, "weighting", data_fname)
if !isdir(weighting_dir)
    @info "$(weighting) is not a directory"
end
background_dir = joinpath(input_dir, "background_fields", data_fname)
if !isdir(background_dir)
    @info "$(background_dir) is not a directory"
end

dy = dx
lonr = replace(lonr, "dx" => string(dx))
latr = replace(latr, "dy" => string(dy))
lonr= eval(Meta.parse(lonr))
latr= eval(Meta.parse(latr))

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

# Figures
figdir = "./resultat/figures/$(savevar)/";
if !isdir(figdir)
    mkpath(figdir)
    mkpath(joinpath(figdir, "test"))
end

# ## Load data big files created by program "syrekartor_data_proc"
#data_fname = "EMODNET_SHARK_ICES"
@time obsval,obslon,obslat,obsdepth,obstime,obsid = loadbigfile(joinpath(input_dir, "$(data_fname).txt"))

# Sätt horisontell uppplösning i grader
# 0.05 motsvarra ca 5km

# Time origin for the NetCDF file
timeorigin = DateTime(1900,1,1,0,0,0);
aspect_ratio = 1/cos(mean(latr) * pi/180);

# ## Extract the bathymetry
# It is used to delimit the domain where the interpolation is performed.
# Modify bathname according to the resolution required.
# Bathymetry needs to be namned "bat"
#bath_file_name = "bat_elevation_Baltic_Sea_masked"
#bath_file_name = "gebco_30sec_4"
bathname = joinpath(input_dir, "$(bath_file_name).nc")
@show("extract bathymetry from $(bathname)")
bathisglobal = true;
@time bx,by,b = DIVAnd.extract_bath(bathname,bathisglobal,lonr,latr);

basin = replace(basin,' '=>'_')

# ## MASK
# #### 1) Create a mask from bathymetry and your selected depth vector
# First set the depth resolotion, this will be the depths used in DIVArun
# Then set the horizontal correlation length (should be twice the resolution)
#depthr = [0.,  10., 20., 25., 30., 35., 40., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 145.,150.,175.,200.,250.,300.];
#lenz_ =  [20., 20., 20., 10., 10., 10., 20., 20., 10., 10., 10., 10.,  10., 10., 10., 10., 10.,  50.,  50.,  50.,  50.,  50.,  50.,  50.,  50.,  50.,  50., 50., 50.,100., 100.,100.];
xmask,ymask,mmask = load_mask(bathname,true,lonr,latr,depthr);

# #### 2) apply floodfill.
# Karin: Varför depthr? borde gå med 0?
# It gives you an array which has the same size has the mask, and contains values equal to 0, 1, 2...
# 1 corresponds to the largest area, label 2 the 2nd largest and so on... fråga karin

# #### 3) Create a new mask
# which means it will have the value true only in the main sea area, not on land or in the small isolated pixels.
# I also attach a plot with the 2 masks (before and after the floodfill). You can see that your problem
# with the pixel in the middle of the island has disappeared. Note also that the northwesternmost part of the
# domain has also been masked, as this area was the 2nd largest sea area.
label = DIVAnd.floodfill(mmask);
new_mask = (label .== 1);

# #### 4) mask for Skagerrak
grid_bx = [i for i in xmask, j in ymask];
grid_by = [j for i in xmask, j in ymask];
mask_edit = copy(new_mask);

#Always exclude Skagerrak
sel_mask1 = (grid_by .>= 57.75) .& (grid_bx .<= 12.2);                              #Skagerrak
sel_mask2 = (grid_by .>= 57.4) .& (grid_by .< 57.75) .& (grid_bx .<= 10.4);         #Skagerrak
sel_mask3 = (grid_by .>= 57.) .& (grid_by .< 57.4) .& (grid_bx .<= 10.);            #Skagerrak

## If basin is....exclude other basins.
if basin == "Gulf_of_Bothnia"
    sel_mask4 = (grid_by .<= 60.2) .& (grid_bx .>= 9.) .& (grid_bx .<= 31.);                                 #Eg. Östersjön & Kattegatt
elseif basin == "Baltic_Proper"
    sel_mask4 = ((grid_by .>= 60.2) .& (grid_bx .>= 15.) .& (grid_bx .<= 25.)) .& (grid_bx .<= 9.);         #Bottniska viken & Kattegatt
elseif basin == "Kattegat"
    sel_mask4 = (grid_by .>= 53.) .& (grid_bx .>= 15.);                                                       #Eg. Östersjön och Bottniska viken
end
new_mask = mask_edit .* .!sel_mask1 .* .!sel_mask2 .* .!sel_mask3 .* .!sel_mask4;

# Ökas slen så minskas avståndet från kusten där len är mindre.
sz = (length(lonr),length(latr),length(depthr));
lenx = fill(lenf,sz)
leny = fill(lenf,sz)
#Horisontell korrelationslängd
lenz =  [lenz_[k] for i = 1:sz[1], j = 1:sz[2], k = 1:sz[3]];
#lenz = fill(10,sz);      # 10m 25 m

len = (lenx, leny, lenz);
lx = lenf
ly = lenf

# To include December from previous year in the analyse
obstime_shifted = copy(obstime)
# Get all dates with dec and nov and add one year to these specific dates.
# Winter period will be jan-feb + nov-dec from previous year.
#obstime_shifted[Dates.month.(obstime) .== 12 .| Dates.month.(obstime) .== 11] .+= Dates.Year(1)
# Iterera över elementen och justera datum om månaden är 11 eller 12
for i in eachindex(obstime)
    month = Dates.month(obstime[i])  # Hämta månaden
    if month == 11 || month == 12
        obstime_shifted[i] += Year(1)  # Lägg till ett år
    end
end

# För att köra cross-validation filtrera bort en del av data
# Filtrera bort data innan viktningen av data så att viktningen görs på data som är kvar
# Det betyder en ny rdiag fil för varje validering.
if cv_mode
    @show "run sensitivity mode"
    cv_settings = JSON.parse(read("cv_settings.json", String))
    @show  cv_settings[scenario_name]
    keep_mask = apply_cv_scenario(
        cv_settings[scenario_name],
        obslon,
        obslat,
        obsdepth,
        obsid
    )
    @show length(keep_mask)
    removed_mask = .!keep_mask
    obsval_removed   = obsval[removed_mask]
    obslon_removed   = obslon[removed_mask]
    obslat_removed   = obslat[removed_mask]
    obsdepth_removed = obsdepth[removed_mask]
    obstime_shifted_removed  = obstime_shifted[removed_mask]
    obsid_removed    = obsid[removed_mask]

    obsval   = obsval[keep_mask]
    obslon   = obslon[keep_mask]
    obslat   = obslat[keep_mask]
    obsdepth = obsdepth[keep_mask]
    obstime_shifted  = obstime_shifted[keep_mask]
    obsid    = obsid[keep_mask]

end

years  = Dates.year.(obstime_shifted)
months = Dates.month.(obstime_shifted)

# Settings for DIVAnd-------------------------------------------------------------------------------
error_thresholds = [("L1", 0.3), ("L2", 0.5)];

# Modify data weight
# Compute the new weights that takes into account close points.
# If the dataset is large, this can take a few minutes.
# The maximal and mean values provide an indication of the spatial proximity between the data.
# If you apply this technique, you need to adapt epsilon2:

# Make a test with different values for depth and days for rdiag.
# 60m och 30 dagar - funkar ej
# 5m och 7 dagar - funkar ej
# 20m och 7 dagar - funkar ej
# 10m och 7 dagar - funkar ej
# 2.5m och 7 dagar - funkar bäst just nu
# 5m och 2 dagar - funkar bäst just nu
# 2.5m och 2 dagar
# 2.5m och 30 dagar


# One metadata set up per season
metadata=Array{DataStructures.OrderedDict{String,Any}}(undef,4) ;

# Ange vilket intervall du vill ha på bakgrundfältet
start_year = 1960
end_year   = 2024

month_list_background = [[1,2,3,4,5,6,7,8,9,10,11,12]];  # 3 whole years climatology
stats_file = joinpath(results_dir, "rdiag_statistics.txt")

for year in year_list
    @show(year)
    #Background file for choosen year
    bkg_filename = "Background_$(varname)_$(year)_All_$(epsilon_background)_$(lenf_background)_$(dx)_$(w_depth)_$(w_days)_$(basin).nc"
    @show(bkg_filename)
    bkg_filepath = joinpath(input_dir, bkg_filename)
    cp(bkg_filepath, joinpath(results_dir, joinpath("DIVArun", bkg_filename)); force=true)
    TSbackground = DIVAnd.TimeSelectorYearListMonthList([year],month_list_background);
    @show(TSbackground)
    diva_background_arg = DIVAnd.backgroundfile(bkg_filepath,varname,TSbackground)
    # Lägg in subsetting av data här så vi bara använder de tre åren i fönstret för bakgrundsfältet

    for monthlist_index in 1:length(month_list)
        season = seasons[monthlist_index]
        @show(season)
        # räkna rdiag för endast den år-säsong vi kör och vikta med endast de data
        # spara rdiag filen, behövs en per år-säsong till istället för en för hela perioden
        year_month_mask = subset_year_month_mask(years, months, year, month_list[monthlist_index])
        obsval_sub   = obsval[year_month_mask]
        obslon_sub   = obslon[year_month_mask]
        obslat_sub   = obslat[year_month_mask]
        obsdepth_sub = obsdepth[year_month_mask]
        obstime_shifted_sub  = obstime_shifted[year_month_mask]
        obsid_sub    = obsid[year_month_mask]
        rdiag_jldfile = joinpath(weighting_dir, "$(data_fname)_$(year)_$(season)_weighted_$(w_depth)_$(w_days).jld")
        if isfile(rdiag_jldfile) && !cv_mode
            @load rdiag_jldfile rdiag
            @info "Loading saved rdiag file with w_depth = $(w_depth) and w_days = $(w_days)"
        else
            @info "Calculating rdiag with w_depth = $(w_depth) and w_days = $(w_days)!"
            @time rdiag = 1 ./ DIVAnd.weight_RtimesOne((obslon_sub,obslat_sub,obsdepth_sub,float.(Dates.dayofyear.(obstime_shifted_sub))),(0.10,0.10,w_depth,w_days));
            # only save for normal runs
            if !cv_mode
                @info "saving new rdiag file"
                @save rdiag_jldfile rdiag
            end
        end
        
        epsilon_weighted = epsilon * rdiag;

        stats_df = rdiag_stats(
            year,
            season,
            obsid_sub,
            obslon_sub,
            obslat_sub,
            obsdepth_sub,
            obstime_shifted_sub,
            rdiag,
            epsilon_weighted
        )
        if isfile(stats_file)
            CSV.write(stats_file, stats_df; append=true)
        else
            CSV.write(stats_file, stats_df)
        end
        
        @info("Creating metadata dicitionary for the season $(season)")
        metadata_season = OrderedDict(
            # set attributes for DVIA run from our settings
            "threshold_list" => "$threshold_list",
            "season" => season,
            "epsilon" => "$epsilon",
            "horizontal correlation length m" => "$lx",
            "start year" => string(year_list[1]),
            "end year" => string(year_list[end]),
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
              
        @info("$(Dates.now()) starting DIVAnd computations for $(year) $(seasons[monthlist_index])")

        # Time selection for the analyse. This was already defined together with yearlist, month_list, seasons
        TS = DIVAnd.TimeSelectorYearListMonthList([year],month_list[monthlist_index:monthlist_index])
        @show(TS)

        # File name based on the variable (but all spaces are replaced by _)
        nc_filename = "$(replace(varname,' '=>'_'))_$(year)_$(season)_$(epsilon)_$(lx)_$(dx)_$(w_depth)_$(w_days)_$(basin)_varcorrlenz.nc"
        nc_filename_res = "$(replace(varname,' '=>'_'))_$(year)_$(season)_residuals.nc"
        nc_filepath = joinpath("$(results_dir)/DIVArun", nc_filename)
        nc_filepath_res = joinpath("$(results_dir)/DIVArun", nc_filename_res)

        if isfile(nc_filepath)
            continue
            # rm(nc_filepath) # delete the previous analysis
        end

        @info("Will write results in $nc_filepath")
        # create attributes for the netcdf file (need an internet connexion) and does sometimes not work. We fixed this by doing our own attributes.
        #ncglobalattrib,ncvarattrib = SDNMetadata(metadata_season,nc_filepath,varname,lonr,latr)
        @show(bkg_filename)
        @time dbinfo = diva3d((lonr, latr, depthr, TS),
                (obslon_sub, obslat_sub, obsdepth_sub, obstime_shifted_sub),
                obsval_sub,
                len,
                epsilon_weighted, # error variance of the observations (normalized by the error variance of the background field)
                nc_filepath,
                varname,
                bathname = bathname,
                bathisglobal = bathisglobal,
                ncvarattrib = ncvarattrib, # dictionary of the netcdf variable attributes
                ncglobalattrib = ncglobalattrib, # dictionary of the netcdf global attributes
                timeorigin = timeorigin,
                # Nedan anges epsilon inte epsilon2: Dvs ngt litet. Det är här Karin och Örjan använt en annan epsilon än den ovan
                # vi har döpt om epsilon2 till epsilon_weighted för att förtydliga vad för ett epsilon det är
                #transform = Anam.loglin(0.01),
                minfield = 0.44662,  # Keep the results within min-max field
                maxfield = 600,
                #transform = Anam.loglin(270.; epsilon = 0.05),
                mask = new_mask,
                solver = :direct,
                niter_e = 1,
                background = diva_background_arg,
                error_thresholds = error_thresholds,
                surfextend = true,
                alphabc = 0,
                #stat_per_timeslice = true,
                MEMTOFIT = 250
        );

        #residuals = dbinfo[:residuals]
        residuals = get(dbinfo, :residuals, 0)
        @show(keys(dbinfo))
        # Residuals with NaNs removed
        sel = .!isnan.(residuals)
        @show(length(residuals))
        @show(length(obsval_sub))

        #@show extrema(res2);
        #@show quantile(res2, [0.01, 0.99]);
        @info("Lowest and highest residuals might be a error sample...>250 and <-250")
        indices = findall(x -> x > 250 || x < -250, residuals)
        println("Number of possible data errors: $length(indices)")
        # Visa både index och värde
        for i in indices
            println("Residual: $(residuals[i]), obsval: $(obsval_sub[i]), DIVAnd: $(obsval_sub[i]-residuals[i]), obsdepth:$(obsdepth_sub[i]), obsid: $(obsid_sub[i])")
        end
        # large_residuals_df = DataFrame(
        #     Residual = residuals[indices],
        #     obsval = obsval[indices],
        #     DIVAnd = obsval_sub[indices]-obsval[indices],
        #     obsdepth = obsdepth_sub[indices],
        #     obsid = obsid_sub[indices],
        #     Residual = residuals[indices],
        #     Residual = residuals[indices],
        #     )
        # if isfile(large_residuals_file)
        #     CSV.write(large_residuals_file, large_residuals_df; append=true)
        # else
        #     CSV.write(large_residuals_file, large_residuals_df)
        # end
    #     @info("Get the residuals...")
    #     selection_per_timeslice = dbinfo[:selection_per_timeslice]
    #     residuals_per_timeslice = dbinfo[:residuals_per_timeslice]
    #     selection_per_timeslice = dbinfo[:selection_per_timeslice]
    #
    #     max_residuals = fill(-Inf,length(residuals_per_timeslice))
    #     for n = 1:length(selection_per_timeslice)
    #         sel = selection_per_timeslice[n]
    #         max_residuals[sel] = max.(max_residuals[sel],residuals_per_timeslice[n])
    #     end
    #     @show(max_residuals)

        # Save the observation metadata in the NetCDF file
        DIVAnd.saveobs(nc_filepath,"Oxygen_data", obsval_sub, (obslon_sub, obslat_sub, obsdepth_sub, obstime_shifted_sub),obsid_sub, used = dbinfo[:used])
        DIVAnd.saveobs(nc_filepath_res,"$(varname)_residual",residuals[sel], (obslon_sub[sel], obslat_sub[sel], obsdepth_sub[sel], obstime_shifted_sub[sel]),obsid_sub[sel])

        if cv_mode
            ds = NCDataset(nc_filepath)
            @show size(ds[varname])
            field = ds[varname]
            nc_filename_removed = "$(replace(varname,' '=>'_'))_$(year)_$(season)_removed.nc"
            nc_filename_residuals_removed = "$(replace(varname,' '=>'_'))_$(year)_$(season)_residuals_removed.nc"
            nc_filepath_removed = joinpath("$(results_dir)/DIVArun", nc_filename_removed)
            nc_filepath_residuals_removed = joinpath("$(results_dir)/DIVArun", nc_filename_residuals_removed)
            # pred_removed = DIVAnd.interp(
            #     field,
            #     (obslon_removed, obslat_removed, obsdepth_removed, obstime_removed)
            # )
            # res_removed = obsval_removed .- pred_removed
            # DIVAnd.saveobs(nc_filepath_residuals_removed, "$(varname)_residual" res_removed, 
            # (obslon_removed, obslat_removed, obsdepth_removed, obstime_shifted_removed), obsid_removed)
            DIVAnd.saveobs(nc_filepath_removed, "Oxygen_data", obsval_removed, 
            (obslon_removed, obslat_removed, obsdepth_removed, obstime_shifted_removed), obsid_removed)
        end

        # Sparar undan alla indata, Divananlays och residualer i varje indata punkt.
        res_filepath = joinpath("$(results_dir)/DIVArun", "$(varname)_$(year)_$(season)_residual.txt")
        diva_result = obsval_sub[sel] .- residuals[sel]
        res_data = [obsval_sub[sel] diva_result residuals[sel] obslon_sub[sel]  obslat_sub[sel]  obsdepth_sub[sel]  obstime_shifted_sub[sel]  obsid_sub[sel]]
        
        # Definiera header som en sträng med tab-separerade kolumnnamn
        #obsval	diva	residual	lat	long	obsdepth	time	id
        header = "obsval\tdiva\tresidual\tlong\tlat\tobsdepth\ttime\tid"
        
        # Skriv till fil med header
        open(res_filepath, "w") do io
            write(io, header * "\n")
            writedlm(io, res_data, '\t')
        end 
        
        # Spara till fil
        #writedlm(res_filepath, res_data, '\t')

    end # end of season loop
end # end of year loop   
