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
freja_location = "path/on/freja"
on_freja = false
if isdir(freja_location)
    location = freja_location
    on_freja = true
end

outputdir = joinpath(location, "data/");
if !isdir(outputdir)
    @info("directory path for data does not exist $(outputdir)")
end
# Figures
figdir = "./resultat/figures/$(savevar)/";
if on_freja
    figdir = joinpath(outputdir, "resultat/figures")
end
if !isdir(figdir)
    mkpath(figdir)
    mkpath(joinpath(figdir, "test"))
end

# ## Load data big files created by program "data_handling"
data_fname = "SHARK_SYKE_IOW_EMODNET_ICES_250325"
#data_fname = "EMODNET_SHARK_ICES_SYKE_241216"
#data_fname = "mat_file_1960_2024_reordered"
@time obsval,obslon,obslat,obsdepth,obstime,obsid = loadbigfile(joinpath(location, "data/$data_fname.txt"));

dx, dy = 0.05, 0.05         #~5km?

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
depthr = Float64.(settings[basin]["depthr"])
lenz_ = Float64.(settings[basin]["lenz_"])
lenf = Float64.(settings[basin]["lenf"])
yearlist_json = settings[basin]["yearlist_background"]
# Konvertera varje par i yearlist till ett intervall (range) i Julia
year_list = [year[1]:year[2] for year in yearlist_json]
years = settings[basin]["years"]
threshold_list = settings[basin]["threshold_list"]

#BACKGROUND

month_list = [ [11,12,1,2], [3,4,5], [6,7,8], [8,9,10]];  # Seasonal climatology
TSbackground = DIVAnd.TimeSelectorYearListMonthList(year_list,month_list);
@show(TSbackground)
seasons=["Winter","Spring","Summer","Autumn"]
months=["(Nov-Feb)","(Mar-May)","(June-Aug)","(Aug-Oct)"];

# Time origin for the NetCDF file
timeorigin = DateTime(1900,1,1,0,0,0);
aspect_ratio = 1/cos(mean(latr) * pi/180);

"BATHYMETRY"
#bathname = joinpath(location, "bathymetry/gebco_30sec_4.nc")
bathname = joinpath(location, "bathymetry/bat_elevation_Baltic_Sea_masked.nc")
bath_file_name = split(bathname,"/")[end]
bathisglobal = true;
bx,by,b = DIVAnd.extract_bath(bathname,bathisglobal,lonr,latr);

# ## MASK
xmask,ymask,mmask = load_mask(bathname,true,lonr,latr,depthr);

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
#lenf = 40_000.
lenx = fill(lenf,sz)   # 200 km
leny = fill(lenf,sz)   # 200 km
#lenz = [min(max(10.,depthr[k]/150),300.) for i = 1:sz[1], j = 1:sz[2], k = 1:sz[3]]
#lenz = fill(10,sz);      # 25 m
lenz =  [lenz_[k] for i = 1:sz[1], j = 1:sz[2], k = 1:sz[3]];
len = (lenx, leny, lenz);

lx = lenf
ly = lenf

epsilon = 0.1;
epsilon = Float64.(settings["Global"]["epsilon_background"])

w_depth = 5.
w_days = 2.

rdiag_jldfile = joinpath(location, "data/$(data_fname)_weighted_$(w_depth)_$(w_days).jld")
@show rdiag_jldfile
if isfile(rdiag_jldfile)
    @load rdiag_jldfile rdiag
    @info "Loading saved rdiag file with w_depth = $(w_depth) and w_days = $(w_days)"
else
    @info "Calculating rdiag with w_depth = $(w_depth) and w_days = $(w_days)!"
    @time rdiag = 1 ./ DIVAnd.weight_RtimesOne((obslon,obslat,obsdepth,float.(Dates.dayofyear.(obstime))),(0.10,0.10,w_depth,w_days));
    @save rdiag_jldfile rdiag
end

@show maximum(rdiag),mean(rdiag)
epsilon_weighted = epsilon * rdiag;


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

# Settings for DIVAnd-------------------------------------------------------------------------------
error_thresholds = [("L1", 0.3), ("L2", 0.5)];

# Filnamn för utmatning
output_file = "obstime_output.txt"

# Skriv till fil
open(output_file, "w") do io
    write(io, "Original obstime and shifted obstime:\n")
    for (orig, shifted) in zip(obstime, obstime_shifted)
        write(io, "Original: $orig, Shifted: $shifted\n")
    end
end

#println("Data written to $output_file")

#filenamebackground = joinpath(outputdir, "$(replace(varname,' '=>'_'))_$(basin)_$(lenf)_$(epsilon)_$(years)_background_weighted_0.05_field_$(bath_file_name)")

#  dbinfo = @time diva3d((lonr,latr,depthr,TSbackground),
#             (obslon,obslat,obsdepth,obstime_shifted), obsval,
#             len, epsilon_weighted,
#             filenamebackground,varname,
#             bathname=bathname,
#             mask = new_mask,
#             fitcorrlen = false,
#             niter_e = 1,
#             solver = :direct,
#             MEMTOFIT = 120,
#         );

#A test to make separate files for each season in the background field to be able to calculate the
#A list of created files

# One metadata set up per season
metadata=Array{DataStructures.OrderedDict{String,Any}}(undef,4) ;
file_list = []

for monthlist_index in 1:length(month_list)
    season = seasons[monthlist_index]

    @info("Creating metadata dicitonary for the season $(season)")
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

    @show(metadata_season)
    @show(ncvarattrib)

    @info("starting DIVAnd computations for $(seasons[monthlist_index])")
    @info(Dates.now())

    # Time selection for the analyse. This was already defined together with yearlist, month_list, seasons
    TS = DIVAnd.TimeSelectorYearListMonthList(year_list,month_list[monthlist_index:monthlist_index])
    @show(TS)
    @info("$(month_list[monthlist_index:monthlist_index])")

    # File name based on the variable (but all spaces are replaced by _)
    #nc_filename = "Background_$(replace(varname,' '=>'_'))_$(years)_$(season)_$(epsilon)_$(lx)_$(dx)_$(w_depth)_$(w_days)_$(bath_file_name)_varcorrlenz.nc"
    nc_filename = "Background_$(replace(varname,' '=>'_'))_$(years)_$(season)_$(epsilon)_$(lx)_$(dx)_$(w_depth)_$(w_days)_$(bath_file_name)"
    #nc_filename_res = "$(replace(varname,' '=>'_'))_$(minimum(year_list))-$(maximum(year_list))_$(season)_residuals.nc"
    nc_filepath = joinpath(outputdir, nc_filename)
    #nc_filepath_res = joinpath("$(results_dir)/DIVArun", nc_filename_res)

    #Append the created files to file_list
    push!(file_list, nc_filename)

    if isfile(nc_filepath)
       rm(nc_filepath) # delete the previous analysis
    end

    @info("Will write results in $nc_filepath")
    # create attributes for the netcdf file (need an internet connexion) and does sometimes not work. We fixed this by doing our own attributes.
    #ncglobalattrib,ncvarattrib = SDNMetadata(metadata_season,nc_filepath,varname,lonr,latr)

    dbinfo = @time diva3d((lonr,latr,depthr,TS),
           (obslon,obslat,obsdepth,obstime_shifted),
           obsval,
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
           fitcorrlen = false, # inte med i vanliga
           niter_e = 1,
           error_thresholds = error_thresholds,
           MEMTOFIT = 120,
       );

    # Save the observation metadata in the NetCDF file
    DIVAnd.saveobs(nc_filepath,"Oxygen_data", obsval, (obslon,obslat,obsdepth,obstime_shifted),obsid, used = dbinfo[:used])

end
