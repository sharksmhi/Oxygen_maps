#=
syrekartor_data_proc:
- Julia version: 
- Author: a001109
- Date: 2023-08-30

Program to handle O2 data before DIVAnd.
Emodnet data and SHARKweb data
Removes duplicates etc.
# ### Lena Viktorsson & Martin Hansson
=#

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
outputdir = "./resultat/nc/$(savevar)/" ;
if !isdir(outputdir)
    mkpath(outputdir)
end
# Figures
figdir = "./resultat/figures/$(savevar)/";
if !isdir(figdir)
    mkpath(figdir)
    mkpath(joinpath(figdir, "test"))
end

# ## Load data big files
# bot_no_header.txt är all data från ICES - 1960 till 2018 den använder vi inte!
# sharkweb_btlctd_O2 är alla data från SHARK 1960-2022
# emodnet BTL + CTD

#fname = "data/bot_no_header.txt"
#@time obsval,obslon,obslat,obsdepth,obstime,obsid = loadbigfile(fname);
#@show(length(obsval));

fname_shark = "C:/Work/DIVAnd/Oxygen_maps/data/sharkweb_btlctd_02.txt"
@time obsval_shark,obslon_shark,obslat_shark,obsdepth_shark,obstime_shark,obsid_shark = loadbigfile(fname_shark);

#plot(obsdepth, obsval, "ko", markersize=0.5);
#sel = (Dates.year.(obstime) .== 2021) .& (Dates.month.(obstime) .>= 8) .& (Dates.month.(obstime) .<= 10) .& (obsdepth .== 125.)
#@show(length(obsval[sel]))
#plot(obsdepth[sel], obsval[sel], "ro", markersize=0.5);

datafile_emod_btl = "C:/Work/DIVAnd/Oxygen_maps/data/emodnet_02_1960_2023_BTL.txt"
@time obsval_emod_btl,obslon_emod_btl,obslat_emod_btl,obsdepth_emod_btl,obstime_emod_btl,obsid_emod_btl = ODVspreadsheet.load(Float64,[datafile_emod_btl],
                           ["Water body dissolved oxygen concentration"]; nametype = :localname );

datafile_emod_ctd = "C:/Work/DIVAnd/Oxygen_maps/data/emodnet_02_1960_2023_CTD.txt"
@time obsval_emod_ctd,obslon_emod_ctd,obslat_emod_ctd,obsdepth_emod_ctd,obstime_emod_ctd,obsid_emod_ctd = ODVspreadsheet.load(Float64,[datafile_emod_ctd],
                           ["Water body dissolved oxygen concentration"]; nametype = :localname );
@show(length(obsval_emod_btl));
@show(length(obsval_emod_ctd));
@show(length(obsval_shark));

# ## Remove low res CTD when BTL is available.
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
newpoints = isempty.(dupl);
@info("Number of new points: $(sum(newpoints)))")

obslon = [obslon_emod; obslon_shark[newpoints]];
obslat = [obslat_emod; obslat_shark[newpoints]];
obsdepth = [obsdepth_emod; obsdepth_shark[newpoints]];
obstime = [obstime_emod; obstime_shark[newpoints]];
obsval = [obsval_emod; obsval_shark[newpoints]];
obsid = [obsid_emod; obsid_shark[newpoints]];

# ## Create a plot showing the additional data points:
figure("Additional-Data")
ax = subplot(1,1,1)
ax.tick_params("both",labelsize=6)
ylim(53, 65);
xlim(9, 27);
#contourf(bx, by, permutedims(Float64.(mask_edit[:,:,1]),[2,1]),
#    levels=[-1e5,0],cmap="binary");
plot(obslon_emod, obslat_emod, "bo", markersize=.2, label="Emodnet")
plot(obslon_shark[newpoints], obslat_shark[newpoints], "go",
   markersize=.2, label="Additional data\nfrom SHARKweb")
legend(loc=3, fontsize=4)
gca().set_aspect(aspectratio)
figname = "additional_data.png"
@show joinpath("$figdir/temp", figname)
PyPlot.savefig(joinpath("$figdir/temp", figname), dpi=300);
PyPlot.close_figs()


# ## Load data from Emodnet Kattegatt
#ncfile_1 = "C:/Work/DIVAnd/Oxygen_maps/data/data_from_Eutrophication_NorthSea_non-nutrient_profiles_2022_unrestricted_1960_2022_02.nc"
#@time obsval_1, obslon_1, obslat_1, obsdepth_1, obstime_1 = NCODV.load(Float64, ncfile_1, "Water body dissolved oxygen concentration");
## ## Load data from Emodnet Eg. Ö (tre filer för ODV kunde inte exportera allt...)
#ncfile_2 = "C:/Work/DIVAnd/Oxygen_maps/data/data_from_Eutrophication_Baltic_profiles_2022_unrestricted_O2_1960_1980.nc"
#@time obsval_2, obslon_2, obslat_2, obsdepth_2, obstime_2 = NCODV.load(Float64, ncfile_2, "Water body dissolved oxygen concentration");
#ncfile_3 = "C:/Work/DIVAnd/Oxygen_maps/data/data_from_Eutrophication_Baltic_profiles_2022_unrestricted_O2_1981_2000.nc"
#@time obsval_3 , obslon_3 , obslat_3, obsdepth_3, obstime_3 = NCODV.load(Float64, ncfile_3, "Water body dissolved oxygen concentration");
#ncfile_4 = "C:/Work/DIVAnd/Oxygen_maps/data/data_from_Eutrophication_Baltic_profiles_2022_unrestricted_O2_2001_2005.nc"
#@time obsval_4 , obslon_4 , obslat_4, obsdepth_4, obstime_4 = NCODV.load(Float64, ncfile_4, "Water body dissolved oxygen concentration");
#ncfile_5 = "C:/Work/DIVAnd/Oxygen_maps/data/data_from_Eutrophication_Baltic_profiles_2022_unrestricted_O2_2006_2008.nc"
#@time obsval_5 , obslon_5 , obslat_5, obsdepth_5, obstime_5 = NCODV.load(Float64, ncfile_5, "Water body dissolved oxygen concentration");

# ## Slå ihop alla edmodnetdata till ett dataset för att sedan kolla dubbletter.
#obsval_emod   = [obsval_1; obsval_2; obsval_3; obsval_4; obsval_5];
#obslon_emod   = [obslon_1; obslon_2; obslon_3; obslon_4; obslon_5];
#obslat_emod   = [obslat_1; obslat_2; obslat_3; obslat_4; obslat_5];
#obsdepth_emod = [obsdepth_1; obsdepth_2; obsdepth_3; obsdepth_4; obsdepth_5];
#obstime_emod  = [obstime_1; obstime_2; obstime_3; obstime_4; obstime_5];
#obsid    = [obsid; obsidns; obsid2];


#Så här kan man lägga ihop olika dataset
#obsval   = [obsval; obsvalns; obsval2];
#obslon   = [obslon; obslonns; obslon2];
#obslat   = [obslat; obslatns; obslat2];
#obsdepth = [obsdepth; obsdepthns; obsdepth2];
#obstime  = [obstime; obstimens; obstime2];
#obsid    = [obsid; obsidns; obsid2];





#Gör ett test på data för att hitta outliers
figure("Data")
ax = subplot(1,1,1)
plot(obslon[sel], obslat[sel], "ko", markersize=.1)
aspectratio = 1/cos(mean(54:0.05:61) * pi/180)
#ax.tick_params("both",labelsize=6)
gca().set_aspect(aspectratio)
figname = "observations.png"
@show joinpath("$figdir/temp", figname)
PyPlot.savefig(joinpath("$figdir/temp", figname), dpi=300);
PyPlot.close_figs()
