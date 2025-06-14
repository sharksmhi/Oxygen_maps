"""
julia julia_code\oxygen_analysis.jl, år, säsong, DIVAsettings
python python_code\calculata_areas.py, år, säsong, DIVASettings
"""
import shutil
import subprocess
import json
from python_code import calculate_areas as calculate_areas
from python_code import plot_result
from pathlib import Path
import numpy as np
from sys import stdout
import datetime as dt

def run_julia_function(args):
    try:
        # Run Julia script using subprocess and pass function name and arguments
        # process = subprocess.Popen(['julia', 'julia_code/my_julia.jl', function_name, arg1, arg2], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # ropa på vår oxygen analysis funktion i oxygen analysis juila scriptet. 
        # Lägg till fler argument efter outputdir med komma mellan
        subprocess.run(args, check=True)

    except FileNotFoundError:
        print("Julia executable not found. Make sure Julia is installed and added to the system PATH.")
        
if __name__ == "__main__":

    # Data input directory should contain files with :
    # - raw data (text file, of the format bigfile that DIVAnd reads)
    # - bathymetry (netcdf file)
    # - backgroundfields (netcdf file either from DIVAnd or modeldata)
    input_dir = "data"
    freja_input_dir = "/path/on/freja/input_data"
    on_freja = False
    if Path(freja_input_dir).is_dir():
        input_dir = freja_input_dir
        on_freja = True

    # Input data filename
    #data_fname = "mat_file_1960_2024_reordered.txt"
    data_fname = "SHARK_SYKE_IOW_EMODNET_ICES_250325.txt"

    # Definiera basins
    #basin = "Kattegat"
    basin = "Baltic_Proper"
    #basin = "Gulf_of_Bothnia"

    # Läs in JSON-filen
    with open('settings.json', 'r') as file:
        settings = json.load(file)
   
    lonr_range = settings[basin]["lonr"]
    latr_range = settings[basin]["latr"]
    dx = json.dumps(settings[basin]["dx"])
    dy = json.dumps(settings[basin]["dy"])
    # Skapa intervall (range) i Julia med angivet dx
    lonr = f"{lonr_range[0]}:{dx}:{lonr_range[1]}"
    latr = f"{latr_range[0]}:{dy}:{latr_range[1]}"
    epsilon = json.dumps(settings[basin]["epsilon"])
    depthr = json.dumps(settings[basin]["depthr"])
    lenz_ = json.dumps(settings[basin]["lenz_"])
    lenf = json.dumps(settings[basin]["lenf"])
    threshold_list = json.dumps(settings[basin]["threshold_list"])
    years = settings[basin]["years"] # utan json.dumps så det passar i bkg_filename strängen.
    yearlist_background = json.dumps(settings[basin]["yearlist_background"])

    epsilon_background = json.dumps(settings["Global"]["epsilon_background"])

    # Visa extracted information
    print("Show the setup for the choosen basin: ", basin)
    print("lonr:", lonr)
    print("latr:", latr)
    print("lenf:", lenf)
    print("epsilon:", epsilon)
    print("dx:", dx)
    print("threshold list:", threshold_list)
    print("depthr list:", depthr)
    print("lenz_ list:", lenz_)

    #Samlar ihop resultaten:
    today = dt.datetime.now().strftime("%Y%m%d_%H%M")
    results_dir = Path(f"resultat/{basin.replace(' ', '_')}/{today}/")
    if on_freja:
        results_dir = Path(f"path/on/freja/resultat/{basin.replace(' ', '_')}/{today}/")
    results_dir.mkdir(parents=True, exist_ok=True)
    Path(results_dir, "figures/").mkdir(parents=True, exist_ok=True)
    Path(results_dir, "DIVArun/").mkdir(parents=True, exist_ok=True)
    Path(results_dir, "processed/").mkdir(parents=True, exist_ok=True)
    
    with open(Path(results_dir, 'settings.json'), 'w') as file:
        json.dump(obj=settings[basin],fp=file)
 
    # Years, month and seasons to be analysed
    year_list = json.dumps([1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978,
    1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997,
    1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
    2017, 2018, 2019, 2020, 2021, 2022])
    year_list = json.dumps([1960, 1965, 1970, 1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020])
    #year_list = json.dumps([2015])
    print(f"calculating for years {year_list}")

    seasons_dict = {
                "Winter": [11,12,1,2],
                "Spring": [3,4,5],
                "Summer": [6,7,8],
                "Autumn": [8,9,10]
               }
    seasons = []
    month_list = []
    for season, months in seasons_dict.items():
        seasons.append(season)
        month_list.append(months)

    seasons = json.dumps(seasons)
    month_list = json.dumps(month_list)

    # Correlation length
    # Vi bör köra med lite längre lenf troligen 80_000km då vi har ca 40nm mellan våra station i eg.Östersjön
    #lenf = json.dumps(80000)    #Km
    # Resolution
    #dx = json.dumps(0.05)       #deg
    # Signal to noise ratio
    # low epsilon means higher noise in data and result is more smoothed
    # high epsilon means lower noise in data and result is less smoothed and each observation is seen more
    #epsilon = json.dumps(0.2)
    #Bathymetry file
    bath_file_name = "bat_elevation_Baltic_Sea_masked"
    varname = "Oxygen"

    print("Bathymetry file: ", bath_file_name)

    # Modify data weight
    w_depth = json.dumps(5.)
    w_days = json.dumps(2.)

    background_filename = f"Background_Oxygen_{years}_{season}_{epsilon_background}_{lenf}_{dx}_{w_depth}_{w_days}_{bath_file_name}.nc"

    bkg_filename = json.dumps(f"{input_dir}/{background_filename}")

    #Set True if you want to save area_data to file (for time-series bar plots)
    save_area_data=True

    print("running DIVAnd in Julia...")
    args = ['julia', 'julia_code/oxygen_analysis.jl', input_dir, results_dir, data_fname, year_list, month_list, seasons, lenf, epsilon, dx, bath_file_name, w_depth, w_days, depthr, lenz_, lonr, latr, basin, threshold_list, bkg_filename, yearlist_background, years, epsilon_background]

    # #Call the function and save a json-file with a file_list containing the results. That we can send to the calculate_areas function.
    try:
        run_julia_function(args)
    except Exception as e:
        # If exception occurs, prompt user
        print(e)
        user_input = input(f"An error occurred. Do you want to remove the path resultat/{basin.replace(' ', '_')}/{today}/? (y/n): ")
        
        if user_input.lower() == 'y':
            path_to_remove = Path(f"resultat/{basin.replace(' ', '_')}/{today}/")
            
            # Check if the path exists before trying to remove it
            if results_dir.exists() and results_dir.is_dir():
                shutil.rmtree(results_dir)
                print(f"Path {results_dir} has been removed.")
                exit()
            else:
                print(f"Path {results_dir} does not exist.")
                exit()
        else:
            print("No path was removed.")
            exit()

    file_list = []

    # Add backgroundfield and analysis to file_list to cal areas and plot
    for season in json.loads(seasons):
        file_list.append(f"Background_Oxygen_{years}_{season}_{epsilon_background}_{lenf}_{dx}_{w_depth}_{w_days}_{bath_file_name}.nc")
        file_list.append(f"Oxygen_{min(json.loads(year_list))}-{max(json.loads(year_list))}_{season}_{json.loads(epsilon)}_{lenf}_{json.loads(dx)}_{w_depth}_{w_days}_{bath_file_name}_varcorrlenz.nc")
    
    # #Calculate areas from DIVA-results and save in a new nc-file. Results in file_list
    print("calculating areas...")
    print(file_list)
    calculate_areas.calculate_areas(results_dir, file_list, json.loads(threshold_list), save_area_data)

    # Read and plot areas in file_list
    print("plotting...")
    plot_result.read_processed_nc(results_dir,file_list,year_list,yearlist_background)

print("DIVAnd is done with its stuff...")
### extract values that are within our limits, save to a new variable and nc-file. ####
# 1 ml/l of O2 is approximately 43.570 µmol/kg
# (assumes a molar volume of O2 of 22.392 l/mole and
# a constant seawater potential density of 1025 kg/m3).
# https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html

