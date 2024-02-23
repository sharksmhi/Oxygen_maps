"""
julia julia_code\oxygen_analysis.jl, år, säsong, DIVAsettings
python python_code\calculata_areas.py, år, säsong, DIVASettings
"""
import subprocess
"""# year: list, seasons: list, DIVAsettings: dict
# location = "//winfs-proj/proj/havgem/DIVA/syrekartor/"
# outputdir = joinpath(location, "resultat/nc/$(savevar)/");
# data_fname = "EMODNET_SHARK_ICES.txt"
# year_list_background = [1960:1969,1970:1979,1980:1989,1990:1999,2000:2009,2010:2021];
# month_list_background = [ [11,12,1,2], [3,4,5], [6,7,8], [8,9,10]];
# filenamebackground
# year_list = []
# month_list = [ [11,12,1,2], [3,4,5], [6,7,8], [8,9,10]];
# seasons=["Winter","Spring","Summer","Autumn"]
# months=["(Nov-Feb)","(Mar-May)","(June-Aug)","(Aug-Oct)"];
# depthr = [0.,  10., 20., 25., 30., 35., 40., 50., 55., 60., 65., 70.,  75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 145.,150.,175.,200.,250.,300.];
"""

def run_julia_function(function_name, outputdir):
    try:
        # Run Julia script using subprocess and pass function name and arguments
        # process = subprocess.Popen(['julia', 'julia_code/my_julia.jl', function_name, arg1, arg2], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # ropa på vår oxygen analysis funktion i oxygen analysis juila scriptet. 
        # Lägg till fler argument efter outputdir med komma mellan
        process = subprocess.Popen(['julia', 'julia_code\oxygen_analysis.jl', 'run_oxygen', outputdir], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            print("Error occurred while executing Julia script:")
            print(stderr)
        else:
            print("Julia script executed successfully.")
            print("Output:")
            print(stdout)
    except FileNotFoundError:
        print("Julia executable not found. Make sure Julia is installed and added to the system PATH.")
        
if __name__ == "__main__":

    # Specify function name and arguments
    # function_name = "my_julia_function"
    # argument1 = "value1"
    # argument2 = "value2"

    function_name = "run_oxygen"
    location = "//winfs-proj/proj/havgem/DIVA/syrekartor/"
    argument2 = "outputdir.txt"

    # Call the function
    run_julia_function(function_name, location)



