# ##Program to convert a ICES-file to a "big file" readable for DIVAnd
# ##ICES data from ices.dk BTL and CTD low res dataset
# ##

import pandas as pd
import numpy as np

# Define file path and data types for each column
#file_path = "//winfs-proj/proj/havgem/DIVA/syrekartor/data/ICES_1960_2022/"
file_path = "C:/Work/DIVAnd/Oxygen_maps/data/all_baltic/"
file_name = "data_from_3dfff340-bf23-4e51-8103-b2da4c6a0803.txt"

dtype_dict = {
    "Cruise": str,
    "Station": str,
    "Latitude [degrees_north]": str,
    "Longitude [degrees_east]": str,
    "mon/day/yr":str,
    "hh:mm":str,
    #"Year": str,
    #"Month": str,
    #"Day": str,
    #"Hour": str,
    #"Minute": str,
    "Depth (ADEPZZ01_ULAA) [m]": float,
    "Temperature (TEMPPR01_UPAA) [degC]": float,
    "Salinity (PSALPR01_UUUU) [dmnless]": float,
    "Oxygen (DOXYZZXX_UMLL) [ml/l]": float,
    "Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]": float
}

# Load data from CSV file into a pandas DataFrame
#df = pd.read_csv(file_path+file_name, dtype=dtype_dict)
#df = pd.read_table(file_path+file_name, delimiter="\t",low_memory=False)
df = pd.read_table(file_path+file_name, delimiter="\t",dtype=dtype_dict)

# Rename the '%Platform' column to 'ID'
#df.rename(columns={"Sampling platform (code)": "ID"}, inplace=True)

# Konvertera datumkolumnen till datetime
df['datum'] = pd.to_datetime(df['mon/day/yr'], format='%m/%d/%Y')
# Extrahera månaden, dag, år
df['Month'] = df['datum'].dt.month.astype(str)
df['Day'] = df['datum'].dt.day.astype(str)
df['Year'] = df['datum'].dt.year.astype(str)

df["Visit event identifier"] = df["Year"] + '-' + df["Cruise"] + '-' + df["Station"]

# Replace invalid hours with '00:00'
df.loc[df["hh:mm"].isin(['-9', '24', '',np.nan]), "hh:mm"] = '00:00'
#df.loc[df["hh:mm"].isin(['-9', '60', '',np.nan]), "hh:mm"] = '00:00'

#df["Time"] = df["hh:mm"][0:1].str.zfill(2) + ":" + df["hh:mm"][3:4].str.zfill(2)

# Combine the date and hour columns into a new column 'date_time'
df['date_time'] = df["Year"] + "-" + df["Month"].str.zfill(2) + "-" + df["Day"].str.zfill(2) + "T" + df["hh:mm"]

# Make a new combined O2 variable.
# Base new variable on O2 bottle and add O2 CTD when O2 bottle is missing.
def get_doxy(row):
    if row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"] and not np.isnan(row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"]):
        return row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"] * -0.04488
        #return row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"] * 0 + 0.01
    if row['Oxygen (DOXYZZXX_UMLL) [ml/l]'] and not np.isnan(row['Oxygen (DOXYZZXX_UMLL) [ml/l]']):
        return row['Oxygen (DOXYZZXX_UMLL) [ml/l]']
    return np.nan

df["DOXY_ml"] = df.apply(get_doxy, axis=1)

# Convert O2 ml/l to µmol/l-> 1ml/l = 44.661 µmol/l
#df['DOXY_umol'] = df['DOXY_ml'].multiply(44.661)
df['DOXY_umol'] = df['DOXY_ml'].apply(lambda x: x*44.661)
print(df.dtypes)
# Filter out rows with missing data
#df_filtered = df.loc[df["DOXY_umol"] != np.nan]
df_filtered = df.dropna(subset=["DOXY_umol"])


#Ta bort alla data som är ryska(90))
# Lista med mönster att ta bort
remove = ['90']

# Filtrera ut rader som matchar något av mönstren
removed_ICES_data = df[df["Cruise"].str.startswith(tuple(remove))]
# Behåll endast rader som inte matchar mönstren
df = df[~df["Cruise"].str.startswith(tuple(remove))]

# Spara de borttagna raderna till en ny CSV-fil (valfritt)
removed_ICES_data.to_csv("removed_data.csv", index=False)

# Define the order of columns for the output files
#column_list = ["Longdeg", "Latdeg", "DOXY(umol/l)", "Pressure(Dbars)", "TEMP", "PSAL", "H2SX(umol/l)", "St_No", "Year", "date_time", "ID", "Mnth", "Dy", "Hr"]
#column_list = ["Sample longitude (DD)", "Sample latitude (DD)", "DOXY_umol", "Sampling depth (m)", "Temperature bottle (C)", "Salinity bottle (o/oo psu)", "Hydrogen sulphide H2S (umol/l)", "Visit event identifier", "Year", "date_time", "Sampling platform (code)", "Month", "Day", "Sampling time (start)"]
column_list = ["Longitude [degrees_east]", "Latitude [degrees_north]", "DOXY_umol", "Depth (ADEPZZ01_ULAA) [m]", "Temperature (TEMPPR01_UPAA) [degC]", "Salinity (PSALPR01_UUUU) [dmnless]", "Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]", "Cruise", "Year", "date_time", "Visit event identifier", "Month", "Day", "hh:mm"]

# Write the filtered data to two output files, one with headers and one without
# df_filtered[column_list].to_csv(file_path+"bot.txt", index=False, sep='\t')
df_filtered[column_list].to_csv(file_path+"ICES_btl_lowres_ctd_02_241107.txt", index=False, header=False, sep='\t')



