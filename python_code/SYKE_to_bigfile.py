# ##Program to convert a ICES-file to a "big file" readable for DIVAnd
# ##ICES data from ices.dk BTL and CTD low res dataset
# ##

import pandas as pd
import numpy as np

# Define file path and data types for each column
#file_path = "//winfs-proj/proj/havgem/DIVA/syrekartor/data/ICES_1960_2022/"
file_path = "C:/Work/DIVAnd/Oxygen_maps/data/"
file_name = "Vattenkvalitet_all_SYKE.csv"

dtype_dict = {
    "dataSource": str,
    "site": str,
    "siteLatitudeWGS84": str,
    "siteLongitudeWGS84": str,
    #"mon/day/yr":str,
    "time":str,
    #"Year": str,
    #"Month": str,
    #"Day": str,
    #"Hour": str,
    #"Minute": str,
    "sampleDepthM": float,
    "analyteName": str,
    #"Salinity (PSALPR01_UUUU) [dmnless]": float,
    "value": float,
    #"Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]": float
}

# Load data from CSV file into a pandas DataFrame
#df = pd.read_csv(file_path+file_name, dtype=dtype_dict)
#df = pd.read_table(file_path+file_name, delimiter="\t",low_memory=False)
df = pd.read_table(file_path+file_name, delimiter=";",dtype=dtype_dict)

# Rename the '%Platform' column to 'ID'
#df.rename(columns={"Sampling platform (code)": "ID"}, inplace=True)

# Konvertera datumkolumnen till datetime
df['Year'] = df['time'].dt.year.astype(str)
# Extrahera månaden, dag, år
df['Month'] = df['time'].dt.month.astype(str)
df['Day'] = df['time'].dt.day.astype(str)
df['hh:mm'] = df['time'].dt.time.astype(str)

df["Visit event identifier"] = df["Year"] + '-' + df["dataSource"] + '-' + df["site"]

# Replace invalid hours with '00:00'
#df.loc[df["hh:mm"].isin(['-9', '24', '',np.nan]), "hh:mm"] = '00:00'
#df.loc[df["hh:mm"].isin(['-9', '60', '',np.nan]), "hh:mm"] = '00:00'
#df["Time"] = df["hh:mm"][0:1].str.zfill(2) + ":" + df["hh:mm"][3:4].str.zfill(2)

# Combine the date and hour columns into a new column 'date_time'
#df['date_time'] = df["Year"] + "-" + df["Month"].str.zfill(2) + "-" + df["Day"].str.zfill(2) + "T" + df["hh:mm"]

if row['Dissolved oxygen']:
    df['DOXY_ml'] = row['Dissolved oxygen'] * 0.7

# Make a new combined O2 variable.
# Base new variable on O2 bottle and add O2 CTD when O2 bottle is missing.
def get_doxy(row):
    if row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"] and not np.isnan(row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"]):
        return row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"] * -0.04488
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

# Define the order of columns for the output files
column_list = ["siteLongitudeWGS84", "siteLatitudeWGS84", "DOXY_umol", "sampleDepthM", "DOXY_umol", "DOXY_umol", "Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]", "Visit event identifier", "Year", "time", "Cruise", "Month", "Day", "hh:mm"]

# Write the filtered data to two output files, one with headers and one without
# df_filtered[column_list].to_csv(file_path+"bot.txt", index=False, sep='\t')
df_filtered[column_list].to_csv(file_path+"ICES_btl_lowres_ctd_02_NEW.txt", index=False, header=False, sep='\t')



