# ##Program to convert a IOW-file to a "big file" readable for DIVAnd
# ##IOW data from ODIN2 BTL and CTD
# ##

import pandas as pd
import numpy as np

# Define file path and data types for each column
file_path = "C:/Work/DIVAnd/Oxygen_maps/data/ODIN2/ODV/"
file_name = "data_from_Spreadsheet_Collection_2025-03-25T12-45-01.txt"

dtype_dict = {
    "Cruise": str,
    "Station": str,
    "Type": str,
    "yyyy-mm-ddThh:mm:ss.sss": str,
    "Longitude [degrees_east]": float,
    "Latitude [degrees_north]": float,
    "depth [m]": float,
    "DOXY4IDD [ml/l]": float,
    "DOXY6TID [ml/l]": float
}

# Load data from CSV file into a pandas DataFrame
df = pd.read_table(file_path+file_name, delimiter="\t",dtype=dtype_dict, comment="*")

df["Year"] = df["yyyy-mm-ddThh:mm:ss.sss"].str.slice(0, 4).astype(str)
df["Month"] = df["yyyy-mm-ddThh:mm:ss.sss"].str.slice(5, 7).astype(str)
df["Day"] = df["yyyy-mm-ddThh:mm:ss.sss"].str.slice(8, 10).astype(str)
df["hh:mm"] = df["yyyy-mm-ddThh:mm:ss.sss"].str.slice(11, 16).astype(str)

#df["datum"] = pd.to_datetime(df["yyyy-mm-ddThh:mm:ss.sss"], format="%Y-%m-%dT%H:%M:%S.%f")

# df["Year"] = pd.to_datetime(df["yyyy-mm-ddThh:mm:ss.sss"], errors="coerce").dt.year
# df["Month"] = pd.to_datetime(df["yyyy-mm-ddThh:mm:ss.sss"], errors="coerce").dt.month
# df["Day"] = pd.to_datetime(df["yyyy-mm-ddThh:mm:ss.sss"], errors="coerce").dt.day

# df['datum'] = pd.to_datetime(df["yyyy-mm-ddThh:mm:ss.sss"], format='%YYYY-%mm-%ddT%hh:%mm:%ss:%sss')
# # Extrahera månaden, dag, år
#df['Month'] = df['datum'].dt.month.astype(str)
#df['Day'] = df['datum'].dt.day.astype(str)
#df['Year'] = df['datum'].dt.year.astype(str)

df["Visit event identifier"] = df["yyyy-mm-ddThh:mm:ss.sss"] + '-' + df["Cruise"] + '-' + df["Station"]

# Make a new combined O2 variable.
# Base new variable on O2 bottle and add O2 CTD when O2 bottle is missing.
def get_doxy(row):
     if row['DOXY6TID [ml/l]'] and not np.isnan(row['DOXY6TID [ml/l]']):
         return row['DOXY6TID [ml/l]']
         #return row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"] * 0 + 0.01
     if row['DOXY4IDD [ml/l]'] and not np.isnan(row['DOXY4IDD [ml/l]']):
         return row['DOXY4IDD [ml/l]']
     return np.nan
#df["DOXY_ml"] = df["DOXY4IDD [ml/l]"].fillna(df["DOXY6TID [ml/l]"])
df['DOXY_ml'] = df.apply(get_doxy, axis=1)

# Convert O2 ml/l to µmol/l-> 1ml/l = 44.661 µmol/l
#df['DOXY_umol'] = df['DOXY_ml'].multiply(44.661)
df['DOXY_umol'] = df["DOXY_ml"].apply(lambda x: x*44.661)
print(df.dtypes)
# Filter out rows with missing data
#df_filtered = df.loc[df["DOXY_umol"] != np.nan]
df_filtered = df.dropna(subset=["DOXY_umol"])


# Define the order of columns for the output files
#column_list = ["Longdeg", "Latdeg", "DOXY(umol/l)", "Pressure(Dbars)", "TEMP", "PSAL", "H2SX(umol/l)", "St_No", "Year", "date_time", "ID", "Mnth", "Dy", "Hr"]
#column_list = ["Sample longitude (DD)", "Sample latitude (DD)", "DOXY_umol", "Sampling depth (m)", "Temperature bottle (C)", "Salinity bottle (o/oo psu)", "Hydrogen sulphide H2S (umol/l)", "Visit event identifier", "Year", "date_time", "Sampling platform (code)", "Month", "Day", "Sampling time (start)"]
column_list = ["Longitude [degrees_east]", "Latitude [degrees_north]", "DOXY_umol", "depth [m]", "depth [m]", "depth [m]", "DOXY_ml", "Cruise", "Year", "yyyy-mm-ddThh:mm:ss.sss", "Visit event identifier", "Month", "Day", "hh:mm"]

# Write the filtered data to two output files, one with headers and one without
# df_filtered[column_list].to_csv(file_path+"bot.txt", index=False, sep='\t')
df_filtered[column_list].to_csv(file_path+"IOW_1959_2024.txt", index=False, header=False, sep='\t')



