# ##Program to convert a sharkwebexport to a "big file" readable for DIVAnd
# ##Export with decoding: UTF-8
# ##

import pandas as pd
import numpy as np

# Define file path and data types for each column
#file_path = "C:/Work/DIVAnd/Oxygen_maps/data/all_baltic/"
file_path = "/home/sm_mahan/DIVAnd/Oxygen_maps/data/all_baltic/Original_data/"
file_name = "sharkweb_data_1960_2024_utf8.txt"

# Load data from CSV file into a pandas DataFrame
#df = pd.read_csv(file_path+file_name, dtype=dtype_dict)
df = pd.read_table(file_path+file_name, delimiter="\t")

# Set flag columns from object to str
df["Q-flag Dissolved oxygen O2 bottle"] = df["Q-flag Dissolved oxygen O2 bottle"].astype(str)
df["Q-flag Dissolved oxygen O2 CTD"] = df["Q-flag Dissolved oxygen O2 CTD"].astype(str)
df["Q-flag Hydrogen sulphide H2S"] = df["Q-flag Hydrogen sulphide H2S"].astype(str)

# Rename the '%Platform' column to 'ID'
#df.rename(columns={"Sampling platform (code)": "ID"}, inplace=True)

# Replace invalid hours with '0000'
df.loc[df["Sampling time (start)"].isin(['-9', '2400', '',np.nan]), "Sampling time (start)"] = '00:00'

# Get day from "Sampling date"
df['Day'] = df["Sampling date"].apply(lambda x: x.split("-")[2])

# Combine the date and hour columns into a new column 'date_time'
df['date_time'] = df["Sampling date"].astype(str) + "T" + df["Sampling time (start)"]
#df['date_time'] = pd.to_datetime(df["Year"] + "-" + df["Mnth"] + "-" + df["Dy"] + " " + df["Hr"].str.slice(0,2)).dt.strftime("%Y-%m-%dT%H:%M")
# df['date'] = df["Year"] + "-" + df["Mnth"] + "-" + df["Dy"] + " " + df["Hr"].str.slice(0,2) #+ ":00"
# # format date column to YYYY-MM-DDTHH:MM
# df['date_time'] = pd.to_datetime(df["date"]).dt.strftime("%Y-%m-%dT%H:%M")

# Make a new combined O2 variable, o2 btl, h2s, o2 ctd that gives only + o2: H2s = 0.01 ml/l o2
# Then make a neg o2 where H2S is calculated to negativ oxygen

def DIVA_oxygen(data: pd.DataFrame):

    valid_btl = np.logical_and(~pd.isna(data["Dissolved oxygen O2 bottle (ml/l)"]), ~data["Q-flag Dissolved oxygen O2 bottle"].str.contains("B|S|<"))
    below_det_btl = np.logical_and(~pd.isna(data['Dissolved oxygen O2 bottle (ml/l)']), data['Q-flag Dissolved oxygen O2 bottle'].str.contains("<"))
    #z_btl = np.logical_and(~pd.isna(data['Dissolved oxygen O2 bottle (ml/l)']), data['Q-flag Dissolved oxygen O2 bottle'].str.contains("Z"))
    
    valid_ctd = np.logical_and(~pd.isna(data['Dissolved oxygen O2 CTD (ml/l)']), ~data['Q-flag Dissolved oxygen O2 CTD'].str.contains("B|S|<"))
    below_det_ctd = np.logical_and(~pd.isna(data['Dissolved oxygen O2 CTD (ml/l)']), data['Q-flag Dissolved oxygen O2 CTD'].str.contains("<"))
    #z_ctd = np.logical_and(~pd.isna(data['Dissolved oxygen O2 CTD (ml/l)']), data['Q-flag Dissolved oxygen O2 CTD'].str.contains("Z"))
    
    valid_h2s = np.logical_and(~pd.isna(data["Hydrogen sulphide H2S (umol/l)"]), ~data["Q-flag Hydrogen sulphide H2S"].str.contains("B|S|Z|<"))
    below_det_h2s = np.logical_and(~pd.isna(data["Hydrogen sulphide H2S (umol/l)"]), data["Q-flag Hydrogen sulphide H2S"].str.contains("<"))
    #above_det_h2s = np.logical_and(~pd.isna(data["Hydrogen sulphide H2S (umol/l)"]), ~data["Q-flag Hydrogen sulphide H2S"].str.contains(">"))
    z_h2s = np.logical_and(~pd.isna(data["Hydrogen sulphide H2S (umol/l)"]), data["Q-flag Hydrogen sulphide H2S"].str.contains("Z"))
    zero_h2s = np.logical_and((data["Hydrogen sulphide H2S (umol/l)"] == 0),~data["Q-flag Hydrogen sulphide H2S"].str.contains("B|S|Z|<|>"))

    # Apply nested np.where for all conditions
    data["DOXY_ml"] = (
        #h2s below_det_h2s -> h2s (NaN)
        np.where(((below_det_h2s)|(z_h2s)) & ((~valid_btl) & (~valid_ctd)),np.nan,
            #h2s == 0 and valid_btl -> valid btl
            np.where((zero_h2s) & (valid_btl), data['Dissolved oxygen O2 bottle (ml/l)'],
                #h2s == 0 and valid_ctd -> valid ctd
                np.where((zero_h2s) & (valid_ctd), data['Dissolved oxygen O2 CTD (ml/l)'],          
                    # h2s == 0 and (btl and ctd not valid) -> 0.01 Vi säker helt enkelt att om noll är uppmätt så är det noll. 
                    np.where((zero_h2s) & ((~valid_btl) & (~valid_ctd)), 0.01,
                        # h2s valid-> h2s default (0.01)
                        np.where(valid_h2s, 0.01,
                            # h2s not valid and o2< gives h2s default (0.01)
                            np.where((~valid_h2s) & (below_det_btl), 0.01,
                                #  o2 BTL is valid gives o2 BTL
                                np.where((valid_btl), data['Dissolved oxygen O2 bottle (ml/l)'],
                                    # O2 CTD exists and O2 CTD is not S gives O2 CTD
                                    np.where((valid_ctd),data['Dissolved oxygen O2 CTD (ml/l)'],
                                        np.where((below_det_ctd),0.01,
                                        # Default case gives NaN
                                        np.nan
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )    
        )
    )    
    return data

def DIVA_neg_oxygen(data: pd.DataFrame):

    valid_btl = np.logical_and(~pd.isna(data["Dissolved oxygen O2 bottle (ml/l)"]), ~data["Q-flag Dissolved oxygen O2 bottle"].str.contains("B|S|<"))
    below_det_btl = np.logical_and(~pd.isna(data['Dissolved oxygen O2 bottle (ml/l)']), data['Q-flag Dissolved oxygen O2 bottle'].str.contains("<"))
    #z_btl = np.logical_and(~pd.isna(data['Dissolved oxygen O2 bottle (ml/l)']), data['Q-flag Dissolved oxygen O2 bottle'].str.contains("Z"))
    
    valid_ctd = np.logical_and(~pd.isna(data['Dissolved oxygen O2 CTD (ml/l)']), ~data['Q-flag Dissolved oxygen O2 CTD'].str.contains("B|S|<"))
    below_det_ctd = np.logical_and(~pd.isna(data['Dissolved oxygen O2 CTD (ml/l)']), data['Q-flag Dissolved oxygen O2 CTD'].str.contains("<"))
    #z_ctd = np.logical_and(~pd.isna(data['Dissolved oxygen O2 CTD (ml/l)']), data['Q-flag Dissolved oxygen O2 CTD'].str.contains("Z"))
    
    valid_h2s = np.logical_and(~pd.isna(data["Hydrogen sulphide H2S (umol/l)"]), ~data["Q-flag Hydrogen sulphide H2S"].str.contains("B|S|Z|<"))
    below_det_h2s = np.logical_and(~pd.isna(data["Hydrogen sulphide H2S (umol/l)"]), data["Q-flag Hydrogen sulphide H2S"].str.contains("<"))
    #above_det_h2s = np.logical_and(~pd.isna(data["Hydrogen sulphide H2S (umol/l)"]), ~data["Q-flag Hydrogen sulphide H2S"].str.contains(">"))
    z_h2s = np.logical_and(~pd.isna(data["Hydrogen sulphide H2S (umol/l)"]), data["Q-flag Hydrogen sulphide H2S"].str.contains("Z"))
    zero_h2s = np.logical_and((data["Hydrogen sulphide H2S (umol/l)"] == 0),~data["Q-flag Hydrogen sulphide H2S"].str.contains("B|S|Z|<|>"))

    # Apply nested np.where for all conditions
    data["NEG_DOXY_ml"] = (
        #h2s below_det_h2s -> h2s (NaN)
        np.where(((below_det_h2s)|(z_h2s)) & ((~valid_btl) & (~valid_ctd)),np.nan,
            #h2s == 0 and valid_btl -> valid btl
            np.where((zero_h2s) & (valid_btl), data['Dissolved oxygen O2 bottle (ml/l)'],
                #h2s == 0 and valid_ctd -> valid ctd
                np.where((zero_h2s) & (valid_ctd), data['Dissolved oxygen O2 CTD (ml/l)'],          
                    # h2s == 0 and (btl and ctd not valid) -> 0.01 Vi säker helt enkelt att om noll är uppmätt så är det noll. 
                    np.where((zero_h2s) & ((~valid_btl) & (~valid_ctd)), 0.01,
                        # h2s valid-> h2s default (0.01)
                        np.where(valid_h2s, data["Hydrogen sulphide H2S (umol/l)"] * -0.04488,
                            # h2s not valid and o2< gives h2s default (0.01)
                            np.where((~valid_h2s) & (below_det_btl), 0.01,
                                #  o2 BTL is valid gives o2 BTL
                                np.where((valid_btl), data['Dissolved oxygen O2 bottle (ml/l)'],
                                    # O2 CTD exists and O2 CTD is not S gives O2 CTD
                                    np.where((valid_ctd),data['Dissolved oxygen O2 CTD (ml/l)'],
                                        np.where((below_det_ctd),0.01,
                                        # Default case gives NaN
                                        np.nan
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )    
        )
    )
    return data

# Call functions for neg o2 and o2
DIVA_neg_oxygen(df)
DIVA_oxygen(df)

# Convert O2 ml/l to µmol/l-> 1ml/l = 44.661 µmol/l
#df['DOXY_umol'] = df['DOXY_ml'].multiply(44.661)
df['NEG_DOXY_umol'] = df['NEG_DOXY_ml'].apply(lambda x: x*44.661)
df['DOXY_umol'] = df['DOXY_ml'].apply(lambda x: x*44.661)

# Filter out rows with missing data
df_filtered = df.dropna(subset=["DOXY_umol"])
df_filtered = df_filtered.replace("nan", "")


df_filtered["SHARK_ID"] = "SHARK-" + df_filtered["Year"].astype(str) +'-'+ df_filtered["Sampling platform (code)"].astype(str) +'-'+ df_filtered["Visit event identifier"].astype(str)

# Define the order of columns for the output files

column_list = ["Sample longitude (DD)", "Sample latitude (DD)", "NEG_DOXY_umol", "Sampling depth (m)", "DOXY_umol", "Salinity bottle (o/oo psu)", "Hydrogen sulphide H2S (umol/l)", "Visit event identifier", "Year", "date_time", "SHARK_ID", "Month", "Day", "Sampling time (start)"]
column_list_test =["Sampling depth (m)","Dissolved oxygen O2 bottle (ml/l)", "Q-flag Dissolved oxygen O2 bottle","Dissolved oxygen O2 CTD (ml/l)","Q-flag Dissolved oxygen O2 CTD","Hydrogen sulphide H2S (umol/l)","Q-flag Hydrogen sulphide H2S","DOXY_ml","DOXY_umol","NEG_DOXY_umol"]

# Write the filtered data to two output files, one with headers and one without
df_filtered[column_list].to_csv(file_path+"sharkweb_btlctd_02_251128.txt", index=False, header=False, sep='\t')
df_filtered[column_list_test].to_csv(file_path+"sharkweb_btlctd_02_251128_test.txt", index=False, header=True, sep='\t')


