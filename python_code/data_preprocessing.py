import pandas as pd
import polars as pl
import numpy as np
from duplicate_check import find_duplicate_candidates, remove_duplicates
from data_reader import load_ices_btl_data, load_ices_ctd_data, load_emodnet_btl_data, load_emodnet_ctd_data
import os
from pathlib import Path

# order of columns in the outputfile have to be:
# ['lon', 'lat', 'value', 'depth', 'date', 'id']
# 0 lom
# 1 lat
# 2 value
# 3 depth
# index 4 to 8 can contain any data, but cannot be left empty
# 9 date 
# 10 id
# 
# date column have to be a string that har the format "%Y-%m-%dT%H:%M"

def range_check(df, column, min_val, max_val):
    """Removes values outside the specified range."""
    return df[(df[column] >= min_val) & (df[column] <= max_val)]

def set_neg_depth_to(df, value):
    "Sets all negative depth to value "
    df["depth"] = df["depth"].astype(float)
    df.loc[df["depth"] < 0, "depth"] = value
    return df

def set_neg_oxyg_to(df, near_zero):
    "Sets all negative oxygen to value"
    df["value"] = df["value"].astype(float)
    df.loc[df["value"] < 0, "value"] = near_zero
    return df

def remove_shallow_low_ox_values(df, column, shallow_ox_threshold = 180, shallow_depth = 20):
    """Removes values shallower than given depth and with values below shallow_ox_threshold"""
    return df[~((df[column] <= shallow_ox_threshold) & (df["depth"] <= shallow_depth))]

# Function to downscale the precision of timestamps to the least precise one
def downscale_to_least_precise(ts1, ts2):
    # Compare the precision (e.g., 'ns' for nanosecond precision)
    if ts1.microsecond == 0 and ts2.microsecond == 0:
        ts1 = ts1.normalize()  # Remove time (downscale to day precision)
        ts2 = ts2.normalize()  # Remove time (downscale to day precision)
    elif ts1.minute == 0 and ts2.minute == 0:
        ts1 = ts1.replace(second=0, microsecond=0)  # Remove seconds and microseconds
        ts2 = ts2.replace(second=0, microsecond=0)  # Remove seconds and microseconds
    return ts1, ts2

def remove_matching_rows(df, bad_data, column, log_file):
    """Removes rows in df that match rows in bad_data based on id column and logs removed rows."""
    removed_rows = df[df['id'].isin(bad_data['id'])]
    cleaned_df = df[~df['id'].isin(bad_data['id'])]

    if not removed_rows.empty:
        with open(log_file, "a") as f:
            removed_rows.to_csv(f, sep="\t", index=False, header=f.tell()==0)
    
    return cleaned_df

def format_output(df):
    """Ensures the output dataset follows the specified 11-column format."""
    df['date'] = pd.to_datetime(df['date']).dt.strftime('%Y-%m-%dT%H:%M')
    filler_cols = pd.DataFrame(np.nan, index=df.index, columns=['dummy_4', 'dummy_5', 'dummy_6', 'dummy_7', 'dummy_8'])
    df = pd.concat([df[['lon', 'lat', 'value', 'depth']], filler_cols, df[['date', 'id']]], axis=1)
 

def load_iow_data(file_path):
    """Loads IOW CTD/BTL dataset with specific format.
        DOXY6TID [ml/l] = BTL (inkl neg o2)
        DOXY4IDD [ml/l] = CTD
    """
    df = pd.read_csv(file_path, sep="\t", skiprows=20)
    df = df[['Longitude [degrees_east]', 'Latitude [degrees_north]', 'DOXY4IDD [ml/l]', 'DOXY6TID [ml/l]', 'depth [m]', 'yyyy-mm-ddThh:mm:ss.sss', 'Cruise', 'Station','Type']]

    df['id'] = "iow-" + df['yyyy-mm-ddThh:mm:ss.sss'].astype(str) + "-" + df['Station'].astype(str)

    def get_doxy(row):
        if row['DOXY6TID [ml/l]'] and not np.isnan(row['DOXY6TID [ml/l]']):
            return row['DOXY6TID [ml/l]']
        if row['DOXY4IDD [ml/l]'] and not np.isnan(row['DOXY4IDD [ml/l]']):
            return row['DOXY4IDD [ml/l]']
        return np.nan

    df['Water body dissolved oxygen concentration [ml/l]'] = df.apply(get_doxy, axis=1)
    # Change unit from ml/l to µmol/l
    df['Water body dissolved oxygen concentration [µmol/l]'] = df['Water body dissolved oxygen concentration [ml/l]'].apply(lambda x: x * 44.661)

    df.rename(columns=
        {
            "Longitude [degrees_east]": "lon",
            "Latitude [degrees_north]": "lat",
            "Water body dissolved oxygen concentration [µmol/l]": "value",
            "depth [m]": "depth",
            "yyyy-mm-ddThh:mm:ss.sss": "date",
        }, inplace=True
    )
    df['date'] = df['date'].str[:16]  # Trim milliseconds

    return df

def load_shark_data(file_path):
    """Loads shark tab-delimited dataset."""
    cols = ['depth', 'Dissolved oxygen O2 bottle (ml/l)', 'Q-flag Dissolved oxygen O2 bottle',	'Dissolved oxygen O2 CTD (ml/l)', 'Q-flag Dissolved oxygen O2 CTD',	'Hydrogen sulphide H2S (umol/l)', 'Q-flag Hydrogen sulphide H2S', 'DOXY_ml', 'value', 'NEG_DOXY_umol']

    df = pd.read_csv(file_path, sep="\t", names=cols)
    df = df[['lon', 'lat', 'value', 'depth', 'date', 'id']]
    return df

def load_big_file_format(file_path):
    """Loads generic tab-delimited dataset."""
    cols = ['lon', 'lat', 'value', 'depth', 'drop1', 'drop2', 'drop3', 'drop4', 'drop5', 'date', 'id']
    df = pd.read_csv(file_path, sep="\t", header=None)
    df.columns = cols + df.columns[len(cols):].tolist()
    df = df[['lon', 'lat', 'value', 'depth', 'date', 'id']]
    
    return df

def process_dataset(file_name, is_emodnet_btl=False, is_emodnet_ctd=False, is_iow=False, is_ices_btl=False, is_ices_ctd=False):
    if is_emodnet_btl:
        df = load_emodnet_btl_data(file_name)
    elif is_emodnet_ctd:
        df = load_emodnet_ctd_data(file_name)
    elif is_iow:
        df = load_iow_data(file_name)
    elif is_ices_btl:
        df = load_ices_btl_data(file_name)
    elif is_ices_ctd:
        df = load_ices_ctd_data(file_name)
    else: 
        df = load_big_file_format(file_name)
    
    df['date'] = pd.to_datetime(df['date'], format="mixed")
    df = range_check(df, 'value', -500, 600)
    """No the below mess up the Kattegatt analysis /MH"""
    #df = remove_shallow_low_ox_values(df, 'value')

    """Lägg till check på att depth inte är neg. Sätt neg till noll. """
    df = set_neg_depth_to(df, 0)

    """Lägg till check på att alla negativa värden sätts till ngt lågt...
    obsval[obsval .<= 0] .= 0.44661"""
    df = set_neg_oxyg_to(df, 0.44661)

    if not df.loc[~df.index.isin(df.dropna(subset=['lon', 'lat', 'value', 'depth','date']).index)].empty:
        print("Warning")
        print("\tFound nan values in one or more of the required columns:")
        print(f"\t{df.loc[~df.index.isin(df.dropna(subset=['lon', 'lat', 'value', 'depth','date']).index)]}")
        print("\tRemoving these from the dataset")
        df.dropna(subset=['lon', 'lat', 'value', 'depth','date'], inplace=True)

    return df

def main():
    data_dir = Path("data")
        
    freja_data_dir = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/data/all_baltic"
    
    print(freja_data_dir)
    on_freja = False
    if Path(freja_data_dir).is_dir():
        data_dir = freja_data_dir
        on_freja = True
    print(f"load data from {data_dir}")
    # output_dir = os.path.join(data_dir, "processed/")
    # os.makedirs(output_dir, exist_ok=True)
    # log_file = os.path.join(output_dir, "removed_rows_log.txt")
        
    # Load and process EMODNET and ICES BTL & CTD datasets
    print("Load EMODnet btl & ctd....")
    emodnet_btl_df = process_dataset(os.path.join(data_dir, "BTL_data_from_EMODnet_Eutrophication_European_2024_unrestricted.txt"), is_emodnet_btl=True)
    emodnet_ctd_df = process_dataset(os.path.join(data_dir, "CTD_data_from_EMODnet_Eutrophication_European_2024_unrestricted.txt"), is_emodnet_ctd=True)

    # print("Load IOW data....")
    # iow_df = process_dataset(os.path.join(data_dir, "IOW_data_from_Spreadsheet_Collection_2025-03-25T12-45-01.txt"), is_iow=True)
    print("Load ICES BTL data....")
    ices_btl_df = process_dataset(os.path.join(data_dir, "ce850252-8a93-49c9-a24a-61d67bef48fe.csv"), is_ices_btl=True)
    # print("Load ICES CTD data....")
    ices_ctd_df = process_dataset(os.path.join(data_dir, "1c9c37d5-2236-42a8-8755-f553c1f4fa0f.csv"), is_ices_ctd=True)
    # print("Load SYKE data...")
    # syke_df = process_dataset(os.path.join(data_dir, "syke_data_no_header_241107.txt"))
    ## SHARK problem, data has no lat, lom, move shark preprocessing from excel_to_tact_sharkweb to this code
    print("Load SHARK data...")
    shark_df = process_dataset(os.path.join(data_dir, "sharkweb_btlctd_02_241107.txt"))

    _, ices_btl_duplicates = find_duplicate_candidates(shark_df, ices_btl_df, suffix="ices_btl")

    with pl.Config(tbl_cols=-1):
        print(ices_btl_duplicates.head())
        print(ices_btl_duplicates.height)
    
    # remove duplicates in ices df
    print("remove shark duplicates from ices btl")
    cleaned_ices_btl = remove_duplicates(shark_df, ices_btl_df)
    print("remove shark duplicates from ices ctd")
    cleaned_ices_ctd = remove_duplicates(shark_df, ices_ctd_df)

    print("remove shark duplicates from emdonet ctd")
    cleaned_emodnet_ctd = remove_duplicates(shark_df, emodnet_ctd_df)
    print("remove shark duplicates from emodnet btl")
    cleaned_emodnet_btl = remove_duplicates(shark_df, emodnet_btl_df)


    _, ices_ices_duplicates = find_duplicate_candidates(cleaned_ices_ctd, cleaned_ices_btl)
    with pl.Config(tbl_cols=-1):
        print(ices_ices_duplicates.head())
        print(f"{ices_ices_duplicates.height} duplicates between ices ctd and btl")

    _, emodnet_emodnet_duplicates = find_duplicate_candidates(cleaned_emodnet_ctd, cleaned_emodnet_btl)
    with pl.Config(tbl_cols=-1):
        print(emodnet_emodnet_duplicates.head())
        print(f"{emodnet_emodnet_duplicates.height} duplicates between emodnet ctd and btl")
    
"""    
    print("\nCheck_duplicates for iow vs emodnet/ices/syke/shark")
    merged_emodnet_ices_syke_shark_iow, emodnet_ices_syke_shark_iow_duplicates = find_duplicate_candidates(iow_df, merged_emodnet_ices_syke_shark, time_threshold=4, vertical_threshold=2, oxygen_value_threshold=50, name1='iow', name2='emodnet_ices_syke_shark')
    merged_emodnet_ices_syke_shark_iow.to_csv(os.path.join(output_dir, "merged_emodnet_ices_syke_shark_iow.txt"), sep="\t", index=False)
    emodnet_ices_syke_shark_iow_duplicates.drop_duplicates(subset=['lon', 'lat', 'date', 'source']).to_csv(os.path.join(output_dir, "emodnet_ices_syke_shark_iow_duplicatecheck_emodneticessykeshark_removed.txt"), sep="\t",index=False)

    formatted_df = format_output(merged_emodnet_ices_syke_shark_iow)
    formatted_df.to_csv(os.path.join(output_dir, "duplicate_checked.txt"), sep="\t", index=False)
    
    # formatted_df = load_data(os.path.join(output_dir, "merged_cleaned.txt"))
    # Load bad data
    bad_data_df = pd.read_csv(os.path.join(data_dir, "bad_data.txt"), sep = '\t')
    bad_data_df['date'] = pd.to_datetime(bad_data_df['obstime'], format="mixed")
    bad_data_df['id'] = bad_data_df['obsid']
    bad_data_df['lon'] = bad_data_df['obslon']
    bad_data_df['lat'] = bad_data_df['obslat']
    bad_data_df['value'] = bad_data_df['obsvalue']
    bad_data_df['depth'] = bad_data_df['obsdepth']

    # Remove bad data after merging and log removed rows
    
    cleaned_df = remove_matching_rows(formatted_df, bad_data_df, ['lon', 'lat', 'depth', 'date', 'value', 'id'], log_file)
    # Format output
    df = format_output(cleaned_df)

    # Save cleaned dataset
    df.to_csv(os.path.join(output_dir, "cleaned_df_using_id_col.txt"), sep="\t", index=False)

    cleaned_df_by_dupl_check, bad_data_removed = find_duplicate_candidates(formatted_df, bad_data_df,
                     suffix ='bad_data')
    
    # Format output
    formatted_df = format_output(cleaned_df_by_dupl_check)

    # Save cleaned dataset
    formatted_df.to_csv(os.path.join(output_dir, "cleaned_df_using_dupl_250619.txt"), sep="\t", index=False)
    bad_data_removed.to_csv(os.path.join(output_dir, "bad_data_removed_using_dupl.txt"), sep="\t", index=False)"""
    
if __name__ == "__main__":
    main()
