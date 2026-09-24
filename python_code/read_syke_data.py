from pathlib import Path
import pandas as pd
import numpy as np

"""The flag explanations are as follows:
C= Confirmed result that is under or exceeds a set “alarm” 
L= Value that is below the detection limit
LC= Confirmed value that is under the detection and “alarm” limit
W= uncertainty is higher than normal (I usually exclude these)

The method codes are as follows:

NA= method not available (only old data)
EL= electrometric or ion-selective determination
ELF= electrometric or ion-selective determination, in field (CTD measurements have this code)
TI= titrimetric, potentiometric
TIJ= titrimetric, in field

"""

def transform_syke_data(path):
    print(f"reading file {path}\n... ... ...")
    df = pd.read_csv(
        path, encoding='utf-8',
        sep=",",
    )
    print("done reading file")
    org_nans = df.value.notna()
    df['value'] = pd.to_numeric(df['value'], errors='coerce')

    new_nan = (
        org_nans &
        df['value'].isna()
    )

    print(f"Antal värden som inte kunde konverteras till numeric: {new_nan.sum()}")

    # Filtrera bort rader där PARAM inte är O2D eller O2S
    data = df.query(
        'parameter_code in ["O2D", "H2SS"]')
    print(f"dropped {len(df)-len(data)} rows with other parameters then O2D and H2SS")
    print(f"size of df decreased by {100*(len(df)-len(data))/len(df):.0f}%")
    print(f"visits in original df: {df.groupby(['time', 'wgs84_lat', 'wgs84_long', 'site_id']).ngroups}")
    print(f"visits left: {data.groupby(['time', 'wgs84_lat', 'wgs84_long', 'site_id']).ngroups}")
    print("visits missing O2D or H2SS:\n"
          f"          {df.groupby(['time', 'wgs84_lat', 'wgs84_long', 'site_id']).ngroups - data.groupby(['time', 'wgs84_lat', 'wgs84_long', 'site_id']).ngroups}")

    # Skapa en ny kolumn för att lagra OXYGEN-beräkningar
    # döp om lite granna
    rename_dict = {
            "site_id": "SERNO",
            "site_depth": "WADEP",
            "wgs84_lat": "LATIT", 
            "wgs84_long": "LONGI",
            "depth_upper": "DEPH",
            "time": "SDATE",
            "site": "STATN",
            "flag": "Q_flag",
        }
    # Byt namn på de kolumner som finns i rename_dict
    data.rename(columns=rename_dict, inplace = True)
    data["ID"] = 'SYKE-' + data['SERNO'].astype(str) + '-' + data['SDATE'].astype(str)

    visit_cols = ['SDATE', 'LATIT', 'LONGI', 'ID']

    # kontrollera om det finns flera mätningar på samma parameter och metod för en visit
    duplicates = (
        data.groupby(
            visit_cols + ['DEPH', 'parameter_code', 'method_code'],
            dropna=False
        )['value']
        .nunique(dropna=False)
    )
    problematic = duplicates[duplicates > 1]

    if not problematic.empty:
        problematic_rows = data.set_index(visit_cols + ['DEPH', 'parameter_code', 'method_code']).index.isin(problematic.index)
        problematic_groups = problematic.groupby(visit_cols + ['DEPH', 'parameter_code', 'method_code'])

        print(f"Number of problematic groups: {problematic_groups.ngroups}")

        if problematic_groups.ngroups <= 50:
            print("\nProblematic index:")
            print(problematic.index)
            for name, group in data[problematic_rows].groupby(visit_cols + ['DEPH', 'parameter_code', 'method_code']):
                print(f"group of {visit_cols + ['DEPH', 'parameter_code', 'method_code']}\n {name}")
                print(group[visit_cols + ['DEPH', 'parameter_code', 'method_code', 'value', 'Q_flag']])

            data = data.loc[~problematic_rows]
            print(f"Removed {problematic_groups.ngroups} depths with >1 measurement.")
        else:
            print("\nProblematic group keys:")
            print(problematic.index.to_frame(index=False).to_string(index=False))
            raise ValueError(
                "Different values found for > 100 groups of"
                "visit/depth/parameter/method combination:\n"
                f"{problematic_rows[visit_cols + ['DEPH', 'parameter_code', 'method_code', 'value']]}"
            )

    group_cols = ['STATN', 'SDATE', 'LATIT', 'LONGI', "DEPH", "ID"] 
    parameter_cols = [
            'value',
            'Q_flag',
            'method_code'
        ]

    wide = (
        data.pivot_table(
            index=group_cols,
            columns='parameter_code',
            values=parameter_cols,
            aggfunc='first'
        )
    )
    wide.columns = [
        f'{parameter}_{column}'
        for column, parameter in wide.columns
    ]

    wide = wide.reset_index()
    data['value'] = data['value'].astype(float)

    # set to negative oxygen equivalents, µg/l till µmol/l (tror vi)
    wide['OXYGEN_H2SS'] = (
        wide['H2SS_value'] * -0.029342
        )

    # change unit from mg/l to umol/l
    wide['OXYGEN_O2D'] = (
        wide['O2D_value'] * 0.700 * 44.661
    )

    # prepare the dataframe with columns for the combined oxygen parameter
    wide['OXYGEN'] = np.nan
    wide['OXYGEN_Qflag'] = ""

    # sätt värden till kolumnerna utefter följand villkor:

    # 1. H2S finns uppmätt och Q_flag för H2S != L. Använd H2S omräknat till neg syre
    condition = (
        (wide['OXYGEN_H2SS'] < 0) &
        (wide['H2SS_Q_flag'] != 'L')
    )

    wide.loc[condition, 'OXYGEN'] = wide.loc[condition, 'OXYGEN_H2SS']

    # 2. O2D method_code == TI alltså alla flaskdata utan flagga. Använd flaskdatan och spara flaggan för denna
    condition = (
        wide['OXYGEN'].isna() &
        (wide['O2D_method_code'] == 'TI')
    )

    wide.loc[condition, 'OXYGEN'] = wide.loc[condition, 'OXYGEN_O2D']
    wide.loc[condition, 'OXYGEN_Q_flag'] = wide.loc[condition, 'O2D_Q_flag']

    # 3. O2D method_code == ELF. alltså alla ctd-data utan flagga. Fyll på med ctd data där det fortfarande saknas värden.
    condition = (
        wide['OXYGEN'].isna() &
        (wide['O2D_method_code'] == 'ELF')
    )

    wide.loc[condition, 'OXYGEN'] = wide.loc[condition, 'OXYGEN_O2D']
    wide.loc[condition, 'OXYGEN_Q_flag'] = wide.loc[condition, 'O2D_Q_flag']

    # 4. O2D method_code == EL. Fyll på med CTD data med metodkod EL där det fortfarande saknas värden
    condition = (
        wide['OXYGEN'].isna() &
        (wide['O2D_method_code'] == 'EL')
    )

    wide.loc[condition, 'OXYGEN'] = wide.loc[condition, 'OXYGEN_O2D']
    wide.loc[condition, 'OXYGEN_Q_flag'] = wide.loc[condition, 'O2D_Q_flag']

    # 5. O2D method_code is NaN. Fyll på med syre med okänd metodkod
    condition = (
        wide['OXYGEN'].isna() &
        wide['O2D_method_code'].isna()
    )

    wide.loc[condition, 'OXYGEN'] = wide.loc[condition, 'OXYGEN_O2D']
    wide.loc[condition, 'OXYGEN_Q_flag'] = wide.loc[condition, 'O2D_Q_flag']

    # 6. Data in H2SS column but not in O2D column, set OXYGEN to zero
    condition = (
            wide['OXYGEN'].isna() &
            wide['OXYGEN_O2D'].isna() &
            wide['OXYGEN_H2SS'].notna() 
        )

    wide.loc[condition, 'OXYGEN'] = 0

    # set all L flagged to zero
    condition = (
            wide['OXYGEN_Q_flag'] == "L"
        )
    wide.loc[condition, 'OXYGEN'] = 0

    # data.drop(data[(data['Q_flag'] != "L")].index, inplace=True)
    wide.drop(wide[(wide['OXYGEN'] < -800)].index, inplace=True)
    wide.drop(wide[(wide['OXYGEN'] > 600)].index, inplace=True)
    # data = data.dropna(subset=['OXYGEN'])
    headers = [
        "SDATE",
        "ID",
        "STATN",
        "LATIT",
        "LONGI",
        "DEPH",
    ]
    wide.sort_values(by=headers, inplace=True)

    # Definiera ordningen på kolumnerna för loadbigfile
    column_order = [
        'LONGI',  # 0
        'LATIT',  # 1
        'OXYGEN',    # 2
        'DEPH',   # 3
        'OXYGEN_H2SS', 'OXYGEN_O2D', 'STATN',  'H2SS_value', 'OXYGEN_Q_flag',  # 4 till 8
        'SDATE',  # 9
        'ID'   # 10
    ]

    # Om kolumner 4-8 är tomma, fyll dem med dummyvärden (här är alla fyllda)
    for col in column_order[4:9]:
        if col not in wide or wide[col].isnull().all():  # Kontrollera om kolumnen saknas eller är helt tom
            wide[col] = ''  # Fyll med en dummy-sträng eller annat värde

    # Ordna kolumnerna enligt ordern vi har specificerat
    # skriv datumsträngen så som loadbigfile vill ha
    wide["SDATE"] = wide["SDATE"].apply(
        lambda row: row.replace(" ", "T") if isinstance(row, str) else row
    )
    oxygen = wide[column_order]  
        
    return oxygen
    

path = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/data/all_baltic/Original_data/SYKE/"
open_sea = transform_syke_data(path=Path(path, "Syke_oxygen_opensea.csv"))
coastal = transform_syke_data(path=Path(path, "Syke_oxygen_coastal.csv"))
data = pd.concat([open_sea, coastal])
data.to_csv(Path("/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/data/all_baltic/", "syke_data_no_header_260924.txt"), sep="\t", index = False, header = False)

# parameter_codes in data
# ['TEMP' 'O2D' 'SAL' 'O2S' 'H2SS']
# order of columns in the outputfile have to be:
# 0 LONGI
# 1 LATIT
# 2 O2D
# 3 DEPH
# index 4 to 8 can contain any data, but cannot be left empty
# 9 SDATE
# 10 SERNO