from pathlib import Path
import pandas as pd
import numpy as np

# Funktion för att applicera logiken för att välja OXYGEN värde i varje grupp
def choose_oxygen_value(group):
    # 1. Om det finns ett negativt OXYGEN värde, välj det
    negative_oxygen = group[group['OXYGEN'] < 0 & (group['Q_flag'] != "L")] 
    if not negative_oxygen.empty:
        return negative_oxygen['OXYGEN'].iloc[0]  # Returnerar första negativa värdet

    # 2. Om det finns ett värde med method_code == "TI" och Q_flag != "L", välj det
    ti_oxygen = group[(group['method_code'] == "TI") & (group['Q_flag'] != "L")]
    if not ti_oxygen.empty:
        return ti_oxygen['OXYGEN'].iloc[0]  # Returnerar första värdet som matchar villkoret

    # 3. Om det finns ett värde med method_code == "ELF", välj det
    elf_oxygen = group[(group['method_code'] == "ELF")  & (group['Q_flag'] != "L")]
    if not elf_oxygen.empty:
        return elf_oxygen['OXYGEN'].iloc[0]

    # 4. Om det finns ett värde med method_code == "EL", välj det
    el_oxygen = group[(group['method_code'] == "EL") & (group['Q_flag'] != "L")]
    if not el_oxygen.empty:
        return el_oxygen['OXYGEN'].iloc[0]

    # 5. I sista hand, om method_code är NaN, välj det värdet
    nan_oxygen = group[group['method_code'].isna()]
    if not nan_oxygen.empty:
        return nan_oxygen['OXYGEN'].iloc[0]

    # Om inget av ovanstående finns, returnera NaN
    return np.nan

def transform_syke_data(path):
    print("reading file {}\n... ... ...".format(path))
    df = pd.read_csv(
        open(path, encoding='utf-8'),
        sep=",",
    )
    print("done reading file")

    # Filtrera bort rader där PARAM inte är O2D eller O2S
    data = df.query(
        'parameter_code in ["O2D", "H2SS"]')

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

    data['OXYGEN'] = np.nan
    data['value'] = data['value'].astype(float)

    # Filtrera för PARAM == "H2SS" och multiplicera VALUE med -0.44
    data.loc[data['parameter_code'] == 'H2SS', 'OXYGEN'] = data['value'] * -0.029342

    # Filtrera för PARAM == "O2D" och multiplicera VALUE med 0.700 * 44.661
    data.loc[data['parameter_code'] == 'O2D', 'OXYGEN'] = data['value'] * 0.700 * 44.661

    data.drop(data[(data['Q_flag'] != "L")].index, inplace=True)
    data.drop(data[(data['OXYGEN'] < -800)].index, inplace=True)
    data.drop(data[(data['OXYGEN'] > 600)].index, inplace=True)
    # data = data.dropna(subset=['OXYGEN'])
    headers = [
        "SDATE",
        "SERNO",
        "STATN",
        "LATIT",
        "LONGI",
        "DEPH",
    ]
    data.sort_values(by=headers, inplace=True)
    
        # Definiera ordningen på kolumnerna för loadbigfile
    column_order = [
        'LONGI',  # 0
        'LATIT',  # 1
        'OXYGEN',    # 2
        'DEPH',   # 3
        'parameter_code', 'Q_flag', 'STATN',  'WADEP', 'random_8',  # 4 till 8
        'SDATE',  # 9
        'SERNO'   # 10
    ]

    # Om kolumner 4-8 är tomma, fyll dem med dummyvärden (här är alla fyllda)
    for col in column_order[4:9]:
        if col not in data or data[col].isnull().all():  # Kontrollera om kolumnen saknas eller är helt tom
            data[col] = ''  # Fyll med en dummy-sträng eller annat värde

    # Ordna kolumnerna enligt ordern vi har specificerat
    # skriv datumsträngen så som loadbigfile vill ha
    data["SDATE"] = data["SDATE"].apply(
        lambda row: row.replace(" ", "T") if isinstance(row, str) else row
    )
    data = data[column_order]    
    return data
    

path = "//winfs-proj/data/proj/havgem/DIVA/syrekartor/data/Syke_oxygen_opensea.csv"
print(Path("//winfs-proj/data/proj/havgem/DIVA/syrekartor/data/", "Syke_oxygen_opensea.txt"))
open_sea = transform_syke_data(path=Path("//winfs-proj/data/proj/havgem/DIVA/syrekartor/data/", "Syke_oxygen_opensea.csv"))
coastal = transform_syke_data(path=Path("//winfs-proj/data/proj/havgem/DIVA/syrekartor/data/", "Syke_oxygen_coastal.csv"))
data = pd.concat([open_sea, coastal])
data.to_csv(Path("//winfs-proj/data/proj/havgem/DIVA/syrekartor/data/", "syke_data_coastal_formatted.txt"), sep="\t", index = False, header = True)
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