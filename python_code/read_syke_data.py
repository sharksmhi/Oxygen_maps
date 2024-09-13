from pathlib import Path
import pandas as pd

def transform_syke_data(path):
    print("reading file {}\n... ... ...".format(path))
    df = pd.read_csv(
        open(path, encoding='utf-8'),
        sep=",",
    )
    print("done reading file")

    # Filtrera bort rader där PARAM inte är O2D eller O2S
    data = df.query(
        'parameter_code in ["O2D", "O2S"]'
    )

    # Kombinera parameterkod och metod så att det går att skilja det vi tror är BTL och CTD
    data.loc[:, "parameter_code"] = data["parameter_code"].fillna('NA').astype(str)
    data.loc[:, 'PARAM_METHOD'] = data['parameter_code'] + "_" + data['method_code']
    
    # konvertera från mg/l till umol/l
    # unit conversion mg/l -> ml/l *0.7, ml/l -> umol/l *44.661:
    data.loc[:, 'value'] = data['value']*0.700*44.661

    # döp om lite granna
    rename_dict = {
            "site_id": "SERNO",
            "site_depth": "WADEP",
            "wgs84_lat": "LATIT", 
            "wgs84_long": "LONGI",
            "depth_upper": "DEPH",
            "time": "SDATE",
            "site": "STATN",
            "PARAM_METHOD": "PARAM",
            "flag": "Q_flag",
            "value": "VALUE",
        }
    # Hämta en lista över kolumner som ska behållas (nycklar i rename_dict)
    columns_to_keep = list(rename_dict.keys())

    # Droppa alla kolumner som inte finns i rename_dict
    data = data[columns_to_keep]

    # Byt namn på de kolumner som finns i rename_dict
    data.rename(columns=rename_dict, inplace = True)

    data["Q_flag"] = data["Q_flag"].fillna('').astype(str)
    
    # skriv datumsträngen så som loadbigfile vill ha
    data["SDATE"] = data["SDATE"].apply(
        lambda row: row.replace(" ", "T") if isinstance(row, str) else row
    )
    
    index_list = [
        "SDATE",
        "SERNO",
        "WADEP",
        "LATIT",
        "LONGI",
        "STATN",
        "DEPH",
    ]
    data.sort_values(by=index_list, inplace=True)
    
    # Steg 1: Hitta alla rader som är dubbletter baserat på index_list + ["PARAM"] + VALUE är samma
    duplicates_to_drop = data.duplicated(subset=index_list + ["PARAM", "VALUE"], keep='first')

    # Steg 2: Ta bort dubbletter där både index_list + ["PARAM"] och VALUE är samma
    data_cleaned = data[~duplicates_to_drop]

    # Steg 3: Hitta duplikat baserat på index_list + ["PARAM"] men där VALUE skiljer sig
    potential_duplicates = data_cleaned.loc[data_cleaned.duplicated(subset=index_list + ["PARAM"], keep=False)]
    # duplicates_with_diff_value = potential_duplicates.groupby(index_list + ["PARAM"]).filter(lambda x: x["VALUE"].nunique() > 1)

    # Printar de rader där VALUE skiljer sig
    print(potential_duplicates)

    # Definiera ordningen på kolumnerna för loadbigfile
    column_order = [
        'LONGI',  # 0
        'LATIT',  # 1
        'VALUE',    # 2
        'DEPH',   # 3
        'PARAM', 'Q_flag', 'STATN',  'WADEP', 'random_8',  # 4 till 8
        'SDATE',  # 9
        'SERNO'   # 10
    ]

    # Om kolumner 4-8 är tomma, fyll dem med dummyvärden (här är alla fyllda)
    for col in column_order[4:9]:
        if col not in data_cleaned or data_cleaned[col].isnull().all():  # Kontrollera om kolumnen saknas eller är helt tom
            data_cleaned[col] = ''  # Fyll med en dummy-sträng eller annat värde

    # Ordna kolumnerna enligt ordern vi har specificerat
    data_cleaned = data_cleaned[column_order]

    # Detta är btl och odefinierad metod gissar vi
    data_cleaned.loc[data_cleaned["PARAM"].isin(["O2D_TI", "O2D_NA"])].to_csv(Path("//winfs-proj/data/proj/havgem/DIVA/syrekartor/data/", "syke_data_O2D_TI_NA.txt"), sep="\t", index = False)
    # Detta är ctd gissar vi, kanske är också EL CTD?
    data_cleaned.loc[data_cleaned["PARAM"] == "O2D_ELF"].to_csv(Path("//winfs-proj/data/proj/havgem/DIVA/syrekartor/data/", "syke_data_O2D_ELF.txt"), sep="\t", index = False)
    

path = "//winfs-proj/data/proj/havgem/DIVA/syrekartor/data/Syke_oxygen_opensea.csv"
print(Path("//winfs-proj/data/proj/havgem/DIVA/syrekartor/data/", "syke_data.txt"))
transform_syke_data(path=path)
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