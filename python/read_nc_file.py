import xarray as xr
import numpy as np
### open netcdf file ###
ds = xr.open_dataset("C:\\Work\\DIVAnd\\Oxygen_maps\\resultat\\nc\\O2\\Oxygen_Winter.4Danl.nc")
#print(df)
### extract values of variable ####

var_name = "Oxygen_L2"
ds["HYPOX"]=xr.where((ds[var_name]<20),ds[var_name]/ds[var_name]*1,ds[var_name],keep_attrs=True)
ds["HYPOX"]=xr.where(((ds["HYPOX"]<50) & (ds["HYPOX"]>=20)),ds["HYPOX"]/ds["HYPOX"]*2,ds["HYPOX"],keep_attrs=True)
ds["HYPOX"]=xr.where(((ds["HYPOX"]>=50)),ds["HYPOX"]*np.nan,ds["HYPOX"],keep_attrs=True)

print(ds[var_name].depth[0])

xr.Dataset.stack(ds["HYPOX"].depth)
#ds["HYPOX_depth"] = xr.where((ds[var_name] < 20), depth, ds[var_name],
 #                                keep_attrs=True)
#for depth in ds["HYPOX_depth"].depth:

#ds["HYPOX_depth"]=xr.where((ds[var_name.depth[0]]<20),ds[var_name]/ds[var_name]*1,ds[var_name],keep_attrs=True)
#ds["HYPOX_depth"]=xr.where(((ds["HYPOX_depth"]<50) & (ds["HYPOX_depth"]>=20)),ds["HYPOX_depth"]/ds["HYPOX_depth"]*2,ds["HYPOX_depth"],keep_attrs=True)
#ds["HYPOX_depth"]=xr.where(((ds["HYPOX_depth"]>=50)),ds["HYPOX_depth"]*np.nan,ds["HYPOX_depth"],keep_attrs=True)

#ds["HYPOX"]=xr.where((ds[var_name]<50),ds[var_name]/ds[var_name]*2,ds[var_name]*np.nan,keep_attrs=True)
#xr.where((ds[var_name]<0.5),ds[var_name]/ds[var_name]*0,ds[var_name]*np.nan,keep_attrs=True)
ds.to_netcdf('modified_test.nc') # rewrite to netcdf


#O2_L2 = df['Oxygen_L2']
#print(df.Oxygen_L2)
#df.Oxygen_L2.sel()
#O2_L2_flat = df['Oxygen_L2'].values.flatten() # 1-d datadf

