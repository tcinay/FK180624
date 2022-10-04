## Import Libraries
import pandas as pd
import numpy as np
import gsw
import PyCO2SYS as pyco2

# Import Parameters and Results
falkor = pd.read_csv('../Python Version/falkor.csv')

# Select Flagged Bad Data
st2_cast5 = falkor.index[(falkor['station']==2)&(falkor['cast']==5)].tolist()
falkor.loc[st2_cast5, 'sal_flag']=0 # Cast was removed as it was potentially influenced by Hurricane Fabio. No significant impact on the further analysis. 
flag_good = np.array([2,6])
idx_flag = np.where((np.isin(falkor['sal_flag'],flag_good))&(np.isin(falkor['O2_flag'],flag_good))&
(np.isin(falkor['NOx_flag'],flag_good))&(np.isin(falkor['phosphate_flag'],flag_good))&
(np.isin(falkor['NH4_flag'],flag_good))&(np.isin(falkor['pH_tot_flag'],flag_good))&
(np.isin(falkor['TA_flag'],flag_good))&(np.isin(falkor['NO2_flag'],flag_good))&(np.isin(falkor['fluor_flag'],flag_good)))

# Calculate Density 
SA =  gsw.SA_from_SP(falkor['sal'], falkor['press'], falkor['lon'], falkor['lat'])
CT = gsw.CT_from_t(SA, falkor['temperature'], falkor['press'])
falkor['sigma0'] = gsw.density.sigma0(SA, CT)
falkor['rho'] = gsw.density.rho(SA, CT, falkor['press'])
 
# Calculate Carbonate System 
Z = pyco2.sys(par1 = falkor['TA'], par2 = falkor['pH_tot'], par1_type = 1, par2_type= 3, salinity= falkor['sal'], temperature=25, temperature_out=falkor['temperature'], pressure= 0, pressure_out=falkor['press'],
total_phosphate= falkor['phosphate'], total_ammonia=falkor['NH4'], opt_pH_scale= 1, opt_k_carbonic= 10, opt_total_borate= 2, opt_k_fluoride= 2)
falkor['DIC'] = Z['dic']
falkor['pH_insitu'] = Z['pH_out']

# Remove Bad Data and Save New Vectors
T = np.array(falkor['temperature'])[idx_flag] # degrees C
S = np.array(falkor['sal'])[idx_flag] # psu
P = np.array(falkor['press'])[idx_flag] # dbar
rho = np.array(falkor['density'])[idx_flag] # kg/m3
sigma0 = np.array(falkor['sigma0'])[idx_flag] # kg/m3 # 
DIC = np.array(falkor['DIC'])[idx_flag] #umol/kg  
DIP = np.array(falkor['phosphate'])[idx_flag]/rho*1000 #umol/kg
NO2 = np.array(falkor['NO2'])[idx_flag]/rho*1000 #umol/kg
NO3 = (np.array(falkor['NOx'])[idx_flag]-np.array(falkor['NO2'])[idx_flag])/rho*1000 #umol/kg
NH4 = np.array(falkor['NH4'])[idx_flag]/rho*1000 #umol/kg
Nstar = (NO2+NO3+NH4)-16*DIP+2.9 # umol/kg  ## Change the DIP coefficient (either 11.4 or 16) for sensitivity analyses. 
TA = np.array(falkor['TA'])[idx_flag] #umol/kg
pH = np.array(falkor['pH_insitu'])[idx_flag] 
O2 = np.array(falkor['O2'])[idx_flag] #umol/kg
station = np.array(falkor['station'])[idx_flag]
cast= np.array(falkor['cast'])[idx_flag]
lat = np.array(falkor['lat'])[idx_flag]
lon = np.array(falkor['lon'])[idx_flag]
pH_tot = np.array(falkor['pH_tot'])[idx_flag]
pH_temp =  np.array(falkor['pH_temp'])[idx_flag]
fluor = np.array(falkor['fluor'])[idx_flag]
falkor_df_data = np.array([station, cast, lat, lon, T, S, P, rho,sigma0, DIC,DIP,NO3, NO2,NH4, Nstar, TA, pH, pH_tot, pH_temp, O2, fluor]).T
falkor_clean= pd.DataFrame(data= falkor_df_data, columns=['Station', 'cast', 'lat', 'lon','T', 'S', 'P', 'rho','sigma0','DIC','DIP','NO3', 'NO2', 'NH4', 'Nstar', 'TA', 'pH', 'pH_tot', 'pH_temp', 'O2',  'fluor'])

# Save the clean version
falkor_clean.to_csv('falkor_clean.csv')