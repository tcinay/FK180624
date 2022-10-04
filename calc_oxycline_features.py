# Import Libraries
import numpy as np
import pandas as pd
import statsmodels.api as sm
import gsw

# Load Data
falkor = pd.read_csv('falkor_clean.csv')
falkor = falkor[['Station', 'lon', 'lat' ,'T', 'S', 'P', 'O2', 'sigma0', 'rho', 'DIC', 'DIP', 'NO3', 'NO2', 'NH4']]
st = np.array((1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))
falkor = falkor.iloc[np.isin(falkor['Station'], st)]

# Calculate AOU 
SA =  gsw.SA_from_SP(falkor['S'], falkor['P'], falkor['lon'], falkor['lat'])
PT = gsw.pt0_from_t(SA, falkor['T'], falkor['P'])
O2sat = gsw.O2sol_SP_pt(falkor['S'], PT)
AOU = O2sat-falkor['O2']
falkor['AOU'] = AOU 
falkor['DIN'] = falkor['NO3']+falkor['NO2']+falkor['NH4']

# Calculate -O2:P, AOU, C:N:P for different O2 ranges (robust fit)
oxycline = falkor.where((falkor['P']>45)&(falkor['P']<95))
oxycline = oxycline.dropna(axis=0)

# RLM
xx = sm.add_constant(oxycline['DIP'])
rO2 = sm.RLM(oxycline['O2'], xx, M = sm.robust.norms.HuberT()).fit()
rAOU = sm.RLM(oxycline['AOU'], xx, M = sm.robust.norms.HuberT()).fit()
rC = sm.RLM(oxycline['DIC'], xx, M = sm.robust.norms.HuberT()).fit()
rN = sm.RLM(oxycline['DIN'], xx, M = sm.robust.norms.HuberT()).fit()

# Calculate respiration quotient and Cox
Cox1 = round(4-4*(rAOU.params[1]-2*rN.params[1])/rC.params[1],2)
Cox2 = round(4-4*(rAOU.params[1])/rC.params[1],2)

# Print Results
text1 = 'RLS C:N:P is {}+-{} : {}+-{} : 1 (se)'.format(round(rC.params[1],3), round(rC.bse[1],3), round(rN.params[1],3), round(rN.bse[1],3))
text2 = 'RLS AOU:P is {}+-{}'.format(round(rAOU.params[1],3), round(rAOU.bse[1],3))
text3 = 'RLS Cox calculated with AOU:P is {}, and {} with no nitrification'.format(Cox1,Cox2)


print(text1,text2, text3, sep='\n')
