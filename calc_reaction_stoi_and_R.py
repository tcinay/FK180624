# Import Libraries
import numpy as np

# Set Directory
output_path = 'output/OM_variations/{}'

# Select Method
method = 'redfield_oxidation_states' # select: "experimental", "anderson", "redfield", "oxidation_states" "redfield_oxidation_states"
# Experimental refers to CNP ratio etc calculated in cacl_oxycline_features
# _oxidation states refers to the sensitivity analysis

# Redfield & Anderson
if method == 'redfield':
    # Parameters
    C = 106
    H = 263
    O = 110
    N = 16
    P = 1
    O2 = 138
    charge = 0
    Cox = 0
    # Stoichiometry calculations
    DNRN_NO3 = round(2*C+0.5*H-O-1.5*N+2.5*P,2)
    DNRN_NO2 = DNRN_NO3
    DNRN_water = round(0.5*H-1.5*N-1.5*P,2)
    denit_NO2 = round(4/3*C+1/3*H-2/3*O-N+5/3*P,2)
    denit_N2 = round(2/3*C+1/6*H-1/3*O-0.5*N+5/6*P,2)
    denit_water = round(2/3*C+2/3*H-1/3*O-2*N-2/3*P,2)
    reactions = ['OM C:H:O:N:P ratio is {}:{}:{}:{}:{}'.format(C,H,O,N,P), 
                'DNRN: OM + {} NO3- + {} H+ -> {} CO2 + {} NH4+ + {} H3PO4 + {} NO2- + {} H20'.format(DNRN_NO3, N, C, N, P, DNRN_NO2, DNRN_water),
                'Denitrification: OM + {} NO2- + {}  H+ -> {} CO2 + {} NH4+ + {} H3PO4 + {} N2 + {} H2O'.format(denit_NO2,denit_NO2+N, C, N, P, denit_N2, denit_water)]
    with open(output_path.format('redfield/reactions.txt'), 'w') as f:
        f.write('\n'.join(reactions))
    # R Matrix
    R_DNRN_NO3 = -DNRN_NO3/C
    R_DNRN_NO2 = DNRN_NO2/C
    R_DNRN_NH4 = N/C
    R_DNRN_Nstar = (-DNRN_NO3+DNRN_NO2+N-16*P)/C
    R_DNRN_TA = (N-P)/C
    R_denit_NO2 = -denit_NO2/C
    R_denit_NH4 = N/C
    R_denit_Nstar = (-denit_NO2+N-16*P)/C
    R_denit_TA = (denit_NO2+N-P)/C
    R = np.array([(R_DNRN_NO3, R_DNRN_NO2,R_DNRN_NH4,R_DNRN_Nstar,R_DNRN_TA,1),
                    (0, R_denit_NO2,R_denit_NH4,R_denit_Nstar,R_denit_TA,1),
                    (2.909,-12,-9.091,-18.182,0.182,-1),
                    (46.296, -46.296,0,0,0,-1),
                    (0,0,0,0,2,1) ], dtype=float).T 
    np.savetxt(output_path.format('redfield/R.txt'), R, delimiter=',')
    print('done')

elif method == 'anderson':
    # Parameters
    C = 106
    H = 175
    O = 42
    N = 16
    P = 1
    O2 = 150
    charge = 0
    Cox = -0.45
    # Stoichiometry calculations
    DNRN_NO3 = round(2*C+0.5*H-O-1.5*N+2.5*P,2)
    DNRN_NO2 = DNRN_NO3
    DNRN_water = round(0.5*H-1.5*N-1.5*P,2)
    denit_NO2 = round(4/3*C+1/3*H-2/3*O-N+5/3*P,2)
    denit_N2 = round(2/3*C+1/6*H-1/3*O-0.5*N+5/6*P,2)
    denit_water = round(2/3*C+2/3*H-1/3*O-2*N-2/3*P,2)
    reactions = ['OM C:H:O:N:P ratio is {}:{}:{}:{}:{}'.format(C,H,O,N,P), 
                'DNRN: OM + {} NO3- + {} H+ -> {} CO2 + {} NH4+ + {} H3PO4 + {} NO2- + {} H20'.format(DNRN_NO3, N, C, N, P, DNRN_NO2, DNRN_water),
                'Denitrification: OM + {} NO2- + {}  H+ -> {} CO2 + {} NH4+ + {} H3PO4 + {} N2 + {} H2O'.format(denit_NO2,denit_NO2+N, C, N, P, denit_N2, denit_water)]
    with open(output_path.format('anderson/reactions.txt'), 'w') as f:
        f.write('\n'.join(reactions))
    # R Matrix
    R_DNRN_NO3 = -DNRN_NO3/C
    R_DNRN_NO2 = DNRN_NO2/C
    R_DNRN_NH4 = N/C
    R_DNRN_Nstar = (-DNRN_NO3+DNRN_NO2+N-16*P)/C
    R_DNRN_TA = (N-P)/C
    R_denit_NO2 = -denit_NO2/C
    R_denit_NH4 = N/C
    R_denit_Nstar = (-denit_NO2+N-16*P)/C
    R_denit_TA = (denit_NO2+N-P)/C
    R = np.array([(R_DNRN_NO3, R_DNRN_NO2,R_DNRN_NH4,R_DNRN_Nstar,R_DNRN_TA,1),
                    (0, R_denit_NO2,R_denit_NH4,R_denit_Nstar,R_denit_TA,1),
                    (2.909,-12,-9.091,-18.182,0.182,-1),
                    (46.296, -46.296,0,0,0,-1),
                    (0,0,0,0,2,1) ], dtype=float).T 
    np.savetxt(output_path.format('anderson/R.txt'), R, delimiter=',')
    print('done')

# Experimental Chanrged OM
elif method == 'experimental':
    # Parameters
    C = 112.5
    N = 11.4
    P = 1
    charge = 0
    Cox = 0.18
    twoO_minus_H = Cox*C-charge-3*N+5*P
    #Stoichiometry calculations
    DNRN_NO3 = round(2*C-0.5*(twoO_minus_H)-1.5*N+2.5*P,2)
    DNRN_NO2 = DNRN_NO3
    DNRN_water = '?'
    denit_NO2 = round(4/3*C-1/3*(twoO_minus_H)-N+5/3*P,2)
    denit_N2 = round(2/3*C-1/6*(twoO_minus_H)-0.5*N+5/6*P,2)
    denit_water = '?'
    reactions = ['OM C:H:O:N:P ratio is {}:a:b:{}:{}'.format(C,N,P), 
                '2b-a = {}'.format(round(twoO_minus_H,2)),
                'DNRN: OM + {} NO3- + {} H+ -> {} CO2 + {} NH4+ + {} H3PO4 + {} NO2- + {} H20'.format(DNRN_NO3, N, C, N, P, DNRN_NO2, DNRN_water),
                'Denitrification: OM + {} NO2- + {}  H+ -> {} CO2 + {} NH4+ + {} H3PO4 + {} N2 + {} H2O'.format(denit_NO2,denit_NO2+N, C, N, P, denit_N2, denit_water)]
    with open(output_path.format('experimental2/reactions.txt'), 'w') as f:
        f.write('\n'.join(reactions))
    # R Matrix
    R_DNRN_NO3 = -DNRN_NO3/C
    R_DNRN_NO2 = DNRN_NO2/C
    R_DNRN_NH4 = N/C
    R_DNRN_Nstar = (-DNRN_NO3+DNRN_NO2+N-N*P)/C
    R_DNRN_TA = (N-P)/C
    R_denit_NO2 = -denit_NO2/C
    R_denit_NH4 = N/C
    R_denit_Nstar = (-denit_NO2+N-N*P)/C
    R_denit_TA = (denit_NO2+N-P)/C
    R = np.array([(R_DNRN_NO3, R_DNRN_NO2,R_DNRN_NH4,R_DNRN_Nstar,R_DNRN_TA,1),
                    (0, R_denit_NO2,R_denit_NH4,R_denit_Nstar,R_denit_TA,1),
                    (2.909,-12,-9.091,-18.182,0.182,-1),
                    (46.296, -46.296,0,0,0,-1),
                    (0,0,0,0,2,1) ], dtype=float).T 
    np.savetxt(output_path.format('experimental2/R.txt'), R, delimiter=',')
    print('done')

# Oxidation State Variations
elif method == 'oxidation_states':
    # Various Carbon Oxidation States
    Cox_list = np.arange(-4.0, 4.5, 0.5)
    charge = 0
    fpath = list()
    for i in Cox_list:
        if i < 0:
            fpath.append('cox_neg{}'.format(round(abs(i),2)))
        elif i == 0:
            fpath.append('cox_zero')
        elif i>0:
            fpath.append('cox_pos{}'.format(round(i,1)))
    
    for i in np.arange(0,len(Cox_list)):
        # Parameters
        C = 112.5
        N = 11.4
        P = 1
        Cox = Cox_list[i]
        twoO_minus_H = Cox*C-charge-3*N+5*P
        #Stoichiometry calculations
        DNRN_NO3 = round(2*C-0.5*(twoO_minus_H)-1.5*N+2.5*P,2)
        DNRN_NO2 = DNRN_NO3
        DNRN_water = '?'
        denit_NO2 = round(4/3*C-1/3*(twoO_minus_H)-N+5/3*P,2)
        denit_N2 = round(2/3*C-1/6*(twoO_minus_H)-0.5*N+5/6*P,2)
        denit_water = '?'
        reactions = ['OM C:H:O:N:P ratio is {}:a:b:{}:{}'.format(C,N,P), 
                '2b-a = {}'.format(round(twoO_minus_H,2)),
                'DNRN: OM + {} NO3- + {} H+ -> {} CO2 + {} NH4+ + {} H3PO4 + {} NO2- + {} H20'.format(DNRN_NO3, N, C, N, P, DNRN_NO2, DNRN_water),
                'Denitrification: OM + {} NO2- + {}  H+ -> {} CO2 + {} NH4+ + {} H3PO4 + {} N2 + {} H2O'.format(denit_NO2,denit_NO2+N, C, N, P, denit_N2, denit_water)]
        with open(output_path.format('{}/reactions.txt'.format(fpath[i])), 'w') as f:
            f.write('\n'.join(reactions))
        # R Matrix
        R_DNRN_NO3 = -DNRN_NO3/C
        R_DNRN_NO2 = DNRN_NO2/C
        R_DNRN_NH4 = N/C
        R_DNRN_Nstar = (-DNRN_NO3+DNRN_NO2+N-N*P)/C
        R_DNRN_TA = (N-P)/C
        R_denit_NO2 = -denit_NO2/C
        R_denit_NH4 = N/C
        R_denit_Nstar = (-denit_NO2+N-N*P)/C
        R_denit_TA = (denit_NO2+N-P)/C
        R = np.array([(R_DNRN_NO3, R_DNRN_NO2,R_DNRN_NH4,R_DNRN_Nstar,R_DNRN_TA,1),
                    (0, R_denit_NO2,R_denit_NH4,R_denit_Nstar,R_denit_TA,1),
                    (2.909,-12,-9.091,-18.182,0.182,-1),
                    (46.296, -46.296,0,0,0,-1),
                    (0,0,0,0,2,1) ], dtype=float).T 
        np.savetxt(output_path.format('{}/R.txt'.format(fpath[i])), R, delimiter=',')   
    
    print('done')

# Redfield Oxidation State Variations
elif method == 'redfield_oxidation_states':
    # Various Carbon Oxidation States
    Cox_list = np.arange(-4.0, 4.5, 0.5)
    charge = 0
    fpath = list()
    for i in Cox_list:
        if i < 0:
            fpath.append('cox_redfield_neg{}'.format(round(abs(i),2)))
        elif i == 0:
            fpath.append('cox_redfield_zero')
        elif i>0:
            fpath.append('cox_redfield_pos{}'.format(round(i,1)))
    
    for i in np.arange(0,len(Cox_list)):
        # Parameters
        C = 106
        N = 16
        P = 1
        Cox = Cox_list[i]
        twoO_minus_H = Cox*C-charge-3*N+5*P
        #Stoichiometry calculations
        DNRN_NO3 = round(2*C-0.5*(twoO_minus_H)-1.5*N+2.5*P,2)
        DNRN_NO2 = DNRN_NO3
        DNRN_water = '?'
        denit_NO2 = round(4/3*C-1/3*(twoO_minus_H)-N+5/3*P,2)
        denit_N2 = round(2/3*C-1/6*(twoO_minus_H)-0.5*N+5/6*P,2)
        denit_water = '?'
        reactions = ['OM C:H:O:N:P ratio is {}:a:b:{}:{}'.format(C,N,P), 
                '2b-a = {}'.format(round(twoO_minus_H,2)),
                'DNRN: OM + {} NO3- + {} H+ -> {} CO2 + {} NH4+ + {} H3PO4 + {} NO2- + {} H20'.format(DNRN_NO3, N, C, N, P, DNRN_NO2, DNRN_water),
                'Denitrification: OM + {} NO2- + {}  H+ -> {} CO2 + {} NH4+ + {} H3PO4 + {} N2 + {} H2O'.format(denit_NO2,denit_NO2+N, C, N, P, denit_N2, denit_water)]
        with open(output_path.format('{}/reactions.txt'.format(fpath[i])), 'w') as f:
            f.write('\n'.join(reactions))
        # R Matrix
        R_DNRN_NO3 = -DNRN_NO3/C
        R_DNRN_NO2 = DNRN_NO2/C
        R_DNRN_NH4 = N/C
        R_DNRN_Nstar = (-DNRN_NO3+DNRN_NO2+N-N*P)/C
        R_DNRN_TA = (N-P)/C
        R_denit_NO2 = -denit_NO2/C
        R_denit_NH4 = N/C
        R_denit_Nstar = (-denit_NO2+N-N*P)/C
        R_denit_TA = (denit_NO2+N-P)/C
        R = np.array([(R_DNRN_NO3, R_DNRN_NO2,R_DNRN_NH4,R_DNRN_Nstar,R_DNRN_TA,1),
                    (0, R_denit_NO2,R_denit_NH4,R_denit_Nstar,R_denit_TA,1),
                    (2.909,-12,-9.091,-18.182,0.182,-1),
                    (46.296, -46.296,0,0,0,-1),
                    (0,0,0,0,2,1) ], dtype=float).T 
        np.savetxt(output_path.format('{}/R.txt'.format(fpath[i])), R, delimiter=',')   
    
    print('done')

else:
    print('Wrong Method')