# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 12:22:31 2020

@author: chad
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


rho_water = 997.00           # density of water in kg/m^3
T_water_in = 298.00          # temperature of water entering condenser in K
T_water_out = 311.00         # temperature of water exiting condenser in K
T_steam_in2 = 319.00         # temperature of steam entering condenser in K
#T_steam_out2 = 332.00        # temperature of steam exiting condenser in K
P_water = 3.17
h_water = 30.00              # heat transfer coefficient of water in W/m^2*K
h_steam2 = 257.25            # heat transfer coefficient of steam in W/m^2*K
mdot_water = 470.00          # mass flow rate of water in kg/s
mdot_steam2 = 7.00           # mass flow rate of steam in kg/s
d_tube2 = 0.1905             # tube diameter in m
w_tube2 = 0.0254             # tube width in m
t_tube2 = 0.000127           # tube thickness in m
h_fin2 = 0.0254              # height of each fin in m
t_fin2 = 0.000254            # thickness of each fin in m
d_fin2 = 0.1651              # depth of each fin in m
s_fin2 = 0.00254             # spacing bewteen each fin in m
perimeter_fin2 = 0.3307      # fin perimeter in m
k_wall2 = 60                 # thermal conductivity of steel in W/m*K
h_wcc = 15.00                # height of condenser in m
w_wcc = 10.00                # width of condenser in m
mu_water = 1.002e-3          # dynamic viscocity of water in kg/m*s
SP2 = range(1,101)           # number of segments in the condenser
UA2 = np.zeros(100)          # heat transfer coefficient of each segment
LMTD2 = np.zeros(100)        # log mean temprature difference of each section
Q2 = np.zeros(100)           # heat transfer of each section in MW

"""
Calculations based on parameters
"""

p_tube2 = 2*h_fin2+w_tube2
num_tube2 = (2*w_wcc)/p_tube2
A_tube_cs2 = np.pi/4*(w_tube2-2*t_tube2)**2+(d_tube2-w_tube2)*(w_tube2-2*t_tube2)
A_tube_in2 = h_wcc*(2*(d_tube2-w_tube2)+np.pi*(w_tube2-2*t_tube2))
D_hy_tube2 = (4*A_tube_cs2)/(2*(d_tube2-w_tube2)+np.pi*(w_tube2-2*t_tube2))
p_fin2 = 2*t_fin2+s_fin2
num_fin2 = h_wcc/p_fin2
A_fin2 = (2*h_fin2*d_fin2)+(2*h_fin2*t_fin2)
A_fin_cs2 = d_fin2*t_fin2
A_fin_base_tube2 = 2*num_fin2*A_fin_cs2
A_fin_tube2 = 2*num_fin2*A_fin2
A_fr2 = 2*h_wcc*w_wcc
A_tube_bare2 = h_wcc*(2*(d_tube2-w_tube2)+np.pi*w_tube2)-A_fin_base_tube2
A_os_tube2 = A_tube_bare2+A_fin_tube2
A_os_total2 = A_os_tube2*num_tube2
m_fin2 = np.sqrt((perimeter_fin2*h_water)/(k_wall2*A_fin_cs2))
n_fin2 = np.tanh(m_fin2*h_fin2)/(m_fin2*h_fin2)
n_total2 = 1-(2*num_fin2*A_fin2*(1-n_fin2)/A_os_tube2)
A_eff2 = A_tube_bare2+n_fin2*A_fin_tube2
num_channel2 = 2*(h_wcc/s_fin2)*((2*w_wcc)/(2*h_fin2-w_tube2))
A_ff2 = num_channel2*s_fin2*h_fin2
Re_water = (rho_water*mdot_water*h_wcc)/mu_water
fv2 = (0.79*np.log(Re_water)-1.64)**-2
Vms = mdot_water/(np.pi*(d_tube2/2)*rho_water)
deltaP = 0.5*fv2*rho_water*Vms**2*h_wcc/d_tube2*1e-3
Power2 = np.pi*(D_hy_tube2/2)**2*Vms*deltaP
T_steam_out2 = n_total2*(T_water_in-T_steam_in2)+T_steam_in2

LMTD2 = ((T_steam_in2-T_water_out)-(T_steam_out2-T_water_in))/(np.log((T_steam_in2-T_water_out)/(T_steam_out2-T_water_in)))
     
for i in range(len(SP2)):
    UA2[i] = 1/((1/(h_steam2*A_tube_in2*1))+(t_tube2/(k_wall2*A_tube_in2*1))+(1/(h_water*A_os_total2*n_total2*1)))
    Q2[i] = UA2[i]*LMTD2
    Q2[i] = Q2[i]/1e3

n2 = {'fin efficiency':pd.Series(n_fin2),
     'total efficiency':pd.Series(n_total2),
     '# of tubes':pd.Series(num_tube2),
     '# of fins':pd.Series(num_fin2)}
ntable2 = pd.DataFrame(n2)

print(ntable2)
print('Mean Velocity =',Vms,' m/s')
print('Pressure Drop =',deltaP,' kPa')
print('Power =',Power2,' kW')

T_water_out_check = n_total2*(T_steam_in2-T_water_in)+T_water_in
n_check2 = (T_water_out-T_water_in)/(T_steam_in2-T_water_in)*100

print('Efficiency Check =',n_check2,' %')