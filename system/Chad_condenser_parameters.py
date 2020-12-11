# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 17:03:36 2020

@author: Chad
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

w_acc = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 15, 17, 18,
         19, 20, 21, 22, 23, 24, 25, 26, 27]
h_acc = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
         24, 25, 26, 27, 28, 29, 30, 31, 32, 33]
d_tube = [0.14, 0.15, 0.16, 0.17, 0.18, 0.1905, 0.2, 0.21, 0.22, 0.23,
          0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33,
          0.34, 0.35, 0.36, 0.37]
w_tube = [0.005, 0.01, 0.015, 0.02, 0.0225, 0.0254, 0.03, 0.035, 0.04, 0.045,
          0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1,
          0.105, 0.11, 0.115]
t_tube = [0.00008, 0.00009, 0.0001, 0.00011, 0.00012, 0.000127, 0.00013, 
          0.00014, 0.00015, 0.00016, 0.00016, 0.00017, 0.00018, 0.00019,
          0.0002, 0.00021, 0.00022, 0.00023, 0.00024, 0.00025, 0.00026,
          0.00027, 0.00028, 0.00029]
h_fin = [0.005, 0.01, 0.015, 0.02, 0.0225, 0.0254, 0.03, 0.035, 0.04, 0.045,
          0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1,
          0.105, 0.11, 0.115]
t_fin = [0.00016, 0.00018, 0.0002, 0.00022, 0.00024, 0.000254, 0.00026, 
          0.00028, 0.0003, 0.00032, 0.00032, 0.00034, 0.00036, 0.00038,
          0.0004, 0.00042, 0.00044, 0.00046, 0.00048, 0.0005, 0.00052,
          0.00054, 0.00056, 0.00058]
d_fin = [0.14, 0.145, 0.15, 0.155, 0.16, 0.1651, 0.17, 0.175, 0.18, 0.185,
         0.19, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,
         0.75, 0.8]
s_fin = [0.0005, 0.001, 0.0015, 0.002, 0.00225, 0.00254, 0.003, 0.0035, 0.004, 0.0045,
          0.005, 0.0055, 0.006, 0.0065, 0.007, 0.0075, 0.008, 0.0085, 0.009, 0.0095, 0.01,
          0.0105, 0.011, 0.0115]
k_wall = 60
rho_air = 1.177
mu_air = 1.846*10**(-5)
T_air_in = np.array([303, 305, 304, 306, 307, 310, 311, 314, 316, 317, 317, 317,
            317, 317, 317, 317, 317, 317, 316, 315, 314, 306, 305, 303])
T_steam_in = np.array([318, 320, 319, 319, 320, 319, 320, 319, 319, 319, 319, 319,
              319, 319, 319, 319, 319, 319, 319, 319, 319, 319, 319, 319])
h_air = 38.36
h_steam = 257.25 
mdot_air = 470
perimeter_fin = [0.28, 0.29, 0.3, 0.31, 0.32, 0.3307, 0.34, 0.35, 0.36, 
                 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 
                 0.46, 0.47, 0.48, 0.49, 0.5, 0.51]
T_air_out = np.array([317.13, 317.4, 314.7, 314.08, 314.7, 315.0, 315.47,
                      316.27, 317.26, 317.79, 317.73, 317.698, 317.7, 317.71,
                      317.711, 317.712, 317.711, 317.71, 317.06, 316.41, 315.75,
                      310.53, 309.85, 308.5])
UA = np.zeros(24)          # heat transfer coefficient of each segment
LMTD = np.zeros(24)        # log mean temprature difference of each section
Q = np.zeros(24)           # heat transfer through each section
p_tube = np.zeros(24)
num_tube = np.zeros(24)
A_tube_cs = np.zeros(24)
A_tube_in = np.zeros(24)
D_hy_tube = np.zeros(24)
p_fin = np.zeros(24)
num_fin = np.zeros(24)
A_fin = np.zeros(24)
A_fin_cs = np.zeros(24)
A_fin_base_tube = np.zeros(24)
A_fin_tube = np.zeros(24)
A_fr = np.zeros(24)
A_tube_bare = np.zeros(24)
A_os_tube = np.zeros(24)
A_os_total = np.zeros(24)
m_fin = np.zeros(24)
n_fin = np.zeros(24)
n_total = np.zeros(24)
A_eff = np.zeros(24)
num_channel = np.zeros(24)
A_ff = np.zeros(24)
T_steam_out = np.zeros(24)
#T_air_out = np.zeros(24)
Ums = np.zeros(24)
Re_air = np.zeros(24)
fv = np.zeros(24)
deltaP = np.zeros(24)
Power = np.zeros(24)
"""
Calculations based on parameters
"""
for i in range(len(w_acc)):
    p_tube[i] = 2*h_fin[i]+w_tube[i]
    num_tube[i] = (2*w_acc[i])/p_tube[i]
    A_tube_cs[i] = np.pi/4*(w_tube[i]-2*t_tube[i])**2+(d_tube[i]-w_tube[i])*(w_tube[i]-2*t_tube[i])
    A_tube_in[i] = h_acc[i]*(2*(d_tube[i]-w_tube[i])+np.pi*(w_tube[i]-2*t_tube[i]))
    D_hy_tube[i] = (4*A_tube_cs[i])/(2*(d_tube[i]-w_tube[i])+np.pi*(w_tube[i]-2*t_tube[i]))
    p_fin[i] = 2*t_fin[i]+s_fin[i]
    num_fin[i] = h_acc[i]/p_fin[i]
    A_fin[i] = (2*h_fin[i]*d_fin[i])+(2*h_fin[i]*t_fin[i])
    A_fin_cs[i] = d_fin[i]*t_fin[i]
    A_fin_base_tube[i] = 2*num_fin[i]*A_fin_cs[i]
    A_fin_tube[i] = 2*num_fin[i]*A_fin[i]
    A_fr[i] = 2*h_acc[i]*w_acc[i]
    A_tube_bare[i] = h_acc[i]*(2*(d_tube[i]-w_tube[i])+np.pi*w_tube[i])-A_fin_base_tube[i]
    A_os_tube[i] = A_tube_bare[i]+A_fin_tube[i]
    A_os_total[i] = A_os_tube[i]*num_tube[i]
    m_fin[i] = np.sqrt((perimeter_fin[i]*h_air)/(k_wall*A_fin_cs[i]))
    n_fin[i] = np.tanh(m_fin[i]*h_fin[i])/(m_fin[i]*h_fin[i])
    n_total[i] = 1-(2*num_fin[i]*A_fin[i]*(1-n_fin[i])/A_os_tube[i])
    A_eff[i] = A_tube_bare[i]+n_fin[i]*A_fin_tube[i]
    num_channel[i] = 2*(h_acc[i]/s_fin[i])*((2*w_acc[i])/(2*h_fin[i]-w_tube[i]))
    A_ff[i] = num_channel[i]*s_fin[i]*h_fin[i]
    T_steam_out[i] = n_total[i]*(T_air_in[i]-T_steam_in[i])+T_steam_in[i]
    #T_air_out[i] = n_total[i]*(T_steam_in[i]-T_air_in[i])+T_air_in[i]
    LMTD[i] = ((T_steam_in[i]-T_air_out[i])-(T_steam_out[i]-T_air_in[i]))/(np.log((T_steam_in[i]-T_air_out[i])/(T_steam_out[i]-T_air_in[i])))
    UA[i] = 1/((1/(h_steam*A_tube_in[i]))+(t_tube[i]/(k_wall*A_tube_in[i]))+(1/(h_air*A_os_total[i]*n_total[i])))
    Q[i] = UA[i]*LMTD[i]
    Q[i] = Q[i]/1e3
    Ums[i] = mdot_air/(np.pi*(d_tube[i]/2)*rho_air)
    Re_air[i] = (rho_air*mdot_air*h_acc[i])/mu_air
    fv[i] = (0.79*np.log(Re_air[i])-1.64)**(-2)
    deltaP[i] = 0.5*fv[i]*rho_air*Ums[i]**2*h_acc[i]/d_tube[i]*1e-3
    Power[i] = np.pi*(D_hy_tube[i]/2)**2*Ums[i]*deltaP[i]
    
n = {'fin efficiency':pd.Series(n_fin),
     'total efficiency':pd.Series(n_total),
     '# of tubes':pd.Series(num_tube),
     '# of fins':pd.Series(num_fin)}
ntable = pd.DataFrame(n)
print(ntable)

x = A_fr
y = n_total
z = deltaP
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
curve1, = ax1.plot(x, y, label="Efficiency", color='r')
curve2, = ax2.plot(x, z, label="Area", color='b')
curves = [curve1, curve2]
ax1.legend(curves, [curve.get_label() for curve in curves])
ax1.set_xlabel("Area of Condenser (m^2)", color='k')
ax1.set_ylabel("Efficiency (%)", color='k')
ax2.set_ylabel("Pressure Drop (kPa)", color='k')
plt.title("Efficiency and Pressure Drop Change with Area")
plt.show()

x = A_fr
y = n_total
z = Power
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
curve1, = ax1.plot(x, y, label="Efficiency", color='r')
curve2, = ax2.plot(x, z, label="Area", color='b')
curves = [curve1, curve2]
ax1.legend(curves, [curve.get_label() for curve in curves])
ax1.set_xlabel("Area of Condenser (m^2)", color='k')
ax1.set_ylabel("Efficiency (%)", color='k')
ax2.set_ylabel("Power (kW)", color='k')
plt.title("Efficiency and Power Change with Area")
plt.show()

x = A_fr
y = n_total
z = Q
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
curve1, = ax1.plot(x, y, label="Efficiency", color='r')
curve2, = ax2.plot(x, z, label="Area", color='b')
curves = [curve1, curve2]
ax1.legend(curves, [curve.get_label() for curve in curves])
ax1.set_xlabel("Area of Condenser (m^2)", color='k')
ax1.set_ylabel("Efficiency (%)", color='k')
ax2.set_ylabel("Heat Transfer (kW)", color='k')
plt.title("Efficiency and Heat Transfer with Area")
plt.show()