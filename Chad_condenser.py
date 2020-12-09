# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


"""
Parameters for the air-cooled condenser (ACC)
"""
w_acc = 10.00               # width of the condenser in m
h_acc = 15.00               # height of condenser in m
beta = 60.00                # angle between base and wall
d_tube = 0.1905             # diameter of each tube in m
w_tube = 0.0254             # width of each tube in m
t_tube = 0.000127           # thickness of each tube in m
h_fin = 0.0254              # height of each fin in m
t_fin = 0.000254            # thickness of each fin in m
d_fin = 0.1651              # depth of each fin in m
s_fin = 0.00254             # spacing bewteen each fin in m
k_wall = 60                 # thermal conductivity of the wall in W/m*K
T_air_in = 303.00           # temperature of outside air in K
T_air_out = 312.00          # temperature of air leaving condenser in K
#T_steam_in = Variations.T_steam_out
T_steam_in = 319.00
h_air = 38.36               # heat transfer coefficient of air in W/m^2*K
h_steam = 257.25            # heat transfer coefficient of steam in W/m^2*K
mdot_steam = 7.00           # mass flow rate of steam in kg/s
mdot_air = 470.00           # mass flow rate of air in kg/s
perimeter_fin = 0.3307      # fin perimeter in m
SP = range(1,101)           # number of segments in the condenser
UA = np.zeros(100)          # heat transfer coefficient of each segment
LMTD = np.zeros(100)        # log mean temprature difference of each section
Q = np.zeros(100)           # heat transfer through each section
    
"""
Calculations based on parameters
"""
p_tube = 2*h_fin+w_tube
num_tube = (2*w_acc)/p_tube
A_tube_cs = np.pi/4*(w_tube-2*t_tube)**2+(d_tube-w_tube)*(w_tube-2*t_tube)
A_tube_in = h_acc*(2*(d_tube-w_tube)+np.pi*(w_tube-2*t_tube))
D_hy_tube = (4*A_tube_cs)/(2*(d_tube-w_tube)+np.pi*(w_tube-2*t_tube))
p_fin = 2*t_fin+s_fin
num_fin = h_acc/p_fin
A_fin = (2*h_fin*d_fin)+(2*h_fin*t_fin)
A_fin_cs = d_fin*t_fin
A_fin_base_tube = 2*num_fin*A_fin_cs
A_fin_tube = 2*num_fin*A_fin
A_fr = 2*h_acc*w_acc
A_tube_bare = h_acc*(2*(d_tube-w_tube)+np.pi*w_tube)-A_fin_base_tube
A_os_tube = A_tube_bare+A_fin_tube
A_os_total = A_os_tube*num_tube
m_fin = np.sqrt((perimeter_fin*h_air)/(k_wall*A_fin_cs))
n_fin = np.tanh(m_fin*h_fin)/(m_fin*h_fin)
n_total = 1-(2*num_fin*A_fin*(1-n_fin)/A_os_tube)
A_eff = A_tube_bare+n_fin*A_fin_tube
num_channel = 2*(h_acc/s_fin)*((2*w_acc)/(2*h_fin-w_tube))
A_ff = num_channel*s_fin*h_fin
#T_air_out = n_total*(T_steam_in-T_air_in)+T_air_in
T_steam_out = n_total*(T_air_in-T_steam_in)+T_steam_in
#LMTD = ((T_steam_in-T_air_out)-(T_steam_out-T_air_in))/(np.log((T_steam_in-T_air_out)/(T_steam_out-T_air_in)))
    
"""
Split the condenser into 100 segments and calculate heat transfer for each section
"""
for i in range(len(SP)):
    LMTD = ((T_steam_in-T_air_out)-(T_steam_out-T_air_in))/(np.log((T_steam_in-T_air_out)/(T_steam_out-T_air_in)))
    UA[i] = 1/((1/(h_steam*A_tube_in*1))+(t_tube/(k_wall*A_tube_in*1))+(1/(h_air*A_os_total*n_total*1)))
    Q[i] = UA[i]*LMTD
    Q[i] = Q[i]/1e3
        
n = {'fin efficiency':pd.Series(n_fin),
     'total efficiency':pd.Series(n_total),
     '# of tubes':pd.Series(num_tube),
     '# of fins':pd.Series(num_fin)}
ntable = pd.DataFrame(n)
    
"""
Power and Pressure drop calculations
"""
    
x = 0.15
rho_air = 1.177
Ums = mdot_air/(np.pi*(d_tube/2)*rho_air)
sigma = A_ff/A_fr
G = mdot_steam/(h_acc*w_acc)
Kc = 0.42*(1-sigma**2)**2
Ki = 0.7
Pr_air = 0.707
Pr_steam = 9.46
rho_steam = 0.598
mu_air = 1.846*10**(-5)
mu_steam = 0.02*10**(-3)
xtt = ((1-x)/x)**(0.9)*(rho_steam/rho_air)**(0.5)*(mu_air/mu_steam)**(0.1)
F_beta = (1+(1-x)**0.2*np.cos(beta-10))/x**0.4
Re_air = (rho_air*mdot_air*h_acc)/mu_air
Nu_air=1.09*Re_air**(0.45)*F_beta**(0.3)*np.sqrt(Pr_air/xtt)
Re_steam = (rho_steam*mdot_steam*d_tube)/mu_steam
Nu_steam = 1.09*Re_steam**(0.45)*F_beta**(0.3)*np.sqrt(Pr_steam/xtt)
fv = (0.79*np.log(Re_air)-1.64)**(-2)
#deltaP1 = Ki*mdot_air**2/3*(h_acc**3-3*x*h_acc+3*x**2*h_acc)*1e-3
#deltaPa = (G)**2*((1/rho_air)-(1/rho_steam))
#deltaPc = 2*fv*h_acc*G**2/(D_hy_tube*rho_steam)
#deltaPi = 0.5*(1-sigma**2-Kc)*G**2/rho_steam
#deltaPe = -0.5*(1-sigma**2-Kc)*G**2/rho_air
#deltaPtotal = deltaPa+deltaPc+deltaPi+deltaPe
deltaPf = 0.5*fv*rho_air*Ums**2*h_acc/d_tube*1e-3
    
Power = np.pi*(D_hy_tube/2)**2*Ums*deltaPf
    
print(ntable)
print('Mean Velocity =',Ums,' m/s')
print('Change in Pressure =', deltaPf, 'kPa')
print('Power =', Power, 'kW')
    
T_air_out_check = n_total*(T_steam_in-T_air_in)+T_air_in
n_check = (T_air_out-T_air_in)/(T_steam_in-T_air_in)*100
    
print('Efficiency Check =',n_check,' %')
    
