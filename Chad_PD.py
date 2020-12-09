# -*- coding: utf-8 -*-

"""

Created on Wed Oct 14 18:51:17 2020



@author: chad

"""



import numpy as np

import math







#Properties of Hot Fluid (Steam)

 

Cph = 1.996 #Specific Heat (kJ/kg*K)

Vh = 0.028e-3 #Viscosity (Pa*s)

Dh = 0.598 #Density (kg/m^3)

TCh = 0.058 # Thermal Conductivity (W/m*K)

Prh = 0.61 # Prandtl Number





#Properties of Cold Fluid (air) at 1 Bar and 25C



Cpc = 1.007 #Specific Heat (kJ/kg*K)

Vc = 18.37e-6 # Viscosity of air (Pa*s)

Dc = 1.184 #Density (kg/m^3)

TCc = 0.02551 # Thermal Conductivity (W/m*K)

Prc = 0.7296 # Prandtl Number

Twsat = 25 # Temperature of saturation (C)





#Mass Flow Rate of Hot Fluid



Mh = 1500.00 # (kg/s)



#Enthalpy Values of Water for different temperatures and positions



h1 = 2801 # Enthalpy of super heated steam at exit of SG 252C (kJ/kg)



h2 = 1832 # Enthalpy of saturated steam at exit of steam generation section 226C (kJ/kg)



h3 = 30.184 # Enthalpy of saturated air at entrance of steam generation section 270C (kJ/kg)



#Temperature of air entering the steam generator



Twi = 30.00 # (C)


#Steam Inlet Temperature

Tsi = 494.00


Thso = 200.00 # (C)


Tso = 60.00


#TUBE PROPERTIES 



L = 15.00 # Lenght of tubes (m)



d_o= 0.1905 #outer tube diameter in m 



t_w= 0.000127 #tube wall thickness 



d_i= d_o-2*t_w #inner tube diameter 



N = 263.00 # Number of tubes



#Energy Balance Set Up



#Part 3 - Steam Cooling



E1 = h1 - h2



E2 = Mh * Cph


E3 = E2*Tsi


#Part 2 Steam Generation



E4 = h3 - h2



#Part 1 Water Heater



E5 = Cpc*(Twsat - Twi)



E6 = E2*Tso


#Set Up 2 equations and solve for unknowns


A = np.array([[E1, -E2, 0], [E4, -E2, E2], [E5, 0, E2]])


B = np.array([-E3, 0, E6])


C = np.linalg.solve(A, B)


#Mass Flow Rate of Cold Fluid (Calcuated) (kg/s)



Mc = round(abs(C[0])) # (kg/s)



print('Mass Flow Rate of Air', Mc, 'kg/s')



#Temperature of Salt at Different Locations



Ts1 = round(C[1]) # Temp at exit of steam heater (C)



print('Temperatue of air at exit of condenser', Ts1, 'C')



Ts2 = round(C[2]) # Temp at exit of steam generation (C)



print('Temperatue of steam at exit of condenser', Ts2, 'C')



#Energy Balance Solutions 



q3 = Mc * E1 #Steam Heater (W/m)



q2 = abs(Mc * E4) # Steam Generation (W/m)



q1 = Mc * E5 # Water Heater (W/m)



qtot = q1 + q2 + q3 # Total Heat transfer per unit lenght (kW/m)



print('Heat Transfer per unit Lenght of System', qtot, 'kW/m')



#Weighted Average of Energies



# Steam Heater



qavg3 = q3/qtot # Percent of energy for steam heater



#Steam Generation



qavg2 = q2/qtot # Percent of energy for steam generation



#Water Heater



qavg1 = q1/qtot # Percent of energy for water heater 



# Log Mean Temperature Difference



#Steam Heater



LMTD3 = ((Tsi - Twsat) - (Ts1 - Thso)) / (np.log((Tsi - Twsat) / (Ts1 - Thso)))



#Steam Generation



LMTD2 = ((Ts1 - Twsat) - (Ts2 - Twsat)) / (np.log((Ts1 - Twsat) / (Ts2 - Twsat)))



# Water Heater



LMTD1 = ((Ts2 - Twi) - (Tso - Twsat)) / (np.log((Ts2 - Twi) / (Tso - Twsat)))



#Weighted Average LMTD 



LMTD = (LMTD3 * qavg3) + (LMTD2 * qavg2) + (LMTD1 * qavg1)



#Surface Area of Tubes



SA = N * (2 * np.pi * L * d_o)



#Solve for Overall Heat Transfer Coefficent U



U = qtot * 1000 / (LMTD * SA)



print('Overall Heat Transfer Coefficent U =', U, 'W/m^2*K')