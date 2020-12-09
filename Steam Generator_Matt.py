import math
import numpy as np
from ThermalEnergyStorageSystem import *
TES = ThermalEnergyStorageSys()
#Mass Flow Rate of Hot Fluid

Mh = TES.massflowrate # (kg/s)

#Mass Flow Rate of Cold Fluid

Mc = 290 # (kg/s)

#Specific Heats 

Cph = 1.52 #Specific Heat of Solar Salt(kJ/kg*K)

Cpl = 4.488 #Specific Heat (kJ/kg*K) liquid water at 70 Bar and 200C

Cps = 2.888 #Specific Heat (kJ/kg*K) superheated steam at 70 Bar and 370C

#Temperature of water entering the steam generator

Twi = 200 # (C)

#Temperature of Solar Salt entering

Tsi = TES.hotT - 273 # (C)

#Temperature of Solar Salt Leaving

Tso = TES.coldT - 273 # (C)

Twsat = 270 # Temperature of saturation (C)

U = 400 # Overall Heat Transfer Coefficent (W/m^2*K)


def Hx1(m_dot_hot, c_p_hot, T_cold_in, m_dot_cold, c_p_cold, T_cold_out, U, A, HE_Type):
    T_cold_in = T_cold_in + 273
    T_cold_out = T_cold_out + 273
    C_hot = m_dot_hot * c_p_hot
    C_cold = m_dot_cold * c_p_cold
    C_min = min(C_hot, C_cold)
    C_max = max(C_hot, C_cold)
    C_r = C_min / C_max
    NTU = U * A / C_min
    epsilon = effectiveness(NTU, C_r, HE_Type)
    Q = C_cold * (T_cold_out - T_cold_in)
    Q_max = Q / epsilon
    T_hot_in1 = T_cold_in + Q_max / C_min
    T_hot_out1 = T_hot_in1 - Q / C_hot
    
    print('Hx1')
    # print(C_hot)
    # print(C_cold)
    # print(C_min)
    # print(C_max)
    # print(C_r)
    # print('NTU =', NTU)
    print('epsilon =', round(epsilon, 2))
    # print(Q)
    # print(Q_max)
    print('T_hot_in1 (k)=', round(T_hot_in1, 2))
    print('T_hot_out1 (k)=', round(T_hot_out1, 2))
    return T_hot_in1, T_hot_out1
    

def Hx2(m_dot_hot, c_p_hot, T_hot_out2, m_dot_cold, c_p_cold, T_cold_in, U, A, HE_Type):
    T_cold_in = T_cold_in + 273
    C_hot = m_dot_hot * c_p_hot
    C_cold = np.inf
    C_min = min(C_hot, C_cold)
    C_max = np.inf
    C_r = C_min / C_max
    NTU = U * A / C_min
    epsilon = effectiveness(NTU, C_r, HE_Type)
    Q = m_dot_cold*(2750 - 1185)
    Q_max = Q / epsilon
    T_hot_in2 = T_cold_in + Q_max / C_min
    
    print(' ')
    print('HX2')
    # print(C_hot)
    # print(C_cold)
    # print(C_min)
    # print(C_max)
    # print(C_r)
    # print('NTU =', NTU)
    print('epsilon =', round(epsilon, 2))
    # print(Q)
    # print(Q_max)
    print('T_hot_in1 (k)=', round(T_hot_in2, 2))
    return T_hot_in2
    
def Hx3(m_dot_hot, c_p_hot, T_hot_out3, m_dot_cold, c_p_cold, T_hot_in3, U, A, HE_Type):
    T_hot_in3 = T_hot_in3 + 273
    C_hot = m_dot_hot * c_p_hot
    C_cold = m_dot_cold * c_p_cold
    C_min = min(C_hot, C_cold)
    C_max = max(C_hot, C_cold)
    C_r = C_min / C_max
    NTU = U * A / C_min
    epsilon = effectiveness(NTU, C_r, HE_Type)
    Q = C_hot * (T_hot_in3 - T_hot_out3)
    Q_max = Q / epsilon
    T_cold_in3 = T_hot_in3 - Q_max / C_min
    T_cold_out3 = T_cold_in3 + Q / C_cold
    
    print(' ')
    print('HX3')
    # print(C_hot)
    # print(C_cold)
    # print(C_min)
    # print(C_max)
    # print(C_r)
    # print('NTU =', NTU)
    print('epsilon =', round(epsilon, 2))
    # print(Q)
    # print(Q_max)
    print('T_cold_in (k)=', round(T_cold_in3, 2))
    print('Final Steam Out Temp (k)=', round(T_cold_out3, 2))
    return T_cold_out3, T_cold_in3

def effectiveness(NTU, C_r, HE_Type):
    if C_r == 0:
        epsilon = 1 - math.exp(-NTU)
        return epsilon

    if HE_Type == 'Parallel Flow':
        epsilon = (1 - math.exp(-NTU * (1 + C_r))) / (1 + C_r)
    elif HE_Type == 'Counter Flow':
        if C_r == 1:
            epsilon = NTU/(1+NTU)
        else:
            epsilon = (1-math.exp(-NTU*(1-C_r)))/(1-C_r*math.exp(-NTU*(1-C_r)))
    elif HE_Type == 'Double Shell Pass':
        N = 2.0
        NTUN = NTU/N
        epsilon1 = 2/(1+C_r+math.sqrt(1+C_r**2)*(1+math.exp(-NTUN*math.sqrt(1+C_r**2)))/(1+math.exp(-NTUN*math.sqrt(1+C_r**2))))
        epsilon = (((1-epsilon1*C_r)/(1-epsilon1))**(N-1)) / (((1-epsilon1*C_r)/(1-epsilon1))*(N-C_r))

    return epsilon


#Hx1(m_dot_hot, c_p_hot, T_cold_in, m_dot_cold, c_p_cold, T_cold_out, U, A, HE_Type):
X, W = Hx1(Mh, Cph, Twi, Mc, Cpl, Twsat, U, 6, 'Counter Flow')

#Hx2(m_dot_hot, c_p_hot, T_hot_out2, m_dot_cold, c_p_cold, T_cold_in, U, A, HE_Type):
Y = Hx2(Mh, Cph, X, Mc, Cpl, Twsat, U, 10.35, 'Counter Flow')

#Hx3(m_dot_hot, c_p_hot, T_hot_out3, m_dot_cold, c_p_cold, T_salt_in, U, A, HE_Type):
Z, V = Hx3(Mh, Cph, Y, Mc, Cps, Tsi, U, 3, 'Counter Flow')


#Check Q values

#Using Salt Delta T

QS = Mh*Cph*((Tsi + 273) - W)

#Using Water Properties

QW = Mc*Cpl*(Twsat-Twi) + Mc*(2760 - 1185) + Mc * Cps * (Z-V)

e = QW/QS

# Q for each section

Q1 = Mc*Cpl*(Twsat-Twi)

Q2 = Mc*(2760 - 1185)

Q3 = Mc * Cps * (Z-V)

# Percent energy of each section

P1 = 100 * Q1/QW

P2 = 100 * Q2/QW

P3 = 100 * Q3/QW

print('')

print('Salt Energy Balance = ', round(QS, 2))

print('')

print('Water/Steam Energy Balance =', round(QW, 2))

print('')

print('Percent Energy for Part 1 =', round(P1, 2), '%')

print('')

print('Percent Energy for Part 2 =', round(P2, 2), '%')

print('')

print('Percent Energy for Part 3 =', round(P3, 2), '%')



