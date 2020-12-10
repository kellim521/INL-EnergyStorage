# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 10:46:06 2020

@author: Riley
"""

import numpy as np
import matplotlib.pyplot as plt
import PowerBlock1_Riley as PB1
import PowerBlock2_Riley as PB2
import PowerBlock3_Riley as PB3
import SteamGenerator_Matt as SG
#LOAD VARIATION ---------------------------------------------------------------

steamFlow = 155
nmax,Qmax,Pmax = PB1.Cycle(1.0,steamFlow)

#MISO DATA
#JAN 1
demandCurve = [10669,10259,9987,9847,9825,9953,10373,10977,11597,11778,11751,11613,11440,11192,11076,11012,11076,11454,12102,12022,11688,11395,11113,10894]

#APRIL 1
#demandCurve = [8493,8303,8226,8250,8407,8848,9559,10079,10092,9882,9729,9708,9674,9582,9405,9222,9184,9389,9606,9741,9885,9601,9170,8820]

#JULY 1
#demandCurve = [10056,9487,9169,9018,9075,9498,10255,11112,11755,12218,12730,13219,13755,14143,14331,14479,14612,14715,14674,14465,14064,13694,12995,11921]

#OCTOBER 1
#demandCurve = [8076,7781,7653,7603,7698,8202,9091,9788,9832,9625,9579,9702,9787,9868,9817,9696,9772,9931,10088,10225,10046,9580,9004,8472]

#CAISO DATA
#JAN 1
#demandCurve = [21533,20681,19933,19439,19198,19209,19808,20409,20122,19487,18677,18294,17910,18027,18488,19550,20927,23138,25060,24907,24482,23719,22683,21236]

#APRIL 1
#demandCurve = [20025,19205,18614,18266,18246,18799,20158,21382,21558,20894,19896,18710,18298,17794,17608,17852,18442,19704,21695,23317,24974,24269,23074,21651]

#JULY 1
#demandCurve = [25092,23695,22671,21988,21717,22111,22840,23740,24518,24971,25084,25185,25420,26039,26862,28017,29436,30990,32162,32440,31616,31141,29344,26985]

#OCTOBER 1
#demandCurve = [27653,25979,24546,23604,23130,23447,24576,26031,27126,28121,29317,31225,33590,36605,39470,41601,42985,43262,42128,40615,38409,35674,32975,30075]

demandPercent = np.zeros(len(demandCurve))
efficiency = np.zeros(len(demandCurve))
hour = np.zeros(len(demandCurve))
requiredPower = np.zeros(len(demandCurve))
requiredHeat  = np.zeros(len(demandCurve))

for i in range(len(demandCurve)):
    demandPercent[i] = demandCurve[i]/np.max(demandCurve)
    nth,Q,P = PB1.Cycle(1.0,steamFlow)
    for j in range(40):
        load = 1 - j/100
        nth,Q,P = PB1.Cycle(load,steamFlow)
        hour[i] = i+1
        if(P >= demandPercent[i]*Pmax):
            requiredPower[i] = P
            requiredHeat[i] = Q
            efficiency[i] = nth
        else:
            break

x = hour
y = efficiency
z = demandCurve
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
curve1, = ax1.plot(x, y, label="Efficiency", color='r')
curve2, = ax2.plot(x, z, label="Demand", color='b')
curves = [curve1, curve2]
ax1.legend(curves, [curve.get_label() for curve in curves])
ax1.set_xlabel("Time of Day (hr)", color='k')
ax1.set_ylabel("Efficiency (%)", color='k')
ax2.set_ylabel("Demand (MW)", color='k')
plt.title("Efficiency with Demand Curve Variation")
plt.show()

#TEMP VARIATION ---------------------------------------------------------------

T_hot_inputs = [590, 580, 570, 550, 530, 490, 500, 540, 510, 560, 600, 620,
              650, 660, 690, 680, 670, 650, 660, 630, 610, 600, 610, 590]

T_steam_out = np.zeros(len(T_hot_inputs))
hourT = np.zeros(len(T_hot_inputs))
efficiencytemp = np.zeros(len(T_hot_inputs))

for i in range(len(T_hot_inputs)):
    Z, V = SG.Hx3(SG.Mh, SG.Cph, SG.Y, SG.Mc, SG.Cps, T_hot_inputs[i], SG.U, 3, 'Counter Flow')
    T_steam_out[i] = Z
    PB2.mainTemp = T_steam_out[i]
    PB2.reheatTemp = T_steam_out[i]
    nth,Q,P = PB2.Cycle(1.0,steamFlow)
    #nth,Q,P = power2.Cycle(1.0,steamFlow,power2.mainTemp,power2.reheatTemp)
    hourT[i] = i+1
    efficiencytemp[i] = nth

x = hourT
y = efficiencytemp
z = T_steam_out
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
curve1, = ax1.plot(x, y, label="Efficiency", color='r')
curve2, = ax2.plot(x, z, label="Temperature", color='b', linestyle='dashed')
curves = [curve1, curve2]
ax1.legend(curves, [curve.get_label() for curve in curves])
ax1.set_xlabel("Time of Day (hr)", color='k')
ax1.set_ylabel("Efficiency (%)", color='k')
ax2.set_ylabel("Temperature (K)", color='k')
plt.title("Efficiency with Temperature Input Variation")
plt.show()

#OTHER PARAMETER VARIATION ----------------------------------------------------------

nmax2,Qmax2,Pmax2, dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8, dn9, dn10, dn11, dn12 = PB3.Cycle(1.0,steamFlow)
demandPercent2 = np.zeros(len(demandCurve))
efficiency2 = np.zeros(len(demandCurve))
hour2 = np.zeros(len(demandCurve))
requiredPower2 = np.zeros(len(demandCurve))
requiredHeat2  = np.zeros(len(demandCurve))
hpturbineinlettemp = np.zeros(len(demandCurve))
ipturbineinlettemp = np.zeros(len(demandCurve))
lpturbineinlettemp = np.zeros(len(demandCurve))
condenseroutlettemp = np.zeros(len(demandCurve))
massflowhpturbineinlet = np.zeros(len(demandCurve))
massflowipturbineinlet = np.zeros(len(demandCurve))
massflowlpturbineinlet = np.zeros(len(demandCurve))
massflowcondenseroutlet = np.zeros(len(demandCurve))
pressurehpturbineinlet = np.zeros(len(demandCurve))
pressureipturbineinlet = np.zeros(len(demandCurve))
pressurelpturbineinlet = np.zeros(len(demandCurve))
pressurecondenseroutlet = np.zeros(len(demandCurve))

for i in range(len(demandCurve)):
    demandPercent2[i] = demandCurve[i]/np.max(demandCurve)
    nth2,Q2,P2,M1T,M2T,M4T,M6T,M1m,M2m,M4m,M6m,M1P,M2P,M4P,M6P = PB3.Cycle(1.0,steamFlow)
    for j in range(40):
        load2 = 1 - j/100
        nth2,Q2,P2,M1T,M2T,M4T,M6T,M1m,M2m,M4m,M6m,M1P,M2P,M4P,M6P = PB3.Cycle(load2,steamFlow)
        hour2[i] = i+1
        if(P2 >= demandPercent2[i]*Pmax2):
            requiredPower2[i] = P2
            requiredHeat2[i] = Q2
            efficiency2[i] = nth2
            hpturbineinlettemp[i] = M1T
            ipturbineinlettemp[i] = M2T
            lpturbineinlettemp[i] = M4T
            condenseroutlettemp[i] = M6T
            massflowhpturbineinlet[i] = M1m
            massflowipturbineinlet[i] = M2m
            massflowlpturbineinlet[i] = M4m
            massflowcondenseroutlet[i] = M6m
            pressurehpturbineinlet[i] = M1P
            pressureipturbineinlet[i] = M2P
            pressurelpturbineinlet[i] = M4P
            pressurecondenseroutlet[i] = M6P
        else:
           break

x = hour2
y1 = hpturbineinlettemp
y2 = ipturbineinlettemp
y3 = lpturbineinlettemp
y4 = condenseroutlettemp
fig, ax1 = plt.subplots()
curve1, = ax1.plot(x, y1, label="HP Inlet", color='r')
curve2, = ax1.plot(x, y2, label="IP Inlet", color='b')
curve3, = ax1.plot(x, y3, label="LP Inlet", color='g')
curve4, = ax1.plot(x, y4, label="Condenser Outlet", color='c')
curves = [curve1, curve2, curve3, curve4]
ax1.legend(curves, [curve.get_label() for curve in curves])
ax1.set_xlabel("Time of Day (hr)", color='k')
ax1.set_ylabel("Temperature (K)", color='k')
plt.title("Temperatures with Demand Curve Variation")
plt.show()

x = hour2
y1 = massflowhpturbineinlet
y2 = massflowipturbineinlet
y3 = massflowlpturbineinlet
y4 = massflowcondenseroutlet
fig, ax1 = plt.subplots()
curve1, = ax1.plot(x, y1, label="HP Inlet", color='r')
curve2, = ax1.plot(x, y2, label="IP Inlet", color='b')
curve3, = ax1.plot(x, y3, label="LP Inlet", color='g')
curve4, = ax1.plot(x, y4, label="Condenser Outlet", color='c')
curves = [curve1, curve2, curve3, curve4]
ax1.legend(curves, [curve.get_label() for curve in curves])
ax1.set_xlabel("Time of Day (hr)", color='k')
ax1.set_ylabel("Mass Flow (kg/s)", color='k')
plt.title("Mass Flow Rates with Demand Curve Variation")
plt.show()

x = hour2
y1 = pressurehpturbineinlet
y2 = pressureipturbineinlet
y3 = pressurelpturbineinlet
y4 = pressurecondenseroutlet
fig, ax1 = plt.subplots()
curve1, = ax1.plot(x, y1, label="HP Inlet", color='r')
curve2, = ax1.plot(x, y2, label="IP Inlet", color='b')
curve3, = ax1.plot(x, y3, label="LP Inlet", color='g')
curve4, = ax1.plot(x, y4, label="Condenser Outlet", color='c')
curves = [curve1, curve2, curve3, curve4]
ax1.legend(curves, [curve.get_label() for curve in curves])
ax1.set_xlabel("Time of Day (hr)", color='k')
ax1.set_ylabel("Pressure (MPa)", color='k')
plt.title("Pressures with Demand Curve Variation")
plt.show()
