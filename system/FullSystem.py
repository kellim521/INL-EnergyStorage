# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 10:46:06 2020

@author: Riley
"""

import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt
import PowerBlock_Riley as PB
import PowerBlock1_Riley as PB1
import PowerBlock2_Riley as PB2
import PowerBlock3_Riley as PB3
import SteamGenerator_Matt as SG
from thermo.chemical import Chemical

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

#From PowerBlock_Riley.py
#This line runs the PowerBlock_Riley.py file to get the values of every temperature.
M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,A,B,C,D,E,F,A2,B2,D2,E2,F2= PB.Cycle(1.0,155)

#These lines convert the temperature values to integers with two decimal places.
M1 = str(round(M1, 2))
M2 = str(round(M2, 2))
M3 = str(round(M3, 2))
M4 = str(round(M4, 2))
M5 = str(round(M5, 2))
M6 = str(round(M6, 2))
M7 = str(round(M7, 2))
M8 = str(round(M8, 2))
M9 = str(round(M9, 2))
M10 = str(round(M10, 2))
M11 = str(round(M11, 2))
M12 = str(round(M12, 2))
M13 = str(round(M13, 2))
M14 = str(round(M14, 2))
A = str(round(A, 2))
B = str(round(B, 2))
C = str(round(C, 2))
D = str(round(D, 2))
E = str(round(E, 2))
F = str(round(F, 2))
A2 = str(round(A2, 2))
B2 = str(round(B2, 2))
D2 = str(round(D2, 2))
E2 = str(round(E2, 2))
F2 = str(round(F2, 2))

#These lines are the code for the GUI. 
#The text2.insert on line 54 is where you can change what is being displayed on the right of the GUI.
root = tk.Tk()
text1 = tk.Text(root, height=50, width=120)
photo = tk.PhotoImage(file='./Power Block.png')
text1.insert(tk.END, '\n')
text1.image_create(tk.END, image=photo)
text1.pack(side=tk.LEFT)
text2 = tk.Text(root, height=50, width=35)
scroll = tk.Scrollbar(root, command=text2.yview)
text2.configure(yscrollcommand=scroll.set)
text2.tag_configure('big', font=('Arial', 15))
text2.insert(tk.END,"\nM1 = " + str(M1) + "\nM2 = " + str(M2) + "\nM3 = " + str(M3) + "\nM4 = " + str(M4) + "\nM5 = " + str(M5) + "\nM6 = " + str(M6) + "\nM7 = " + str(M7) + "\nM8 = " + str(M8) + "\nM9 = " + str(M9) + "\nM10 = " + str(M10) + "\nM11 = " + str(M11) + "\nM12 = " + str(M12) + "\nM13 = " + str(M13) + "\nM14 = " + str(M14) + "\nA = " + str(A) + "\nB = " + str(B) + "\nC = " + str(C) + "\nD = " + str(D) + "\nE = " + str(E) + "\nF = " + str(F) + "\nA2 = " + str(A2) + "\nB2 = " + str(B2) + "\nD2 = " + str(D2) + "\nE2 = " + str(E2) + "\nF2 = " + str(F2) + "\n", 'big')
text2.pack(side=tk.LEFT)
scroll.pack(side=tk.RIGHT, fill=tk.Y)
root.mainloop()
        