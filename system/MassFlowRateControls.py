# -*- coding: utf-8 -*-
"""
Mass Flow Rate Control

"""
import csv
import numpy as np 

# From Steam Generator, into cold side of TES (kg/s)
mdot_1 = 1500

# From cold side to solar
mdot_2 = 0

# From solar to Hot side of TES
mdot_3 = 0

# From cold side of TES to reactor
mdot_4 = 0

# From reactor to hot side of TES
mdot_5 = 0

# From hot side to Steam Generator
mdot_6 = 1500

# import mass flow rates
solarData = np.empty(0)
reactorData = np.empty(0)

with open('Solar.csv', newline='') as Solar:
    reader = csv.DictReader(Solar)
    for row in reader:
        solarData = np.append(solarData, row)
        
with open('Reactor.csv', newline='') as Reactor:
    reader = csv.DictReader(Reactor)
    for row in reader:
        reactorData = np.append(reactorData, row)

# array to save tank volumes
hotTankVolumes = np.zeros(145)
coldtankVolumes = np.zeros(145)

solar = solarData[0]        
mdot_3 = float(solar['mdot_s'])
mdot_2 = mdot_3

reactor = reactorData[0]
mdot_5 = float(reactor['mdot_r'])
mdot_4 = mdot_5


# loops through a day in 10 min intervals based on solar data to determine 
# ideal reactor mass flow rate

for i in range(0,145):
    idealReactorData = [0,0]

    solar = solarData[i]        
    mdot_3 = float(solar['mdot_s'])
    mdot_2 = mdot_3
    
    idealReactorData[1] = mdot_6 - mdot_3
    idealReactorData[0] = i+1
    with open('IdealReactor.csv', 'a', newline='') as Reactor:
        writer = csv.writer(Reactor)
        writer.writerow(idealReactorData)
    
            





