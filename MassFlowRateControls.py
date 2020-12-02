# -*- coding: utf-8 -*-
"""
Mass Flow Rate Control

"""
import csv
import numpy as np
import math

class tank:
        
    def __init__(self, currentVolume, totalVolume):
        #density of solar salt (kg/m^3)
        self.rhoSalt = 1804
        self.volume_salt = currentVolume
        self.volume_tank = totalVolume
        self.timeOverflow = math.inf
        self.timeEmpty = math.inf
    
    def VolChange(self, mOut, mIn): 
        vOut = mOut/self.rhoSalt #m^3/s
        vIn = mIn/self.rhoSalt #m^3/s
        vStore = vIn - vOut #m^3/s
        
        if vStore > 0:
            self.timeOverflow = (self.volume_tank - self.volume_salt)/vStore

            # add the volume of salt that will accumulate in 10 mins
            self.volume_salt += vStore*600
        
        elif vStore < 0:
            self.timeEmpty = self.volume_salt/-vStore
            
            # remove the volume of salt that will accumulate in 10 mins
            self.volume_salt += vStore*600
 

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
mdot_6 = 0

# storage volumes for hot and cold side (m^3)
# currently just test values
hotTank = tank(1000, 6000)
coldTank = tank(3000, 6000)

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

# loops through a day in 10 min intervals based on solar and reactor data
for i in range(145):
    solar = solarData[i]        
    mdot_3 = float(solar['mdot_s'])
    mdot_2 = mdot_3
    
    reactor = reactorData[i]
    mdot_5 = float(reactor['mdot_r'])
    mdot_4 = mdot_5
    
    mdot_1 = 1500
    mdot_6 = mdot_3 + mdot_5
            
    if mdot_6 != mdot_1:
        
        # adjust flow out of hot side to match desired flow for steam generator
        mdot_6 = mdot_1
        
        # change volume of hot side
        hotTank.VolChange(mdot_6, mdot_3+mdot_5)
            
        # change volume of cold side
        coldTank.VolChange(mdot_2+mdot_4, mdot_1)
        
        hotTankVolumes[i] = hotTank.volume_salt
        coldtankVolumes[i] = coldTank.volume_salt
    
            





