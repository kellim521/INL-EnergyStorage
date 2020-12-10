# -*- coding: utf-8 -*-
"""
Thermal Energy Storage

"""

import numpy as np
    
class storage():
    def __init__(self, hotT, coldT, mass_flowrateHOT, mass_flowrateCOLD, numTanks, diameter, wallThickness, tankHeight=10, thermocline_position=5):
        """
        
        Parameters
        ----------
        hotT : Float
            Hot Storage Temp (K)
        coldT : Float
            Cold Storage Temp (K)
        mass_flowrateHOT : Float
            Mass flow rate for hot storage (kg/s)
        mass_flowrateCOLD : Float
            Mass flow rate for cold storage (kg/s)
        numTanks : int
            Number of tanks for system
        diameter : float
            Tank diameter (m)
        wallThickness : float
            Tank wall thickness (m)
        tankHeight : float, optional
            Tank height (m)
            The defualt is 10 m
        thermocline_position : float, optional
            Thermocline original position (m) 
            The default is 5 (thermocline in middle of tank)

        Returns
        -------
        None.

        """
        Nt = numTanks   # Number of tanks in system
        Od = diameter   # Outer Diameter of tank (m)
        Wall = wallThickness   # Wall thickness of tank (m)
        Id = Od - 2 * (Wall)   # Inner Diameter of tank (m)
        A = np.pi * Id**2 /4   # Area of one tank end (m^2)
        At = A * Nt   # Area of all tanks in system (m^2)
        L = tankHeight   # Height of tank (m)
        V = At * L   # Total volume of tanks (m^3)   
        t = 600   # 10 minutes (s)
        self.Area=At 
        
        Rhos = 1804  # Density of salt 'kg/m3'                                                                                                                                     #Power out 'MW'
        
        # Volumetric flow rates (m^3/s)
        volumetric_flowrateHOT = mass_flowrateHOT/Rhos
        volumetric_flowrateCOLD = mass_flowrateCOLD/Rhos
        
        # Thermocline
        tf = 0.05   # Size of Thermocline as fraction of tank height
        Thh = L * tf    # Height of the Thermocline (m)
        Thop = thermocline_position   # Thermocline original position from bottom of tank (m)
        Thtc = (hotT - coldT) / Thh   # Rate of change of temperature in thermocline 'K/m"
        
        # Thermocline Velocity (m/s)
        # Thermocline velocity is positive if moving 'up', negative if moving 'down'
        # More mass being stored in hot side than cold side, thermocline moves down                                                  
        thermocline_velocity = (volumetric_flowrateCOLD - volumetric_flowrateHOT) / At
            
        Thl = Thop + thermocline_velocity * t   #Location of Thermocline 'm'        

        self.thermocline_height = Thh
        self.thermocline_location = Thl
        
        # Cold/Hot Storage Volumes
        V_cold = Thl*At
        V_hot = V-V_cold
        
        # Energy
        K = 12                                                                         #Thermal Conductivity of Solid     
        Vtt = V * .64                                                                  #Total volume of filler
        r = .002                                                                       #Radius of filler
        Vpb = 4/3 * 3.14159265 * r**3                                                  #Volume of 1 ball
        Sapb = 4 * 3.14159265 * r**2                                                   #Surface area of 1 ball
        Nb = Vtt / Vpb                                                                 #Number of balls
        Sat = Nb * Sapb                                                                #Total Surface area
        Saic = Sat * tf                                                                #Surface area in contact with thermocline
        U = 1/(1/13+1/K)                                                               #Overall heat transfer rate
    
        # Material properties
        Rhof = 3943                                                                    #Density of alumina 'kg/m3'
        Cp = 1165                                                                      #Specific heat of alumina 'J/kgK'
        Vf = .44                                                                       #Void Fraction
        Cps = 1520  # specific heat of salt 'J/kg'
        
        Ed = Rhof * Cp                                                                 #Energy density of filler medium 'J/m3K'
        Eds = Rhos * Cps                                                               #Energy density of salt 'J/m3K'
        Edsy = (Ed * (1 - Vf)) + (Eds * (Vf))                                          #Combined energy density 'J/m3K'
        
        Edsymwh = Edsy / 3.6e9 * (hotT - coldT)                                        #Combined energy density 'MWh/m3'
        
        Vfill = At * (L - Thl)                                                         #Volume of tank currently full 'm3'
        
        Se = Edsymwh * Vfill                                                           #Energy currently stored in tank 'MWh'
        self.energy_stored = Se
        
        self.Q_salt = Cps * (hotT - coldT) * (mass_flowrateCOLD-mass_flowrateHOT) / 1000 / 1000
            
                
        # Heat Transfer    
        
        Aet = np.pi * Od**2 /4                                                         #Area of top of one tank 'm3'
        Aett = Aet * Nt / .9068996821                                                  #Area of top of hexagon 'm3'
        a = np.sqrt((2 * Aett) / (3 * np.sqrt(3)))                                     #Side length of the hexagon 'm'
        Aw = L * a                                                                     #Area of side wall 'm3'
        Tea = Aett + (Aw * 6)
        
        ks = 22.6                                                                      #Heat capacity of AISI 304 steel 'W/mK'
        Thos = (Od - Id) / 2                                                           #Thickness of steel wall 'm'
        Rdps = Thos / ks                                                               #Thermal resistance of steel tank 'm2K/W'
        
        ki = .1226                                                                     #Heat capacity of insulation
        Thoi = .1508                                                                   #Thickness of insulation 'm'
        Rdpi = Thoi / ki                                                               #Thermal resistance of insulation 'm2K/W'
        
        h = 50                                                                         #Convection Coefficient of air 'W/m2K'
        Rdpa = 1 / h                                                                   #Thermal resistance of air 'm2K/W'
        
        sigma = 5.67e-8                                                                #Sigma
        eps = .9                                                                       #Epsilon
        Tsurr = 298                                                                    #Tsurrounding 'K'
        Tsurf = 321.12                                                                 #Tsurface 'K'
        
        hr = sigma * eps * (Tsurf - Tsurr) * (Tsurf**2 + Tsurr**2)                     #hr radiation 'W/m2K'
        Rdpr = 1 / hr                                                                  #Thermal Resistance of radiation 'm2K/W'
        
        Rdpt = Rdps + Rdpi + ((( 1 / Rdpa) + (1 / Rdpr))**(-1))                        #Total thermal resistance 'm2K/W'
        
        qdp = (hotT - 298) / Rdpt                                                       #Heat flux 'W/m2'
        
        qt = qdp * Tea / 1000 / 1000                                                   #Heat transfer 'MW'
        
        Thcs = 13
        Thca = 11
        
        Tr = 1/11 + 1/13
        
        U = 1 / Tr
        Hlr = (hotT - coldT) / Se  #Heat loss rate 'K/hr'  