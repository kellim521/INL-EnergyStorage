"""
full cycle as a function
"""

import components
import SteamGenerator_Matt as SG
import Condenser
from iapws import IAPWS97 as steam

class state():
    T = 0.0 # inlet temperature (K)
    P = 0.0 # inlet pressure (MPa)
    m = 0.0 # mass flow rate (kg/s)
    h = 0.0 # specific enthalpy (kJ/kg)
    s = 0.0 # specific entropy (kJ/kgK)

def Cycle(load,massFlow):
    #--------------------------------------------------------------------------
    # Initial Design Conditions
    #--------------------------------------------------------------------------
    
    #load = 1.0                  # load as percentage of designed mass flow
    #massFlow = 450              # mass flow rate at 100% load (kg/s)
    
    mainTemp, V = SG.Hx3(SG.Mh, SG.Cph, SG.Y, SG.Mc, SG.Cps, SG.Tsi, SG.U, 3, 'Counter Flow')      # temperature at inlet of HP turbine (K)
    reheatTemp = mainTemp
    #mainTemp = 585 + 273
    #reheatTemp = 585 + 273      # temperature at inlet of IP turbine (K)
    
    highPressure = 19           # pressure at inlet of HP turbine (MPa)
    intPressure = 2             # pressure at inlet of IP turbine (MPa)
    lowPressure = 0.8           # pressure at inlet of LP turbine (MPa)
    condPressure = 0.01         # pressure at inlet of condensor  (MPa)
    
    pE_A = 10/100               # percent of mass flow extracted at A
    pE_B = 10/100               # percent of mass flow extracted at B
    pE_C = 5/100                # percent of mass flow extracted at C
    #pE_D = 000                 # percent of mass flow extracted at D
    #pE_E = 000                 # percent of mass flow extracted at E
    #pE_F = 000                 # percent of mass flow extracted at F
    
    TTD = 5/1.8                 # terminal temperature difference (K)
    DCA = 10/1.8                # drain cooler approach (K)
    deaeratorInletPressure = 0.2
    
    #--------------------------------------------------------------------------
    # Power Cycle
    #--------------------------------------------------------------------------
    
    # set inital state at the high pressure turbine inlet
    M1 = state()
    M1.T = mainTemp
    M1.P = highPressure
    M1.m = massFlow * load
    M1.h = steam(T=M1.T, P=M1.P).h
    M1.s = steam(T=M1.T, P=M1.P).s

    # set deaeration conditions for determining extraction flows
    deaerationConditions = components.preheatConditions()
    deaerationConditions.condPressure = condPressure
    deaerationConditions.P = deaeratorInletPressure
    deaerationConditions.T = steam(P=deaerationConditions.P,x=0).T
    deaerationConditions.h = steam(P=deaerationConditions.P,x=0).h
    deaerationConditions.s = steam(P=deaerationConditions.P,x=0).s
    deaerationConditions.TTD = TTD
    deaerationConditions.DCA = DCA
    
    deaeratorCondensate = components.preheatConditions()
    deaeratorCondensate.condPressure = 0.02
    deaeratorCondensate.P = M1.P
    deaeratorCondensate.TTD = -5/1.8
    deaeratorCondensate.DCA = 10/1.8
    
    # high pressure turbine states and function
    M2 = state()
    A = state()
    PR_HP = intPressure / M1.P   # pressure ratio of the high pressure turbine
    PowerHP = components.Turbine(M1,M2,A,None,None,pE_A,None,None,None,PR_HP,load)
    deaeratorCondensate.T = steam(P=A.P,x=0).T
    
    # reheat after high pressure turbine
    M3 = state()
    Qin_reheat = components.Reheat(M2, M3, reheatTemp)
    
    # intermediate pressure turbine states and function
    M4 = state()
    B = state()
    C = state()
    PR_IP = lowPressure / M3.P    # pressure ratio of the intermediate pressure turbine
    PowerIP = components.Turbine(M3,M4,B,C,None,pE_B,pE_C,None,None,PR_IP,load)
    
    # low pressure turbine states and function
    M5 = state()
    D = state()
    E = state()
    F = state()
    PR_LP = condPressure / M4.P    # pressure ratio of the low pressure turbine
    PowerLP = components.Turbine(M4,M5,D,E,F,None,None,None,deaerationConditions,PR_LP,load)
    
    # extracted steam through closed feedwater heaters
    A2 = state()
    B2 = state()
    D2 = state()
    E2 = state()
    F2 = state()
    components.FW_extracted(A,A2,B,B2,None,None,deaeratorCondensate)
    components.FW_extracted(D,D2,E,E2,F,F2,deaerationConditions)
    
    # [main + extracted steam from LP] though condenser
    M6 = state()
    Qout = components.Condenser(F2, M5, M6)
    M6.T = Condenser.condenser(M5.T)
    
    M7 = state()
    M7.P = deaerationConditions.P
    PowerPump1 = components.Pump(M6, M7, 0.85)
    
    # feedwater heaters for LP section
    M8 = state()
    M9 = state()
    M10 = state()
    M10.m = M9.m = M8.m = M7.m
    components.FW_main(E2,F,F2,M7,M8)
    components.FW_main(D2,E,E2,M8,M9)
    components.FW_main(None,D,D2,M9,M10)
    
    # deaerator combining extracted steam at C, 
    # combined feedwater from HP & IP, and main steam
    M11 = state()
    components.Deaerator(B2,C,M10,M11,deaeratorCondensate)
    
    # pump 2 to reach high pressure for the HP turbine
    M12 = state()
    M12.P = M1.P
    PowerPump2 = components.Pump(M11, M12, 0.85)
    
    # feedwater heating from HP and IP turbines
    M14 = state()
    M13 = state()
    M13.m=M14.m = M12.m
    components.FW_main(A2,B,B2,M12,M13)
    components.FW_main(None,A,A2,M13,M14)
    
    # heat in from the steam generator
    Qin_main = components.SteamGenerator(M14, M1, 0.95)
    
    # efficiency calculations
    Power = PowerHP+PowerIP+PowerLP-PowerPump1-PowerPump2
    Qin = Qin_reheat + Qin_main
    nth = Power/Qin
    nc = 1 - M6.T / M1.T

    
    return nth, Power/1000, Qin/1000   
    #values = Cycle(load,massFlow)
    #print(values)
    
def CycleQ(massFlow, load):
    #--------------------------------------------------------------------------
    # Initial Design Conditions
    #--------------------------------------------------------------------------
    
    #load = 1.0                  # load as percentage of designed mass flow
    #massFlow = 450              # mass flow rate at 100% load (kg/s)
    
    mainTemp =   585 + 273      # temperature at inlet of HP turbine (K)
    reheatTemp = 585 + 273      # temperature at inlet of IP turbine (K)
    
    highPressure = 19           # pressure at inlet of HP turbine (MPa)
    intPressure = 2             # pressure at inlet of IP turbine (MPa)
    lowPressure = 0.8           # pressure at inlet of LP turbine (MPa)
    condPressure = 0.01         # pressure at inlet of condensor  (MPa)
    
    pE_A = 10/100               # percent of mass flow extracted at A
    pE_B = 10/100               # percent of mass flow extracted at B
    pE_C = 5/100                # percent of mass flow extracted at C
    #pE_D = 000                 # percent of mass flow extracted at D
    #pE_E = 000                 # percent of mass flow extracted at E
    #pE_F = 000                 # percent of mass flow extracted at F
    
    TTD = 5/1.8                 # terminal temperature difference (K)
    DCA = 10/1.8                # drain cooler approach (K)
    deaeratorInletPressure = 0.2
    
    #--------------------------------------------------------------------------
    # Power Cycle
    #--------------------------------------------------------------------------
    
    # set inital state at the high pressure turbine inlet
    M1 = state()
    M1.T = mainTemp
    M1.P = highPressure
    M1.m = massFlow * load
    M1.h = steam(T=M1.T, P=M1.P).h
    M1.s = steam(T=M1.T, P=M1.P).s

    # set deaeration conditions for determining extraction flows
    deaerationConditions = components.preheatConditions()
    deaerationConditions.condPressure = condPressure
    deaerationConditions.P = deaeratorInletPressure
    deaerationConditions.T = steam(P=deaerationConditions.P,x=0).T
    deaerationConditions.h = steam(P=deaerationConditions.P,x=0).h
    deaerationConditions.s = steam(P=deaerationConditions.P,x=0).s
    deaerationConditions.TTD = TTD
    deaerationConditions.DCA = DCA
    
    deaeratorCondensate = components.preheatConditions()
    deaeratorCondensate.condPressure = 0.02
    deaeratorCondensate.P = M1.P
    deaeratorCondensate.TTD = -5/1.8
    deaeratorCondensate.DCA = 10/1.8
    
    # high pressure turbine states and function
    M2 = state()
    A = state()
    PR_HP = intPressure / M1.P   # pressure ratio of the high pressure turbine
    PowerHP = components.Turbine(M1,M2,A,None,None,pE_A,None,None,None,PR_HP,load)
    deaeratorCondensate.T = steam(P=A.P,x=0).T
    
    # reheat after high pressure turbine
    M3 = state()
    Qin_reheat = components.Reheat(M2, M3, reheatTemp)
    
    # intermediate pressure turbine states and function
    M4 = state()
    B = state()
    C = state()
    PR_IP = lowPressure / M3.P    # pressure ratio of the intermediate pressure turbine
    PowerIP = components.Turbine(M3,M4,B,C,None,pE_B,pE_C,None,None,PR_IP,load)
    
    # low pressure turbine states and function
    M5 = state()
    D = state()
    E = state()
    F = state()
    PR_LP = condPressure / M4.P    # pressure ratio of the low pressure turbine
    PowerLP = components.Turbine(M4,M5,D,E,F,None,None,None,deaerationConditions,PR_LP,load)
    
    # extracted steam through closed feedwater heaters
    A2 = state()
    B2 = state()
    D2 = state()
    E2 = state()
    F2 = state()
    components.FW_extracted(A,A2,B,B2,None,None,deaeratorCondensate)
    components.FW_extracted(D,D2,E,E2,F,F2,deaerationConditions)
    
    # [main + extracted steam from LP] though condenser
    M6 = state()
    Qout = components.Condenser(F2, M5, M6)
    
    
    M7 = state()
    M7.P = deaerationConditions.P
    PowerPump1 = components.Pump(M6, M7, 0.85)
    
    # feedwater heaters for LP section
    M8 = state()
    M9 = state()
    M10 = state()
    M10.m = M9.m = M8.m = M7.m
    components.FW_main(E2,F,F2,M7,M8)
    components.FW_main(D2,E,E2,M8,M9)
    components.FW_main(None,D,D2,M9,M10)
    
    # deaerator combining extracted steam at C, 
    # combined feedwater from HP & IP, and main steam
    M11 = state()
    components.Deaerator(B2,C,M10,M11,deaeratorCondensate)
    
    # pump 2 to reach high pressure for the HP turbine
    M12 = state()
    M12.P = M1.P
    PowerPump2 = components.Pump(M11, M12, 0.85)
    
    # feedwater heating from HP and IP turbines
    M14 = state()
    M13 = state()
    M13.m=M14.m = M12.m
    components.FW_main(A2,B,B2,M12,M13)
    components.FW_main(None,A,A2,M13,M14)
    
    # heat in from the steam generator
    Qin_main = components.SteamGenerator(M14, M1, 0.95)
    
    # efficiency calculations
    Power = PowerHP+PowerIP+PowerLP-PowerPump1-PowerPump2
    Qin = Qin_reheat + Qin_main
    nth = Power/Qin
    nc = 1 - M6.T / M1.T

    
    return Qin/1000
 
def CycleHP(load,massFlow):
    #--------------------------------------------------------------------------
    # Initial Design Conditions
    #--------------------------------------------------------------------------
    
    #load = 1.0                  # load as percentage of designed mass flow
    #massFlow = 450              # mass flow rate at 100% load (kg/s)
    
    mainTemp =   585 + 273      # temperature at inlet of HP turbine (K)
    reheatTemp = 585 + 273      # temperature at inlet of IP turbine (K)
    
    highPressure = 19           # pressure at inlet of HP turbine (MPa)
    intPressure = 2             # pressure at inlet of IP turbine (MPa)
    lowPressure = 0.8           # pressure at inlet of LP turbine (MPa)
    condPressure = 0.01         # pressure at inlet of condensor  (MPa)
    
    pE_A = 10/100               # percent of mass flow extracted at A
    pE_B = 10/100               # percent of mass flow extracted at B
    pE_C = 5/100                # percent of mass flow extracted at C
    #pE_D = 000                 # percent of mass flow extracted at D
    #pE_E = 000                 # percent of mass flow extracted at E
    #pE_F = 000                 # percent of mass flow extracted at F
    
    TTD = 5/1.8                 # terminal temperature difference (K)
    DCA = 10/1.8                # drain cooler approach (K)
    deaeratorInletPressure = 0.2
    
    #--------------------------------------------------------------------------
    # Power Cycle
    #--------------------------------------------------------------------------
    
    # set inital state at the high pressure turbine inlet
    M1 = state()
    M1.T = mainTemp
    M1.P = highPressure
    M1.m = massFlow * load
    M1.h = steam(T=M1.T, P=M1.P).h
    M1.s = steam(T=M1.T, P=M1.P).s

    # set deaeration conditions for determining extraction flows
    deaerationConditions = components.preheatConditions()
    deaerationConditions.condPressure = condPressure
    deaerationConditions.P = deaeratorInletPressure
    deaerationConditions.T = steam(P=deaerationConditions.P,x=0).T
    deaerationConditions.h = steam(P=deaerationConditions.P,x=0).h
    deaerationConditions.s = steam(P=deaerationConditions.P,x=0).s
    deaerationConditions.TTD = TTD
    deaerationConditions.DCA = DCA
    
    deaeratorCondensate = components.preheatConditions()
    deaeratorCondensate.condPressure = 0.02
    deaeratorCondensate.P = M1.P
    deaeratorCondensate.TTD = -5/1.8
    deaeratorCondensate.DCA = 10/1.8
    
    # high pressure turbine states and function
    M2 = state()
    A = state()
    PR_HP = intPressure / M1.P   # pressure ratio of the high pressure turbine
    PowerHP = components.Turbine(M1,M2,A,None,None,pE_A,None,None,None,PR_HP,load)
    deaeratorCondensate.T = steam(P=A.P,x=0).T
    
    # reheat after high pressure turbine
    M3 = state()
    Qin_reheat = components.Reheat(M2, M3, reheatTemp)
    
    # intermediate pressure turbine states and function
    M4 = state()
    B = state()
    C = state()
    PR_IP = lowPressure / M3.P    # pressure ratio of the intermediate pressure turbine
    PowerIP = components.Turbine(M3,M4,B,C,None,pE_B,pE_C,None,None,PR_IP,load)
    
    # low pressure turbine states and function
    M5 = state()
    D = state()
    E = state()
    F = state()
    PR_LP = condPressure / M4.P    # pressure ratio of the low pressure turbine
    PowerLP = components.Turbine(M4,M5,D,E,F,None,None,None,deaerationConditions,PR_LP,load)
    
    # extracted steam through closed feedwater heaters
    A2 = state()
    B2 = state()
    D2 = state()
    E2 = state()
    F2 = state()
    components.FW_extracted(A,A2,B,B2,None,None,deaeratorCondensate)
    components.FW_extracted(D,D2,E,E2,F,F2,deaerationConditions)
    
    # [main + extracted steam from LP] though condenser
    M6 = state()
    Qout = components.Condenser(F2, M5, M6)
    
    
    M7 = state()
    M7.P = deaerationConditions.P
    PowerPump1 = components.Pump(M6, M7, 0.85)
    
    # feedwater heaters for LP section
    M8 = state()
    M9 = state()
    M10 = state()
    M10.m = M9.m = M8.m = M7.m
    components.FW_main(E2,F,F2,M7,M8)
    components.FW_main(D2,E,E2,M8,M9)
    components.FW_main(None,D,D2,M9,M10)
    
    # deaerator combining extracted steam at C, 
    # combined feedwater from HP & IP, and main steam
    M11 = state()
    components.Deaerator(B2,C,M10,M11,deaeratorCondensate)
    
    # pump 2 to reach high pressure for the HP turbine
    M12 = state()
    M12.P = M1.P
    PowerPump2 = components.Pump(M11, M12, 0.85)
    
    # feedwater heating from HP and IP turbines
    M14 = state()
    M13 = state()
    M13.m=M14.m = M12.m
    components.FW_main(A2,B,B2,M12,M13)
    components.FW_main(None,A,A2,M13,M14)
    
    # heat in from the steam generator
    Qin_main = components.SteamGenerator(M14, M1, 0.95)
    
    # efficiency calculations
    Power = PowerHP+PowerIP+PowerLP-PowerPump1-PowerPump2
    Qin = Qin_reheat + Qin_main
    nth = Power/Qin
    nc = 1 - M6.T / M1.T
    
    return PowerHP/1000







