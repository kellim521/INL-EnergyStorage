# -*- coding: utf-8 -*-
"""
function definitions for steam cycle components
"""

from iapws import IAPWS97 as steam

class preheatConditions():
    T = 0.0                 # temperature (K)
    P = 0.0                 # pressure (MPa)
    m = 0.0                 # mass flow rate at condensor inlet(kg/s)
    h = 0.0                 # specific enthalpy (kJ/kg)
    s = 0.0                 # specific entropy (kJ/kgK)
    condPressure = 0.0      # pressure at inlet of condensor (MPa)
    pumpEfficiency = 0.85   # pump isentropic efficiency
    TTD = 5/1.8             # terminal temperature difference
    DCA = 10/1.8            # drain cooler approach

# -----------------------------------------------------------------------------
def Turbine(Mi, Me, A, B, C, pE_A, pE_B, pE_C, preheatConditions, PR, load):
    
    efficiency = (0.86-0.72)/(1.0-0.6)*load + 0.51
    
    # find exit properties
    ses = Mi.s # isentropic exit entropy (kJ/kg/K)
    Me.P = Mi.P * PR
    # set rest of state with steam tables 
    hes = steam(P=Me.P, s=ses).h
    Me.h = Mi.h - efficiency * (Mi.h - hes)
    Me.T = steam(h=Me.h, P=Me.P).T
    Me.s = steam(h=Me.h, P=Me.P).s
    
    numExtract = 0
    if(A):
        numExtract += 1
    if(B):
        numExtract += 1
    if(C):
        numExtract += 1
    
    if(preheatConditions):
        if(numExtract == 3):
            
            P0 = preheatConditions.condPressure
            h0 = steam(P=P0, x=0).h
            s0 = steam(P=P0, x=0).s

            m2 = m3 = Mi.m            
            P1 = P2 = P3 = preheatConditions.P
            s1s = s0
            h1s = steam(s=s1s, P=P1).h
            h1 = (h1s - h0)/preheatConditions.pumpEfficiency + h0
            T1 = steam(P=P1, h=h1).T
            
            T4 = preheatConditions.T
            h4 = preheatConditions.h
            
            T2 = 1*(T4-T1)/3 + T1
            T3 = 2*(T4-T1)/3 + T1
            h2 = steam(T=T2,P=P2).h
            h3 = steam(T=T3,P=P3).h
            
            A.P = PA2 = steam(T=(T4+preheatConditions.TTD),x=0).P
            TA2 = T3 + preheatConditions.DCA
            hA2 = steam(T=TA2,P=PA2).h
            hAs = steam(P=A.P, s=ses).h
            A.h = Mi.h - efficiency * (Mi.h - hAs)
            A.T = steam(P=A.P, h=A.h).T
            A.s = steam(P=A.P, h=A.h).s
            A.m = mA2 = m3*(h4-h3)/(A.h-hA2)
            
            B.P = PB2 = steam(T=(T3+preheatConditions.TTD),x=0).P
            TB2 = T2 + preheatConditions.DCA
            hB2 = steam(T=TB2,P=PB2).h
            hBs = steam(P=B.P, s=ses).h
            B.h = Mi.h - efficiency * (Mi.h - hBs)
            B.T = steam(P=B.P, h=B.h).T
            B.s = steam(P=B.P, h=B.h).s
            B.m = (m2*(h3-h2)-mA2*(hA2-hB2))/(B.h-hB2)
            mB2 = B.m + mA2
            
            C.P = PC2 = steam(T=(T2+preheatConditions.TTD),x=0).P
            TC2 = T1 + preheatConditions.DCA
            hC2 = steam(T=TC2,P=PC2).h
            hCs = steam(P=C.P, s=ses).h
            C.h = Mi.h - efficiency * (Mi.h - hCs)
            C.T = steam(P=C.P, h=C.h).T
            C.s = steam(P=C.P, h=C.h).s
            C.m = (m2*(h2-h1)-mB2*(hB2-hC2))/(C.h-hC2)
                
            Me.m = Mi.m - A.m - B.m - C.m
            
            Power = Mi.m*(Mi.h-A.h) + (Mi.m-A.m)*(A.h-B.h) + (Mi.m-A.m-B.m)*(B.h-Me.h)
        else:
            print("ERROR: Turbine()")
    else:
        if(numExtract == 1):
            A.P = Me.P
            A.h = Me.h
            A.T = Me.T
            A.s = Me.s
            A.m = Mi.m * pE_A
            Me.m = Mi.m - A.m
            Power = Mi.m*(Mi.h-Me.h)
        elif(numExtract == 2):
            TH = steam(P=Mi.P,x=1).T
            TL = steam(P=Me.P,x=1).T
            TA = (TH-TL)/(numExtract) + TL
            A.P = steam(T=TA, x=1).P
            hAs = steam(P=A.P, s=ses).h
            A.h = Mi.h - efficiency * (Mi.h - hAs)
            A.T = steam(P=A.P, h=A.h).T
            A.s = steam(P=A.P, h=A.h).s
            A.m = Mi.m * pE_A
            B.P = Me.P
            B.h = Me.h
            B.T = Me.T
            B.s = Me.s
            B.m = Mi.m * pE_B
            Me.m = Mi.m - A.m - B.m
            Power = Mi.m*(Mi.h-A.h) + (Mi.m-A.m)*(A.h-Me.h)
        elif(numExtract == 3):
            TH = steam(P=Mi.P,x=1).T
            TL = steam(P=Me.P,x=1).T
            TA = 2*(TH-TL)/(numExtract) + TL
            A.P = steam(T=TA, x=1).P
            hAs = steam(P=A.P, s=ses).h
            A.h = Mi.h - efficiency * (Mi.h - hAs)
            A.T = steam(P=A.P, h=A.h).T
            A.s = steam(P=A.P, h=A.h).s
            A.m = Mi.m * pE_A
            TB = (TH-TL)/(numExtract) + TL
            B.P = steam(T=TB, x=1).P
            hBs = steam(P=B.P, s=ses).h
            B.h = Mi.h - efficiency * (Mi.h - hBs)
            B.T = steam(P=B.P, h=B.h).T
            B.s = steam(P=B.P, h=B.h).s
            B.m = Mi.m * pE_B
            C.P = Me.P
            C.h = Me.h
            C.T = Me.T
            C.s = Me.s
            C.m = Mi.m * pE_C
            Me.m = Mi.m - A.m - B.m - C.m
            Power = Mi.m*(Mi.h-A.h) + (Mi.m-A.m)*(A.h-B.h) + (Mi.m-A.m-B.m)*(B.h-Me.h)
        else:
            Me.m = Mi.m
            Power = Mi.m*(Mi.h - Me.h)
    
    return Power
# -----------------------------------------------------------------------------
def Reheat(Mi, Me, temp):
    
    Me.T = temp
    Me.P = Mi.P
    Me.m = Mi.m
    Me.h = steam(T=Me.T, P=Me.P).h
    Me.s = steam(T=Me.T, P=Me.P).s
    
    Q = Me.m * (Me.h - Mi.h)
    return Q
# -----------------------------------------------------------------------------
def FW_extracted(A, A2, B, B2, C, C2, preheatConditions):
    
    numExtract = 0
    if(A):
        numExtract += 1
    if(B):
        numExtract += 1
    if(C):
        numExtract += 1
    
    if(preheatConditions):
        P0 = preheatConditions.condPressure
        h0 = steam(P=P0, x=0).h
        s0 = steam(P=P0, x=0).s
        
        P1 = preheatConditions.P
        s1s = s0
        h1s = steam(s=s1s, P=P1).h
        h1 = (h1s - h0)/preheatConditions.pumpEfficiency + h0
        T1 = steam(P=P1, h=h1).T
        
        if(numExtract == 1):
            T2 = preheatConditions.T
            
            A2.P = steam(T=(T2+preheatConditions.TTD),x=0).P
            A2.T = T1 + preheatConditions.DCA
            A2.h = steam(T=A2.T,P=A2.P).h
            A2.s = steam(T=A2.T,P=A2.P).s
            A2.m = A.m
            
        elif(numExtract == 2):
            T3 = preheatConditions.T
            T2 = (T3-T1)/2 + T1
            
            A2.P = steam(T=(T3+preheatConditions.TTD),x=0).P
            A2.T = T2 + preheatConditions.DCA
            A2.h = steam(T=A2.T,P=A2.P).h
            A2.s = steam(T=A2.T,P=A2.P).s
            A2.m = A.m
            
            B2.P = steam(T=(T2+preheatConditions.TTD),x=0).P
            B2.T = T1 + preheatConditions.DCA
            B2.h = steam(T=B2.T,P=B2.P).h
            B2.s = steam(T=B2.T,P=B2.P).s
            B2.m = B.m + A2.m
            
        elif(numExtract == 3):
            T4 = preheatConditions.T
            T2 = 1*(T4-T1)/3 + T1
            T3 = 2*(T4-T1)/3 + T1
            
            A2.P = steam(T=(T4+preheatConditions.TTD),x=0).P
            A2.T = T3 + preheatConditions.DCA
            A2.h = steam(T=A2.T,P=A2.P).h
            A2.s = steam(T=A2.T,P=A2.P).s
            A2.m = A.m
            
            B2.P = steam(T=(T3+preheatConditions.TTD),x=0).P
            B2.T = T2 + preheatConditions.DCA
            B2.h = steam(T=B2.T,P=B2.P).h
            B2.s = steam(T=B2.T,P=B2.P).s
            B2.m = B.m + A2.m
            
            C2.P = steam(T=(T2+preheatConditions.TTD),x=0).P
            C2.T = T1 + preheatConditions.DCA
            C2.h = steam(T=C2.T,P=C2.P).h
            C2.s = steam(T=C2.T,P=C2.P).s
            C2.m = C.m + B2.m
    else:
        print("ERROR: no preheat conditions")
    
# -----------------------------------------------------------------------------    
def Condenser(F2, Mi, Me):
    
    if(F2):
        mMF = Me.m = Mi.m + F2.m
        hMF = (Mi.m*Mi.h + F2.m*F2.h)/mMF
    else:
        mMF = Me.m = Mi.m
        hMF = Mi.h
    
    Me.P = Mi.P
    Me.h = steam(P=Me.P, x=0).h
    Me.s = steam(P=Me.P, x=0).s
    Me.T = steam(P=Me.P, x=0).T
    
    Q = mMF*(hMF - Me.h)
    return Q
# -----------------------------------------------------------------------------    
def Pump(Mi, Me, efficiency):
    
    Me.m = Mi.m
    ses = Mi.s
    hes = steam(s=ses, P=Me.P).h
    Me.h = (hes - Mi.h)/efficiency + Mi.h
    Me.T = steam(P=Me.P, h=Me.h).T
    Me.s = steam(P=Me.P, h=Me.h).s
    
    Power = Me.m*(Me.h - Mi.h)
    return Power
# -----------------------------------------------------------------------------    
def FW_main(A2, B, B2, Mi, Me):
    
    if(A2):
        Me.h = Mi.h + (B.m*(B.h-B2.h)+A2.m*(A2.h-B2.h))/Mi.m
    else:
        Me.h = Mi.h + B.m*(B.h-B2.h)/Mi.m
    
    Me.P = Mi.P
    Me.T = steam(P=Me.P, h=Me.h).T
    Me.s = steam(P=Me.P, h=Me.h).s
# -----------------------------------------------------------------------------    
def Deaerator(A2, B, Mi, Me, deaerationCondensate):
    
    if(A2):
        Me.m = A2.m + B.m + Mi.m
        Me.h = (Mi.h*Mi.m + A2.h*A2.m + B.h*B.m)/Me.m
    else:
        Me.m = B.m + Mi.m
        Me.h = (Mi.h*Mi.m + B.h*B.m)/Me.m
    
    Me.P = deaerationCondensate.condPressure
    Me.T = deaerationCondensate.T = steam(P=Me.P, h=Me.h).T
    Me.s = deaerationCondensate.s = steam(P=Me.P, h=Me.h).s
    deaerationCondensate.h = steam(P=Me.P, h=Me.h).h
# -----------------------------------------------------------------------------        
def SteamGenerator(Mi, Me, efficiency):
    
    if(Me.m - 0.0001 <= Mi.m <= Me.m + 0.00001):
        Q = Mi.m * (Me.h - Mi.h)
        Q += Q*(1-efficiency)
    
    return Q
# -----------------------------------------------------------------------------







