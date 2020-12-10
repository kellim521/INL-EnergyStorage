# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:19:07 2020

@author: Yugi
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import ode
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math
import pandas as pd
import itertools
#from RandomMs import SolarM
'''
    
    Parameters
    ----------
    params : 
    alpha= #average neutron life time
    lamb= #decay constant
    beta= #delayed neutron fraction
    c_pf= #specific heat of moderator
    c_pc= #specific heat of coolant
    m_f=#mass of fuel
    m_c=#mass of coolant
    mdot_h=#mass flow rate
    T_cine=#Temperature in
    a_f=#neutrons to thermal factor
    n_e= #scaling for power at equillibrium
    alpha_f= #Change in reactivity based on temp of fuel
    alpha_c=#Change in reactivity based on temp for moderator
    h=#heat transfer coefficient and total area of fuel
    Rf =#fouling factor
    Returns
    -------
    derivs : change in neutron density, precursor density, fuel temp, coolant temp
'''
"""specific outputs of the reactor and heatexchanger"""
Power = [] #MWth, power of reactor over transient period
Powere =[] #MWe, eletric power of reactor generated from 45% efficiency
InletTemp = [] # Kelvin, Variation of inlet temperature
finaltemps = []# Kelvin, Variation of outlet temperature
MdotO = []# m/s, Variation of velocity
Htransfer = [] # W/K, Variation of the heat transfer coefficient
Q = [] # W, Variation of heat transfer accross the heat exchanger
NeutronD = [] #normalized neutron density
Vprofile=[] # m/s profile of velocity over 12 sections in 1m increments
Delp = [] #Pa Pressure Drop over the reactor with respect to cycles
FuelT =[]# K, Equilibrium Temperature of the fuel
V = []# velocity at each time step in the reactor

"""initialization parameters of the reactor active core region based off of 
steady state calculations without thermohydraulic considerations"""
IDfuelContainer = .35 #meter of the inner diameter fuel tank container
ODfuelcontainer =1.05 # meter of the outerdiameter fuel container
HydD = (4*(math.pi/4)*(ODfuelcontainer**2 - IDfuelContainer**2))/(math.pi*(IDfuelContainer+ODfuelcontainer)) # meter the hydraulic diamter of the reactor
ACore =(math.pi/4)*(ODfuelcontainer**2 - IDfuelContainer**2) # m^2 the area of the core
Aout = .25 # m^2 the outlet arwa of the core
Npebbles = 441000
Hreflector = [0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0] #meters Hieght of the total reactor core for profiles
PebbleD = 0.03 # meters diameter of the fuel pebble
SurfaceP = 4*math.pi*(PebbleD/2)**2
TotSurf = SurfaceP*Npebbles
h = 16000*TotSurf# W/K heat transfer area considering pebble surface area, initial estimate
F1 = 0.05 # friction factor of the pebble bed reactor with a porosity of approximately .3-.4
K =1700# constant value when solving out the frictional pressure drop integration solved from annular system of equations
g = 9.81 #m/s^2 as the graviational constant
u0 = 0.65 # initialized inlet velocity of the core
porosity = .4 # constant porosity of the reactor
DelP = 66300 # initialized pressure drop over reactor core

T_cine = 600+273 # kelvin, initialized inlet temperature
T_h_in = 750+273# kelvin, initialized outlet temperature found from steady state conditions
mdot_h = 1319 # kg/s, initialized mass flow rate from stead state conditions
b = 1/T_cine # K^-1, thermal expansion coefficent based on the ilet temperature of the core
rho_h= 2518-0.406*(T_cine) #density in kg/m^3 in core
mu_h= 0.000116*np.exp(3755/T_cine) #viscosity in Pa*s in core
c_ph= 2415.78  #specific heat in J/kg*K in core
k_h= 0.629697+0.0005*(T_cine) #thermal conductivity in W/mK in core


def dydt(t,y, params):
    x, y, z_f, z_c=y
    alpha,lamb,beta,c_pf,c_pc,m_f,m_c,mdot_h,T_cine,a_f,n_e,alpha_f,alpha_c,h=params

    T_fe=T_cine+(1/(2*mdot_h*c_pc)+(1/h))*a_f*n_e #equillibrium of fuel temp
    T_ce=T_cine+(a_f*n_e/(2*mdot_h*c_pc)) #equillibrium of coolant temp
    u=(T_cine-T_cine)/T_cine
    w=(1300-mdot_h)/mdot_h
    Power = 1 #percentage
    p_c=(Power-(x))*10
    p=p_c+alpha_c*T_ce*z_c+alpha_f*T_fe*z_f
    
    
    dydt1 = -(beta*x/alpha)+(beta*y/alpha)+(p/alpha)+(p*x/alpha)
    dydt2 = (x-y)*lamb
    dydt3 = ((a_f*n_e*x)/(m_f*c_pf*T_fe))-(h*z_f/(m_f*c_pf))+(h*T_ce*z_c/(m_f*c_pf*T_fe))
    dydt4 = (h*T_fe*z_f/(m_c*c_pc*T_ce))-((2*c_pc*mdot_h+h)*z_c/(m_c*c_pc))+((2*mdot_h*T_cine*u)/(m_c*T_ce))
    -(2*mdot_h*w*(T_ce-T_cine)/(m_c*T_ce))-(2*mdot_h*w*z_c/m_c)+(2*mdot_h*T_cine*u*w/(m_c*T_ce))
    
    derivs=[dydt1, dydt2, dydt3, dydt4]

    return derivs, p

def tempoutput(params2):
    c_pc,mdot_h,T_cine,a_f,h,finaltempchange = params2
    T_fe = T_cine+(1/(2*mdot_h*c_pc)+(1/h))*a_f*n_e #equillibrium of fuel temp
    T_ce = T_cine+(a_f*n_e/(2*mdot_h*c_pc)) #equillibrium of coolant temp
    Tout = (T_fe-T_ce)/(Rf*mdot_h*c_pc) + T_ce + finaltempchange
    power = mdot_h*c_pc*(Tout-T_cine)
    return Tout, power, T_cine, T_fe, T_ce

"""takes data from the solar power data file for the mwe being produced by the solar power,
    the corresponding mass flow rate of the reactors heat exchanger on the cold side,
    and the reactivity table for ranging values of the power which change according to
    the amount of power being produced by the solar power data"""
#    
MyData = pd.read_excel("C:/Users/Yugi/Documents/Senior Design 2/Solar_Power_Data.xlsx", sheet_name = "CloudTransient(SolP)") #extracting values of solar power generated in 10 min steps
#MyData1 = pd.read_excel("C:/Users/Yugi/Documents/Senior Design 2/Solar_Power_Data.xlsx", sheet_name = "Sheet4")#extracting values of reactivity inserted for corresponding solar power 10 min steps
MyData2 =pd.read_excel("C:/Users/Yugi/Documents/Senior Design 2/Solar_Power_Data.xlsx", sheet_name = "CloudTransient(ReactorM)")# extracting mass flow rate values for each 10 min cycle on the cold side of the heat exchanger
MyData3 =pd.read_excel("C:/Users/Yugi/Documents/Senior Design 2/Solar_Power_Data.xlsx", sheet_name = "CloudTransient(SolM)")
#
df = MyData.values.tolist()
#df1 = MyData1.values.tolist()
df2 = MyData2.values.tolist()
df3 = MyData3.values.tolist()
#
mergedf = list(itertools.chain.from_iterable(df))
#mergedf1 =list(itertools.chain.from_iterable(df1))
mergedf2 = list(itertools.chain.from_iterable(df2))
mergedf3 = list(itertools.chain.from_iterable(df3))
#
SolarPower = mergedf[0:90]
#Reactivity =mergedf1[3:80:5]
mdot_cA =mergedf2[0:90]
mdot_solar = mergedf3[0:90]

for j in range(90):

#    if SolarPower[j] < 50:
#       n_e = 230.05 
#    if 50 <= SolarPower[j] <100:    
#       deln_e = 172.25*(Reactivity[0]**2) +204.19*Reactivity[0]+0.0517
#       n_e =200 + deln_e
#    if 100 <=SolarPower[j] <150:
#       deln_e = 172.25*(Reactivity[1]**2) +204.19*Reactivity[1]+0.0517
#       n_e =200 + deln_e       
#    if 150 <= SolarPower[j] <200:
#       deln_e = 172.25*(Reactivity[2]**2) +204.19*Reactivity[2]+0.0517
#       n_e =200 + deln_e        
#    if 200 <= SolarPower[j] <250:
#       deln_e = 172.25*(Reactivity[3]**2) +204.19*Reactivity[3]+0.0517
#       n_e =200 + deln_e        
#    if 250 <= SolarPower[j] < 300:
#       deln_e = 172.25*(Reactivity[4]**2) +204.19*Reactivity[4]+0.0517
#       n_e =200 + deln_e        
#    if 300 <= SolarPower[j] <350:
#       deln_e = 172.25*(Reactivity[5]**2) +204.19*Reactivity[5]+0.0517
#       n_e =200 + deln_e    
#    if 350 <= SolarPower[j] <400:
#       deln_e = 172.25*(Reactivity[6]**2) +204.19*Reactivity[6]+0.0517
#       n_e =200 + deln_e    
#    if 400<= SolarPower[j] < 450:
#       deln_e = 172.25*(Reactivity[7]**2) +204.19*Reactivity[7]+0.0517
#       n_e =200 + deln_e        
#    if 450<= SolarPower[j] < 500:
#       deln_e = 172.25*(Reactivity[8]**2) +204.19*Reactivity[8]+0.0517
#       n_e =200 + deln_e        
#    if 500 <= SolarPower[j] < 550:
#       deln_e = 172.25*(Reactivity[9]**2) +204.19*Reactivity[9]+0.0517
#       n_e =200 + deln_e        
#    if 550<= SolarPower[j] < 600:
#       deln_e = 172.25*(Reactivity[10]**2) +204.19*Reactivity[10]+0.0517
#       n_e =200 + deln_e        
#    if 600 <= SolarPower[j] <=700:
#       deln_e = 172.25*(Reactivity[11]**2) +204.19*Reactivity[11]+0.0517
#       n_e =200 + deln_e            
        
    if mdot_solar[j] < 90:
        deln_e = -117
        n_e = 200 +deln_e
    if 90 <= mdot_solar[j] <195:    
       deln_e = -127
       n_e =200 + deln_e
    if 195 <= mdot_solar[j] <295: 
       deln_e = -135
       n_e =200 + deln_e 
    if 295 <=mdot_solar[j] <395:
       deln_e = -142
       n_e =200 + deln_e       
    if 395<= mdot_solar[j] <495:
       deln_e = -155
       n_e =200 + deln_e        
    if 495 <= mdot_solar[j]<595:
       deln_e = -160
       n_e =200 + deln_e        
    if 595 <= mdot_solar[j] < 695:
       deln_e = -167
       n_e =200 + deln_e        
    if 695 <= mdot_solar[j] <795:
       deln_e = -175
       n_e =200 + deln_e    
    if 795 <= mdot_solar[j] <850:
       deln_e = -175
       n_e =200 + deln_e    
    if 850 <= mdot_solar[j] <895:
       deln_e = -190
       n_e =200 + deln_e    
    if 895 <= mdot_solar[j]<995:
       deln_e = -191
       n_e =200 + deln_e        
    '''Initial parameters of the core'''
    PebbleN = 441000 # number of pebbles in core
    alpha=0.001
    lamb=0.1
    beta=7.5*10**-3
    c_pf=717 #specific heat of graphite moderator
    c_pc=2414.7 #specific heat of FliBE
    m_f=PebbleN*(1.5/1000) #mass of u235 in 470,000 pellets
    m_c=90830.8 #mass of coolant
    a_f=7.0e6
    alpha_f=-5.4e-6 #Change in reactivity based on temp of fuel
    alpha_c=-1.8e-5 #Change in reactivity based on temp for moderator
    Rf = .0005 #fouling factor    
    
    params=[alpha,lamb,beta,c_pf,c_pc,m_f,m_c,mdot_h,T_cine,a_f,n_e,alpha_f,alpha_c,h]
    x0=0.0 #starting neutron pop
    y0=0.0 #starting precursors
    z_f0=1.0 #starting fuel temp try changing to 0
    z_c0=1.0 #starting moderator temp
    y0=[x0,y0,z_f0,z_c0]
    t0=0
    A=[]
    
    # Solver
    r = ode(dydt).set_integrator('dopri5', method='nsteps')
    r.set_initial_value(y0, t0).set_f_params(params)
    '''run each cycle for an arbitrary volume of fluid over 1 minute through 1 cycle'''
    t1 =1.0
    dt = 0.01
    T=[]
    while r.successful() and r.t < t1:
        r.integrate(r.t+dt)
        T=np.append(T,r.t)
        A=np.append(A,r.y)
    #print np.size(A)
    B= A.reshape(np.size(T),4)
    
    finaltempchange = sum(B[:,3])
    """ these append the values of needed variables after each ten minute cycle"""    
    params2 = [c_pc,mdot_h,T_cine,a_f,h,finaltempchange]    
    Power.append(tempoutput(params2)[1]/1e6)
    Powere.append(tempoutput(params2)[1]/(1e6*2.222))
    FuelT.append(tempoutput(params2)[3])
    InletTemp.append(tempoutput(params2)[2])
    finaltemps.append(tempoutput(params2)[0])
    MdotO.append(mdot_h)
    Htransfer.append(h)
    NeutronD.append(B[:,0])
    '''Solving for the velocity, Total Heat transfer ocefficinet in the core, 
       and mass flow rate of the reactor'''   
    

    Vavg = math.sqrt((g*b*Hreflector[9]**2*(finaltemps[j]-T_cine)*rho_h*ACore*u0)/(Hreflector[9]*(1+K+(F1*Hreflector[9]/PebbleD)))) # average velocity through core solving continuity equation w/ assumptions
    
    Recore = (Vavg*PebbleD*rho_h)/(mu_h) #Re number through porous bed
    Prcore = (mu_h*c_ph)/k_h # the prandtl number of the coolant
    Nucore = 2+1.1*(Prcore**(1/3))*(Recore**(.6)) # nusselt value for Re>50 por=.35-.4
    h_FLiBe = (6*(1-porosity)*k_h*Nucore*TotSurf)/PebbleD # W/(K)Thermal convection Coefficient
    h=h_FLiBe #setting H value in ODE to the h value determined from fluid mechanics
    mdot_h= Vavg*rho_h*Aout #mass flow rate in kg/s 
    
    V.append(Vavg)
    
    """Setting the input and output temperatures of the cold side heat exchanger"""
    T_h_in= finaltemps[j] #ALL TEMPERATURES LISTED IN KELVIN, K 
    T_c_in= 450
    T_c_out= 700
    
    
    #FLUID PROPERTIES -- FLiBe (shellside hot)
    
    rho_h= 2518-0.406*(T_h_in+T_c_in)/2 #density in kg/m^3
    mu_h= 0.000116*np.exp(3755/((T_h_in+T_c_in)/2)) #viscosity in Pa*s
    c_ph= 2415.78  #specific heat in J/kg*K
    k_h= 0.629697+0.0005*((T_h_in+T_c_in)/2) #thermal conductivity in W/mK
    b = 1/T_cine
    
    ''' Gives the pressure drop and velocity profile over each cycle'''
    DelP = (rho_h/2)*(Vavg**2 - u0**2) + rho_h*g*Hreflector[9] - rho_h*(1-b*(finaltemps[j]-T_cine))
    F1 = (DelP/Hreflector[9])*(PebbleD**2/Vavg**2)*(porosity**3/(1-porosity)**2)
    Delp.append(DelP)
#    for i in range(9):
#        Vh = math.sqrt((g*b*Hreflector[i]**2*(finaltemps[j]-T_cine)*rho_h*ACore*u0)/(Hreflector[i]*(1+K+((DelP/Hreflector[i])*(PebbleD**2/Vavg**2)*(porosity**3/(1-porosity)**2)*Hreflector[i]/PebbleD)))) # average velocity through core solving continuity equation w/ assumptions
#        Vprofile.append(Vh)
    
    CHOICE= 1
    #FLUID PROPERTIES -- Solar Salt (tubeside cold)
        
    rho_c= 1804
    mu_c= 0.00169
    c_pc= 1520 
    k_c= 0.53
    Pr_c= 4.85
        
    #TUBE PROPERTIES 
        
    d_o= 0.02 #outer tube diameter in m 
    t_w= 0.001 #tube wall thickness 
    d_i= d_o-2*t_w #inner tube diameter 
        
    #GUESSES 
        
    U= 100
    U_guess= 200 #Overall HT Coefficient in W/m^2*K
    v_tube_guess= 1.5  #Tube velocity in m/s 
    #Energy Balance 
    """ setting values of the cold side heat exchanger mass flow rate"""
    #mass flow rate in kg/s 
    mdot_c=mdot_cA[j]
    Qdot= (mdot_cA[j])*c_pc*(T_c_out-T_c_in)
    T_h_out= T_h_in-mdot_c*c_pc*(T_c_out-T_c_in)/(mdot_h*c_ph)
    T_cine = T_h_out
    Q.append(Qdot)

#Vprof= np.reshape(Vprofile, (145,9)).T
#Vprof1 = Vprof[0:120:1]
#Vprof2 = Vprof[0:120:1]
#Vprof3 = Vprof[0:120:1]
#Vprof4 = Vprof[0:120:1]
#Vprof5 = Vprof[0:120:1]
#Vprof6 = Vprof[0:120:1]
#Vprof7 = Vprof[0:120:1]
#Vprof8 = Vprof[0:120:1]
#Vprof9 = Vprof[0:120:1]
#Vprof10 = Vprof[0:120:1]
#Vprof11 = Vprof[0:120:1]
#Vprof12 = Vprof[0:120:1]
#Vprof13 = Vprof[0:120:1]

np.savetxt("MassFlowRateSolar(ColdSide).csv", mdot_solar) 
np.savetxt("MassFlowRateReactor(ColdSide).csv", mdot_cA) 

'''plotting the inlet and outlet temperature of the reactor according to changes in power and mass flow rate'''
plt.figure(1)    
plt.plot(InletTemp, label = 'Inlet')  
plt.plot(finaltemps, label = 'Outlet')
plt.xlabel('Cycles')
plt.ylabel('Temperature(Kelvin)')
plt.title('Inlet & Outlet Temperature of Reactor')
ax = plt.subplot(111) 
ax.legend()
'''plotting the power change over each cycle of 20 minutes'''
plt.figure(2)  
plt.plot(Power, label= ' Reactor(MWth)')
plt.plot(SolarPower, label= ' Solar(MWe)')
plt.xlabel('Cycles=(10min/cycle)')
plt.ylabel('Power(MW)')
plt.title('Power')
ax = plt.subplot(111) 
ax.legend()
"""plots the average mass flow rate in the core"""
plt.figure(3)  
plt.plot(MdotO, label= 'Reactor ')
plt.plot(mdot_solar, label='Solar')
plt.plot(mdot_cA, label = 'Reactor IHE' )
plt.xlabel('Cycles')
plt.ylabel('mdot(kg/s)')
plt.title('Mass Flow rate')
ax = plt.subplot(111) 
ax.legend()
"""plots the heat transfer coefficient in the core"""
plt.figure(4)  
plt.plot(Htransfer, label= 'Reactor H')
plt.xlabel('Cycles')
plt.ylabel('h(W/K)')
plt.title('Heat transfer coeff. for Total Surface Area')
ax = plt.subplot(111) 
ax.legend()
#"""plots the normalized neutron density in the core"""
#plt.figure(5)  
#plt.plot(NeutronD, label= 'Normalized neutron Density')
#plt.xlabel('Cycles')
#plt.ylabel('X0')
#plt.title('Neutron Density')
#ax = plt.subplot(111) 
#ax.legend()
""" plots the heat transfer of the intermediate heat exchanger attatched to thermal energy storage system"""
plt.figure(6)  
plt.plot(Q, label= 'IHE heat tansfer')
plt.xlabel('Cycles')
plt.ylabel('Q(W)')
plt.title('Heat Transferred from IHE into TES')
ax = plt.subplot(111) 
ax.legend()
"""plots the velocity profile of the core for each cycle"""
plt.figure(7)  
plt.plot(V, label= 'Average Velocity')
plt.xlabel('24 hours(1 min interval)')
plt.ylabel('V(m/s)')
plt.title('Velocity inside Reactor')
ax = plt.subplot(111) 
ax.legend()

#plt.figure(8)  
#plt.plot(Vprof)
#plt.xlabel('axial profile over H=3(m)')
#plt.ylabel('V(m/s)')
#plt.title('Velocity in Reactor over 24 hours')
#ax = plt.subplot(111) 
#ax.legend()
"""plots the pressure drop along the reactor core"""
plt.figure(9)  
plt.plot(Delp, label= 'pressure drop')
plt.xlabel('cycles')
plt.ylabel('DelP(Pa)')
plt.title('Pressure drop inside of reactor Core')
ax = plt.subplot(111) 
ax.legend()

"""plots the equilibrium fuel temperature of the core"""
plt.figure(10)  
plt.plot(FuelT, label= 'Fuel Temp.')
plt.xlabel('cycles')
plt.ylabel('Temperature(kelvin)')
plt.title('Equilibrium Fuel Temp. of Reactor')
ax = plt.subplot(111) 
ax.legend()