# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 15:08:13 2020

@author: Yugi
"""

#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import math as math
import scipy.optimize as op

def project_univ(*args):
    qww=args[0]
    kfact=args[1]
    Circ=args[2]
    if(len(args)>3):
        Freal=args[3]
        Nf=args[4]
        p_case=args[5]
#    if(len(args)>1):
#        kfact=args(1)
#        if(len(args)>2):
#            Circ=args(2)
#            if(len(args)>3):
#                Freal=args(3)
#                if(len(args)>4):
#                    Nf=args(4)
#                    if(len(args)>5):
#                        p_case=args(5)

    # Universal Constants
    g=9.80665                   # gravity (m/s**2)
    
    # Loop Properties
    #D=0.0625
    D=0.7
    Dhx=D/1.88
                        # Diameter (m)
    H_core=3              # Core height (m)
    H_hot=3               # Hot leg height (m)
    L_hot=3                # Hot leg length (m)
    H_hex=H_core                # Steam generator height (m)
    H_cold=H_hot                # Cold leg height (m)
    L_cold=L_hot                # Cold leg (m)
    H_th=H_cold                 # Height between the center of the core and steam generator (m)
    L=L_hot+L_cold+H_core+H_hex # Loop length (m)
    Perim=np.pi*Dhx               # Wetted perimeter (m)
    Area=np.pi*Dhx**2/4.0         # Cross sectional area (m^2)

    # Reactor Properties
    qw=qww/(math.pi*D*H_core)       # Core power (W)

    # Coolant Properties 
    # Anonymous Functions
    # These equations are used in other codes
#    F_rho = lambda T: 10740.0-1.38*T                    # Density (kg/m^3)
#    F_mu = lambda T: 4.9E-4*math.exp(760.1/(T+273.15))  # Dynamic Viscosity (Pa-s)
#    F_kappa = lambda T: 11.16+0.01023*(T-150.0)         # Thermal Conductivity (W/m-K)
#    F_beta = lambda T: 1.38/(10740.0-1.38*T)            # Thermal Expansion Coefficient (1/K)
#    F_Cp = lambda T: 148.1                              # Specific Heat (J/kg-K)

    # These equations come from the Handbook of LBE
    #F_rho = lambda T: 11096-1.3236*(T+273.15)                           # Density (kg/m^3) 
    #F_mu = lambda T: 4.94E-4*math.exp(754.1/(T+273.15))                 # Dynamic Viscosity (Pa-s)
    #F_kappa = lambda T: 3.61+1.517E-2*(T+273.15)-1.741E-6*(T+273.15)**2 # Thermal Conductivity (W/m-K)
    #F_beta = lambda T: 1.0/(8383.2-(T+273.15))                          # Thermal Expansion Coefficient (1/K)
    #F_Cp = lambda T: 159-2.72E-2*(T+273.15)+7.12E-6*(T+273.15)**2       # Specific Heat (J/kg-K)
    #F_omega = lambda omega, M1r, M1i, FF, L, rho: np.abs(complex(M1r,M1i)-FF-L*rho*complex(omega[0],omega[1]))
    # F_omega is a function used to find the complex frequency, omega, that minimizes |F(omega)|
    ###GALLIUM WITH TEMPERATURE DEPENDENCE###
    F_rho = lambda T: 2518.0 - 0.406*(T)                           # Density (kg/m^3) 
    #F_mu = lambda T: 0.01207-5.754E-5*(T+273.15)+7.891E-8*(T+273.15)**2                 # Dynamic Viscosity (Pa-s)
    #F_kappa = lambda T: -7.448+0.1256*(T+273.15) # Thermal Conductivity (W/m-K)
    #F_beta = lambda T: 118.7E-6*(T+273.15)                        # Thermal Expansion Coefficient (1/K)
    #F_Cp = lambda T: 400    # Specific Heat (J/kg-K)
    #F_omega = lambda omega, M1r, M1i, FF, L, rho: np.abs(complex(M1r,M1i)-FF-L*rho*complex(omega[0],omega[1]))
    ###GALLIUM WITH NO TEMPERATURE DEPENDENCE###
    F_mu = lambda T: 0.000116*np.exp(3755/(T))                # Dynamic Viscosity (Pa-s)
    F_kappa = lambda T: 0.629697+0.0005*(T) # Thermal Conductivity (W/m-K)
    F_beta = lambda T: (.001145)                        # Thermal Expansion Coefficient (1/K)
    F_Cp = lambda T: 2415.78    # Specific Heat (J/kg-K)
    F_omega = lambda omega, M1r, M1i, FF, L, rho: np.abs(complex(M1r,M1i)-FF-L*rho*complex(omega[0],omega[1]))

    #Tw=108.9        # Tsat (C)
    Tw=1200
    Tref=Tw         # Temperature (C) 
    rho=F_rho(Tref) # Density (kg/m^3)
    Cp=F_Cp(Tref)   # Specific heat capacity (J/kg-K)
    u=.65          # Velocity (m/s)

    err=1.0
    tol=1E-5
    count=0.0
    while(err>tol and count<=100.0):
        count=count+1.0
        u_old=u
        
        # Lead Bismuth Eutectic
        rho=F_rho(Tref)             # Density (kg/m**3)                   
        mu=F_mu(Tref)               # Dynamic viscosity 
        kappa=F_kappa(Tref)*kfact   # Coefficient of axial heat conduction (W/m-C)
        beta=F_beta(Tref)           # VolumeTrefic thermal expansion coefficient of the fluid (1/C)
        Cp=F_Cp(Tref)               # Specific heat capacity (J/kg-C)
        Re=(rho*u*Dhx)/mu             # Reynolds number
        if(Re<2100.0):
            f=64.0/4.0*Re**(-1.0)   # Darcy friction factor for the laminar flow for the reactor core 
        elif(Re<30000.0):
            f=0.316/4.0*Re**(-0.25) # Darcy friction factor for the turbulent flow (Blasius) for the reactor core 
        else:
            f=0.184/4.0*Re**(-0.20) # Darcy friction factor for the turbulent flow (McAdams) for the reactor core 
        Pr=mu*Cp/(kappa/kfact)      # Prandt number
        Nu=0.625*(Re**0.4*Pr**0.4)  # Nusselt number 
        h=((Nu*kappa/kfact) / Dhx)    # Convective heat teansfer coefficient (W/m^2-K)
        Km=678                      # Form loss coefficient

        # Constants from the paper 
        N1=rho*Cp*u/kappa
        N2=qw/(rho*Cp*u)*4.0/D
        L1=0.5-math.sqrt(0.25+h*kappa*Perim/Area/(rho*Cp*u)**2.0) 
        L2=0.5+math.sqrt(0.25+h*kappa*Perim/Area/(rho*Cp*u)**2.0)
        
        # If A is singular (det==0) the solution falls apart and matrix inversion will not work.
        # This happens at around qw=1.0E7
        A=np.array([[1.0, 1.0, -math.exp(-N1*L_hot), -1.0, 0.0, 0.0, 0.0, 0.0],
           [N1, 0.0, -N1*math.exp(-N1*L_hot), 0.0, 0.0, 0.0, 0.0, 0.0],
           [0.0, 0.0, 1.0, 1.0, -1.0, -math.exp(-N1*L2*H_hex), 0.0, 0.0],
           [0.0, 0.0, N1, 0.0, -N1*L1, -N1*L2*math.exp(-N1*L2*H_hex), 0.0, 0.0],
           [0.0, 0.0, 0.0, 0.0, math.exp(N1*L1*H_hex), 1.0, -math.exp(-N1*L_cold), -1.0],
           [0.0, 0.0, 0.0, 0.0, N1*L1*math.exp(N1*L1*H_hex), N1*L2, -N1*math.exp(-N1*L_cold), 0.0],
           [-math.exp(-N1*H_core), -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
           [-N1*math.exp(-N1*H_core), 0.0, 0.0, 0.0, 0.0, 0.0, N1, 0.0]])
       
        B=np.array([[-N2*H_core],
           [-N2],
           [Tw],
           [0.0],
           [-Tw],
           [0.0],
           [0.0],
           [N2]],dtype=float)

        X=np.linalg.solve(A, B)

        if(Circ):   # Forward Circulation
            M=rho*beta*g*(N2*H_core**2/2+X[0]/N1+X[1]*H_core-X[0]*math.exp(-N1*H_core)/N1+X[2]*math.exp(-N1*(L_hot-H_hot))/N1 
                +X[3]*H_hot-X[2]*math.exp(-N1*L_hot)/N1-X[4]*math.exp(N1*L1*H_hex)/(N1*L1)-X[5]/(N1*L2) 
                -Tw*H_hex+X[4]/(N1*L1)+X[5]*math.exp(-N1*L2*H_hex)/(N1*L2)-X[6]*math.exp(-N1*(L_cold-H_cold))/N1 
                -X[7]*H_cold+X[6]*math.exp(-N1*L_cold)/N1)
        else:       # Backward Circulation
            M=rho*beta*g*(-N2*H_core**2/2-X[0]/N1-X[1]*H_core+X[0]*math.exp(-N1*H_core)/N1+X[2]*math.exp(-N1*(L_hot-L_cold))/N1 
                -X[2]*math.exp(-N1*(L_hot-L_cold+H_cold))/N1+X[3]*H_cold+X[4]*math.exp(N1*L1*H_hex)/(N1*L1) 
                +X[5]/(N1*L2)+Tw*H_hex-X[4]/(N1*L1)-X[5]*math.exp(-N1*L2*H_hex)/(N1*L1) 
                -X[6]*math.exp(-N1*(L_cold-L_hot))/N1-X[7]*L_hot+X[6]*math.exp(-N1*(L_cold-L_hot+H_hot))/N1+X[7]*(L_hot-H_hot))
        
        # Update the reference temperature
        # Tref=(T_core_out+T_SG_out)/2=(max(T)+min(T))/2
        Tref=((N2*H_core+X[0]+X[1])+(Tw+X[4]*math.exp(N1*L1*H_hex)+X[5]))/2.0
        u=math.sqrt(M/((2*(f*Km)/D)*rho))
        print(u)

        # Avoid divide by zero errors if u=0
        if(u!=0.0):
            err=abs(1.0-u_old/u)

        # Path length for each section
        s1=np.linspace(0.0,H_core,100)
        s2=np.linspace(H_core,H_core+L_hot,100)
        s3=np.linspace(H_core+L_hot,H_core+L_hot+H_hex,100)
        s4=np.linspace(H_core+L_hot+H_hex,H_core+L_hot+H_hex+L_cold,100)
        s=np.concatenate((s1,s2,s3,s4))

        # Temperature for each section
        T_core=N2*s1+X[0]*np.exp(-N1*(H_core-s1))+X[1]
        T_hot=X[2]*np.exp(-N1*(L_hot-(s2-H_core)))+X[3]
        T_SG=Tw+X[4]*np.exp(N1*L1*(s3-(L_hot+H_core)))+X[5]*np.exp(-N1*L2*(H_hex-(s3-(L_hot+H_core))))
        T_cold=X[6]*np.exp(-N1*(L_cold-(s4-(H_hex+L_hot+H_core))))+X[7]
        if(Circ):   # Forward Circulation
            T=np.concatenate((T_core,T_hot,T_SG,T_cold))
        else:       # Backward Circulation
            T=np.concatenate((T_core[::-1],T_cold[::-1],T_SG[::-1],T_hot[::-1]))

    # Final Steady State Parameters
    Pe=rho*Cp*u*D/kappa                 # Peclet Number
    Pe_mod=Pe*H_core/D                  # Modified Peclet Number
    dT=np.max(T)-np.min(T)              # Largest Temperature Difference (K)
    Gr=rho**2*g*beta*dT*H_th**3/mu**2   # Grashoff Number

    if(len(args)<4):
        return s, u, T, dT, Re, Gr


    ###################################
    # INSTABILITY ANALYSIS BEGINS HERE#
    ###################################
    # Define the B matrix to match equation variable names in the paper
    B=X

    # Declare the array of imaginary components based on the desired plot case
    if(p_case==1):      # Nyquist Plot
        Fimag=np.linspace(-2.0,2.0,Nf)
    elif(p_case==2):    # Coarse Surface Plot
        Fimag=np.linspace(-5.0,5.0,Nf)
    elif(p_case==3):    # Fine Surface Plot
        Fimag=np.linspace(-0.5,0.5,Nf)
    else:               # Optimization Plot
        Fimag=np.linspace(-1.0,1.0,Nf)

    # Alocate sufficient memory for F(omega)
    Fw=np.zeros(Nf, complex)

    # For each complex frequency calculate F(omega)
    for i in range(Nf):
        omega=complex(Freal, Fimag[i])
        alpha=kappa/(rho*Cp)            # Thermal Diffusivity (m^2/s)

        # Constants from the paper
        Y1=N1*L1
        Y2=N1*L2
        Y3=(u+np.sqrt(u**2+4.0*alpha*omega))/(2.0*alpha)
        Y4=(u-np.sqrt(u**2+4.0*alpha*omega))/(2.0*alpha)
        Y5=(u+np.sqrt(u**2+4.0*alpha*(omega+h*Perim/(Area*rho*Cp))))/(2.0*alpha)
        Y6=(u-np.sqrt(u**2+4.0*alpha*(omega+h*Perim/(Area*rho*Cp))))/(2.0*alpha)
        Z=h*Perim/(Area*rho*Cp)
        Z1=-alpha*N1**2+u*N1
        Z2=-alpha*Y1**2+u*Y1+Z
        Z3=-alpha*Y2**2+u*Y2+Z

        # Decompose the B matrix just to make entering the equations easier
        B1=B[0]
        B2=B[1]
        B3=B[2]
        B4=B[3]
        B5=B[4]
        B6=B[5]
        B7=B[6]
        B8=B[7]

        A=np.array([[-1.0, -np.exp(Y4[0]*H_core), np.exp(-Y3[0]*L_hot), 1.0, 0.0, 0.0, 0.0, 0.0],
                    [-Y3[0], -Y4[0]*np.exp(Y4[0]*H_core), Y3[0]*np.exp(-Y3[0]*L_hot), Y4[0], 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, -1.0, -np.exp(Y4[0]*L_hot), np.exp(-Y5[0]*H_hex), 1.0, 0.0, 0.0],
                    [0.0, 0.0, -Y3[0], -Y4[0]*np.exp(Y4[0]*L_hot), Y5[0]*np.exp(-Y5[0]*H_hex), Y6[0], 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, -1.0, -np.exp(Y6[0]*H_hex), np.exp(-Y3[0]*L_cold), 1.0],
                    [0.0, 0.0, 0.0, 0.0, -Y5[0], -Y6[0]*np.exp(Y6[0]*H_hex), Y3[0]*np.exp(-Y3[0]*L_cold), Y4[0]],
                    [np.exp(-Y3[0]*H_core), 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -np.exp(Y4[0]*L_cold)],
                    [Y3[0]*np.exp(-Y3[0]*H_core), Y4[0], 0.0, 0.0, 0.0, 0.0, -Y3[0], -Y4[0]*np.exp(Y4[0]*L_cold)]])


# This usage throws an error "cannot convert complex to float.
# It's an error in the NumPy module.
# To work around this pass the complex variable as an array, ie: Y[0]
#        A=np.array([[-1.0, -np.exp(Y4*H_core), np.exp(-Y3*L_hot), 1.0, 0.0, 0.0, 0.0, 0.0],
#                    [-Y3, -Y4*np.exp(Y4*H_core), Y3*np.exp(-Y3*L_hot), Y4, 0.0, 0.0, 0.0, 0.0],
#                    [0.0, 0.0, -1.0, -np.exp(Y4*L_hot), np.exp(-Y5*H_hex), 1.0, 0.0, 0.0],
#                    [0.0, 0.0, -Y3, -Y4*np.exp(Y4*L_hot), Y5*np.exp(-Y5*H_hex), Y6, 0.0, 0.0],
#                    [0.0, 0.0, 0.0, 0.0, -1.0, -np.exp(Y6*H_hex), np.exp(-Y3*L_cold), 1.0],
#                    [0.0, 0.0, 0.0, 0.0, -Y5, -Y6*np.exp(Y6*H_hex), Y3*np.exp(-Y3*L_cold), Y4],
#                    [np.exp(-Y3*H_core), 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -np.exp(Y4*L_cold)],
#                    [Y3*np.exp(-Y3*H_core), Y4, 0.0, 0.0, 0.0, 0.0, -Y3, -Y4*np.exp(Y4*L_cold)]])

        C=np.ones(8)

        C=np.array([[N1*B3*np.exp(-N1*L_hot)/(omega+Z1)-N2/omega-N1*B1/(omega+Z1)],
                    [N1**2*B3/(omega+Z1)*np.exp(-N1*L_hot)-N1**2*B1/(omega+Z1)],
                    [B5*Y1/(omega+Z2)+B6*Y2/(omega+Z3)*np.exp(-Y2*H_hex)-N1*B3/(omega+Z1)],
                    [B5*Y1**2/(omega+Z2)+B6*Y2**2/(omega+Z3)*np.exp(-Y2*H_hex)-N1**2*B3/(omega+Z1)],
                    [N1*B7/(omega+Z1)*np.exp(-N1*L_cold)-B5*Y1/(omega+Z2)*np.exp(Y1*H_hex)-B6*Y2/(omega+Z3)],
                    [N1**2*B7*np.exp(-N1*L_cold)/(omega+Z1)-B5*Y1**2*np.exp(Y1*H_hex)/(omega+Z2)-B6*Y2**2/(omega+Z3)],
                    [N2/omega+N1*B1/(omega+Z1)*np.exp(-N1*H_core)-N1*B7/(omega+Z1)],
                    [N1**2*B1*np.exp(-N1*H_core)/(omega+Z1)-N1**2*B7/(omega+Z1)]])

        X=np.dot(np.linalg.inv(A),C.reshape(8,1))

        # Decompose the C matrix to make entering the equations easier
        C1=X[0]
        C2=X[1]
        C3=X[2]
        C4=X[3]
        C5=X[4]
        C6=X[5]
        C7=X[6]
        C8=X[7]

        # Determine the friction factor and FF based on Reynold's number
        if(Re<2100.0):      # Laminar Flow
            f=(64.0/4.0)*Re**(-1.0)  
            FF=rho*u*f*(4.0*Km/D)*(2.0-1.0)/2.0
        elif(Re<30000.0):   # Turbulent (Blasius)
            f=(0.316/4.0)*Re**(-0.25)
            FF=rho*u*f*(4.0*Km/D)*(2.0-0.25)/2.0
        else:               # Turbulent (McAdams)
            f=(0.184/4.0)*Re**(-0.20)  
            FF=rho*u*f*(4.0*Km/D)*(2.0-0.20)/2.0

        # Evaluate the momentum integral
        if(Circ):   # Forward Circulation
            M1=(rho*beta*g)*( (C1/Y3)+(C2*np.exp(Y4*H_core)/Y4)-(N2*H_core/omega)-(B1/(Z1+omega))
                -(C1*np.exp(-Y3*H_core)/Y3)-(C2/Y4)+(B1*np.exp(-N1*H_core)/(Z1+omega))
                +(C3*np.exp(-Y3*(L_hot-H_hot))/Y3)+(C4*np.exp(Y4*H_hot)/Y4)
                -(B3*np.exp(-N1*(L_hot-H_hot))/(Z1+omega))-(C3*np.exp(-Y3*L_hot)/Y3)-(C4/Y4)
                +(B3*np.exp(-N1*L_hot)/(Z1+omega))-(C5/Y5)-(C6*np.exp(Y6*H_hex)/Y6)
                +(B5*np.exp(Y1*H_hex)/(Z2+omega))+(B6/(Z3+omega))+(C5*np.exp(-Y5*H_hex)/Y5)
                +(C6/Y6)-(B5/(Z2+omega))-(B6*np.exp(-Y2*H_hex)/(Z3+omega))
                -(C7*np.exp(-Y3*(L_cold-H_cold))/Y3)-(C8*np.exp(Y4*H_cold)/Y4)
                +(B7*np.exp(-N1*(L_cold-H_cold))/(Z1+omega))+(C7*np.exp(-Y3*L_cold)/Y3)
                +(C8/Y4)-(B7*np.exp(-N1*L_cold)/(Z1+omega)))

        else:       # Backward Circulation
            M1=(rho*beta*g)*((-C1/Y3)+(-C2*np.exp(Y4*H_core)/Y4)+(N2*H_core/omega)
                +(B1/(Z1+omega))+(C1*np.exp(-Y3*H_core)/Y3)+(C2/Y4)-(B1*np.exp(-N1*H_core)/(Z1+omega))
                +(C3*np.exp(-Y3*(L_hot-L_cold))/Y3)+(C4*np.exp(Y4*L_cold)/Y4)
                -(B3*np.exp(-N1*(L_hot-L_cold))/(Z1+omega))-(C3*np.exp(-Y3*(L_hot-L_cold+H_cold))/Y3)
                -(C4*np.exp(Y4*(L_cold-H_cold))/Y4)+(B3*np.exp(-N1*(L_hot-L_cold+H_cold))/(Z1+omega))
                +(C5/Y5)+(C6*np.exp(Y6*H_hex)/Y6)-(B5*np.exp(Y1*H_hex)/(Z2+omega))-(B6/(Z3+omega))
                -(C5*np.exp(-Y5*H_hex)/Y5)-(C6/Y6)+(B5/(Z2+omega))+(B6*np.exp(-Y2*H_hex)/(Z3+omega))
                -(C7*np.exp(-Y3*(L_cold-L_hot))/Y3)-(C8*np.exp(Y4*L_hot)/Y4)
                +(B7*np.exp(-N1*(L_cold-L_hot))/(Z1+omega))+(C7*np.exp(-Y3*(L_cold-L_hot+H_hot))/Y3)
                +(C8*np.exp(Y4*(L_hot-H_hot))/Y4)-(B7*np.exp(-N1*(L_cold-L_hot+H_hot))/(Z1+omega)))


        # The same issue as above with complex numbers is circumvented 
        # by using the complex variables as a slice from an array
        temp=M1[0]-FF-Km*rho*omega
        Fw[i]=temp[0]

        if(p_case>3):
            # This bit isn't working so well and I really just don't know why 
            #    and I don't have time to look into it further.
            # Find the optimum frequency to minimize F(omega)
#            omega_min=op.fmin(F_omega,np.array([omega.real, omega.imag]),
#                              args=(M1.real,M1.imag,FF,L,rho),
#                              disp=False,
#                              maxiter=1000,
#                              ftol=1.0E-6)

            omega_m=op.minimize(F_omega,np.array([omega.real, omega.imag]),args=(M1.real,M1.imag,FF,L,rho))
            omega_min=omega_m.x
        else:
            omega_min=np.zeros(2)

    temp=np.array([[Y1, Y2, Y3, Y4, Y5, Y6]])
#    temp=np.linalg.eigvals(A) # This doesn't work because the matrix is non-linear 
                               # (the eigenvalues are in the exponent of most of the non-zero terms)
    return s, u, T, dT, Re, Gr, Fw, omega_min


























# from my_fun import *
# project_uniq(33418.338,1000.0,True,0.0625) # controlling diameter
