# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 15:08:12 2020

@author: Yugi
"""

#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import math as math
from McCulloch_Richard import *

# 3D Plot Modules
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter



# :16,324s/^/#/ 
# :16,324s/^#//
########################
# Steady State Analysis#
########################
# 
# Usage: 
# s, u, T, dT, Re, Gr = project_univ(qw, Kfact, Circ)
#
#    Input     Definition     (s=scalar, b=boolean)
# =========================================================================
#     qw (s) = Wall heat flux in the core (W/m^2)
#  Kfact (s) = Parameter to artificially increase the Peclet number
#   Circ (b) = Circulation (True=Forward, False=Backward)
#
#   Output     Definition     (a=array, s=scalar)
# =================================================
#    s (a)  =  path length
#    u (s)  =  velocity 
#    T (a)  =  temperature 
#   dT (s)  =  maximum temperature difference
#   Re (s)  =  Reynold's number
#   Gr (s)  =  Grashoff number
#

qww=1e6

# FORWARD CIRCULATION TEMPERATURE PLOTS
s1, u1, T1, dT, Re, Gr = project_univ(qww, 1, True) # Forward Circulation, Pe=100
#s2, u2, T2, dT, Re, Gr = project_univ(qww, 10.0, True) # Forward Circulation, Pe=10
#s3, u3, T3, dT, Re, Gr = project_univ(qww, 100.0, True) # Forward Circulation, Pe=1

# Plot Settings
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Plot the Temperature Profile
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1) # rows, columns, axis number
ax1.plot(s1,T1,'b',linewidth=2.0) #,s2,T2,'k',s3,T3,':r'
ax1.plot([s1[100], s1[100]],[30.0, 300.0],'--k')
ax1.plot([s1[200], s1[200]],[30.0, 300.0],'--k')
ax1.plot([s1[300], s1[300]],[30.0, 300.0],'--k')
ax1.set_xlim([0.0, 1.44])
ax1.set_ylim([30.0, 220])
ax1.set_xlabel('Path Length (m)')
ax1.set_ylabel('Temperature ($^o$C)')
#ax1.legend(['$Pe_m$=100','$Pe_m$=10','$Pe_m$=1'], 
 #          loc=0, shadow=True, fontsize='x-large').get_frame().set_facecolor('#F2F2F2')
fig1.show()



## BACKWARD CIRCULATION TEMPERATURE PLOTS
#s1, u1, T1, dT, Re, Gr = project_univ(qww,   1.0, False) # Backward Circulation, Pe=100
#s2, u2, T2, dT, Re, Gr = project_univ(qww,   10.0, False) # Backward Circulation, Pe=10 #9461
#s3, u3, T3, dT, Re, Gr = project_univ(qww,   100.0, False) # Backward Circulation, Pe=1
#
## Plot the Temperature Profile
#fig2 = plt.figure()
#ax2 = fig2.add_subplot(1, 1, 1) # rows, columns, axis number
#ax2.plot(s1,T1,'--b',s2,T2,'k',s3,T3,':r',linewidth=2.0)
#ax2.plot([s2[100], s2[100]],[50.0, 300.0],'--k')
#ax2.plot([s2[200], s2[200]],[50.0, 300.0],'--k')
#ax2.plot([s2[300], s2[300]],[50.0, 300.0],'--k')
#ax2.set_xlim([0.0, 1.44])
#ax2.set_ylim([50.0, 500.0])
#ax2.set_xlabel('Path Length (m)')
#ax2.set_ylabel('Temperature ($^o$C)')
#ax2.set_title('Backward Circulation, $q_w$=1E+06')
#ax2.legend(['$Pe_m$=100','$Pe_m$=10','$Pe_m$=1'], 
#           loc=0, shadow=True, fontsize='x-large').get_frame().set_facecolor('#F2F2F2')
#fig2.show()




# VELOCITY, TEMPERATURE CHANGE and Gr VS Re PLOTS
N=100
uf=np.zeros(N)
#ub=np.zeros(N)
dTf=np.zeros(N)
#dTb=np.zeros(N)
Ref=np.zeros(N)
#Reb=np.zeros(N)
Grf=np.zeros(N)
#Grb=np.zeros(N)
##q=np.linspace(10.0,1.0E6,N)
q=np.linspace(10,10E3,N)
D=0.7
H_core=12
#qw=q/(np.pi*D*H_core)

for i in range(0,N):
    s, uf[i], T, dTf[i], Ref[i], Grf[i] = project_univ(q[i], 1.0, True)
    #s, ub[i], T, dTb[i], Reb[i], Grb[i] = project_univ(qw[i],  981.0, False)


# Plot the Velocity Profile
fig3 = plt.figure()
ax3 = fig3.add_subplot(1, 1, 1) # rows, columns, axis number
ax3.plot(q,uf,'k',linewidth=2.0)
#ax3.plot(q,ub,'--b',linewidth=2.0)
ax3.set_xlim([0.0, 7.0E3])
ax3.set_ylim([0.0, 0.01])
ax3.set_xlabel('Power (W)')
ax3.set_ylabel('Velocity (m/s)')
#ax3.legend(['Forward Circulation','Backward Circulation'], 
 #          loc=0, shadow=True, fontsize='x-large').get_frame().set_facecolor('#F2F2F2')
fig3.show()
#
#Plot the Temperature Change Profile
#fig4 = plt.figure()
#ax4 = fig4.add_subplot(1, 1, 1) # rows, columns, axis number
#ax4.plot(qww,dTf,'k',linewidth=2.0)
#ax4.plot(qww,dTb,'--b',linewidth=2.0)
#ax4.set_xlim([0.0, 7.0E3])
#ax4.set_ylim([0.0, 200.0])
#ax4.set_xlabel('Power (W)')
#ax4.set_ylabel('Maximum Temperature Change ($^o$C)')
#ax3.set_title('$Pe_m$=1,000')
#ax4.legend(['Forward Circulation','Backward Circulation'], 
#           loc=0, shadow=True, fontsize='x-large').get_frame().set_facecolor('#F2F2F2')
#fig4.show()
#
## Plot Gr vs Re
#fig5 = plt.figure()
#ax5 = fig5.add_subplot(1, 1, 1) # rows, columns, axis number
#ax5.loglog(Ref,Grf,'k',linewidth=2.0)
#ax5.loglog(Reb,Grb,'--b',linewidth=2.0)
#ax5.set_xlim([1.0E0, 1E8])
#ax5.set_ylim([1.0E4, 1.0E12])
#ax5.set_xlabel('Re')
#ax5.set_ylabel('Gr')
#ax5.grid(True, which="both")
##ax5.grid(True, which="both", ls="-")
#ax5.legend(['Forward Circulation','Backward Circulation'], 
#           loc=0, shadow=True, fontsize='x-large').get_frame().set_facecolor('#F2F2F2')
#fig5.show()




######################
## Stability Analysis#
######################
## 
## Usage: 
## s, u, T, dT, Re, Gr, Fw, omega_min = project_univ(qw, Kfact, Circ, omega_real, Nf, P_case)
##
##      Input       Definition     (s=scalar, b=boolean)
## =========================================================================
##         qw (s) = Wall heat flux in the core (W/m^2)
##      Kfact (s) = Parameter to artificially increase the Peclet number
##       Circ (b) = Circulation (True=Forward, False=Backward)
## omega_real (s) = Real component of the complex frequency
##         Nf (s) = Number of imaginary nodes (ie: linspace(-2,2,Nf)
##     P_case (s) = Plotting case (1=Nyquist, 2=Coarse Surface, 3=Fine Surface, 4=Optimum omega)
##
##      Output     Definition     (a=array, s=scalar, c=complex)
## =========================================================================
##         s (a) = path length
##         u (s) = velocity 
##         T (a) = temperature 
##        dT (s) = maximum temperature difference
##        Re (s) = Reynold's number
##        Gr (s) = Grashoff number
##        Fw (a) = F(omega) from the paper
## omega_min (c) = complex frequency corresponding to the minimum of F(omega)
##
#
## NYQUIST PLOTS
s1, u1, T1, dT, Re, Gr, Fw1, omega_min1 = project_univ(qww, 1.0,  True, 0, 1000, 1)  #  Forward Circulation, Pe=100
#s2, u2, T2, dT, Re, Gr, Fw2, omega_min2 = project_univ(3.3418E6, 100.0, False, 0, 1000, 1)  # Backward Circulation, Pe=100
##s1, u1, T1, dT, Re, Gr, Fw1, omega_min1 = project_univ(1E6, 100.0,  True, 0, 4000, 1)  #  Forward Circulation, Pe=100
##s2, u2, T2, dT, Re, Gr, Fw2, omega_min2 = project_univ(1E6, 100.0, False, 0, 4000, 1)  # Backward Circulation, Pe=100
#
#
## Forward Circulation
fig6 = plt.figure()
ax6 = fig6.add_subplot(1, 1, 1) # rows, columns, axis number
ax6.plot(Fw1[0:len(Fw1)/2].real,Fw1[0:len(Fw1)/2].imag,'k',linewidth=2.0)
ax6.plot(Fw1[len(Fw1)/2:].real,Fw1[len(Fw1)/2:].imag,'--b',linewidth=2.0)
##ax6.set_xlim([-8.0E4, 0.0])
##ax6.set_ylim([-2.5E5, 2.5E5])
ax6.set_xlabel('$F_{real}(\omega)$')
ax6.set_ylabel('$F_{imag}(\omega)$')
ax6.grid(True, which="both")
fig6.show()
#
#
## Backward Circulation
#fig7 = plt.figure()
#ax7 = fig7.add_subplot(1, 1, 1) # rows, columns, axis number
#ax7.plot(Fw2[0:len(Fw2)/2].real,Fw2[0:len(Fw2)/2].imag,'k',linewidth=2.0)
#ax7.plot(Fw2[len(Fw2)/2:].real,Fw2[len(Fw2)/2:].imag,'--b',linewidth=2.0)
##ax7.set_xlim([-8.0E4, 0.0])
##ax7.set_ylim([-2.5E5, 2.5E5])
#ax7.set_xlabel('$F_{real}(\omega)$')
#ax7.set_ylabel('$F_{imag}(\omega)$')
#ax7.grid(True, which="both")
#fig7.show()
#
#
#
## F(omega) Surface Plots
## Coarser Plot
#Nf=100;
#Freal=np.linspace(-5.0,5.0,Nf)
#Fcoord=np.zeros((Nf,Nf))
#Fmag=np.zeros((Nf,Nf))
#for i in range(Nf):
#    s, u, T, dT, Re, Gr, Fw, omega_min = project_univ(1.0E6, 100.0, True, Freal[i], Nf, 2)  #  Forward Circulation, Pe=100
#    Fmag[i,:]=np.abs(Fw)
#    if(i%(Nf/10)==0):
#        print "Completed loop {} out of {}".format(i,Nf)
#    
#[xg,yg]=np.meshgrid(Freal,Freal)
#np.savetxt('Fmagc.txt', Fmag.reshape(Nf*Nf))
#np.savetxt('Xmagc.txt',   xg.reshape(Nf*Nf))
#np.savetxt('Ymagc.txt',   yg.reshape(Nf*Nf))
#
#fig0c = plt.figure()
#ax0c = fig0c.gca(projection='3d')
#surf = ax0c.plot_surface(xg, yg, Fmag, cmap=cm.jet, linewidth=0, antialiased=False, shade=True, cstride=1, rstride=1)
#ax0c.view_init(azim=-150,elev=43)
#fig0c.show()
#
#ca,cb = surf.get_clim()
#
#
## Finer Plot
#Nf=100;
#Freal=np.linspace(-0.5,0.5,Nf)
#Fcoord=np.zeros((Nf,Nf))
#Fmag=np.zeros((Nf,Nf))
#for i in range(Nf):
#    s, u, T, dT, Re, Gr, Fw, omega_min = project_univ(1.0E6, 100.0, True, Freal[i], Nf, 3)  #  Forward Circulation, Pe=100
#    Fmag[i,:]=np.abs(Fw)
#    if(i%(Nf/10)==0):
#        print "Completed loop {} out of {}".format(i,Nf)
#    
#np.savetxt('Fmagf.txt', Fmag.reshape(Nf*Nf))
#np.savetxt('Xmagf.txt',   xg.reshape(Nf*Nf))
#np.savetxt('Ymagf.txt',   yg.reshape(Nf*Nf))
#    
#Fmag[Fmag>=1.0E5]=1.0E5
#
#[xg,yg] = np.meshgrid(Freal,Freal)
#fig0f = plt.figure()
#ax0f = fig0f.gca(projection='3d')
#ax0f.set_zlim(-1E4,10E4)
#surf = ax0f.plot_surface(xg, yg, Fmag, cmap=cm.jet, linewidth=0, antialiased=False, shade=True, cstride=1, rstride=1)
#surf.set_clim(ca, cb)
#ax0f.view_init(azim=-150,elev=43)
#fig0f.show()
#
#
#
#
#
#
#
#
#
## VELOCITY, TEMPERATURE CHANGE and Gr VS Re PLOTS
#N=100
#Ref=np.zeros(N)
#Reb=np.zeros(N)
#omega_minf=np.zeros((N,2))
#omega_minb=np.zeros((N,2))
#q=np.linspace(10.0,7.0E3,N)
#D=0.0525
#H_core=0.36
#qw=q/(np.pi*D*H_core)
#
#for i in range(0,N):
#    s, u, T, dT, Ref[i], Gr, Fw, omega_minf[i,:] = project_univ(qw[i], 1000.0,  True, 0, 1, 4)
#    s, u, T, dT, Reb[i], Gr, Fw, omega_minb[i,:] = project_univ(qw[i],  981.0, False, 0, 1, 4)
#
#
## Plot the Stability Minimums 
## Real Part
#fig9 = plt.figure()
#ax9 = fig9.add_subplot(1, 1, 1) # rows, columns, axis number
#ax9.plot(Ref,omega_minf[:,0],'k',linewidth=2.0)
#ax9.plot(Reb,omega_minb[:,0],'--b',linewidth=2.0)
##ax9.set_xlim([0.0, 1.0E6])
##ax9.set_ylim([-0.005, 0.005])
#ax9.set_xlabel('Re')
#ax9.set_ylabel('$F_{real}(\omega)$')
#ax9.set_title('$Pe_m$=1,000')
#ax9.legend(['Forward Circulation','Backward Circulation'], 
#           loc=0, shadow=True, fontsize='x-large').get_frame().set_facecolor('#F2F2F2')
#fig9.show()
#
## Complex Part
#fig10 = plt.figure()
#ax10 = fig10.add_subplot(1, 1, 1) # rows, columns, axis number
#ax10.semilogx(Ref,omega_minf[:,1],'k',linewidth=2.0)
#ax10.semilogx(Reb,omega_minb[:,1],'--b',linewidth=2.0)
##ax10.plot(Ref,omega_minf[:,1],'k',linewidth=2.0)
##ax10.plot(Reb,omega_minb[:,1],'--b',linewidth=2.0)
##ax10.set_xlim([0.0, 1.0E6])
##ax10.set_ylim([0.01, 1.0])
#ax10.set_xlabel('Re')
#ax10.set_ylabel('$F_{imag}(\omega)$')
#ax10.set_title('$Pe_m$=1,000')
#ax10.legend(['Forward Circulation','Backward Circulation'], 
#           loc=0, shadow=True, fontsize='x-large').get_frame().set_facecolor('#F2F2F2')
#fig10.show()
#
#
#
#
#
#

#########################
##  Stability  Analysis #
#########################
## 
## Usage: 
## E = project_uniq(qw, Kfact, Circ, H_hot)
##
##    Input     Definition     (s=scalar, b=boolean)
## =========================================================================
##     qw (s) = Wall heat flux in the core (W/m^2)
##  Kfact (s) = Parameter to artificially increase the Peclet number
##   Circ (b) = Circulation (True=Forward, False=Backward)
##
##   Output     Definition     (a=array, s=scalar)
## =================================================
##    s (a)  =  path length
##    u (s)  =  velocity 
##    T (a)  =  temperature 
##   dT (s)  =  maximum temperature difference
##   Re (s)  =  Reynold's number
##   Gr (s)  =  Grashoff number
##
#
## FIND STABILITY LIMIT FOR CHANGING LENGTHS
#N=8
#omega_minf=np.zeros((N,2))-1.0
#omega_minb=np.zeros((N,2))-1.0
#Ef=np.zeros((N,2))
#Eb=np.zeros((N,2))
#q=np.linspace(10.0,1.0E6,N)
#D=0.0625
#H_core=1.5240
#qw=q/(np.pi*D*H_core)
#
#for i in range(0,N):
#    s, u, T, dT, Ref[i], Gr, Fw, omega_minf[i,:], Ef[i,:] = project_univ(qw[i], 1000.0,  True, 0, 4000, 5)
#    s, u, T, dT, Reb[i], Gr, Fw, omega_minb[i,:], Eb[i,:] = project_univ(qw[i],  981.0, False, 0, 4000, 5)
#
#
## Plot the Stability Minimums 
## Real Part
#fig11 = plt.figure()
#ax11 = fig11.add_subplot(1, 1, 1) # rows, columns, axis number
#ax11.plot(Ref,omega_minf[:,0],'k',linewidth=2.0)
#ax11.plot(Reb,omega_minb[:,0],'--b',linewidth=2.0)
##ax11.set_xlim([0.0, 1.0E6])
##ax11.set_ylim([-0.005, 0.005])
#ax11.set_xlabel('Re')
#ax11.set_ylabel('$F_{real}(\omega)$')
#ax11.set_title('$Pe_m$=1,000')
#ax11.legend(['Forward Circulation','Backward Circulation'], 
#           loc=0, shadow=True, fontsize='x-large').get_frame().set_facecolor('#F2F2F2')
#fig11.show()














