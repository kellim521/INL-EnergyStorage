# -*- coding: utf-8 -*-
"""
Thermal Energy Storage System
"""

from storage import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import csv
import numpy as np
import os
from PIL import Image
import glob

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
        
m_SteamGen = 1500  
energy = np.array([])
thermocline_position = 5 
k=0
for j in range(7):
    for i in range(144):
        solar = solarData[i]        
        m_Solar = float(solar['mdot_s'])*0.7
        
        reactor = reactorData[i]        
        m_Reactor = float(reactor['mdot_r'])
        
        if 40 < i < 100:
            m_HOT = (m_Solar + m_Reactor)*0.92
            
        else:
            m_HOT = m_Solar + m_Reactor
        
        m_COLD = m_SteamGen
        
        TankStorage = storage(863, 563, m_HOT, m_COLD, 70, 3.5, 0.026, thermocline_position=thermocline_position)
        
        if k % 50 == 0:
            width = 0.35
            fig = plt.figure()
            ax = fig.add_axes([0,0,1,1])
            ax.bar(0, TankStorage.thermocline_location, 0.4, color='b')
            ax.bar(0, TankStorage.thermocline_height, 0.4, bottom=TankStorage.thermocline_location, color='rebeccapurple')
            ax.bar(0, 10-(TankStorage.thermocline_location+TankStorage.thermocline_height), 0.4,bottom=TankStorage.thermocline_location+TankStorage.thermocline_height, color='r')
            ax.set_ylabel('Hot and Cold Storage')
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(10))
            label1 = 'Day: '+str(j+1)
            if i < 6:
                label2 = 'Time: 12:'+str(i%6).zfill(2)+' AM'
            elif 6 <= i < 72:
                label2 = 'Time: '+str(int(i/6)).zfill(2)+':'+str(i%6).zfill(2)+' AM'
            elif 72 <= i < 78:
                label2 = 'Time: '+str(int(i/6)).zfill(2)+':'+str(i%6).zfill(2)+' AM'
            else:
                label2 = 'Time: '+str(int((i/6))-12).zfill(2)+':'+str(i%6).zfill(2)+' PM'
            label3 = 'Stored Energy: ' + str(round(TankStorage.energy_stored)) +' MWh'
            ax.set_xlabel(label1+'\n'+label2+'\n'+label3)
            ax.set_title('Thermal Energy Storage Tank')
            ax.set_xticks([])
            ax.set_yticks([0,1,2,3,4,5,6,7,8,9,10])
            filename = os.getcwd()+'\\graphs\\tank'+str(k)+'.png'
            plt.savefig(filename, dpi=250, bbox_inches='tight')
            plt.close()
        
        energy = np.append(energy, TankStorage.energy_stored)
        thermocline_position = TankStorage.thermocline_location
        area = TankStorage.Area
        k+=1

# Create .gif
frames = []
imgs = sorted(glob.glob('**/*.png', recursive=True), key=os.path.getmtime)
for i in imgs:
    nf = Image.open(i)
    frames.append(nf)

frames[0].save('tank.gif', format='GIF', append_images=frames[1:],save_all=True, duration=200, loop=0)

# Show stored energy graph        
fig2 = plt.figure(dpi=300)
ax = fig2.add_axes([0,0,1,1])
ax.plot(np.arange(1008), energy)
ax.set_xlabel('Time (hr)')
ax.set_ylabel('Energy Stored (MWh)')
ax.set_title('Energy Stored over 1 week')
ax.set_xlim(left=0,right=1009)
plt.show()


