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

class ThermalEnergyStorageSys:
    def __init__(self):
        solarData = np.empty(0)
        reactorData = np.empty(0)
        
        hotT = 863
        coldT = 563
        
        # import solar and reactor mass flow rate data
        with open('Solar.csv', newline='') as Solar:
            reader = csv.DictReader(Solar)
            for row in reader:
                solarData = np.append(solarData, row)
                
        with open('Reactor.csv', newline='') as Reactor:
            reader = csv.DictReader(Reactor)
            for row in reader:
                reactorData = np.append(reactorData, row)
        
        # set steam generator mass flow rate
        m_SteamGen = 1500
        
        # creat array to hold stored energy amounts
        energy = np.array([])
        
        # assume thermocline is initially in average position
        thermocline_position = 5
        
        # loop through 7 days, with 1440 1-min intervals per day
        k=0
        for j in range(7):
            for i in range(1440):
                solar = solarData[i]        
                m_Solar = float(solar['Mass flow rate'])
                
                reactor = reactorData[i]        
                m_Reactor = float(reactor['Mass flow rate'])
                
                m_HOT = m_Solar + m_Reactor
                
                m_COLD = m_SteamGen
                #print(m_HOT)
                
                TankStorage = storage(hotT, coldT, m_HOT, m_COLD, 30, 2.5, 0.026, thermocline_position=thermocline_position)
                
                # show plot at specified times, in minutes
                # so k % 480 means when k is a multiple of 480 minutes (8 hours)
                # the plot will display. This number can be adjusted for how
                # often you wish to see the plot.
                if k % 480 == 0:
                    fig = plt.figure()
                    ax = fig.add_axes([0,0,1,1])
                    ax.bar(0, TankStorage.thermocline_location, 0.4, color='b')
                    ax.bar(0, TankStorage.thermocline_height, 0.4, bottom=TankStorage.thermocline_location, color='rebeccapurple')
                    ax.bar(0, 10-(TankStorage.thermocline_location+TankStorage.thermocline_height), 0.4,bottom=TankStorage.thermocline_location+TankStorage.thermocline_height, color='r')
                    ax.set_ylabel('Hot and Cold Storage')
                    ax.yaxis.set_major_formatter(mtick.PercentFormatter(10))
                    label1 = 'Day: '+str(j+1)
                    if i < 60:
                        label2 = 'Time: 12:'+str(i%60).zfill(2)+' AM'
                    elif 60 <= i < 720:
                        label2 = 'Time: '+str(int(i/60)).zfill(2)+':'+str(i%60).zfill(2)+' AM'
                    elif 720 <= i < 780:
                        label2 = 'Time: '+str(int(i/60)).zfill(2)+':'+str(i%60).zfill(2)+' AM'
                    else:
                        label2 = 'Time: '+str(int((i/60))-12).zfill(2)+':'+str(i%60).zfill(2)+' PM'
                    label3 = 'Stored Energy: ' + str(round(TankStorage.energy_stored)) +' MWh'
                    ax.set_xlabel(label1+'\n'+label2+'\n'+label3)
                    ax.set_title('Thermal Energy Storage Tank')
                    ax.set_xticks([])
                    ax.set_yticks([0,1,2,3,4,5,6,7,8,9,10])
                    filename = os.getcwd()+'\\graphs\\tank'+str(k)+'.png'
                    # Uncomment the following line to save the plots
                ### MUST be uncommented to create .gif file
                    #plt.savefig(filename, dpi=250, bbox_inches='tight')
                    plt.show()
                
                energy = np.append(energy, TankStorage.energy_stored)
                thermocline_position = TankStorage.thermocline_location
                area = TankStorage.Area
                k+=1
        
        # Uncomment the lines following 6 lines to save the TES animation as a .gif
        # frames = []
        # imgs = sorted(glob.glob('**/*.png', recursive=True), key=os.path.getmtime)
        # for i in imgs:
        #     nf = Image.open(i)
        #     frames.append(nf)
        
        # frames[0].save('tank.gif', format='GIF', append_images=frames[1:],save_all=True, duration=200, loop=0)
        
        # Show stored energy graph        
        fig2 = plt.figure()
        ax = fig2.add_axes([0,0,1,1])
        ax.plot(np.arange(10080), energy)
        ax.set_xlabel('Time (hr)')
        ax.set_ylabel('Energy Stored (MWh)')
        ax.set_title('Energy Stored over 1 week')
        ax.set_xlim(left=0,right=10090)
        plt.savefig('exergy-1week.png', dpi=250, bbox_inches='tight')
        plt.show()

        # Energy stored over 1 week:
        self.StoredEnergy = energy
        
        # Parameters for steam generator
        self.hotT = hotT
        self.coldT = coldT
        self.massflowrate = m_SteamGen

if __name__ == '__main__':
    TES = ThermalEnergyStorageSys()