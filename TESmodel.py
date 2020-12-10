# -*- coding: utf-8 -*-
"""
TES Exergy Model

"""

import numpy as np
import math
from thermo.chemical import Chemical


class tmodel:
    def __init__(self, cod, peit, hotT):
        print('started tmodel')
        
        h = 0.001  # grid spacing
        k = 0.0001  # time step
        ka = 0.2  # Thermal capacity ratio(Gas/Solid)
        p = int(1/h)

        # Initializing Arrays
        th = np.zeros(p)  # Fluid temperature

        ths = np.zeros(p)  # Storage bed temperature

        pei = peit  # conjugate heat transfer-stanton number
        pd = np.zeros(p)  # diffusion part
        pdr = np.zeros(p)  # diffusion part-recovery
        thr = np.ones(p)
        tdr = np.zeros(p)
        td = np.zeros(p)
        dth = np.zeros(p)
        thrs = np.ones(p)
        thss = np.zeros(90000)
        thrss = np.zeros(90000)

        # storage
        c = 0
        while min(th) < 0.999:
            
            for i in range(p):
                if i == 0:
                    th[0] = 1  # inlet temperature of hot gas
                    td[i] = th[i] * hotT + 298
                    ths[i] = ths[i] + pei * k * (th[i] - ths[i]) * ka
                if 0 < i < (p - 1):
                    td[i] = th[i] * hotT + 298
                    pd[i] = cod * (k / h ** 2) * (th[i + 1] - 2 * th[i] + th[i - 1])
                    th[i] = th[i] - (k / h) * (th[i] - th[i - 1]) - pei * k * (th[i] - ths[i]) + pd[i] - (10**-4)*(th[i]-th[0])
                    dth[i] = pei * k * (th[i] - ths[i]) * ka
                    ths[i] = ths[i] + dth[i]
                if i == (p - 1):
                    th[i] = th[i - 1]
                    td[i] = th[i] * hotT + 298
                    ths[i] = ths[i] + pei * k * (th[i] - ths[i]) * ka
                    thss[c] = ths[i]
            c+=1
            print(c)
        
        # recovery
        u = 0
        while max(thr) > 0.001:
            
            for i in range(p):
                if i == 0:
                    thr[0] = 0  # Inlet temperature of cold gas
                    tdr[i] = thr[i] * hotT + 298
                    thrs[i] = thrs[i] + pei * k * (thr[i] - thrs[i]) * ka
                if 0 < i < (p - 1):
                    pdr[i] = cod * (k / h ** 2) * (thr[i + 1] - 2 * thr[i] + thr[i - 1])
                    thr[i] = thr[i] - (k / h) * (thr[i] - thr[i - 1]) - pei * k * (thr[i] - thrs[i]) + pdr[i] - (10**-4)*(th[i]-th[0])
                    tdr[i] = thr[i] * hotT + 298
                    thrs[i] = thrs[i] + pei * k * (thr[i] - thrs[i]) * ka
                if i == (p - 1):
                    thr[i] = thr[i - 1]
                    tdr[i] = thr[i] * hotT + 298
                    thrs[i] = thrs[i] + pei * k * (thr[i] - thrs[i]) * ka
                    thrss[u] = thr[i]
            u+=1
            print(u)
        self.th = th
        self.ths = ths
        
        self.thr = thr
        self.tdr = tdr
        self.td = td
        self.dth = dth
        self.thrs = thrs

        self.stin = thss
        self.stout = thrss
        self.sttime = c
        self.retime = u
        print('finished tmodel')
        
class exergy:    
    def __init__(self, ratio, TES):
        """
        Parameters
        ----------
        ratio : ratio of hot fluid temperature entering 
            heat exchanger over ambient temp

        TES: tmodel object

        """
        x = 0
        k = 0.01
    
        for i in range(TES.retime):
            x = x + (TES.stout[i] + (1 / ratio) * math.log((TES.stout[i]) * ratio + 1)) * k  # exergy integration
    
        y = x / (1.5 * k)  # fractional exergy recovery
    
        exev = y / TES.retime
    
        self.exergy_efficiency = exev-0.06
    
        
        
if __name__ == "__main__":
    
    hotT = 863

    alumina = Chemical('PubChem=14769', T = hotT)
    sodiumNitrate = Chemical('sodium nitrate', T = hotT)

    k = alumina.kl*.65 + sodiumNitrate.kl*.35
    rho = alumina.rho*.65 + sodiumNitrate.rho*.35
    Cp = alumina.Cp*.65 + sodiumNitrate.Cp*.35

    cod = (k/(rho*Cp))
    
    Th_in = 863
    ratio = Th_in/298
    
    peit = 2000
    
    TES = tmodel(cod, peit, hotT)
    exergy = exergy(ratio, TES)
    
    print('Exergy Efficiency: ', exergy.exergy_efficiency)
    
    st_in = TES.stin
    st_out = TES.stout
    
    thr = TES.thr
    tdr = TES.tdr
    td = TES.td
    dth = TES.dth
    thrs = TES.thrs
    
    fluid_temp = TES.th
    bed_temp = TES.ths