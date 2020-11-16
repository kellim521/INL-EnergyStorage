# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 19:15:31 2020

@author: kelli
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
        p = 100
        q = 108388
        o = q

        # Initializing Arrays
        th = np.zeros(p)  # Fluid temperature

        ths = np.zeros(p)  # Storage bed temperature

        pei = peit  # conjugate heat transfer-stanton number
        pe = 0.0003  # Peclet number
        pd = np.zeros(p)  # diffusion part
        pdr = np.zeros(p)  # diffusion part-recovery
        thr = np.ones(p)
        tdr = np.zeros(p)
        td = np.zeros(p)
        dth = np.zeros(p)
        thrs = np.ones(p)
        thss = np.zeros(70000)
        thrss = np.zeros(70000)

        """ dp = 0.002
        L = 5
        muf = 0.000036
        x = 0
        xp = 0
        xpr = 0 """
        epsi = 0.35
        """ fac1 = 25.4 * ((1 - epsi) ** 1.33) / ((1 - (1 - epsi) ** 0.33) * (1 - (1 - epsi) ** 0.66) ** 2)
        fac2 = (1 - (1 - epsi) ** 0.66) ** 2 """
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
                    th[i] = th[i] - (k / h) * (th[i] - th[i - 1]) - pei * k * (th[i] - ths[i]) + pd[i]
                    dth[i] = pei * k * (th[i] - ths[i]) * ka
                    ths[i] = ths[i] + dth[i]
                if i == (p - 1):
                    th[i] = th[i - 1]
                    td[i] = th[i] * hotT + 298
                    ths[i] = ths[i] + pei * k * (th[i] - ths[i]) * ka
                    thss[c] = ths[i]
            c = c + 1

        u = 0
        print('c = ', c)
        while max(thr) > 0.001:
            for i in range(p):
                if i == 0:
                    thr[0] = 0  # Inlet temperature of cold gas
                    tdr[i] = thr[i] * hotT + 298
                    thrs[i] = thrs[i] + pei * k * (thr[i] - thrs[i]) * ka
                if 0 < i < (p - 1):
                    pdr[i] = cod * (k / h ** 2) * (thr[i + 1] - 2 * thr[i] + thr[i - 1])
                    thr[i] = thr[i] - (k / h) * (thr[i] - thr[i - 1]) - pei * k * (thr[i] - thrs[i]) + pdr[i]
                    tdr[i] = thr[i] * hotT + 298
                    thrs[i] = thrs[i] + pei * k * (thr[i] - thrs[i]) * ka
                if i == (p - 1):
                    thr[i] = thr[i - 1]
                    tdr[i] = thr[i] * hotT + 298
                    thrs[i] = thrs[i] + pei * k * (thr[i] - thrs[i]) * ka
                    thrss[u] = thr[i]
            u = u + 1

        self.stin = thss
        self.stout = thrss
        self.sttime = c
        self.retime = u
        print('u = ', u)
        print('finished tmodel')
        
if __name__ == "__main__":
    
    hotT = 500

    alumina = Chemical('PubChem=14769', T = hotT)
    sodiumNitrate = Chemical('sodium nitrate', T = hotT)

    k = alumina.kl*.65 + sodiumNitrate.kl*.35
    rho = alumina.rho*.65 + sodiumNitrate.rho*.35
    Cp = alumina.Cp*.65 + sodiumNitrate.Cp*.35

    cod = (k/(rho*Cp))    #*1000
    
    peit = 2000
    tmodel(cod, peit, hotT)