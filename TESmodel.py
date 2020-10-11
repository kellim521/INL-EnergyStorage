import math
import numpy as np

cod = 0.0002  # Thermal diffusivity
peit = 2000
hotT = 500


class P2S_HX:
    def __init__(self, m_dot_hot, c_p_hot, T_hot_in, m_dot_cold, c_p_cold, T_hot_out, U, A, HE_Type):
        C_hot = m_dot_hot * c_p_hot
        C_cold = m_dot_cold * c_p_cold
        C_min = min(C_hot, C_cold)
        C_max = max(C_hot, C_cold)
        C_r = C_min / C_max
        NTU = U * A / C_min
        epsilon = effectiveness(NTU, C_r, HE_Type)
        Q = C_hot * (T_hot_in - T_hot_out)
        Q_max = Q / epsilon
        T_cold_in = T_hot_in - Q_max / C_min
        T_cold_out = T_cold_in + Q / C_cold

        self.epsilon = epsilon
        self.T_cold_in = T_cold_in
        self.T_cold_out = T_cold_out


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

        # COD: Thermal Dispersion factor

        pei = peit  # conjugate heat transfer-stanton number
        # peo = 0.0  # Wall heat loss stanton number
        pe = 0.0003  # Peclet number

        pd = np.zeros(p)  # diffusion part

        pdr = np.zeros(p)  # diffusion part-recovery

        # pdp = np.zeros(p)

        # pdpr = np.zeros(p)

        # rhod = np.zeros(p)

        # red = np.zeros(p)

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


def effectiveness(NTU, C_r, HE_Type):
    if C_r == 0:
        epsilon = 1 - math.exp(-NTU)
        return epsilon

    if HE_Type == 'Parallel Flow':
        epsilon = (1 - math.exp(-NTU * (1 + C_r))) / (1 + C_r)
    elif HE_Type == 'Counter Flow':
        if C_r == 1:
            epsilon = NTU/(1+NTU)
        else:
            epsilon = (1-math.exp(-NTU*(1-C_r)))/(1-C_r*math.exp(-NTU*(1-C_r)))
    elif HE_Type == 'Double Shell Pass':
        N = 2.0
        NTUN = NTU/N
        epsilon1 = 2/(1+C_r+math.sqrt(1+C_r**2)*(1+math.exp(-NTUN*math.sqrt(1+C_r**2)))/(1+math.exp(-NTUN*math.sqrt(1+C_r**2))))
        epsilon = (((1-epsilon1*C_r)/(1-epsilon1))**(N-1)) / (((1-epsilon1*C_r)/(1-epsilon1))*(N-C_r))

    return epsilon


def exer(ratio, Ts_in, Th_in, material):

    thsm = tmodel(cod, peit, hotT)
    x = 0
    k = 0.0001

    for i in range(thsm.retime):
        x = x + (thsm.stout[i] + (1 / ratio) * math.log((thsm.stout[i]) * ratio + 1)) * k  # exergy integration

    y = x / (1.5 * k)  # fractional exergy recovery

    if material == 'Alumina':
        exev = y / thsm.retime
    elif material == 'Molten_salt':
        exev = 0.9 * Ts_in / Th_in

    exep = exev-0.06

    return exep


def eco(exer_value, material, block_eff, Tc_out, Tc_in):
    pass


# Test case
cp_h = 1  # KJ/kg-K
m_h = 20  # kg/s
Th_in = 550+273  # K
Th_out = 450+273  # K
cp_c = 1  # KJ/kg-K
m_c = 15  # kg/s
U = 2  # kW/m^2-K
A = 20  # m2
ratio = Th_in/293  # Ratio of T_hot/T_ambient in K
material = 'Molten_salt'

p2s_hx = P2S_HX(m_h, cp_h, Th_in, m_c, cp_c, Th_out, U, A, 'Double Shell Pass')
ep = p2s_hx.epsilon
Tc_in = p2s_hx.T_cold_in
Tc_out = p2s_hx.T_cold_out

Ts_in = Tc_out
exer_value = exer(ratio, Ts_in, Th_in, material)

print("Tc_in =  ", Tc_in, " Tc_out = ", Tc_out)
# st_cost = eco(exer_value, material, 0.4, Tc_out, Tc_in)  # Storage cost in $/KWhr
