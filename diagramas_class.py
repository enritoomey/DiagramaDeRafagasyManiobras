import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

class Diagramas_de_maniobras_y_rafagas(object):
    def __init__(self, datos, w, h, den, units='SI'):
        self.CAM = datos["CAM"]
        self.sw = datos["sw"]
        self.a3D = datos["a3D"]
        self.MTOW = datos["MTOW"]
        self.MLW = datos["MLW"]
        self.W0 = datos["W0"]
        self.MZFW = datos["MZFW"]
        self.Vc = datos["Vc"]
        self.clmax = datos["clmax"]
        self.clmax_flap = datos["clmax_flap"]
        self.clmin = datos["clmin"]
        self.Zmo = datos["Zmo"]
        self.W = w
        self.h = h
        self.den = den
        self.units = units

        self.carga_alar = {}
        self.H = {}
        self.Vs1 = {}
        self.Vs0 = {}
        self.Vsf = {}
        self.Vd = {}
        self.Va = {}
        self.Vf = {}
        self.Vf_n2 = {}
        self.Vb = {}
        self.Uref = {}
        self.Uds = self.U = {}
        self.Ude_25fts = {}
        self.Ude_50fts = {}
        self.Ude_60fts = {}

        self.vel_label = {'IM': 'ft/s', 'SI': 'm/s'}

        # constantes fijas:
        self.ft2m = 0.3048
        self.lb2kg = 0.453592
        self.slugcuft2kgm3 = 515.379

        self.cte_fgz = {'IM': 250000}
        self.cte_fgz['SI'] = self.cte_fgz['IM'] * self.ft2m
        self.s = {'IM': 100.015}
        self.s['SI'] = self.s['IM'] * self.ft2m
        self.gravedad = {'SI': 9.81}
        self.gravedad['IM'] = self.gravedad['SI'] * self.ft2m / self.lb2kg
        self.cte_nmax_1 = {'IM': 24000}
        self.cte_nmax_1['SI'] = self.cte_nmax_1['IM'] * self.lb2kg
        self.cte_nmax_2 = {'IM': 10000}
        self.cte_nmax_2['SI'] = self.cte_nmax_1['IM'] * self.lb2kg

        self.cte_Uref_h1 = {'IM': 15000}
        self.cte_Uref_h1['SI'] = self.cte_Uref_h1['IM'] * self.ft2m
        self.cte_Uref_h2 = {'IM': 50000}
        self.cte_Uref_h2['SI'] = self.cte_Uref_h2['IM'] * self.ft2m
        self.cte_Uref_v1 = {'IM': 56}
        self.cte_Uref_v1['SI'] = self.cte_Uref_v1['IM'] * self.ft2m
        self.cte_Uref_v2 = {'IM': 56}
        self.cte_Uref_v2['SI'] = self.cte_Uref_v2['IM'] * self.ft2m
        self.cte_Uref_v3 = {'IM': 26}
        self.cte_Uref_v3['SI'] = self.cte_Uref_v3['IM'] * self.ft2m

        # Esta constante esta porque hay que usar la pendiente a_cn = dCn/dalpha, y no a_cl = dCl/dalpha, pero no se de donde sale el valor
        self.ad_CN = 0.59248
        self.cte_Vb = {'IM': 498.0}  # lb/s**2
        self.cte_Vb['SI'] = self.cte_Vb['IM'] * self.ft2m ** 4 / self.lb2kg

        # Velocidad de rafadas
        self.cte_Ude_h1 = {'IM': 20000}
        self.cte_Ude_h1['SI'] = self.cte_Ude_h1['IM'] * self.ft2m
        self.cte_Ude_h2 = {'IM': 50000}
        self.cte_Ude_h2['SI'] = self.cte_Ude_h2['IM'] * self.ft2m
        self.cte_25fts_v1 = {'IM': 25}
        self.cte_25fts_v1['SI'] = self.cte_25fts_v1['IM'] * self.ft2m
        self.cte_25fts_v2 = {'IM': 33.34}
        self.cte_25fts_v2['SI'] = self.cte_25fts_v2['IM'] * self.ft2m
        self.cte_25fts_m2 = 0.000417
        self.cte_25fts_v3 = {'IM': 12.5}
        self.cte_25fts_v3['SI'] = self.cte_25fts_v3['IM'] * self.ft2m
        self.cte_50fts_v1 = {'IM': 50}
        self.cte_50fts_v1['SI'] = self.cte_50fts_v1['IM'] * self.ft2m
        self.cte_50fts_v2 = {'IM': 66.77}
        self.cte_50fts_v2['SI'] = self.cte_50fts_v2['IM'] * self.ft2m
        self.cte_50fts_m2 = 0.0008933
        self.cte_50fts_v3 = {'IM': 25}
        self.cte_50fts_v3['SI'] = self.cte_50fts_v3['IM'] * self.ft2m
        self.cte_60fts_v1 = {'IM': 60}
        self.cte_60fts_v1['SI'] = self.cte_60fts_v1['IM'] * self.ft2m
        self.cte_60fts_v2 = {'IM': 60}
        self.cte_60fts_v2['SI'] = self.cte_60fts_v2['IM'] * self.ft2m
        self.cte_60fts_m2 = {'IM': 18}
        self.cte_60fts_m2['SI'] = self.cte_60fts_m2['IM'] * self.ft2m
        self.cte_60fts_v3 = {'IM': 38}
        self.cte_60fts_v3['SI'] = self.cte_60fts_v3['IM'] * self.ft2m

        # constantes relacionadas con el diagrama de ráfagas
        self.R1 = self.MLW[units] / self.MTOW[units]
        self.R2 = self.MZFW[units] / self.MTOW[units]
        self.fgm = np.sqrt(self.R2 * np.tan(np.pi * self.R1 / 4.0))
        self.fgz = 1 - self.Zmo[units] / self.cte_fgz[units]
        self.fg = 0.5 * (self.fgz + self.fgm)

    def calculos(self):
        self.carga_alar[self.units] = self.W[self.units] / self.sw[self.units]
        self.mu_g = 2 * self.carga_alar[self.units] / (self.den[self.units] * self.CAM[self.units] * self.a3D)  # *gravedad[units])
        self.Kg = 0.88 * (self.mu_g / (5.3 + self.mu_g))
        self.Vs1[self.units] = np.sqrt(self.carga_alar[self.units] / (0.5 * self.den[self.units] * self.clmax))
        self.Vs0[self.units] = np.sqrt(-self.carga_alar[self.units] / (0.5 * self.den[self.units] * self.clmin))
        self.Vsf[self.units] = np.sqrt(self.carga_alar[self.units] / (0.5 * self.den[self.units] * self.clmax_flap))

        # Calculo de n_max
        self.n_max = 2.1 + self.cte_nmax_1[self.units] / (self.MTOW[self.units] + self.cte_nmax_2[self.units])
        if self.n_max < 2.5:
            self.n_max = 2.5
        elif self.n_max > 3.8:
            self.n_max = 3.8

        self.Va[self.units] = self.Vs1[self.units] * np.sqrt(self.n_max)
        if self.Va[self.units] > self.Vc[self.units]:
            self.Va[self.units] = self.Vc[self.units]
        self.Vd[self.units] = self.Vc[self.units] / 0.85
        self.Vf[self.units] = max(self.Vs1[self.units] * 1.6, self.Vsf[self.units] * 1.8)

        if self.h[self.units] < self.cte_Uref_h1[self.units]:
            self.Uref[self.units] = self.cte_Uref_v1[self.units] - 12.0 * self.h[self.units] / self.cte_Uref_h1[self.units]
        elif self.h[self.units] < self.cte_Uref_h2[self.units]:
            self.Uref[self.units] = self.cte_Uref_v2[self.units] - 18.0 * (self.h[self.units] - self.cte_Uref_h1[self.units]) / \
                                               (self.cte_Uref_h2[self.units] - self.cte_Uref_h1[self.units])
        else:
            self.Uref[self.units] = self.cte_Uref_v3[self.units]

        self.Vb[self.units] = min(self.Vc[self.units], self.Vs1[self.units] * np.sqrt(1 + self.Kg * self.Uref[self.units] * self.Vc[self.units] *
                                            self.a3D * self.ad_CN / (self.cte_Vb[self.units] * self.carga_alar[self.units])))

        if self.h[self.units] < self.cte_Ude_h1[self.units]:
            self.Ude_25fts[self.units] = self.cte_25fts_v1[self.units]
            self.Ude_50fts[self.units] = self.cte_50fts_v1[self.units]
            self.Ude_60fts[self.units] = self.cte_60fts_v1[self.units]
        elif self.h[self.units] < self.cte_Ude_h2[self.units]:
            self.Ude_25fts[self.units] = self.cte_25fts_v2[self.units] - self.cte_25fts_m2 * self.h[self.units]
            self.Ude_50fts[self.units] = self.cte_50fts_v2[self.units] - self.cte_50fts_m2 * self.h[self.units]
            self.Ude_60fts[self.units] = self.cte_60fts_v2[self.units] - self.cte_60fts_m2[self.units] * \
                                                               (self.h[self.units] - self.cte_Ude_h1[self.units]) \
                                                               /(self.cte_Ude_h2[self.units] - self.cte_Ude_h2[self.units])
        else:
            self.Ude_25fts[self.units] = self.cte_25fts_v3[self.units]
            self.Ude_50fts[self.units] = self.cte_50fts_v3[self.units]
            self.Ude_60fts[self.units] = self.cte_60fts_v3[self.units]

        self.Vf_n2[self.units] = np.sqrt(2 * self.W[self.units] / (0.5 * self.den[self.units] * self.clmax_flap * self.sw[self.units]))

    def n_25fts(self, vel):
        return self.fg * self.Ude_25fts[self.units] * self.a3D * self.ad_CN * vel / (self.cte_Vb[self.units] * self.carga_alar[self.units])

    def n_50fts(self, vel):
        return self.Kg * self.Ude_50fts[self.units] * self.a3D * self.ad_CN * vel / (self.cte_Vb[self.units] * self.carga_alar[self.units])

    def n_60fts(self, vel):
        return self.Kg * self.Ude_60fts[self.units] * self.a3D * self.ad_CN * vel / (self.cte_Vb[self.units] * self.carga_alar[self.units])

    def n_gust_pos(self, vel):
        if 0 <= vel <= self.Vb[self.units]:
            return 1 + self.n_60fts(vel)
        elif vel <= self.Vc[self.units]:
            m = (self.n_50fts(self.Vc[self.units]) - self.n_60fts(self.Vb[self.units])) / (self.Vc[self.units] - self.Vb[self.units])
            b = self.n_50fts(self.Vc[self.units]) - m * self.Vc[self.units]
            return 1 + m * vel + b
        elif vel <= self.Vd[self.units]:
            m = (self.n_25fts(self.Vd[self.units]) - self.n_50fts(self.Vc[self.units])) / (self.Vd[self.units] - self.Vc[self.units])
            b = self.n_25fts(self.Vd[self.units]) - m * self.Vd[self.units]
            return 1 + m * vel + b
        return None

    def n_gust_neg(self, vel):
        if 0 <= vel <= self.Vb[self.units]:
            return 1 - self.n_60fts(vel)
        elif vel <= self.Vc[self.units]:
            m = (self.n_50fts(self.Vc[self.units]) - self.n_60fts(self.Vb[self.units])) / (self.Vc[self.units] - self.Vb[self.units])
            b = self.n_50fts(self.Vc[self.units]) - m * self.Vc[self.units]
            return 1 - m * vel + b
        elif vel <= self.Vd[self.units]:
            m = (self.n_25fts(self.Vd[self.units]) - self.n_50fts(self.Vc[self.units])) / (self.Vd[self.units] - self.Vc[self.units])
            b = self.n_25fts(self.Vd[self.units]) - m * self.Vd[self.units]
            return 1 - m * vel + b
        return None

    def n_stall_pos(self, vel):
        return 0.5 * self.den[self.units] * vel**2 * self.sw[self.units] * self.clmax / self.W[self.units]


    def n_stall_neg(self, vel):
        return 0.5 * self.den[self.units] * vel**2 * self.sw[self.units] * self.clmin / self.W[self.units]


    def n_stall_flap(self, vel):
        return 0.5 * self.den[self.units] * vel**2 * self.sw[self.units] * self.clmax_flap / self.W[self.units]

    def n_manoeuvre_pos(self, vel):
        if 0 <= vel <= self.Va[self.units]:
            return self.n_stall_pos(vel)
        elif vel <= self.Vd[self.units]:
            return self.n_max

        return None

    def n_manoeuvre_neg(self, vel):
        if 0 <= vel <= self.Vs0[self.units]:
            return self.n_stall_neg(vel)
        elif vel <= self.Vc[self.units]:
            return -1.0
        elif vel <= self.Vd[self.units]:
            return -1 + 1 / (self.Vd[self.units] - self.Vc[self.units]) * (vel - self.Vc[self.units])
        return None

    def plot_diagrama_de_rafagas(self, ax, dv):
        ax.plot(np.arange(0, self.Vb[self.units], dv), [1 + self.n_60fts(vel) for vel in np.arange(0, self.Vb[self.units], dv)], color='r')
        ax.plot(np.arange(0, self.Vb[self.units], dv), [1 - self.n_60fts(vel) for vel in np.arange(0, self.Vb[self.units], dv)], color='r')

        ax.plot(np.arange(0, self.Vc[self.units], dv), [1 + self.n_50fts(vel) for vel in np.arange(0, self.Vc[self.units], dv)], color='b')
        ax.plot(np.arange(0, self.Vc[self.units], dv), [1 - self.n_50fts(vel) for vel in np.arange(0, self.Vc[self.units], dv)], color='b')

        ax.plot(np.arange(0, self.Vd[self.units], dv), [1 + self.n_25fts(vel) for vel in np.arange(0, self.Vd[self.units], dv)], color='g')
        ax.plot(np.arange(0, self.Vd[self.units], dv), [1 - self.n_25fts(vel) for vel in np.arange(0, self.Vd[self.units], dv)], color='g')

        ax.plot([self.Vb[self.units], self.Vc[self.units]], [1 + self.n_60fts(self.Vb[self.units]), 1 + self.n_50fts(self.Vc[self.units])], color='m')
        ax.plot([self.Vb[self.units], self.Vc[self.units]], [1 - self.n_60fts(self.Vb[self.units]), 1 - self.n_50fts(self.Vc[self.units])], color='m')

        ax.plot([self.Vc[self.units], self.Vd[self.units]], [1 + self.n_50fts(self.Vc[self.units]), 1 + self.n_25fts(self.Vd[self.units])], color='m')
        ax.plot([self.Vc[self.units], self.Vd[self.units]], [1 - self.n_50fts(self.Vc[self.units]), 1 - self.n_25fts(self.Vd[self.units])], color='m')

        ax.plot([self.Vd[self.units], self.Vd[self.units]], [1 + self.n_25fts(self.Vd[self.units]), 1 - self.n_25fts(self.Vd[self.units])], color='m')
        ax.set_xlabel("Speed [{}]".format(self.vel_label[self.units]))
        ax.set_ylabel("n")
        ax.set_title("Gust Diagram")

    def plot_diagrama_de_maniobras(self, ax, dv):
        ax.plot(np.arange(0, self.Vs1[self.units], dv), [self.n_stall_pos(vel) for vel in np.arange(0, self.Vs1[self.units], dv)], color='m',
                linestyle='--')
        ax.plot([self.Vs1[self.units], self.Vs1[self.units]], [0, self.n_stall_pos(self.Vs1[self.units])], color='m')
        ax.plot(np.arange(self.Vs1[self.units], self.Va[self.units], dv),
                [self.n_stall_pos(vel) for vel in np.arange(self.Vs1[self.units], self.Va[self.units], dv)],
                color='m', linestyle='-')
        ax.plot(np.arange(0, self.Vs0[self.units] + dv, dv), [self.n_stall_neg(vel) for vel in np.arange(0, self.Vs0[self.units] + dv, dv)],
                color='m', linestyle='--')
        ax.plot([self.Vs0[self.units], self.Vs0[self.units]], [0, -1.0], color='m')
        ax.plot([self.Vs1[self.units], self.Vs0[self.units]], [0.0, 0.0], color='m')
        ax.plot([self.Va[self.units], self.Vd[self.units]], [self.n_max, self.n_max], color='m')
        ax.plot([self.Vd[self.units], self.Vd[self.units]], [self.n_max, 0], color='m')
        ax.plot([self.Vs0[self.units], self.Vc[self.units]], [-1.0, -1.0], color='m')
        ax.plot([self.Vc[self.units], self.Vd[self.units]], [-1.0, 0.0], color='m')
        ax.set_xlabel("Speed [{}]".format(self.vel_label[self.units]))
        ax.set_ylabel("n")
        ax.set_title("Manoeuvre Diagram")

    def plot_diagrama_de_maniobras_con_flap(self, ax, dv):
        ax.plot(np.arange(0, self.Vsf[self.units] + dv, dv), [self.n_stall_flap(vel) for vel in np.arange(0, self.Vsf[self.units] + dv, dv)],
                color='b', linestyle='--')
        ax.plot(np.arange(self.Vsf[self.units], self.Vf_n2[self.units] + dv, dv),
                [self.n_stall_flap(vel) for vel in np.arange(self.Vsf[self.units], self.Vf_n2[self.units] + dv, dv)],
                color='b', linestyle='-')
        ax.plot([self.Vsf[self.units], self.Vsf[self.units]], [0.0, self.n_stall_flap(self.Vsf[self.units])], color='b', linestyle='-')
        ax.plot([self.Vf_n2[self.units], self.Vf[self.units]], [2.0, 2.0], color='b', linestyle='-')
        ax.plot([self.Vf[self.units], self.Vf[self.units]], [0.0, 2.0], color='b', linestyle='-')
        ax.plot([self.Vsf[self.units], self.Vf[self.units]], [0.0, 0.0], color='b', linestyle='-')
        ax.set_xlabel("Speed [{}]".format(self.vel_label[self.units]))
        ax.set_ylabel("n")
        ax.set_title("Manoeuvre Diagram")

    def plot_diagrama_de_maniobras_y_rafagas(self, ax, dv):
        # Calculo de las intersecciones:
        if self.n_gust_pos(self.Va[self.units]) > self.n_max:
            # extender stall hasta interseccion con gust y arrancar la comparacion desde ese punto
            def func1(vel):
                return self.n_gust_pos(vel) - self.n_stall_pos(vel)

            v_intersec_pos = fsolve(func1, self.Va[self.units])[0]
        else:
            v_intersec_pos = self.Va[self.units]

        if self.n_gust_pos(self.Vs0[self.units]) < -1.0:
            # extender stall hasta interseccion con gust y arrancar la comparacion desde ese punto
            def func2(vel):
                return self.n_gust_neg(vel) - self.n_stall_neg(vel)

            v_intersec_neg = fsolve(func2, self.Vs0[self.units])[0]
        else:
            v_intersec_neg = self.Vs0[self.units]

        ax.fill_between(np.arange(0, v_intersec_pos, dv), 0,
                        [self.n_stall_pos(vel) for vel in np.arange(0, v_intersec_pos, dv)],
                        color='m', alpha=0.2)
        ax.fill_between(np.arange(v_intersec_pos, self.Vd[self.units], dv), 0, [max(self.n_gust_pos(vel), self.n_manoeuvre_pos(vel))
                                                                      for vel in
                                                                      np.arange(v_intersec_pos, self.Vd[self.units], dv)],
                        color='m', alpha=0.2)
        ax.fill_between(np.arange(0, v_intersec_neg, dv), 0,
                        [self.n_stall_neg(vel) for vel in np.arange(0, v_intersec_neg, dv)],
                        color='m', alpha=0.2)
        ax.fill_between(np.arange(v_intersec_neg, self.Vd[self.units], dv), 0, [min(self.n_gust_neg(vel), self.n_manoeuvre_neg(vel))
                                                                      for vel in
                                                                      np.arange(v_intersec_neg, self.Vd[self.units], dv)],
                        color='m', alpha=0.2)
        ax.fill_between([self.Vd[self.units], self.Vd[self.units]], 0, [max(self.n_manoeuvre_pos(self.Vd[self.units]), self.n_gust_pos(self.Vd[self.units])),
                                                    min(self.n_manoeuvre_neg(self.Vd[self.units]), self.n_gust_neg(self.Vd[self.units]))], color='m',
                        alpha=0.2)
        ax.set_xlabel("Speed [{}]".format(self.vel_label[self.units]))
        ax.set_ylabel("n")
        ax.set_title("Combined Gust & Manoeuvre Diagram")


if __name__ == "__main__":
    ft2m = 0.3048
    lb2kg = 0.453592
    slugcuft2kgm3 = 515.379

    # Input Data:
    CAM = {'SI': 2.461}
    CAM['IM'] = CAM['SI'] / ft2m
    sw = {'SI': 60}
    sw['IM'] = sw['SI'] / ft2m / ft2m
    a3D = 5.0037  # 1/rad
    MTOW = {'SI': 23000}
    MTOW['IM'] = MTOW['SI'] / lb2kg
    MLW = {'SI': 23000}
    MLW['IM'] = MLW['SI'] / lb2kg
    W0 = {'SI': 13766.0}
    W0['IM'] = W0['SI'] / lb2kg
    MZFW = {'SI': 16376.0}
    MZFW['IM'] = MZFW['SI'] / lb2kg
    Vc = {'SI': 151.93}
    Vc['IM'] = Vc['SI'] / ft2m
    clmax = 1.2463
    clmax_flap = 1.499
    clmin = -0.75 * clmax
    Zmo = {'SI': 9999.2}
    Zmo['IM'] = Zmo['SI'] / ft2m

    # Variables
    W = {'SI': 20000}
    W['IM'] = W['SI'] / lb2kg
    h = {'SI': 5000}
    h['IM'] = h['SI'] / ft2m
    den = {'SI': 0.125}
    den['IM'] = den['SI'] / lb2kg * ft2m ** 3

    datos = {
        'CAM': CAM,
        'sw': sw,
        'a3D': a3D,
        'MTOW': MTOW,
        'MLW': MLW,
        'W0': W0,
        'MZFW': MZFW,
        'Vc': Vc,
        'clmax': clmax,
        'clmax_flap': clmax,
        'clmin': clmin,
        'Zmo': Zmo
    }

    diagrama = Diagramas_de_maniobras_y_rafagas(datos, W, h, den, units='SI')
    diagrama.calculos()
    
    fig, ax1 = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, squeeze=True)
    diagrama.plot_diagrama_de_rafagas(ax1, 0.5)

    fig, ax2 = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, squeeze=True)
    diagrama.plot_diagrama_de_maniobras(ax2, 0.5)

    fig, ax3 = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, squeeze=True)
    diagrama.plot_diagrama_de_maniobras_con_flap(ax3, 0.5)

    fig, ax4 = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, squeeze=True)
    diagrama.plot_diagrama_de_maniobras_y_rafagas(ax4, 0.5)

    plt.grid(True)
    plt.show()
    