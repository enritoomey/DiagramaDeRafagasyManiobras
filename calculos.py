# TODO: turn this into a set of functions

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def plot_diagrama_de_rafagas(ax, Vb, Vc, Vd, n_25fts, n_50fts, n_60fts, dv, units, vel_label):
    ax.plot(np.arange(0, Vb[units], dv), [1 + n_60fts(vel) for vel in np.arange(0, Vb[units], dv)], color='r')
    ax.plot(np.arange(0, Vb[units], dv), [1 - n_60fts(vel) for vel in np.arange(0, Vb[units], dv)], color='r')

    ax.plot(np.arange(0, Vc[units], dv), [1 + n_50fts(vel) for vel in np.arange(0, Vc[units], dv)], color='b')
    ax.plot(np.arange(0, Vc[units], dv), [1 - n_50fts(vel) for vel in np.arange(0, Vc[units], dv)], color='b')

    ax.plot(np.arange(0, Vd[units], dv), [1 + n_25fts(vel) for vel in np.arange(0, Vd[units], dv)], color='g')
    ax.plot(np.arange(0, Vd[units], dv), [1 - n_25fts(vel) for vel in np.arange(0, Vd[units], dv)], color='g')

    ax.plot([Vb[units], Vc[units]], [1 + n_60fts(Vb[units]), 1 + n_50fts(Vc[units])], color='m')
    ax.plot([Vb[units], Vc[units]], [1 - n_60fts(Vb[units]), 1 - n_50fts(Vc[units])], color='m')

    ax.plot([Vc[units], Vd[units]], [1 + n_50fts(Vc[units]), 1 + n_25fts(Vd[units])], color='m')
    ax.plot([Vc[units], Vd[units]], [1 - n_50fts(Vc[units]), 1 - n_25fts(Vd[units])], color='m')

    ax.plot([Vd[units], Vd[units]], [1 + n_25fts(Vd[units]), 1 - n_25fts(Vd[units])], color='m')
    ax.set_xlabel("Speed [{}]".format(vel_label[units]))
    ax.set_ylabel("n")
    ax.set_title("Gust Diagram")


def plot_diagrama_de_maniobras(ax, n_stall_pos, n_stall_neg, n_max, Vs1, Vs0, Va, dv):
    ax.plot(np.arange(0, Vs1[units], dv), [n_stall_pos(vel) for vel in np.arange(0, Vs1[units], dv)], color='m',
            linestyle='--')
    ax.plot([Vs1[units], Vs1[units]], [0, n_stall_pos(Vs1[units])], color='m')
    ax.plot(np.arange(Vs1[units], Va[units], dv), [n_stall_pos(vel) for vel in np.arange(Vs1[units], Va[units], dv)],
            color='m', linestyle='-')
    ax.plot(np.arange(0, Vs0[units] + dv, dv), [n_stall_neg(vel) for vel in np.arange(0, Vs0[units] + dv, dv)],
            color='m', linestyle='--')
    ax.plot([Vs0[units], Vs0[units]], [0, -1.0], color='m')
    ax.plot([Vs1[units], Vs0[units]], [0.0, 0.0], color='m')
    ax.plot([Va[units], Vd[units]], [n_max, n_max], color='m')
    ax.plot([Vd[units], Vd[units]], [n_max, 0], color='m')
    ax.plot([Vs0[units], Vc[units]], [-1.0, -1.0], color='m')
    ax.plot([Vc[units], Vd[units]], [-1.0, 0.0], color='m')


def plot_diagrama_de_maniobras_con_flap(ax, n_stall_flap, Vsf, Vf_n2, Vf, dv, units, vel_label):
    ax.plot(np.arange(0, Vsf[units] + dv, dv), [n_stall_flap(vel) for vel in np.arange(0, Vsf[units] + dv, dv)],
            color='b', linestyle='--')
    ax.plot(np.arange(Vsf[units], Vf_n2 + dv, dv), [n_stall_flap(vel) for vel in np.arange(Vsf[units], Vf_n2 + dv, dv)],
            color='b', linestyle='-')
    ax.plot([Vsf[units], Vsf[units]], [0.0, n_stall_flap(Vsf[units])], color='b', linestyle='-')
    ax.plot([Vf_n2, Vf[units]], [2.0, 2.0], color='b', linestyle='-')
    ax.plot([Vf[units], Vf[units]], [0.0, 2.0], color='b', linestyle='-')
    ax.plot([Vsf[units], Vf[units]], [0.0, 0.0], color='b', linestyle='-')
    ax.set_xlabel("Speed [{}]".format(vel_label[units]))
    ax.set_ylabel("n")
    ax.set_title("Manoeuvre Diagram")

def plot_diagrama_de_maniobras_y_rafagas(ax, n_stall_pos, n_stall_neg, n_gust_pos, n_gust_neg, n_manoeuvre_pos,
                                         n_manoeuvre_neg, v_intersec_pos, v_intersec_neg, Vd, dv, units, vel_label):
    ax.fill_between(np.arange(0, v_intersec_pos, dv), 0, [n_stall_pos(vel) for vel in np.arange(0, v_intersec_pos, dv)],
                    color='m', alpha=0.2)
    ax.fill_between(np.arange(v_intersec_pos, Vd[units], dv), 0, [max(n_gust_pos(vel), n_manoeuvre_pos(vel))
                                                                  for vel in np.arange(v_intersec_pos, Vd[units], dv)],
                    color='m', alpha=0.2)
    ax.fill_between(np.arange(0, v_intersec_neg, dv), 0, [n_stall_neg(vel) for vel in np.arange(0, v_intersec_neg, dv)],
                    color='m', alpha=0.2)
    ax.fill_between(np.arange(v_intersec_neg, Vd[units], dv), 0, [min(n_gust_neg(vel), n_manoeuvre_neg(vel))
                                                                  for vel in np.arange(v_intersec_neg, Vd[units], dv)],
                    color='m', alpha=0.2)
    ax.fill_between([Vd[units], Vd[units]], 0, [max(n_manoeuvre_pos(Vd[units]), n_gust_pos(Vd[units])),
                                                min(n_manoeuvre_neg(Vd[units]), n_gust_neg(Vd[units]))], color='m',
                    alpha=0.2)
    ax.set_xlabel("Speed [{}]".format(vel_label[units]))
    ax.set_ylabel("n")
    ax.set_title("Combined Gust & Manoeuvre Diagram")

#import ipdb
if __name__ == '__main__':
    # Units
    units = 'IM' # 'IM'
    ft2m = 0.3048
    lb2kg = 0.453592
    slugcuft2kgm3 = 515.379
    vel_label = {'IM': 'ft/s', 'SI': 'm/s'}

    # Input Data:
    CAM = {'SI': 2.461}
    CAM['IM'] = CAM['SI']/ft2m
    sw = {'SI': 60}
    sw['IM'] = sw['SI']/ft2m/ft2m
    a3D = 5.0037 #1/rad
    MTOW = {'SI': 23000}
    MTOW['IM'] = MTOW['SI']/lb2kg
    MLW = {'SI': 23000}
    MLW['IM'] = MLW['SI']/lb2kg
    W0 = {'SI': 13766.0}
    W0['IM'] = W0['SI']/lb2kg
    MZFW = {'SI': 16376.0}
    MZFW['IM'] = MZFW['SI']/lb2kg
    Vc = {'SI': 151.93}
    Vc['IM'] = Vc['SI']/ft2m
    clmax = 1.2463
    clmax_flap = 1.499
    clmin = -0.75*clmax
    Zmo = {'SI': 9999.2}
    Zmo['IM'] = Zmo['SI']/ft2m


    # Variables
    W = {'SI': 20000}
    W['IM'] = W['SI']/lb2kg

    h = {'SI': 5000}
    h['IM'] = h['SI']/ft2m

    den = {'SI': 0.125}
    den['IM'] = den['SI']/lb2kg*ft2m**3
    print(den)
    # constantes
    cte_fgz = {'IM': 250000}
    cte_fgz['SI'] = cte_fgz['IM']*ft2m
    s = {'IM': 100.015}
    s['SI'] = s['IM']*ft2m
    gravedad = {'SI': 9.81}
    gravedad['IM'] = gravedad['SI']*ft2m/lb2kg
    cte_nmax_1 = {'IM': 24000}
    cte_nmax_1['SI'] = cte_nmax_1['IM']*lb2kg
    cte_nmax_2 = {'IM': 10000}
    cte_nmax_2['SI'] = cte_nmax_1['IM']*lb2kg


    # Constants depending from input data
    carga_alar = {}
    H = {}
    Vs1 = {}
    Vs0 = {}
    Vsf = {}
    Vd = {}
    Va = {}
    Vf = {}
    Vb = {}
    Uref = {}
    Uds = U = {}
    Ude_25fts = {}
    Ude_50fts = {}
    Ude_60fts = {}
    carga_alar[units] = W[units]/sw[units]
    mu_g = 2*carga_alar[units]/(den[units]*CAM[units]*a3D)#*gravedad[units])
    print(mu_g)
    Kg = 0.88*(mu_g/(5.3+mu_g))
    Vs1[units] = np.sqrt(carga_alar[units]/(0.5*den[units]*clmax))
    Vs0[units] = np.sqrt(-carga_alar[units]/(0.5*den[units]*clmin))
    Vsf[units] = np.sqrt(carga_alar[units]/(0.5*den[units]*clmax_flap))

    # Calculo de n_max
    n_max = 2.1+cte_nmax_1[units]/(MTOW[units]+cte_nmax_2[units])
    if n_max < 2.5:
        n_max = 2.5
    elif n_max > 3.8:
        n_max = 3.8

    Va[units] = Vs1[units]*np.sqrt(n_max)
    if Va[units] > Vc[units]:
        Va[units] = Vc[units]
    Vd[units] = Vc[units]/0.85
    Vf[units] = max(Vs1[units]*1.6, Vsf[units]*1.8)

    cte_Uref_h1 = {'IM': 15000}
    cte_Uref_h1['SI'] = cte_Uref_h1['IM']*ft2m
    cte_Uref_h2 = {'IM': 50000}
    cte_Uref_h2['SI'] = cte_Uref_h2['IM']*ft2m
    cte_Uref_v1 = {'IM': 56}
    cte_Uref_v1['SI'] = cte_Uref_v1['IM']*ft2m
    cte_Uref_v2 = {'IM': 56}
    cte_Uref_v2['SI'] = cte_Uref_v2['IM']*ft2m
    cte_Uref_v3 = {'IM': 26}
    cte_Uref_v3['SI'] = cte_Uref_v3['IM']*ft2m

    #ipdb.set_trace()
    if h[units] < cte_Uref_h1[units]:
        Uref[units] = cte_Uref_v1[units]-12.0*h[units]/cte_Uref_h1[units]
    elif h[units] < cte_Uref_h2[units]:
        Uref[units] = cte_Uref_v2[units]-18.0*(h[units]-cte_Uref_h1[units])/(cte_Uref_h2[units]-cte_Uref_h1[units])
    else:
        Uref[units] = cte_Uref_v3[units]

    # Esta constante esta porque hay que usar la pendiente a_cn = dCn/dalpha, y no a_cl = dCl/dalpha, pero no se de donde sale el valor
    ad_CN = 0.59248
    cte_Vb = {'IM':498.0} # lb/s**2
    cte_Vb['SI'] = cte_Vb['IM']*ft2m**4/lb2kg
    Vb[units] = min(Vc[units], Vs1[units]*np.sqrt(1+Kg*Uref[units]*Vc[units]*a3D*ad_CN/(cte_Vb[units]*carga_alar[units])))

    # Velocidad de rafadas
    cte_Ude_h1 = {'IM': 20000}
    cte_Ude_h1['SI'] = cte_Ude_h1['IM']*ft2m
    cte_Ude_h2 = {'IM': 50000}
    cte_Ude_h2['SI'] = cte_Ude_h2['IM']*ft2m
    cte_25fts_v1 = {'IM': 25}
    cte_25fts_v1['SI'] = cte_25fts_v1['IM']*ft2m
    cte_25fts_v2 = {'IM': 33.34}
    cte_25fts_v2['SI'] = cte_25fts_v2['IM']*ft2m
    cte_25fts_m2 = 0.000417
    cte_25fts_v3 = {'IM': 12.5}
    cte_25fts_v3['SI'] = cte_25fts_v3['IM']*ft2m
    cte_50fts_v1 = {'IM': 50}
    cte_50fts_v1['SI'] = cte_50fts_v1['IM']*ft2m
    cte_50fts_v2 = {'IM': 66.77}
    cte_50fts_v2['SI'] = cte_50fts_v2['IM']*ft2m
    cte_50fts_m2 = 0.0008933
    cte_50fts_v3 = {'IM': 25}
    cte_50fts_v3['SI'] = cte_50fts_v3['IM']*ft2m
    cte_60fts_v1 = {'IM': 60}
    cte_60fts_v1['SI'] = cte_60fts_v1['IM']*ft2m
    cte_60fts_v2 = {'IM': 60}
    cte_60fts_v2['SI'] = cte_60fts_v2['IM']*ft2m
    cte_60fts_m2 = {'IM': 18}
    cte_60fts_m2['SI'] = cte_60fts_m2['IM']*ft2m
    cte_60fts_v3 = {'IM': 38}
    cte_60fts_v3['SI'] = cte_60fts_v3['IM']*ft2m


    if h[units] < cte_Ude_h1[units]:
        Ude_25fts[units] = cte_25fts_v1[units]
        Ude_50fts[units] = cte_50fts_v1[units]
        Ude_60fts[units] = cte_60fts_v1[units]
    elif h[units] < cte_Ude_h2[units]:
        Ude_25fts[units] = cte_25fts_v2[units] - cte_25fts_m2 * h[units]
        Ude_50fts[units] = cte_50fts_v2[units] - cte_50fts_m2 * h[units]
        Ude_60fts[units] = cte_60fts_v2[units] - cte_60fts_m2[units]*(h[units]-cte_Ude_h1[units])/(cte_Ude_h2[units]-cte_Ude_h2[units])
    else:
        Ude_25fts[units] = cte_25fts_v3[units]
        Ude_50fts[units] = cte_50fts_v3[units]
        Ude_60fts[units] = cte_60fts_v3[units]

    def n_25fts(vel):
        return fg*Ude_25fts[units]*a3D*ad_CN*vel/(cte_Vb[units]*carga_alar[units])

    def n_50fts(vel):
        return Kg*Ude_50fts[units]*a3D*ad_CN*vel/(cte_Vb[units]*carga_alar[units])

    def n_60fts(vel):
        return Kg*Ude_60fts[units]*a3D*ad_CN*vel/(cte_Vb[units]*carga_alar[units])

    def n_gust_pos(vel):
        if 0 <= vel <= Vb[units]:
            return 1+n_60fts(vel)
        elif vel <= Vc[units]:
            m = (n_50fts(Vc[units])-n_60fts(Vb[units]))/(Vc[units]-Vb[units])
            b = n_50fts(Vc[units])-m*Vc[units]
            return 1+m*vel+b
        elif vel <= Vd[units]:
            m = (n_25fts(Vd[units])-n_50fts(Vc[units]))/(Vd[units]-Vc[units])
            b = n_25fts(Vd[units])-m*Vd[units]
            return 1+m*vel+b
        return None

    def n_gust_neg(vel):
        if 0 <= vel <= Vb[units]:
            return 1-n_60fts(vel)
        elif vel <= Vc[units]:
            m = (n_50fts(Vc[units])-n_60fts(Vb[units]))/(Vc[units]-Vb[units])
            b = n_50fts(Vc[units])-m*Vc[units]
            return 1-m*vel+b
        elif vel <= Vd[units]:
            m = (n_25fts(Vd[units])-n_50fts(Vc[units]))/(Vd[units]-Vc[units])
            b = n_25fts(Vd[units])-m*Vd[units]
            return 1-m*vel+b
        return None

    # Variables definidas pero no utilizadas
    # H[units] = 12.5*CAM[units]
    R1 = MLW[units]/ MTOW[units]
    R2 = MZFW[units]/MTOW[units]
    fgm =  np.sqrt(R2*np.tan(np.pi*R1/4.0))
    fgz = 1 - Zmo[units]/cte_fgz[units]
    fg = 0.5*(fgz+fgm)
    # cte_Uds = {'IM':350}
    # cte_Uds['SI'] = cte_Uds['IM']*ft2m
    # Uds[units] = Uref[units]*fg*(H[units]/cte_Uds[units])**(1.0/6.0)
    # U[units] = 0.5*Uds[units]*(1-np.cos(np.pi*s[units]/H[units]))
    print("Kg = {}, Vb = {}, Vc = {}, Vd = {}".format(Kg,Vb, Vc, Vd))

    dv = 1.0

    # Diagrama de Rafagas
    fig, ax = plt.subplots(nrows=1, ncols= 1, sharex=True, sharey=True, squeeze=True)
    plot_diagrama_de_rafagas(ax, Vb, Vc, Vd, n_25fts, n_50fts, n_60fts, dv, units, vel_label)
    plt.grid(True)

    plt.show()

    def n_stall_pos(vel):
        return 0.5*den[units]*vel**2*sw[units]*clmax/W[units]


    def n_stall_neg(vel):
        return 0.5*den[units]*vel**2*sw[units]*clmin/W[units]


    def n_stall_flap(vel):
        return 0.5*den[units]*vel**2*sw[units]*clmax_flap/W[units]


    def n_manoeuvre_pos(vel):
        if 0 <= vel <= Va[units]:
            return n_stall_pos(vel)
        elif vel <= Vd[units]:
            return n_max

        return None

    def n_manoeuvre_neg(vel):
        if 0 <= vel <= Vs0[units]:
            return n_stall_neg(vel)
        elif vel <= Vc[units]:
            return -1.0
        elif vel <= Vd[units]:
            return -1+1/(Vd[units]-Vc[units])*(vel-Vc[units])

        return None

    # Diagrama de maniobras:
    fig, ax = plt.subplots(nrows=1, ncols= 1, sharex=True, sharey=True, squeeze=True)
    plot_diagrama_de_maniobras(ax, n_stall_pos, n_stall_neg, n_max, Vs1, Vs0, Va, dv)
    #plt.show()

    # Diagrama de maniobras c/flap:
    #fig, ax = plt.subplots(nrows=1, ncols= 1, sharex=True, sharey=True, squeeze=True)
    Vf_n2 = np.sqrt(2*W[units]/(0.5*den[units]*clmax_flap*sw[units]))
    plot_diagrama_de_maniobras_con_flap(ax, n_stall_flap, Vsf, Vf_n2, Vf, dv, units, vel_label)
    plt.grid(True)
    plt.show()

    # Calculo de las intersecciones:
    if n_gust_pos(Va[units]) > n_max:
        # extender stall hasta interseccion con gust y arrancar la comparacion desde ese punto
        def func1(vel):
            return n_gust_pos(vel) - n_stall_pos(vel)
        v_intersec_pos = fsolve(func1, Va[units])[0]
    else:
        v_intersec_pos = Va[units]

    if n_gust_pos(Vs0[units]) < -1.0:
        # extender stall hasta interseccion con gust y arrancar la comparacion desde ese punto
        def func2(vel):
            return n_gust_neg(vel) - n_stall_neg(vel)
        v_intersec_neg = fsolve(func2, Vs0[units])[0]
    else:
        v_intersec_neg = Vs0[units]

    # Plot intersection
    # fig = plt.figure(facecolor='white')
    # axescolor = '#f6f6f6'  # the axes background color
    # ax = fig.add_axes([0, 1, 0, 1], axisbg=axescolor)

    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, squeeze=True)
    plot_diagrama_de_maniobras_y_rafagas(ax, n_stall_pos, n_stall_neg, n_gust_pos, n_gust_neg, n_manoeuvre_pos,
                                             n_manoeuvre_neg, v_intersec_pos, v_intersec_neg, Vd, dv, units, vel_label)
    plt.grid(True)
    plt.show()

