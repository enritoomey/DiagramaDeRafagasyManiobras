# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 00:12:28 2015

@author: Enriquito
"""
import numpy as np
import argparse

def atmosfera_estandar(calculo,atmosfera):

    b_suth = 1.458*10**(-6) #[Kg/(m*s*sqrt(K))]
    s_suth = 110.4 #[K]
    g  =9.81 #m/s**2
    R = 287.0 #J/(kg*K)
    k = 1.4
    Z = [0.0, 11000.0, 20100.0, 32200.0, 52400.0, 61600.0, 80000.0, 95000.0]#[m]
    T = [288.0, 216.5, 215.6, 228.5, 270.5, 252.5, 180.5, 180.5]#[K]
    P = [101325.0, 22628.3619, 5386.81, 840.76, 52.62, 15.822, 0.8455, 0.04951]#Pa

    if calculo == 'altura':
        h = atmosfera['h']
        deltaT = atmosfera['deltaT']
        n = 1
        c = 0
        while c == 0 | n < 7:
            if h < Z[n]:
                c = 1
            else:
                n += 1
        
        if T[n] == T[n-1]:
            t = T[n]+deltaT
            p = P[n-1]*np.e**(g*(h-Z[n-1])/R/T[n])
            rho = p/t/R
        else:
            m = (T[n]-T[n-1])/(Z[n]-Z[n-1])
            temp = m*(h-Z[n-1])+T[n-1]
            t = temp+deltaT
            p = P[n-1]*(temp/T[n-1])**(-g/m/R)
            rho = p/t/R
    
    elif calculo == 'presion':
        p = atmosfera['p']
        deltaT = atmosfera['deltaT']
        n = 1
        c = 0
        while c == 0 | n < 7:
            if p > P[n]:
                c = 1
            else:
                n += 1
            
        if T[n] == T[n-1]:
            t = T[n]+deltaT
            h = R*T[n]/g*np.log(p/P[n-1])+Z[n-1]
            rho = p/t/R
        else:
            m = (T[n]-T[n-1])/(Z[n]-Z[n-1])
            h = ((p/P[n-1])**(-m*R/g)*T[n-1]-T[n-1])/m+Z[n-1]
            t = m*(h-Z[n-1])+T[n-1]+deltaT
            rho = p/t/R                        

    else:
        raise("Option {} yet not implementes".format(calculo))

    mu = b_suth*np.sqrt(t)/(1+s_suth/t)
    vson = np.sqrt(k*R*t)

    atmosfera['h'] = h
    atmosfera['deltaT'] = deltaT
    atmosfera['p'] = p
    atmosfera['t'] = t
    atmosfera['rho'] = rho
    atmosfera['mu']= mu
    atmosfera['Vson'] = vson
    
    return atmosfera


def parse():
    parser = argparse.ArgumentParser(description="Funcion para calcular datos de la atmofera"
                                                 " estandard a partir de altura, presion, temperatura"
                                                 "o densidad, mas un deltaT en grados centigrados")
    parser.add_argument("input_type", type=str, default="altura",
                        choices=('altura', 'presion', 'temperatura', 'densidad'),
                        help="Seleccionar el tipo de dato de entrada")
    parser.add_argument("input1", type=float, help="ingrese el primer dato de entrada en unidades SI")
    parser.add_argument("input2", type=float, help="ingrese el deltaT en C")
    return parser.parse_args()

if __name__ == '__main__':
    args = parse()
    atmosfera = {
        'h': 0,
        'deltaT': 0,
        'p': 0,
        't': 0,
        'rho': 0,
        'mu': 0,
        'vson': 0
    }
    if args.input_type == 'altura':
        atmosfera['h'] = args.input1
    elif args.input_type == 'presion':
        atmosfera['p'] = args.input1
    elif args.input_type == 'temperatura':
        atmosfera['t'] = args.input1
    elif args.input_type == 'densidad':
        atmosfera['rho'] = args.input1

    atmosfera['deltaT'] = args.input2

    atmosfera = atmosfera_estandar(args.input_type, atmosfera)
    print(atmosfera)