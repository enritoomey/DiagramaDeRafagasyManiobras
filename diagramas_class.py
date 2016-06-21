import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import
class Diagramas_de_maniobras_y_rafagas(object):
    def __init__(self, datos, w, h):
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
        self.w = w
        self.h = h
        self.den = atmosfera_estandar(self.h)
