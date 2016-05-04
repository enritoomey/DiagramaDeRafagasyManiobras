# -*- coding: utf-8 -*-
"""

@author: Enriquito
"""

import sys
import matplotlib  # Para los graficos
from PySide.QtCore import *
from PySide.QtGui import *  # importo todas las funciones de pyside
from matplotlib.backends.qt4_editor.formlayout import QDialog

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
# Estas lineas son un poco misteriosas, pero son las que me permite
# vincular los Widget de Qt/Pyside con matplotlib

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
# from excepciones import NumeroNegativoError,MayorAUnoError,TemperaturaIncompatibleError
import numpy as np # Para las cuentas
import layout_DiagramaDeRafagasyManiobras # importo las clases creadas con Qt y pyside
# import GUI_atmosfera_estandar # todavia no voy a usar esta clase
__appName__ = 'Diagrama de Rafagas y Maniobras'

# Creo la clase principal, llamada ""Main Dialog"
class DiagramaDeRafagasyManiobrasDialog(QFrame, layout_DiagramaDeRafagasyManiobras.Ui_Form):

    # Estas lineas son medias magicas, pero siempre van:
    def __init__(self,parent=None):
        super(DiagramaDeRafagasyManiobrasDialog, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle(__appName__)

        # Define input unit's label for SI and IM
        self.length_label = {'IM':'ft', 'SI':'m'}
        self.area_label = {'IM':'ft^2', 'SI':'m^2'}
        self.weight_label = {'IM':'lb', 'SI':'kg'}
        self.speed_label = {'IM':'ft/s', 'SI':'m/s'}
        self.den_label =  {'IM':'slug/ft^2', 'SI':'kg/m^3'}

        ft2m = 0.3048
        lb2kg = 0.453592
        slugcuft2kgm3 = 515.379
        # Creo algunas variables que generales que luego voy a usar
        self.CAM = {'SI':[],'IM':[]}
        self.sw = {'SI':[],'IM':[]}
        self.MTOW = {'SI':[],'IM':[]}
        self.MLW = {'SI':[],'IM':[]}
        self.W0 = {'SI':[],'IM':[]}
        self.MZFW = {'SI':[],'IM':[]}
        self.Vc = {'SI':[],'IM':[]}
        self.Zmo = {'SI':[],'IM':[]}


        self.CAM['IM'] = self.CAM['SI']/ft2m
        self.sw['IM']= self.sw['SI']/ft2m/ft2m
        self.MTOW['IM'] = self.MTOW['SI']/lb2kg
        self.MLW['IM'] = self.MLW['SI']/lb2kg
        self.W0['IM'] = self.W0['SI']/lb2kg
        self.MZFW['IM'] = self.MZFW['SI']/lb2kg
        self.Vc['IM'] = self.Vc['SI']/ft2m
        self.a3D = 5.0037 # 1/rad
        self.clmax = 1.2463
        self.clmax_flap = 1.499
        self.clmin = -0.75*self.clmax
        self.Zmo['IM'] = self.Zmo['SI']/ft2m

        # Variables
        self.W = {'SI':20000}
        self.h = {'SI':5000}
        self.den = {'SI':0.125}

        self.W['IM'] = self.W['SI']/lb2kg
        self.h['IM'] = self.h['SI']/ft2m
        self.den['IM'] = self.den['SI']/lb2kg*ft2m**3

        # constantes
        cte_fgz = {'IM':250000}
        cte_fgz['SI']=cte_fgz['IM']*ft2m
        s = {'IM':100.015}
        s['SI'] = s['IM']*ft2m
        gravedad = {'SI':9.81}
        gravedad['IM'] = gravedad['SI']*ft2m/lb2kg
        cte_nmax_1= {'IM':24000}
        cte_nmax_1['SI'] = cte_nmax_1['IM']*lb2kg
        cte_nmax_2= {'IM':10000}
        cte_nmax_2['SI'] = cte_nmax_1['IM']*lb2kg

        # Actualizo los labels:
        self.update_units()

        # SIGNALS & SLOTS
        self.connect(self.IM_radioButton,SIGNAL("clicked()"), self.update_units)
        self.connect(self.SI_radioButton,SIGNAL("clicked()"), self.update_units)

    def update_units(self):
        if self.IM_radioButton.isChecked():
            self.units = "IM"
            print "IM"
        elif self.SI_radioButton.isChecked():
            self.units = "SI"
            print "SI"
        else: return -1;

        print(self.units)

        # Actualizo unitlabels
        self.CAM_unitlabel.setText(self.length_label[self.units])
        self.sw_unitlabel.setText(self.area_label[self.units])
        self.MTOW_unitlabel.setText(self.weight_label[self.units])
        self.MLW_unitlabel.setText(self.weight_label[self.units])
        self.MZFW_unitlabel.setText(self.weight_label[self.units])
        self.W0_unitlabel.setText(self.weight_label[self.units])
        self.W_unitlabel.setText(self.weight_label[self.units])
        self.Zmo_unitlabel.setText(self.length_label[self.units])
        self.h_unitlabel.setText(self.length_label[self.units])
        self.den_unitlabel.setText(self.den_label[self.units])
        self.Vc_unitlabel.setText(self.speed_label[self.units])

        # Actualizo los valores
        self.CAM_lineEdit.setText(self.CAM[self.units])
        self.sw_lineEdit.setText(self.sw[self.units])
        self.MTOW_lineEdit.setText(self.MROW[self.units])
        self.MLW_lineEdit.setText(self.MLW[self.units])
        self.MZFW_lineEdit.setText(self.MZFW[self.units])
        self.W0_lineEdit.setText(self.W0[self.units])
        self.W_lineEdit.setText(self.W[self.units])
        self.Zmo_lineEdit.setText(self.Zmo[self.units])
        self.h_lineEdit.setText(self.h[self.units])
        self.den_lineEdit.setText(self.den[self.units])
        self.Vc_lineEdit.setText(self.Vc[self.units])

        # Actualizo Lineedits


    # self.Q = []
    # self.gamma_aire = 1.4
    # self.R_aire = 287
    # #Defino las formulas quimicas para selecionar en formulas_comboBox
    # self.formulas = {"AvGas":{'c':7,'h':16,'o':0,'s':0.00},"Diesel":{'c':12,'h':21,'o':0,'s':0.00}}
    # # Utilizo las 'keys' de self.formulas para cargar los items en la combobox
    # for key in self.formulas.keys():
    #     self.formula_comboBox.addItem(key)
    #
    # self.actualizarFormula()

app = QApplication(sys.argv)
form = DiagramaDeRafagasyManiobrasDialog()
form.show()
app.exec_()