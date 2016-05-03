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
        # Creo algunas variables que generales que luego voy a usar

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