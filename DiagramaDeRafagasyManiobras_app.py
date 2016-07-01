# -*- coding: utf-8 -*-
"""

@author: Enriquito
"""

import sys
sys.path.append('./atmosfera_estandar')

import matplotlib  # Para los graficos
from PySide.QtCore import *
from PySide.QtGui import *  # importo todas las funciones de pyside
# from matplotlib.backends.qt4_editor.formlayout import QDialog

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'
# Estas lineas son un poco misteriosas, pero son las que me permite
# vincular los Widget de Qt/Pyside con matplotlib

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
# from excepciones import NumeroNegativoError,MayorAUnoError,TemperaturaIncompatibleError
import numpy as np # Para las cuentas
import layout_DiagramaDeRafagasyManiobras # importo las clases creadas con Qt y pyside
import GUI_atmosfera_estandar # todavia no voy a usar esta clase
__appName__ = 'Diagrama de Rafagas y Maniobras'

# Creo la clase principal, llamada ""Main Dialog"
class DiagramaDeRafagasyManiobrasDialog(QFrame, layout_DiagramaDeRafagasyManiobras.Ui_Form):

    # Estas lineas son medias magicas, pero siempre van:
    def __init__(self, parent=None):
        super(DiagramaDeRafagasyManiobrasDialog, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle(__appName__)

        # Define input unit's label for SI and IM
        self.length_label = {'IM': 'ft', 'SI': 'm'}
        self.area_label = {'IM': 'ft^2', 'SI': 'm^2'}
        self.weight_label = {'IM': 'lb', 'SI': 'kg'}
        self.speed_label = {'IM': 'ft/s', 'SI': 'm/s'}
        self.den_label = {'IM': 'slug/ft^2', 'SI': 'kg/m^3'}

        self.ft2m = 0.3048
        self.lb2kg = 0.453592
        self.slugcuft2kgm3 = 515.379
        # Creo algunas variables generales que luego voy a usar
        self.CAM = {'SI': 0, 'IM': 0}
        self.sw = {'SI': 0, 'IM': 0}
        self.MTOW = {'SI': 0, 'IM': 0}
        self.MLW = {'SI': 0, 'IM': 0}
        self.W0 = {'SI': 0, 'IM': 0}
        self.MZFW = {'SI': 0, 'IM': 0}
        self.Vc = {'SI': 0, 'IM': 0}
        self.Zmo = {'SI': 0, 'IM': 0}
        self.W = {'SI': 0, 'IM': 0}
        self.h = {'SI': 0, 'IM': 0}
        self.den = {'SI': 0, 'IM': 0}

        self.a3D = 5.0037# 1/rad
        self.clmax = 1.2463
        self.clmax_flap = 1.499
        self.clmin = -0.75*self.clmax

        # input_data = {self.CAM:self.CAM_lineEdit, self.sw:self.sw_lineEdit}
        # constantes
        cte_fgz = {'IM': 250000}
        cte_fgz['SI'] = cte_fgz['IM']*self.ft2m
        s = {'IM': 100.015}
        s['SI'] = s['IM']*self.ft2m
        gravedad = {'SI': 9.81}
        gravedad['IM'] = gravedad['SI']*self.ft2m/self.lb2kg
        cte_nmax_1 = {'IM': 24000}
        cte_nmax_1['SI'] = cte_nmax_1['IM']*self.lb2kg
        cte_nmax_2 = {'IM': 10000}
        cte_nmax_2['SI'] = cte_nmax_1['IM']*self.lb2kg

        # Actualizo los labels:
        self.update_units()

        # SIGNALS & SLOTS
        self.connect(self.IM_radioButton, SIGNAL("clicked()"), self.update_units)
        self.connect(self.SI_radioButton, SIGNAL("clicked()"), self.update_units)
        self.connect(self.Altura_Button, SIGNAL("clicked()"), self.seleccionAltura)

        # self.CAM_lineEdit.textChanged()
        # self.connect(self.CAM_lineEdit,SIGNAL("textEdited(const dict&, const QString&)"),
        #              self.lecturadatos)
        self.CAM_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.CAM, float(self.CAM_lineEdit.text()), self.ft2m))
        self.sw_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.sw, float(self.sw_lineEdit.text()), self.ft2m**2))
        self.MTOW_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.MTOW, float(self.MTOW_lineEdit.text()), self.lb2kg))
        self.MLW_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.MLW, float(self.MLW_lineEdit.text()), self.lb2kg))
        self.MZFW_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.MZFW, float(self.MZFW_lineEdit.text()), self.lb2kg))
        self.W0_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.W0, float(self.W0_lineEdit.text()), self.lb2kg))

        self.a3D_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.a3D, float(self.a3D_lineEdit.text()), []))
        self.clmax_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.clmax, float(self.clmax_lineEdit.text())))
        self.clmax_flap_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.clmax_flap, float(self.clmax_flap_lineEdit.text())))

        self.Zmo_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.Zmo, float(self.Zmo_lineEdit.text()), self.ft2m))
        self.Vc_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.Vc, float(self.Vc_lineEdit.text()), self.ft2m))

        self.W_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.W, float(self.W_lineEdit.text()), self.lb2kg))
        self.h_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.h, float(self.h_lineEdit.text()), self.ft2m))
        self.den_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.den, float(self.den_lineEdit.text()), self.lb2kg*self.ft2m**3))

        self.grafiacar_pushButton.clicked.connect(self.Calculos)

    def lecturadatos(self, variable, value, unitConverter = []):
        print(value)
        if not unitConverter:
            variable = value
        elif self.units == 'IM':
            variable['IM'] = value
            variable['SI'] = value*unitConverter
        else:
            variable['SI'] = value
            variable['IM'] = value/unitConverter

    def update_units(self):
        if self.IM_radioButton.isChecked():
            self.units = "IM"
        elif self.SI_radioButton.isChecked():
            self.units = "SI"
        else:
            return -1

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
        self.CAM_lineEdit.setText(str(self.CAM[self.units]))
        self.sw_lineEdit.setText(str(self.sw[self.units]))
        self.MTOW_lineEdit.setText(str(self.MTOW[self.units]))
        self.MLW_lineEdit.setText(str(self.MLW[self.units]))
        self.MZFW_lineEdit.setText(str(self.MZFW[self.units]))
        self.W0_lineEdit.setText(str(self.W0[self.units]))
        self.W_lineEdit.setText(str(self.W[self.units]))
        self.Zmo_lineEdit.setText(str(self.Zmo[self.units]))
        self.h_lineEdit.setText(str(self.h[self.units]))
        self.den_lineEdit.setText(str(self.den[self.units]))
        self.Vc_lineEdit.setText(str(self.Vc[self.units]))

    def Calculos(self):
        print("CAM = {}".format(self.CAM[self.units]))
        print("Sw = {}".format(self.sw[self.units]))
        print("MTOW = {}".format(self.MTOW[self.units]))
        print("MLW = {}".format(self.MLW[self.units]))
        print("MZFW = {}".format(self.MZFW[self.units]))
        print("W0 = {}".format(self.W0[self.units]))
        print("a3D = {}".format(self.a3D))
        print("clmax = {}".format(self.clmax))
        print("clmax_flap = {}".format(self.clmax_flap))
        print("clmin = {}".format(self.clmin))
        print("Zmo = {}".format(self.Zmo[self.units]))
        print("Vc = {}".format(self.Vc[self.units]))

    def seleccionAltura(self):
        dialogo = GUI_atmosfera_estandar.AtmosferaEstandarDialog(unit=self.units)
        if dialogo.exec_():
            self.h[self.units] = dialogo.atmosfera[self.units]['h']
            self.den[self.units] = dialogo.atmosfera[self.units]['rho']
            self.h_lineEdit.setText(str(self.h[self.units]))
            self.den_lineEdit.setText(str(self.den[self.units]))


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