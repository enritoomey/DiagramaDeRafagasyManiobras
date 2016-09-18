# -*- coding: utf-8 -*-
"""

@author: Enriquito
"""
import logging
import sys
import matplotlib.pyplot as plt
sys.path.append('./atmosfera_estandar')

import matplotlib  # Para los graficos
from PySide.QtCore import *
from PySide.QtGui import *  # importo todas las funciones de pyside
# from matplotlib.backends.qt4_editor.formlayout import QDialog
from diagramas_class import Diagramas
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'
# Estas lineas son un poco misteriosas, pero son las que me permite
# vincular los Widget de Qt/Pyside con matplotlib

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from diagramas_class import Diagramas

# from excepciones import NumeroNegativoError,MayorAUnoError,TemperaturaIncompatibleError
import layout_DiagramaDeRafagasyManiobras  # importo las clases creadas con Qt y pyside
from atmosfera_estandar import GUI_atmosfera_estandar

__appName__ = 'Diagrama de Rafagas y Maniobras'

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# Creo la clase principal, llamada ""Main Dialog"
class DiagramaDeRafagasyManiobrasDialog(QFrame, layout_DiagramaDeRafagasyManiobras.Ui_Dialog):

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

        datos = {
            'CAM': {'SI': 2.461, 'IM': 2.461/self.ft2m},
            'sw': {'SI': 60, 'IM': 60 / self.ft2m / self.ft2m},
            'a3D': 5.0037,
            'MTOW': {'SI': 23000, 'IM': 23000 / self.lb2kg},
            'MLW': {'SI': 23000, 'IM': 23000 / self.lb2kg},
            'W0': {'SI': 13766.0, 'IM': 13766 / self.lb2kg},
            'MZFW': {'SI': 16376.0, 'IM': 16376.0 / self.lb2kg},
            'Vc': {'SI': 151.93, 'IM': 151.93 / self.ft2m},
            'clmax': 1.2463,
            'clmax_flap': 1.499,
            'clmin': -0.75*1.2463,
            'Zmo': {'SI': 9999.2, 'IM': 9999.2 / self.ft2m}
        }
        w = {'SI': 20000, 'IM': 20000 / self.lb2kg}
        h = {'SI': 5000, 'IM': 5000 / self.ft2m}
        den = {'SI': 1.225, 'IM': 1.225 / self.lb2kg * self.ft2m**3}
        self.diagramas = Diagramas(datos, w, h, den, units='SI')

        self.update_units()

        # Generamos dos figuras, cada una luego asociada a un canvas, que a su vez tiene como padre una de las pestañas
        # self.tab -> contiene la pestaña titulada "Diagrama P-S"
        # self.tab_2 -> contiene la pestaña titulada "Diagrama T-S"
        self.fig1 = Figure((5.0, 3.0), dpi=72, facecolor=(1, 1, 1), edgecolor=(0, 0, 0))
        self.fig1.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.25)
        self.axes1 = self.fig1.add_subplot(111)
        self.axes1.set_ylabel('n')
        self.axes1.set_xlabel(self.speed_label[self.units])
        self.axes1.set_title('Diagrama de Maniobras')
        self.axes1.ticklabel_format(style="sci", scilimits=(0, 0), axis="both")  # , useOffset=True,useLocale=True)
        self.axes1.tick_params(axis="both", direction='in', length=6, width=2, labelsize="medium")

        self.fig2 = Figure(dpi=72, facecolor=(1, 1, 1), edgecolor=(0, 0, 0))
        self.fig2.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.25)
        self.axes2 = self.fig2.add_subplot(111)
        self.axes2.set_ylabel('n')
        self.axes2.ticklabel_format(style='sci', scilimits=(0, 0), axis="both")
        self.axes2.set_xlabel(self.speed_label[self.units])
        self.axes2.set_title('Diagrama de Rafagas')

        self.fig3 = Figure(dpi=72, facecolor=(1, 1, 1), edgecolor=(0, 0, 0))
        self.fig3.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.25)
        self.axes3 = self.fig3.add_subplot(111)
        self.axes3.set_ylabel('n')
        self.axes3.ticklabel_format(style='sci', scilimits=(0, 0), axis="both")
        self.axes3.set_xlabel(self.speed_label[self.units])
        self.axes3.set_title('Diagrama de Rafagas y maniobras')

        # generate the canvas to display the plot
        self.canvas1 = FigureCanvas(self.fig1)
        self.canvas1.setParent(self.manoeuvre_tab)
        self.canvas1.show()
        self.canvas2 = FigureCanvas(self.fig2)
        self.canvas2.setParent(self.gust_tab)
        self.canvas2.show()
        self.canvas3 = FigureCanvas(self.fig3)
        self.canvas3.setParent(self.combined_tab)
        self.canvas3.show()

        # SIGNALS & SLOTS
        self.connect(self.IM_radioButton, SIGNAL("clicked()"), self.update_units)
        self.connect(self.SI_radioButton, SIGNAL("clicked()"), self.update_units)
        self.connect(self.Altura_Button, SIGNAL("clicked()"), self.seleccionAltura)
        self.connect(self.grafiacar_pushButton, SIGNAL("clicked()"), self.Calculos())

        self.CAM_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.CAM, float(self.CAM_lineEdit.text()), self.ft2m))
        self.sw_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.sw, float(self.sw_lineEdit.text()), self.ft2m**2))
        self.MTOW_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.MTOW, float(self.MTOW_lineEdit.text()), self.lb2kg))
        self.MLW_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.MLW, float(self.MLW_lineEdit.text()), self.lb2kg))
        self.MZFW_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.MZFW, float(self.MZFW_lineEdit.text()), self.lb2kg))
        self.W0_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.W0, float(self.W0_lineEdit.text()), self.lb2kg))

        self.a3D_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.a3D, float(self.a3D_lineEdit.text())))
        self.clmax_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.clmax, float(self.clmax_lineEdit.text())))
        self.clmax_flap_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.clmax_flap, float(self.clmax_flap_lineEdit.text())))

        self.Zmo_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.Zmo, float(self.Zmo_lineEdit.text()), self.ft2m))
        self.Vc_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.Vc, float(self.Vc_lineEdit.text()), self.ft2m))

        self.W_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.W, float(self.W_lineEdit.text()), self.lb2kg))
        self.h_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.h, float(self.h_lineEdit.text()), self.ft2m))
        self.den_lineEdit.editingFinished.connect(lambda: self.lecturadatos(self.diagramas.den, float(self.den_lineEdit.text()), self.lb2kg / self.ft2m**3))

        self.grafiacar_pushButton.clicked.connect(self.Calculos)

    def resizeEvent(self, event):
        self.canvas1.setGeometry(self.PlotArea.rect())
        self.canvas2.setGeometry(self.PlotArea.rect())
        self.canvas3.setGeometry(self.PlotArea.rect())

    def lecturadatos(self, variable, value, unitConverter = 1.0):
        logger.debug(value)
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
            self.diagramas.units = "IM"
        elif self.SI_radioButton.isChecked():
            self.units = "SI"
            self.diagramas.units = "SI"
        else:
            return -1
        self.update_unitLabels()
        self.write_lineEdits()

    def update_unitLabels(self):
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

    def write_lineEdits(self):
        self.CAM_lineEdit.setText(str(self.diagramas.CAM[self.units]))
        self.sw_lineEdit.setText(str(self.diagramas.sw[self.units]))
        self.MTOW_lineEdit.setText(str(self.diagramas.MTOW[self.units]))
        self.MLW_lineEdit.setText(str(self.diagramas.MLW[self.units]))
        self.MZFW_lineEdit.setText(str(self.diagramas.MZFW[self.units]))
        self.W0_lineEdit.setText(str(self.diagramas.W0[self.units]))
        self.W_lineEdit.setText(str(self.diagramas.W[self.units]))
        self.Zmo_lineEdit.setText(str(self.diagramas.Zmo[self.units]))
        self.h_lineEdit.setText(str(self.diagramas.h[self.units]))
        self.den_lineEdit.setText(str(self.diagramas.den[self.units]))
        self.Vc_lineEdit.setText(str(self.diagramas.Vc[self.units]))
        self.a3D_lineEdit.setText(str(self.diagramas.a3D))
        self.clmax_lineEdit.setText(str(self.diagramas.clmax))
        self.clmax_flap_lineEdit.setText(str(self.diagramas.clmax_flap))
        self.clmin_lineEdit.setText(str(self.diagramas.clmin))


    def Calculos(self):
        logger.info("CAM = {}".format(self.diagramas.CAM[self.units]))
        logger.info("Sw = {}".format(self.diagramas.sw[self.units]))
        logger.info("MTOW = {}".format(self.diagramas.MTOW[self.units]))
        logger.info("MLW = {}".format(self.diagramas.MLW[self.units]))
        logger.info("MZFW = {}".format(self.diagramas.MZFW[self.units]))
        logger.info("W0 = {}".format(self.diagramas.W0[self.units]))
        logger.info("a3D = {}".format(self.diagramas.a3D))
        logger.info("clmax = {}".format(self.diagramas.clmax))
        logger.info("clmax_flap = {}".format(self.diagramas.clmax_flap))
        logger.info("clmin = {}".format(self.diagramas.clmin))
        logger.info("Zmo = {}".format(self.diagramas.Zmo[self.units]))
        logger.info("Vc = {}".format(self.diagramas.Vc[self.units]))

        self.diagramas.calculos()

        self.axes1.clear()
        self.diagramas.plot_diagrama_de_maniobras(self.axes1, 0.5)
        self.diagramas.plot_diagrama_de_maniobras_con_flap(self.axes1, 0.5)
        self.canvas1.draw()

        self.axes2.clear()
        self.diagramas.plot_diagrama_de_rafagas(self.axes2, 0.5)
        self.canvas2.draw()

        self.axes3.clear()
        self.diagramas.plot_diagrama_de_maniobras_y_rafagas(self.axes3, 0.5)
        self.canvas3.draw()

    def seleccionAltura(self):
        dialogo = GUI_atmosfera_estandar.AtmosferaEstandarDialog(unit=self.units)
        if dialogo.exec_():
            self.lecturadatos(self.diagramas.h, dialogo.atmosfera[self.units]['h'], unitConverter=self.ft2m)
            self.lecturadatos(self.diagramas.den, dialogo.atmosfera[self.units]['rho'], unitConverter=self.slugcuft2kgm3)
            self.write_lineEdits()

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