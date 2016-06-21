__author__ = 'Enriquito'
import sys
import atmosfera_estandar
import layout_atmosfera_estandar as layout
from PySide.QtCore import *
from PySide.QtGui import *

class AtmosferaEstandarDialog(QDialog, layout.Ui_Dialog):

    def __init__(self, presion=101325.0, temperatura=293.0, parent=None):
        super(AtmosferaEstandarDialog, self).__init__(parent)
        self.setupUi(self)
        self.atmosfera = {'h': {'SI': 0, 'IM': 0},
                          'deltaT': {'SI': 0, 'IM': 0},
                          'p': {'SI': 0, 'IM': 0},
                          't': {'SI': 0, 'IM': 0},
                          'rho': {'SI': 0, 'IM': 0},
                          'mu': {'SI': 0, 'IM': 0},
                          'Vson': {'SI': 0, 'IM': 0},
                          'calculo': {'SI': 0, 'IM': 0}}

        # Define input unit's label for SI and IM
        self.length_label = {'IM': 'ft', 'SI': 'm'}
        self.temp_label = {'IM': 'K', 'SI': 'R'}
        self.speed_label = {'IM': 'ft/s', 'SI': 'm/s'}
        self.den_label = {'IM': 'slug/ft^2', 'SI': 'kg/m^3'}
        self.pressure_label = {'IM': 'Pa', 'SI': 'slug/m^2'}

        self.ft2m = 0.3048
        self.lb2kg = 0.453592
        self.slugcuft2kgm3 = 515.379

        self.units = 'SI'
        self.R = 287.00
        self.gamma = 1.4
        self.atmosfera['h'][self.units] = float(self.lineEdit_h.text())
        self.atmosfera['deltaT'][self.units] = float(self.lineEdit_deltaT.text())
        self.atmosfera['p'][self.units] = float(self.lineEdit_p.text())
        self.atmosfera['t'][self.units] = float(self.lineEdit_t.text())
        self.atmosfera['rho'][self.units] = float(self.lineEdit_rho.text())
        self.atmosfera['mu'][self.units] = float(self.lineEdit_mu.text())
        self.atmosfera['Vson'][self.units] = float(self.lineEdit_Vson.text())

        self.lineEdit_p.setText(str(presion))
        self.actualizarP()
        self.lineEdit_t.setText(str(temperatura))
        self.actualizarT()

        self.connect(self.lineEdit_h, SIGNAL("editingFinished()"), self.actualizarH)
        self.connect(self.lineEdit_deltaT, SIGNAL("editingFinished()"), self.actualizarH)
        self.connect(self.lineEdit_t, SIGNAL("editingFinished()"), self.actualizarT)
        self.connect(self.lineEdit_p, SIGNAL("editingFinished()"), self.actualizarP)
        self.connect(self.lineEdit_rho, SIGNAL("editingFinished()"), self.actualizarRho)
        self.connect(self.acceptButton, SIGNAL("clicked()"), self, SLOT("accept()"))

    def actualizarH(self):
        # TODO: actualizar todos los datos, tanto en SI como en IM
        self.update_values()
        calculo = 'altura'
        results = atmosfera_estandar.atmosfera_estandar(calculo, self.atmosfera['h']['SI'],
                                                        self.atmosfera['deltaT']['SI'])
        # TODO: mejor aun, poner todo esto en una funcion
        # def actualizar_resultados(sefl, results):

        self.atmosfera['h']['SI'] = results[0]
        self.atmosfera['deltaT']['SI'] = results[1]
        self.atmosfera['p']['SI'] = results[2]
        self.atmosfera['t']['SI'] = results[3]
        self.atmosfera['rho']['SI'] = results[4]
        self.atmosfera['mu']['SI'] = results[5]
        self.atmosfera['Vson']['SI'] = results[6]
        self.si2im()
        self.update_lineEdits()

    def actualizarP(self):
        calculo = 'presion'
        results = atmosfera_estandar.atmosfera_estandar(calculo, float(self.lineEdit_p.text()),float(self.lineEdit_deltaT.text()))
        self.atmosfera['h'] = results[0]
        self.lineEdit_h.setText(str(round(results[0], 1)))
        self.atmosfera['deltaT'] = results[1]
        self.lineEdit_deltaT.setText(str(round(results[1], 2)))
        self.atmosfera['p'] = results[2]
        self.lineEdit_p.setText(str(round(results[2], 1)))
        self.atmosfera['t'] = results[3]
        self.lineEdit_t.setText(str(round(results[3], 2)))
        self.atmosfera['rho'] = results[4]
        self.lineEdit_rho.setText(str(round(results[4], 3)))
        self.atmosfera['mu'] = results[5]
        self.lineEdit_mu.setText(str(round(results[5], 7)))
        self.atmosfera['Vson'] = results[6]
        self.lineEdit_Vson.setText(str(round(results[6], 2)))

    def actualizarT(self):
        calculo = 'temperatura'
        temp = float(self.lineEdit_t.text())
        self.atmosfera['deltaT'] =  self.atmosfera['deltaT']+temp - self.atmosfera['t']
        self.lineEdit_deltaT.setText(str(round(self.atmosfera['deltaT'], 2)))
        self.atmosfera['t'] = temp
        self.atmosfera['rho'] = self.atmosfera['p'] / self.atmosfera['t']/self.R
        self.lineEdit_rho.setText(str(round(self.atmosfera['rho'], 3)))

    def actualizarRho(self):
        calculo = 'densidad'
        self.atmosfera['rho'] = float(round(self.lineEdit_rho.text(), 3))
        temp = self.atmosfera['p'] / self.atmosfera['rho']/self.R
        self.atmosfera['deltaT'] = self.atmosfera['deltaT'] + temp - self.atmosfera['t']
        self.lineEdit_deltaT.setText(str(round(self.atmosfera['deltaT'], 2)))
        self.atmosfera['t'] = temp
        self.lineEdit_t.setText(str(round(temp, 2)))

    def update_units(self):
        if self.IM_radioButton.isChecked():
            self.units = "IM"
        elif self.SI_radioButton.isChecked():
            self.units = "SI"
        else:
            return -1
        # Actualizo unitlabels
        self.unitlabel_h.setText(self.length_label[self.units])
        self.unitlabel_deltaT.setText(self.temp_label[self.units])
        self.unitlabel_p.setText(self.pressure_label[self.units])
        self.unitlabel_t.setText(self.temp_label[self.units])
        self.unitlabel_rho.setText(self.den_label[self.units])
        self.unitlabel_mu.setText(self.length_label[self.units])
        self.unitlabel_Vson.setText(self.speed_label[self.units])
        self.update_lineEdits()

    def update_lineEdits(self):
        self.lineEdit_h.setText(str(self.atmosfera['h'][self.units]))
        self.lineEdit_deltaT.setText(str(self.atmosfera['deltaT'][self.units]))
        self.lineEdit_p.setText(str(self.atmosfera['p'][self.units]))
        self.lineEdit_t.setText(str(self.atmosfera['t'][self.units]))
        self.lineEdit_rho.setText(str(self.atmosfera['rho'][self.units]))
        self.lineEdit_mu.setText(str(self.atmosfera['mu'][self.units]))
        self.lineEdit_Vson.setText(str(self.atmosfera['Vson'][self.units]))


    def update_atmosfera(self):
        self.read_from_lineEdits()
        if self.units == 'SI':
            self.si2im()
        elif self.units == 'IM':
            self.im2si()
        else:
            raise

    def read_from_lineEdits(self):
        pass # TODO

    def si2im(self):
        self.atmosfera['h']['IM'] = self.atmosfera['h']['SI']/self.ft2m
        self.atmosfera['deltaT']['IM'] = self.kelvin2rankine(self.atmosfera['deltaT']['SI'])
        #TODO: repetir para el resto de las variables

    def im2si(self):
        self.atmosfera['h']['SI'] = self.atmosfera['h']['IM']*self.ft2m
        self.atmosfera['deltaT']['SI'] = self.rankine2kelvin(self.atmosfera['deltaT']['IM'])
        #TODO: repetir para el resto de las variables

    @staticmethod
    def kelvin2rankine(temp):
        return temp*9/5.0

    @staticmethod
    def rakine2kelvin(temp):
        return temp*5/9.0
if __name__=='__main__':
    app = QApplication(sys.argv)
    dialogo = AtmosferaEstandarDialog()
    dialogo.show()
    app.exec_()