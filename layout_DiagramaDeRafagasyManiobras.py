# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'layout_DiagramaDeRafagasyManiobras.ui'
#
# Created: Tue Jun 21 01:09:53 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.4
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(716, 545)
        self.PlotArea = QtGui.QTabWidget(Form)
        self.PlotArea.setGeometry(QtCore.QRect(250, 30, 431, 361))
        self.PlotArea.setObjectName("PlotArea")
        self.gust_tab = QtGui.QWidget()
        self.gust_tab.setObjectName("gust_tab")
        self.PlotArea.addTab(self.gust_tab, "")
        self.manoeuvre_tab = QtGui.QWidget()
        self.manoeuvre_tab.setObjectName("manoeuvre_tab")
        self.PlotArea.addTab(self.manoeuvre_tab, "")
        self.combined_tab = QtGui.QWidget()
        self.combined_tab.setObjectName("combined_tab")
        self.PlotArea.addTab(self.combined_tab, "")
        self.InputData_label = QtGui.QLabel(Form)
        self.InputData_label.setGeometry(QtCore.QRect(20, 10, 181, 31))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setWeight(75)
        font.setUnderline(True)
        font.setBold(True)
        self.InputData_label.setFont(font)
        self.InputData_label.setObjectName("InputData_label")
        self.layoutWidget = QtGui.QWidget(Form)
        self.layoutWidget.setGeometry(QtCore.QRect(10, 140, 194, 151))
        self.layoutWidget.setObjectName("layoutWidget")
        self.W_group_layout = QtGui.QVBoxLayout(self.layoutWidget)
        self.W_group_layout.setContentsMargins(0, 0, 0, 0)
        self.W_group_layout.setObjectName("W_group_layout")
        self.MTOW_layout = QtGui.QHBoxLayout()
        self.MTOW_layout.setObjectName("MTOW_layout")
        self.MTOW_label = QtGui.QLabel(self.layoutWidget)
        self.MTOW_label.setObjectName("MTOW_label")
        self.MTOW_layout.addWidget(self.MTOW_label)
        self.MTOW_lineEdit = QtGui.QLineEdit(self.layoutWidget)
        self.MTOW_lineEdit.setObjectName("MTOW_lineEdit")
        self.MTOW_layout.addWidget(self.MTOW_lineEdit)
        self.MTOW_unitlabel = QtGui.QLabel(self.layoutWidget)
        self.MTOW_unitlabel.setObjectName("MTOW_unitlabel")
        self.MTOW_layout.addWidget(self.MTOW_unitlabel)
        self.W_group_layout.addLayout(self.MTOW_layout)
        self.MLW_layout = QtGui.QHBoxLayout()
        self.MLW_layout.setObjectName("MLW_layout")
        self.MLW_label = QtGui.QLabel(self.layoutWidget)
        self.MLW_label.setObjectName("MLW_label")
        self.MLW_layout.addWidget(self.MLW_label)
        self.MLW_lineEdit = QtGui.QLineEdit(self.layoutWidget)
        self.MLW_lineEdit.setObjectName("MLW_lineEdit")
        self.MLW_layout.addWidget(self.MLW_lineEdit)
        self.MLW_unitlabel = QtGui.QLabel(self.layoutWidget)
        self.MLW_unitlabel.setObjectName("MLW_unitlabel")
        self.MLW_layout.addWidget(self.MLW_unitlabel)
        self.W_group_layout.addLayout(self.MLW_layout)
        self.MZFW_layout = QtGui.QHBoxLayout()
        self.MZFW_layout.setObjectName("MZFW_layout")
        self.MZFW_label = QtGui.QLabel(self.layoutWidget)
        self.MZFW_label.setObjectName("MZFW_label")
        self.MZFW_layout.addWidget(self.MZFW_label)
        self.MZFW_lineEdit = QtGui.QLineEdit(self.layoutWidget)
        self.MZFW_lineEdit.setObjectName("MZFW_lineEdit")
        self.MZFW_layout.addWidget(self.MZFW_lineEdit)
        self.MZFW_unitlabel = QtGui.QLabel(self.layoutWidget)
        self.MZFW_unitlabel.setObjectName("MZFW_unitlabel")
        self.MZFW_layout.addWidget(self.MZFW_unitlabel)
        self.W_group_layout.addLayout(self.MZFW_layout)
        self.W0_layout = QtGui.QHBoxLayout()
        self.W0_layout.setObjectName("W0_layout")
        self.W0_label = QtGui.QLabel(self.layoutWidget)
        self.W0_label.setObjectName("W0_label")
        self.W0_layout.addWidget(self.W0_label)
        self.W0_lineEdit = QtGui.QLineEdit(self.layoutWidget)
        self.W0_lineEdit.setObjectName("W0_lineEdit")
        self.W0_layout.addWidget(self.W0_lineEdit)
        self.W0_unitlabel = QtGui.QLabel(self.layoutWidget)
        self.W0_unitlabel.setObjectName("W0_unitlabel")
        self.W0_layout.addWidget(self.W0_unitlabel)
        self.W_group_layout.addLayout(self.W0_layout)
        self.layoutWidget1 = QtGui.QWidget(Form)
        self.layoutWidget1.setGeometry(QtCore.QRect(10, 310, 200, 136))
        self.layoutWidget1.setObjectName("layoutWidget1")
        self.Cl_group_layout = QtGui.QVBoxLayout(self.layoutWidget1)
        self.Cl_group_layout.setContentsMargins(0, 0, 0, 0)
        self.Cl_group_layout.setObjectName("Cl_group_layout")
        self.a3D_layout = QtGui.QHBoxLayout()
        self.a3D_layout.setObjectName("a3D_layout")
        self.a3D_label = QtGui.QLabel(self.layoutWidget1)
        self.a3D_label.setObjectName("a3D_label")
        self.a3D_layout.addWidget(self.a3D_label)
        self.a3D_lineEdit = QtGui.QLineEdit(self.layoutWidget1)
        self.a3D_lineEdit.setObjectName("a3D_lineEdit")
        self.a3D_layout.addWidget(self.a3D_lineEdit)
        self.Cl_group_layout.addLayout(self.a3D_layout)
        self.clmax_layout = QtGui.QHBoxLayout()
        self.clmax_layout.setObjectName("clmax_layout")
        self.clmax_label = QtGui.QLabel(self.layoutWidget1)
        self.clmax_label.setObjectName("clmax_label")
        self.clmax_layout.addWidget(self.clmax_label)
        self.clmax_lineEdit = QtGui.QLineEdit(self.layoutWidget1)
        self.clmax_lineEdit.setObjectName("clmax_lineEdit")
        self.clmax_layout.addWidget(self.clmax_lineEdit)
        self.Cl_group_layout.addLayout(self.clmax_layout)
        self.clmax_flap_layout = QtGui.QHBoxLayout()
        self.clmax_flap_layout.setObjectName("clmax_flap_layout")
        self.clmax_flap_label = QtGui.QLabel(self.layoutWidget1)
        self.clmax_flap_label.setObjectName("clmax_flap_label")
        self.clmax_flap_layout.addWidget(self.clmax_flap_label)
        self.clmax_flap_lineEdit = QtGui.QLineEdit(self.layoutWidget1)
        self.clmax_flap_lineEdit.setObjectName("clmax_flap_lineEdit")
        self.clmax_flap_layout.addWidget(self.clmax_flap_lineEdit)
        self.Cl_group_layout.addLayout(self.clmax_flap_layout)
        self.clmin_layout = QtGui.QHBoxLayout()
        self.clmin_layout.setObjectName("clmin_layout")
        self.clmin_label = QtGui.QLabel(self.layoutWidget1)
        self.clmin_label.setObjectName("clmin_label")
        self.clmin_layout.addWidget(self.clmin_label)
        self.clmin_lineEdit = QtGui.QLineEdit(self.layoutWidget1)
        self.clmin_lineEdit.setObjectName("clmin_lineEdit")
        self.clmin_layout.addWidget(self.clmin_lineEdit)
        self.Cl_group_layout.addLayout(self.clmin_layout)
        self.layoutWidget2 = QtGui.QWidget(Form)
        self.layoutWidget2.setGeometry(QtCore.QRect(250, 400, 196, 123))
        self.layoutWidget2.setObjectName("layoutWidget2")
        self.verticalLayout = QtGui.QVBoxLayout(self.layoutWidget2)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.W_layout = QtGui.QHBoxLayout()
        self.W_layout.setObjectName("W_layout")
        self.W_label = QtGui.QLabel(self.layoutWidget2)
        self.W_label.setObjectName("W_label")
        self.W_layout.addWidget(self.W_label)
        self.W_lineEdit = QtGui.QLineEdit(self.layoutWidget2)
        self.W_lineEdit.setObjectName("W_lineEdit")
        self.W_layout.addWidget(self.W_lineEdit)
        self.W_unitlabel = QtGui.QLabel(self.layoutWidget2)
        self.W_unitlabel.setObjectName("W_unitlabel")
        self.W_layout.addWidget(self.W_unitlabel)
        self.verticalLayout.addLayout(self.W_layout)
        self.Altura_Button = QtGui.QPushButton(self.layoutWidget2)
        self.Altura_Button.setObjectName("Altura_Button")
        self.verticalLayout.addWidget(self.Altura_Button)
        self.h_layout = QtGui.QHBoxLayout()
        self.h_layout.setObjectName("h_layout")
        self.h_label = QtGui.QLabel(self.layoutWidget2)
        self.h_label.setObjectName("h_label")
        self.h_layout.addWidget(self.h_label)
        self.h_lineEdit = QtGui.QLineEdit(self.layoutWidget2)
        self.h_lineEdit.setReadOnly(True)
        self.h_lineEdit.setObjectName("h_lineEdit")
        self.h_layout.addWidget(self.h_lineEdit)
        self.h_unitlabel = QtGui.QLabel(self.layoutWidget2)
        self.h_unitlabel.setObjectName("h_unitlabel")
        self.h_layout.addWidget(self.h_unitlabel)
        self.verticalLayout.addLayout(self.h_layout)
        self.den_layout = QtGui.QHBoxLayout()
        self.den_layout.setObjectName("den_layout")
        self.den_label = QtGui.QLabel(self.layoutWidget2)
        self.den_label.setObjectName("den_label")
        self.den_layout.addWidget(self.den_label)
        self.den_lineEdit = QtGui.QLineEdit(self.layoutWidget2)
        self.den_lineEdit.setReadOnly(True)
        self.den_lineEdit.setObjectName("den_lineEdit")
        self.den_layout.addWidget(self.den_lineEdit)
        self.den_unitlabel = QtGui.QLabel(self.layoutWidget2)
        self.den_unitlabel.setObjectName("den_unitlabel")
        self.den_layout.addWidget(self.den_unitlabel)
        self.verticalLayout.addLayout(self.den_layout)
        self.layoutWidget3 = QtGui.QWidget(Form)
        self.layoutWidget3.setGeometry(QtCore.QRect(10, 60, 191, 66))
        self.layoutWidget3.setObjectName("layoutWidget3")
        self.geometry_group_layout = QtGui.QVBoxLayout(self.layoutWidget3)
        self.geometry_group_layout.setContentsMargins(0, 0, 0, 0)
        self.geometry_group_layout.setObjectName("geometry_group_layout")
        self.CAM_layout = QtGui.QHBoxLayout()
        self.CAM_layout.setObjectName("CAM_layout")
        self.CAM_label = QtGui.QLabel(self.layoutWidget3)
        self.CAM_label.setObjectName("CAM_label")
        self.CAM_layout.addWidget(self.CAM_label)
        self.CAM_lineEdit = QtGui.QLineEdit(self.layoutWidget3)
        self.CAM_lineEdit.setObjectName("CAM_lineEdit")
        self.CAM_layout.addWidget(self.CAM_lineEdit)
        self.CAM_unitlabel = QtGui.QLabel(self.layoutWidget3)
        self.CAM_unitlabel.setObjectName("CAM_unitlabel")
        self.CAM_layout.addWidget(self.CAM_unitlabel)
        self.geometry_group_layout.addLayout(self.CAM_layout)
        self.sw_layout = QtGui.QHBoxLayout()
        self.sw_layout.setObjectName("sw_layout")
        self.sw_label = QtGui.QLabel(self.layoutWidget3)
        self.sw_label.setObjectName("sw_label")
        self.sw_layout.addWidget(self.sw_label)
        self.sw_lineEdit = QtGui.QLineEdit(self.layoutWidget3)
        self.sw_lineEdit.setObjectName("sw_lineEdit")
        self.sw_layout.addWidget(self.sw_lineEdit)
        self.sw_unitlabel = QtGui.QLabel(self.layoutWidget3)
        self.sw_unitlabel.setObjectName("sw_unitlabel")
        self.sw_layout.addWidget(self.sw_unitlabel)
        self.geometry_group_layout.addLayout(self.sw_layout)
        self.layoutWidget4 = QtGui.QWidget(Form)
        self.layoutWidget4.setGeometry(QtCore.QRect(10, 460, 191, 66))
        self.layoutWidget4.setObjectName("layoutWidget4")
        self.performace_group_layout = QtGui.QVBoxLayout(self.layoutWidget4)
        self.performace_group_layout.setContentsMargins(0, 0, 0, 0)
        self.performace_group_layout.setObjectName("performace_group_layout")
        self.Zmo_layout = QtGui.QHBoxLayout()
        self.Zmo_layout.setObjectName("Zmo_layout")
        self.Zmo_label = QtGui.QLabel(self.layoutWidget4)
        self.Zmo_label.setObjectName("Zmo_label")
        self.Zmo_layout.addWidget(self.Zmo_label)
        self.Zmo_lineEdit = QtGui.QLineEdit(self.layoutWidget4)
        self.Zmo_lineEdit.setObjectName("Zmo_lineEdit")
        self.Zmo_layout.addWidget(self.Zmo_lineEdit)
        self.Zmo_unitlabel = QtGui.QLabel(self.layoutWidget4)
        self.Zmo_unitlabel.setObjectName("Zmo_unitlabel")
        self.Zmo_layout.addWidget(self.Zmo_unitlabel)
        self.performace_group_layout.addLayout(self.Zmo_layout)
        self.Vc_layout = QtGui.QHBoxLayout()
        self.Vc_layout.setObjectName("Vc_layout")
        self.Vc_label = QtGui.QLabel(self.layoutWidget4)
        self.Vc_label.setObjectName("Vc_label")
        self.Vc_layout.addWidget(self.Vc_label)
        self.Vc_lineEdit = QtGui.QLineEdit(self.layoutWidget4)
        self.Vc_lineEdit.setObjectName("Vc_lineEdit")
        self.Vc_layout.addWidget(self.Vc_lineEdit)
        self.Vc_unitlabel = QtGui.QLabel(self.layoutWidget4)
        self.Vc_unitlabel.setObjectName("Vc_unitlabel")
        self.Vc_layout.addWidget(self.Vc_unitlabel)
        self.performace_group_layout.addLayout(self.Vc_layout)
        self.grafiacar_pushButton = QtGui.QPushButton(Form)
        self.grafiacar_pushButton.setGeometry(QtCore.QRect(580, 410, 111, 51))
        self.grafiacar_pushButton.setObjectName("grafiacar_pushButton")
        self.units_groupBox = QtGui.QGroupBox(Form)
        self.units_groupBox.setGeometry(QtCore.QRect(460, 409, 121, 60))
        self.units_groupBox.setCheckable(False)
        self.units_groupBox.setObjectName("units_groupBox")
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.units_groupBox)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.SI_radioButton = QtGui.QRadioButton(self.units_groupBox)
        self.SI_radioButton.setEnabled(True)
        self.SI_radioButton.setChecked(True)
        self.SI_radioButton.setObjectName("SI_radioButton")
        self.horizontalLayout_2.addWidget(self.SI_radioButton)
        self.IM_radioButton = QtGui.QRadioButton(self.units_groupBox)
        self.IM_radioButton.setObjectName("IM_radioButton")
        self.horizontalLayout_2.addWidget(self.IM_radioButton)

        self.retranslateUi(Form)
        self.PlotArea.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.PlotArea.setTabText(self.PlotArea.indexOf(self.gust_tab), QtGui.QApplication.translate("Form", "Gust Diagram", None, QtGui.QApplication.UnicodeUTF8))
        self.PlotArea.setTabText(self.PlotArea.indexOf(self.manoeuvre_tab), QtGui.QApplication.translate("Form", "Manoeuvre Diagram", None, QtGui.QApplication.UnicodeUTF8))
        self.PlotArea.setTabText(self.PlotArea.indexOf(self.combined_tab), QtGui.QApplication.translate("Form", "Combined Diagram", None, QtGui.QApplication.UnicodeUTF8))
        self.InputData_label.setText(QtGui.QApplication.translate("Form", "Input Data", None, QtGui.QApplication.UnicodeUTF8))
        self.MTOW_label.setText(QtGui.QApplication.translate("Form", "MTOW =", None, QtGui.QApplication.UnicodeUTF8))
        self.MTOW_unitlabel.setText(QtGui.QApplication.translate("Form", "Kg", None, QtGui.QApplication.UnicodeUTF8))
        self.MLW_label.setText(QtGui.QApplication.translate("Form", "MLW=", None, QtGui.QApplication.UnicodeUTF8))
        self.MLW_unitlabel.setText(QtGui.QApplication.translate("Form", "Kg", None, QtGui.QApplication.UnicodeUTF8))
        self.MZFW_label.setText(QtGui.QApplication.translate("Form", "MZFW=", None, QtGui.QApplication.UnicodeUTF8))
        self.MZFW_unitlabel.setText(QtGui.QApplication.translate("Form", "Kg", None, QtGui.QApplication.UnicodeUTF8))
        self.W0_label.setText(QtGui.QApplication.translate("Form", "W0=", None, QtGui.QApplication.UnicodeUTF8))
        self.W0_unitlabel.setText(QtGui.QApplication.translate("Form", "Kg", None, QtGui.QApplication.UnicodeUTF8))
        self.a3D_label.setText(QtGui.QApplication.translate("Form", "a3D =", None, QtGui.QApplication.UnicodeUTF8))
        self.clmax_label.setText(QtGui.QApplication.translate("Form", "Cl_max =", None, QtGui.QApplication.UnicodeUTF8))
        self.clmax_flap_label.setText(QtGui.QApplication.translate("Form", "Cl_max_flap =", None, QtGui.QApplication.UnicodeUTF8))
        self.clmin_label.setText(QtGui.QApplication.translate("Form", "Cl_min =", None, QtGui.QApplication.UnicodeUTF8))
        self.W_label.setText(QtGui.QApplication.translate("Form", "W =", None, QtGui.QApplication.UnicodeUTF8))
        self.W_unitlabel.setText(QtGui.QApplication.translate("Form", "Kg", None, QtGui.QApplication.UnicodeUTF8))
        self.Altura_Button.setText(QtGui.QApplication.translate("Form", "Seleccionar Altura", None, QtGui.QApplication.UnicodeUTF8))
        self.h_label.setText(QtGui.QApplication.translate("Form", "h =", None, QtGui.QApplication.UnicodeUTF8))
        self.h_unitlabel.setText(QtGui.QApplication.translate("Form", "m", None, QtGui.QApplication.UnicodeUTF8))
        self.den_label.setText(QtGui.QApplication.translate("Form", "den =", None, QtGui.QApplication.UnicodeUTF8))
        self.den_unitlabel.setText(QtGui.QApplication.translate("Form", "kg/m³", None, QtGui.QApplication.UnicodeUTF8))
        self.CAM_label.setText(QtGui.QApplication.translate("Form", "CAM =", None, QtGui.QApplication.UnicodeUTF8))
        self.CAM_unitlabel.setText(QtGui.QApplication.translate("Form", "m", None, QtGui.QApplication.UnicodeUTF8))
        self.sw_label.setText(QtGui.QApplication.translate("Form", "sw =", None, QtGui.QApplication.UnicodeUTF8))
        self.sw_unitlabel.setText(QtGui.QApplication.translate("Form", "m^2", None, QtGui.QApplication.UnicodeUTF8))
        self.Zmo_label.setText(QtGui.QApplication.translate("Form", "Zmo=", None, QtGui.QApplication.UnicodeUTF8))
        self.Zmo_unitlabel.setText(QtGui.QApplication.translate("Form", "m", None, QtGui.QApplication.UnicodeUTF8))
        self.Vc_label.setText(QtGui.QApplication.translate("Form", "Vc =", None, QtGui.QApplication.UnicodeUTF8))
        self.Vc_unitlabel.setText(QtGui.QApplication.translate("Form", "m/s", None, QtGui.QApplication.UnicodeUTF8))
        self.grafiacar_pushButton.setText(QtGui.QApplication.translate("Form", "Graficar", None, QtGui.QApplication.UnicodeUTF8))
        self.units_groupBox.setTitle(QtGui.QApplication.translate("Form", "Units", None, QtGui.QApplication.UnicodeUTF8))
        self.SI_radioButton.setText(QtGui.QApplication.translate("Form", "SI", None, QtGui.QApplication.UnicodeUTF8))
        self.IM_radioButton.setText(QtGui.QApplication.translate("Form", "IM", None, QtGui.QApplication.UnicodeUTF8))

