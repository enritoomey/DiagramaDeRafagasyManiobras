# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'layout_DiagramaDeRafagasyManiobras.ui'
#
# Created: Sun Sep 18 18:03:53 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.4
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(611, 390)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
        self.horizontalLayout_3.setContentsMargins(-1, -1, 11, -1)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.left_Layout = QtGui.QVBoxLayout()
        self.left_Layout.setObjectName("left_Layout")
        self.InputData_label = QtGui.QLabel(Dialog)
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setWeight(75)
        font.setUnderline(True)
        font.setBold(True)
        self.InputData_label.setFont(font)
        self.InputData_label.setObjectName("InputData_label")
        self.left_Layout.addWidget(self.InputData_label)
        self.geometry_group_layout = QtGui.QVBoxLayout()
        self.geometry_group_layout.setObjectName("geometry_group_layout")
        self.CAM_layout = QtGui.QHBoxLayout()
        self.CAM_layout.setObjectName("CAM_layout")
        self.CAM_label = QtGui.QLabel(Dialog)
        self.CAM_label.setObjectName("CAM_label")
        self.CAM_layout.addWidget(self.CAM_label)
        self.CAM_lineEdit = QtGui.QLineEdit(Dialog)
        self.CAM_lineEdit.setObjectName("CAM_lineEdit")
        self.CAM_layout.addWidget(self.CAM_lineEdit)
        self.CAM_unitlabel = QtGui.QLabel(Dialog)
        self.CAM_unitlabel.setObjectName("CAM_unitlabel")
        self.CAM_layout.addWidget(self.CAM_unitlabel)
        self.geometry_group_layout.addLayout(self.CAM_layout)
        self.sw_layout = QtGui.QHBoxLayout()
        self.sw_layout.setObjectName("sw_layout")
        self.sw_label = QtGui.QLabel(Dialog)
        self.sw_label.setObjectName("sw_label")
        self.sw_layout.addWidget(self.sw_label)
        self.sw_lineEdit = QtGui.QLineEdit(Dialog)
        self.sw_lineEdit.setObjectName("sw_lineEdit")
        self.sw_layout.addWidget(self.sw_lineEdit)
        self.sw_unitlabel = QtGui.QLabel(Dialog)
        self.sw_unitlabel.setObjectName("sw_unitlabel")
        self.sw_layout.addWidget(self.sw_unitlabel)
        self.geometry_group_layout.addLayout(self.sw_layout)
        self.left_Layout.addLayout(self.geometry_group_layout)
        spacerItem = QtGui.QSpacerItem(20, 5, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.left_Layout.addItem(spacerItem)
        self.W_group_layout = QtGui.QVBoxLayout()
        self.W_group_layout.setObjectName("W_group_layout")
        self.MTOW_layout = QtGui.QHBoxLayout()
        self.MTOW_layout.setObjectName("MTOW_layout")
        self.MTOW_label = QtGui.QLabel(Dialog)
        self.MTOW_label.setObjectName("MTOW_label")
        self.MTOW_layout.addWidget(self.MTOW_label)
        self.MTOW_lineEdit = QtGui.QLineEdit(Dialog)
        self.MTOW_lineEdit.setObjectName("MTOW_lineEdit")
        self.MTOW_layout.addWidget(self.MTOW_lineEdit)
        self.MTOW_unitlabel = QtGui.QLabel(Dialog)
        self.MTOW_unitlabel.setObjectName("MTOW_unitlabel")
        self.MTOW_layout.addWidget(self.MTOW_unitlabel)
        self.W_group_layout.addLayout(self.MTOW_layout)
        self.MLW_layout = QtGui.QHBoxLayout()
        self.MLW_layout.setObjectName("MLW_layout")
        self.MLW_label = QtGui.QLabel(Dialog)
        self.MLW_label.setObjectName("MLW_label")
        self.MLW_layout.addWidget(self.MLW_label)
        self.MLW_lineEdit = QtGui.QLineEdit(Dialog)
        self.MLW_lineEdit.setObjectName("MLW_lineEdit")
        self.MLW_layout.addWidget(self.MLW_lineEdit)
        self.MLW_unitlabel = QtGui.QLabel(Dialog)
        self.MLW_unitlabel.setObjectName("MLW_unitlabel")
        self.MLW_layout.addWidget(self.MLW_unitlabel)
        self.W_group_layout.addLayout(self.MLW_layout)
        self.MZFW_layout = QtGui.QHBoxLayout()
        self.MZFW_layout.setObjectName("MZFW_layout")
        self.MZFW_label = QtGui.QLabel(Dialog)
        self.MZFW_label.setObjectName("MZFW_label")
        self.MZFW_layout.addWidget(self.MZFW_label)
        self.MZFW_lineEdit = QtGui.QLineEdit(Dialog)
        self.MZFW_lineEdit.setObjectName("MZFW_lineEdit")
        self.MZFW_layout.addWidget(self.MZFW_lineEdit)
        self.MZFW_unitlabel = QtGui.QLabel(Dialog)
        self.MZFW_unitlabel.setObjectName("MZFW_unitlabel")
        self.MZFW_layout.addWidget(self.MZFW_unitlabel)
        self.W_group_layout.addLayout(self.MZFW_layout)
        self.W0_layout = QtGui.QHBoxLayout()
        self.W0_layout.setSizeConstraint(QtGui.QLayout.SetMaximumSize)
        self.W0_layout.setObjectName("W0_layout")
        self.W0_label = QtGui.QLabel(Dialog)
        self.W0_label.setObjectName("W0_label")
        self.W0_layout.addWidget(self.W0_label)
        self.W0_lineEdit = QtGui.QLineEdit(Dialog)
        self.W0_lineEdit.setObjectName("W0_lineEdit")
        self.W0_layout.addWidget(self.W0_lineEdit)
        self.W0_unitlabel = QtGui.QLabel(Dialog)
        self.W0_unitlabel.setObjectName("W0_unitlabel")
        self.W0_layout.addWidget(self.W0_unitlabel)
        self.W_group_layout.addLayout(self.W0_layout)
        self.left_Layout.addLayout(self.W_group_layout)
        spacerItem1 = QtGui.QSpacerItem(20, 5, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.left_Layout.addItem(spacerItem1)
        self.Cl_group_layout = QtGui.QVBoxLayout()
        self.Cl_group_layout.setObjectName("Cl_group_layout")
        self.a3D_layout = QtGui.QHBoxLayout()
        self.a3D_layout.setObjectName("a3D_layout")
        self.a3D_label = QtGui.QLabel(Dialog)
        self.a3D_label.setObjectName("a3D_label")
        self.a3D_layout.addWidget(self.a3D_label)
        self.a3D_lineEdit = QtGui.QLineEdit(Dialog)
        self.a3D_lineEdit.setObjectName("a3D_lineEdit")
        self.a3D_layout.addWidget(self.a3D_lineEdit)
        self.Cl_group_layout.addLayout(self.a3D_layout)
        self.clmax_layout = QtGui.QHBoxLayout()
        self.clmax_layout.setObjectName("clmax_layout")
        self.clmax_label = QtGui.QLabel(Dialog)
        self.clmax_label.setObjectName("clmax_label")
        self.clmax_layout.addWidget(self.clmax_label)
        self.clmax_lineEdit = QtGui.QLineEdit(Dialog)
        self.clmax_lineEdit.setObjectName("clmax_lineEdit")
        self.clmax_layout.addWidget(self.clmax_lineEdit)
        self.Cl_group_layout.addLayout(self.clmax_layout)
        self.clmax_flap_layout = QtGui.QHBoxLayout()
        self.clmax_flap_layout.setObjectName("clmax_flap_layout")
        self.clmax_flap_label = QtGui.QLabel(Dialog)
        self.clmax_flap_label.setObjectName("clmax_flap_label")
        self.clmax_flap_layout.addWidget(self.clmax_flap_label)
        self.clmax_flap_lineEdit = QtGui.QLineEdit(Dialog)
        self.clmax_flap_lineEdit.setObjectName("clmax_flap_lineEdit")
        self.clmax_flap_layout.addWidget(self.clmax_flap_lineEdit)
        self.Cl_group_layout.addLayout(self.clmax_flap_layout)
        self.clmin_layout = QtGui.QHBoxLayout()
        self.clmin_layout.setObjectName("clmin_layout")
        self.clmin_label = QtGui.QLabel(Dialog)
        self.clmin_label.setObjectName("clmin_label")
        self.clmin_layout.addWidget(self.clmin_label)
        self.clmin_lineEdit = QtGui.QLineEdit(Dialog)
        self.clmin_lineEdit.setObjectName("clmin_lineEdit")
        self.clmin_layout.addWidget(self.clmin_lineEdit)
        self.Cl_group_layout.addLayout(self.clmin_layout)
        self.left_Layout.addLayout(self.Cl_group_layout)
        self.horizontalLayout_3.addLayout(self.left_Layout)
        self.right_Layout = QtGui.QVBoxLayout()
        self.right_Layout.setObjectName("right_Layout")
        self.PlotArea = QtGui.QTabWidget(Dialog)
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
        self.right_Layout.addWidget(self.PlotArea)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.performace_group_layout = QtGui.QVBoxLayout()
        self.performace_group_layout.setObjectName("performace_group_layout")
        self.W_layout = QtGui.QHBoxLayout()
        self.W_layout.setObjectName("W_layout")
        self.W_label = QtGui.QLabel(Dialog)
        self.W_label.setObjectName("W_label")
        self.W_layout.addWidget(self.W_label)
        self.W_lineEdit = QtGui.QLineEdit(Dialog)
        self.W_lineEdit.setObjectName("W_lineEdit")
        self.W_layout.addWidget(self.W_lineEdit)
        self.W_unitlabel = QtGui.QLabel(Dialog)
        self.W_unitlabel.setObjectName("W_unitlabel")
        self.W_layout.addWidget(self.W_unitlabel)
        self.performace_group_layout.addLayout(self.W_layout)
        self.Zmo_layout = QtGui.QHBoxLayout()
        self.Zmo_layout.setObjectName("Zmo_layout")
        self.Zmo_label = QtGui.QLabel(Dialog)
        self.Zmo_label.setObjectName("Zmo_label")
        self.Zmo_layout.addWidget(self.Zmo_label)
        self.Zmo_lineEdit = QtGui.QLineEdit(Dialog)
        self.Zmo_lineEdit.setObjectName("Zmo_lineEdit")
        self.Zmo_layout.addWidget(self.Zmo_lineEdit)
        self.Zmo_unitlabel = QtGui.QLabel(Dialog)
        self.Zmo_unitlabel.setObjectName("Zmo_unitlabel")
        self.Zmo_layout.addWidget(self.Zmo_unitlabel)
        self.performace_group_layout.addLayout(self.Zmo_layout)
        self.Vc_layout = QtGui.QHBoxLayout()
        self.Vc_layout.setObjectName("Vc_layout")
        self.Vc_label = QtGui.QLabel(Dialog)
        self.Vc_label.setObjectName("Vc_label")
        self.Vc_layout.addWidget(self.Vc_label)
        self.Vc_lineEdit = QtGui.QLineEdit(Dialog)
        self.Vc_lineEdit.setObjectName("Vc_lineEdit")
        self.Vc_layout.addWidget(self.Vc_lineEdit)
        self.Vc_unitlabel = QtGui.QLabel(Dialog)
        self.Vc_unitlabel.setObjectName("Vc_unitlabel")
        self.Vc_layout.addWidget(self.Vc_unitlabel)
        self.performace_group_layout.addLayout(self.Vc_layout)
        self.horizontalLayout.addLayout(self.performace_group_layout)
        self.Altura_Layout = QtGui.QVBoxLayout()
        self.Altura_Layout.setObjectName("Altura_Layout")
        self.Altura_Button = QtGui.QPushButton(Dialog)
        self.Altura_Button.setObjectName("Altura_Button")
        self.Altura_Layout.addWidget(self.Altura_Button)
        self.h_layout = QtGui.QHBoxLayout()
        self.h_layout.setObjectName("h_layout")
        self.h_label = QtGui.QLabel(Dialog)
        self.h_label.setObjectName("h_label")
        self.h_layout.addWidget(self.h_label)
        self.h_lineEdit = QtGui.QLineEdit(Dialog)
        self.h_lineEdit.setReadOnly(False)
        self.h_lineEdit.setObjectName("h_lineEdit")
        self.h_layout.addWidget(self.h_lineEdit)
        self.h_unitlabel = QtGui.QLabel(Dialog)
        self.h_unitlabel.setObjectName("h_unitlabel")
        self.h_layout.addWidget(self.h_unitlabel)
        self.Altura_Layout.addLayout(self.h_layout)
        self.den_layout = QtGui.QHBoxLayout()
        self.den_layout.setObjectName("den_layout")
        self.den_label = QtGui.QLabel(Dialog)
        self.den_label.setObjectName("den_label")
        self.den_layout.addWidget(self.den_label)
        self.den_lineEdit = QtGui.QLineEdit(Dialog)
        self.den_lineEdit.setReadOnly(False)
        self.den_lineEdit.setObjectName("den_lineEdit")
        self.den_layout.addWidget(self.den_lineEdit)
        self.den_unitlabel = QtGui.QLabel(Dialog)
        self.den_unitlabel.setObjectName("den_unitlabel")
        self.den_layout.addWidget(self.den_unitlabel)
        self.Altura_Layout.addLayout(self.den_layout)
        self.horizontalLayout.addLayout(self.Altura_Layout)
        self.units_Layout = QtGui.QVBoxLayout()
        self.units_Layout.setObjectName("units_Layout")
        self.units_groupBox = QtGui.QGroupBox(Dialog)
        self.units_groupBox.setCheckable(False)
        self.units_groupBox.setObjectName("units_groupBox")
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.units_groupBox)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.IM_radioButton = QtGui.QRadioButton(self.units_groupBox)
        self.IM_radioButton.setObjectName("IM_radioButton")
        self.horizontalLayout_2.addWidget(self.IM_radioButton)
        self.SI_radioButton = QtGui.QRadioButton(self.units_groupBox)
        self.SI_radioButton.setEnabled(True)
        self.SI_radioButton.setChecked(True)
        self.SI_radioButton.setObjectName("SI_radioButton")
        self.horizontalLayout_2.addWidget(self.SI_radioButton)
        self.units_Layout.addWidget(self.units_groupBox)
        self.grafiacar_pushButton = QtGui.QPushButton(Dialog)
        self.grafiacar_pushButton.setObjectName("grafiacar_pushButton")
        self.units_Layout.addWidget(self.grafiacar_pushButton)
        self.horizontalLayout.addLayout(self.units_Layout)
        self.right_Layout.addLayout(self.horizontalLayout)
        self.horizontalLayout_3.addLayout(self.right_Layout)
        self.verticalLayout.addLayout(self.horizontalLayout_3)

        self.retranslateUi(Dialog)
        self.PlotArea.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))
        self.InputData_label.setText(QtGui.QApplication.translate("Dialog", "Input Data", None, QtGui.QApplication.UnicodeUTF8))
        self.CAM_label.setText(QtGui.QApplication.translate("Dialog", "CAM =", None, QtGui.QApplication.UnicodeUTF8))
        self.CAM_unitlabel.setText(QtGui.QApplication.translate("Dialog", "m", None, QtGui.QApplication.UnicodeUTF8))
        self.sw_label.setText(QtGui.QApplication.translate("Dialog", "sw =", None, QtGui.QApplication.UnicodeUTF8))
        self.sw_unitlabel.setText(QtGui.QApplication.translate("Dialog", "m^2", None, QtGui.QApplication.UnicodeUTF8))
        self.MTOW_label.setText(QtGui.QApplication.translate("Dialog", "MTOW =", None, QtGui.QApplication.UnicodeUTF8))
        self.MTOW_unitlabel.setText(QtGui.QApplication.translate("Dialog", "Kg", None, QtGui.QApplication.UnicodeUTF8))
        self.MLW_label.setText(QtGui.QApplication.translate("Dialog", "MLW=", None, QtGui.QApplication.UnicodeUTF8))
        self.MLW_unitlabel.setText(QtGui.QApplication.translate("Dialog", "Kg", None, QtGui.QApplication.UnicodeUTF8))
        self.MZFW_label.setText(QtGui.QApplication.translate("Dialog", "MZFW=", None, QtGui.QApplication.UnicodeUTF8))
        self.MZFW_unitlabel.setText(QtGui.QApplication.translate("Dialog", "Kg", None, QtGui.QApplication.UnicodeUTF8))
        self.W0_label.setText(QtGui.QApplication.translate("Dialog", "W0=", None, QtGui.QApplication.UnicodeUTF8))
        self.W0_unitlabel.setText(QtGui.QApplication.translate("Dialog", "Kg", None, QtGui.QApplication.UnicodeUTF8))
        self.a3D_label.setText(QtGui.QApplication.translate("Dialog", "a3D =", None, QtGui.QApplication.UnicodeUTF8))
        self.clmax_label.setText(QtGui.QApplication.translate("Dialog", "Cl_max =", None, QtGui.QApplication.UnicodeUTF8))
        self.clmax_flap_label.setText(QtGui.QApplication.translate("Dialog", "Cl_max_flap =", None, QtGui.QApplication.UnicodeUTF8))
        self.clmin_label.setText(QtGui.QApplication.translate("Dialog", "Cl_min =", None, QtGui.QApplication.UnicodeUTF8))
        self.PlotArea.setTabText(self.PlotArea.indexOf(self.gust_tab), QtGui.QApplication.translate("Dialog", "Gust Diagram", None, QtGui.QApplication.UnicodeUTF8))
        self.PlotArea.setTabText(self.PlotArea.indexOf(self.manoeuvre_tab), QtGui.QApplication.translate("Dialog", "Manoeuvre Diagram", None, QtGui.QApplication.UnicodeUTF8))
        self.PlotArea.setTabText(self.PlotArea.indexOf(self.combined_tab), QtGui.QApplication.translate("Dialog", "Combined Diagram", None, QtGui.QApplication.UnicodeUTF8))
        self.W_label.setText(QtGui.QApplication.translate("Dialog", "W =", None, QtGui.QApplication.UnicodeUTF8))
        self.W_unitlabel.setText(QtGui.QApplication.translate("Dialog", "Kg", None, QtGui.QApplication.UnicodeUTF8))
        self.Zmo_label.setText(QtGui.QApplication.translate("Dialog", "Zmo=", None, QtGui.QApplication.UnicodeUTF8))
        self.Zmo_unitlabel.setText(QtGui.QApplication.translate("Dialog", "m", None, QtGui.QApplication.UnicodeUTF8))
        self.Vc_label.setText(QtGui.QApplication.translate("Dialog", "Vc =", None, QtGui.QApplication.UnicodeUTF8))
        self.Vc_unitlabel.setText(QtGui.QApplication.translate("Dialog", "m/s", None, QtGui.QApplication.UnicodeUTF8))
        self.Altura_Button.setText(QtGui.QApplication.translate("Dialog", "Seleccionar Altura", None, QtGui.QApplication.UnicodeUTF8))
        self.h_label.setText(QtGui.QApplication.translate("Dialog", "h =", None, QtGui.QApplication.UnicodeUTF8))
        self.h_unitlabel.setText(QtGui.QApplication.translate("Dialog", "m", None, QtGui.QApplication.UnicodeUTF8))
        self.den_label.setText(QtGui.QApplication.translate("Dialog", "den =", None, QtGui.QApplication.UnicodeUTF8))
        self.den_unitlabel.setText(QtGui.QApplication.translate("Dialog", "kg/m³", None, QtGui.QApplication.UnicodeUTF8))
        self.units_groupBox.setTitle(QtGui.QApplication.translate("Dialog", "Units", None, QtGui.QApplication.UnicodeUTF8))
        self.IM_radioButton.setText(QtGui.QApplication.translate("Dialog", "IM", None, QtGui.QApplication.UnicodeUTF8))
        self.SI_radioButton.setText(QtGui.QApplication.translate("Dialog", "SI", None, QtGui.QApplication.UnicodeUTF8))
        self.grafiacar_pushButton.setText(QtGui.QApplication.translate("Dialog", "Graficar", None, QtGui.QApplication.UnicodeUTF8))

