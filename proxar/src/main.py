'''
For proxar design tool
Created on Jan 24, 2013

@author: jifeng
'''
#!/usr/bin/python -d

import os
import time
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QCompleter,  QStringListModel
from Bio import SeqIO
import primer_design

from gui import Ui_Form
import advance
import Proxar_GUI3

import configure

from xlwt import Workbook
class MyForm(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.ui.progressBar.hide()
        self.ui.pushButton_2.hide()
        self.ui.label_2.hide()
        QtCore.QObject.connect(self.ui.pushButton, QtCore.SIGNAL("clicked()"), self.filebrower) #for input file
        QtCore.QObject.connect(self.ui.buttonBox, QtCore.SIGNAL(("accepted()")), self.submit) #submit
        QtCore.QObject.connect(self.ui.buttonBox, QtCore.SIGNAL(("rejected()")), self.reset) #submit
        QtCore.QObject.connect(self.ui.pushButton_6, QtCore.SIGNAL("clicked()"), self.advance) #advanced setting
        QtCore.QObject.connect(self.ui.pushButton_4, QtCore.SIGNAL("clicked()"), self.add_taxid1) #add taxid to include
        QtCore.QObject.connect(self.ui.pushButton_5, QtCore.SIGNAL("clicked()"), self.add_taxid2) #add taxid to exclude
        QtCore.QObject.connect(self.ui.pushButton_2, QtCore.SIGNAL("clicked()"), self.savefilebrower_all) #for input file
        self.ui.lineEdit_3 = QtGui.QLineEdit(self.ui.groupBox_2)
        self.ui.lineEdit_3.setGeometry(QtCore.QRect(70, 20, 270, 21))
        self.ui.lineEdit_3.setObjectName("lineEdit_3")
        self.ui.lineEdit_3.setText(QtGui.QApplication.translate("Form", "", None, QtGui.QApplication.UnicodeUTF8))
        completer = QCompleter()
        self.ui.lineEdit_3.setCompleter(completer)
        model = QStringListModel()
        completer.setModel(model)
        completer.setModelSorting(QCompleter.CaseInsensitivelySortedModel)
        Proxar_GUI3.get_data(model)
        self.ui.lineEdit_13 = QtGui.QLineEdit(self.ui.groupBox_4)
        self.ui.lineEdit_13.setGeometry(QtCore.QRect(70, 20, 270, 21))
        self.ui.lineEdit_13.setObjectName("lineEdit_13")
        self.ui.lineEdit_13.setText(QtGui.QApplication.translate("Form", "", None, QtGui.QApplication.UnicodeUTF8))
        completer = QCompleter()
        self.ui.lineEdit_13.setCompleter(completer)
        model = QStringListModel()
        completer.setModel(model)
        completer.setModelSorting(QCompleter.CaseInsensitivelySortedModel)
        Proxar_GUI3.get_data(model)
        #self.ui.lineEdit_5.setText('Mycobacterium avium (taxid:1764) OR Mycobacterium bovis (taxid:1765)')# include list
        #self.ui.lineEdit_15.setText('homo sapiens (taxid:9606)')
        self.ui.lineEdit.setText('short_example206.txt')

    def filebrower(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file', './')
        self.ui.lineEdit.setText(fname)

    def add_taxid1(self):
        mid=" OR "
        if  self.ui.lineEdit_5.text()=="":
            self.ui.lineEdit_5.setText(self.ui.lineEdit_3.text())
        else:
            self.ui.lineEdit_5.setText(self.ui.lineEdit_5.text()+mid+self.ui.lineEdit_3.text())
        self.ui.lineEdit_3.setText("")

    def add_taxid2(self):
        mid=" OR "
        if  self.ui.lineEdit_15.text()=="":
            self.ui.lineEdit_15.setText(self.ui.lineEdit_13.text())
        else:
            self.ui.lineEdit_15.setText(self.ui.lineEdit_15.text()+mid+self.ui.lineEdit_13.text())
        self.ui.lineEdit_13.setText("")

    def savefilebrower_all(self):
        fName = QtGui.QFileDialog.getSaveFileName(self, "Save as xls file", "Save as new file", self.tr("Excel Files (*.xls)"))
        if fName.isEmpty() == False:
            self.write_xls(fName)
            QtGui.QMessageBox.about(self, "Saved", "%s is generated!" % (fName))
            Proxar_GUI3.delete_file(configure.tempname)
            self.close()

    def write_xls(self,filename):
        book = Workbook()
        sheet1 = book.add_sheet('Whole PT/PX Combinations',cell_overwrite_ok=True)
        sheet2 = book.add_sheet('Desired thermo',cell_overwrite_ok=True)
        sheet3 = book.add_sheet('All submitted templates',cell_overwrite_ok=True)

        row11 = sheet1.row(0)
        row11.write(0,'PT_name')
        row11.write(1,'PT_sequence')
        row11.write(2,'PX_name')
        row11.write(3,'PX_sequence')

        i=1
        for record in SeqIO.parse(seq_file, "fasta") :
            if i%2==1:
                row11 = sheet1.row(i)
                row11.write(0, str(record.id))
                row11.write(1, str(record.seq))
            else:
                row11.write(2, str(record.id))
                row11.write(3, str(record.seq))
                i=i+1



        row21 = sheet2.row(0)
        row21.write(0,'name')
        row21.write(1,'sequence')
        row21.write(2,'bayes_class')
        row21.write(3,'pwm_class')
        row21.write(4,'p90 score')
        row21.write(5,'diff score')
        row21.write(6,'tri-temp Tm')
        row21.write(7,'temp-temp Tm')
        row21.write(8,'bonds')


        row31 = sheet3.row(0)
        row31.write(0,'name')
        row31.write(1,'sequence')
        row31.write(2,'bayes_class')
        row31.write(3,'pwm_class')
        row31.write(4,'p90 score')
        row31.write(5,'diff score')
        row31.write(6,'tri-temp Tm')
        row31.write(7,'temp-temp Tm')
        row31.write(8,'bonds')

        book.save(filename)

    def newgui5(self):
        self.ui.progressBar.show()
        self.ui.progressBar.reset()
        self.ui.progressBar.setMinimum(0)
        self.ui.progressBar.setMaximum(10)
        while configure.progress_simbol>9:
            time.sleep(1)
            self.ui.progressBar.setValue(configure.progress_simbol)
            self.ui.progressBar.update()
            if configure.progress_simbol>9:
                break
        self.ui.pushButton_2.show()
        self.ui.progressBar.hide()

    def submit(self):
        #Proxar_GUI3.main_page(target_file=self.ui.lineEdit.text(),include_line=self.ui.lineEdit_5.text(), exclude_line=self.ui.lineEdit_15.text(), hairpin_file='hairpin.txt',extra_base_file='extra_base.txt', templates_file='good_template.txt',exclude_GAGTC=self.ui.checkBox.isChecked(),two_direction=self.ui.checkBox_2.isChecked(),PT_foot_min=self.ui.lineEdit_6.text(),PT_foot_max=self.ui.lineEdit_9.text(),PT_Tm_above=self.ui.lineEdit_10.text(),PX_foot_min=self.ui.lineEdit_8.text(),PX_foot_max=self.ui.lineEdit_11.text() ,PX_Tm_above=self.ui.lineEdit_12.text(),max_gap = self.ui.lineEdit_7.text(),pt_longer_num=self.ui.lineEdit_14.text(),px_longer_num=self.ui.lineEdit_16.text())
        whole_PT2='good_template.txt_wholePT.txttemp_pass_bind.txt'
        whole_PX1='px_foot_hair_file2.txttemp_pass_bind.txt'
        exclude_line='homo sapiens (taxid:9606)'
        pt_longer_num=19
        px_longer_num=15
        #self.ui.progressBar.show()
        #self.ui.progressBar.reset()
        self.ui.buttonBox.hide()
        self.ui.pushButton_6.hide()
        print "ok"
        time.sleep(5)
        self.ui.pushButton_2.show()                                                                
        primer_design.main_page(self.ui.lineEdit.text(),self.ui.lineEdit_5.text(),self.ui.lineEdit_15.text(),\
                                                            'hairpin.txt','extra_base.txt','good_template.txt',self.ui.checkBox.isChecked(),\
                                                            self.ui.checkBox_2.isChecked(),\
                                                            self.ui.lineEdit_6.text(),self.ui.lineEdit_9.text(),configure.pt_foot_tm,\
                                                            self.ui.lineEdit_8.text(),\
                                                            self.ui.lineEdit_11.text(),configure.px_foot_tm,self.ui.lineEdit_7.text(),\
                                                            self.ui.lineEdit_14.text(),\
                                                            self.ui.lineEdit_16.text(),\
                                                            configure.pt_foot_conc,\
                                                            configure.px_foot_conc, configure.pt_px_foot_match_max,\
                                                            configure.hair_tar_max,\
                                                            configure.overhang_pt_tar_max,\
                                                            configure.overhang_px_tar_max)
       
          #main_page(target_file,include_line, exclude_list, 
          #hairpin_file,extra_base_file, templates_file,exclude_GAGTC,
          #two_direction,
          #PT_foot_min,PT_foot_max,PT_Tm_above,
          #PX_foot_min,
          #PX_foot_max,PX_Tm_above,max_gap,
          ##PT_discard_longer_tahn_default 19)
          #PX_discard_longer_tahn_default 19
          #...)
    def reset(self):
        self.ui.lineEdit.setText('J:\PROXAR\2-23-2012\short_example.txt') #file
        self.ui.lineEdit_6.setText('23')# PT foot min
        self.ui.lineEdit_9.setText('29') # PT foot max
        self.ui.lineEdit_10.setText('4') # PX mismatch
        self.ui.lineEdit_4.setText('4') # PT mismatch
        self.ui.lineEdit_8.setText('21') # PX foot min
        self.ui.lineEdit_11.setText('23') # PX foot max
        self.ui.lineEdit_12.setText('0') # PT above Tm
        self.ui.lineEdit_7.setText('2') # gap
        self.ui.lineEdit_5.setText('Mycobacterium avium (taxid:1764) OR Mycobacterium bovis (taxid:1765)')# include list
        self.ui.lineEdit_15.setText('')# exclude list
        self.ui.lineEdit_14.setText('19')
        self.ui.lineEdit_16.setText('15')
        self.checkBox.setChecked(True)
        self.checkBox_4.setChecked(True)
        self.checkBox_3.setChecked(True)
        self.checkBox_2.setChecked(True)

    def advance(self):
        self.dialog_advance = QtGui.QDialog()
        self.dialog_advance.ui = advance.Ui_Dialog()
        self.dialog_advance.ui.setupUi(self.dialog_advance)
        self.dialog_advance.connect(self.dialog_advance.ui.buttonBox, QtCore.SIGNAL("accepted()"), self.saveAdvancedSettings)
        self.dialog_advance.exec_()


    def saveAdvancedSettings(self):
        #configure.progress_simbol=
        configure.pt_foot_conc=float(self.dialog_advance.ui.lineEdit_2.text())
        configure.px_foot_conc=float(self.dialog_advance.ui.lineEdit_4.text())
        configure.pt_foot_tm=float(self.dialog_advance.ui.lineEdit.text())
        configure.px_foot_tm=float(self.dialog_advance.ui.lineEdit_3.text())
        configure.pt_px_foot_match_max=int(self.dialog_advance.ui.lineEdit_5.text())
        configure.hair_tar_max=int(self.dialog_advance.ui.lineEdit_6.text())
        configure.whole_PT_tm=float(self.dialog_advance.ui.lineEdit_8.text())
        configure.whole_PX_tm=float(self.dialog_advance.ui.lineEdit_9.text())
        configure.overhang_pt_tar_max=int(self.dialog_advance.ui.lineEdit_7.text())
        configure.overhang_px_tar_max=int(self.dialog_advance.ui.lineEdit_10.text())


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Form = QtGui.QWidget()
    if os.system('perl -v >perl.version'):
                QtGui.QMessageBox.critical(Form,"Error", "you haven't properly installed perl yet!")
                exit()
    Newui = MyForm()
    Newui.show()
    sys.exit(app.exec_())
