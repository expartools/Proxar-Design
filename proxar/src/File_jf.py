'''
Created on Jan 25, 2013
this is a package for file and checking stuff
@author: jifeng
'''
import os
from Bio import SeqIO
import re
from PyQt4 import QtGui
def check_file(target_file, hairpin_file, extra_base_file,templates_file):
    if (os.path.isfile(target_file) and os.path.isfile(hairpin_file) and os.path.isfile(extra_base_file) and os.path.isfile(templates_file)):
        print os.path.isfile(target_file),os.path.isfile(hairpin_file),os.path.isfile(extra_base_file),os.path.isfile(templates_file)
        print "all necessary file exist!"
        three_file_check_symbol(target_file,hairpin_file,extra_base_file)
        input_content=open(target_file).read()
        seq_num=[i - 1 for i in range(len(input_content)) if input_content.startswith('>', i - 1)]
        if len(seq_num)!=1:
            print "input should be one and only one sequence in fasta format!"
            return 0
        else:
            for cur_record in SeqIO.parse(open(target_file), "fasta"):
                input_seq=str(cur_record.seq).upper()
        return input_seq
    else:
        QtGui.QMessageBox.critical(None, 'Error', 'critical file is missing!') 
        exit()

def change_symbol(filename):
    p = re.compile( '-')
    (new_file_content,count)=p.subn( '~', open(filename).read())
    f=open(filename,'w')
    if count >0 :
        QtGui.QMessageBox.warning(None,"Warning", filename+" has "+str(count)+'\'-\' in it, chaanged to \'~\'')
        print >>f, new_file_content
    else:
        print >>f, new_file_content
    f.close()

def three_file_check_symbol(target_file,hairpin_file,extra_base_file):
    change_symbol(target_file)
    change_symbol(hairpin_file)
    change_symbol(extra_base_file)