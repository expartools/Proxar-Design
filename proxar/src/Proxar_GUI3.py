import configure

import os
from string import maketrans,join
import alignment
import re
from PyQt4 import QtGui
import shutil
import File_jf
from Bio import SeqIO

def revcomp(dna):
    comp = dna.translate(maketrans("AGCTagct", "TCGAtcga"))
    lcomp = list(comp)
    lcomp.reverse()
    return  join(lcomp,  "")

def get_data(model):
    f=open('taxid9_short.txt').read().split('\n')
    model.setStringList(f)
    #model.setStringList(['above','below','citrus'])

def main_page(target_file,include_line, exclude_line, hairpin_file,extra_base_file, templates_file,\
              exclude_GAGTC,two_direction,PT_foot_min,PT_foot_max,PT_Tm_above,PX_foot_min,PX_foot_max,PX_Tm_above,\
              max_gap,pt_longer_num,px_longer_num,pt_foot_conc,px_foot_conc,pt_px_foot_match_max,hair_tar_max,overhang_pt_tar_max,overhang_px_tar_max, target_mismatch):
    print "enter",target_file, hairpin_file, extra_base_file, templates_file
    if (os.path.isfile(target_file) and os.path.isfile(hairpin_file) and os.path.isfile(extra_base_file) and os.path.isfile(templates_file)):
        print os.path.isfile(target_file),os.path.isfile(hairpin_file),os.path.isfile(extra_base_file),os.path.isfile(templates_file)
        print "all necessary file exist!"
        File_jf.three_file_check_symbol(target_file,hairpin_file,extra_base_file)
        input_content=open(target_file).read()
        seq_num=[i - 1 for i in range(len(input_content)) if input_content.startswith('>', i - 1)]
        if len(seq_num)!=1:
            print "input should be one and only one sequence in fasta format!"
            return 0
        else:
            for cur_record in SeqIO.parse(open(target_file), "fasta"):
                input_seq=str(cur_record.seq).upper()
    else:
        #Form = QtGui.QWidget()
        #QtGui.QMessageBox.critical(Form,"Error", "you haven't properly installed perl yet!")
        print "error, necessary file is not there"
        exit()

        return 0
    print "what's wrong with this?working on obtain pt and px foot"
    configure.progress_simbol=1
    check_file_size(1)
    ptpx_foot_list=ptpx_foot(input_seq,PT_foot_min,PT_foot_max,PX_foot_min,PX_foot_max,e,'p')
    if two_direction :
        ptpx_foot_list=ptpx_foot_list|ptpx_foot(revcomp(input_seq),PT_foot_min,PT_foot_max,PX_foot_min,PX_foot_max,e,'n') ###################get the pt and px foot sequence
    f_p=open(configure.temper_print,'w')
    f_p2=open('del.txt','w')
    print >>f_p2, ptpx_foot_list
    print "done with working on obtain pt and px foot"
    write_ptpx_foot_file_check_foot(ptpx_foot_list,input_seq,pt_px_foot_match_max,configure.pt_foot_file1, configure.px_foot_file1,configure.tar_file1)###################get the pt and px foot sequence stored in files
    print "done with writing pt/px/target files"
    configure.progress_simbol=2
    #check_file_size(2)
    print PT_Tm_above,PX_Tm_above
    Tm_pfoot(configure.pt_foot_file1,configure.px_foot_file1,configure.tar_file1, PT_Tm_above,PX_Tm_above,configure.pt_foot_file2,configure.px_foot_file2,configure.tar_file2, configure.pt_bind_filename, configure.px_bind_filename, pt_foot_conc,px_foot_conc)
    print configure.pt_foot_file2,configure.px_foot_file2,configure.tar_file2,"are generated!"
    configure.progress_simbol=3
    check_file_size(3)
    if include_line!='': # if no include organism
        combine_check_include(include_line,target_mismatch,configure.tar_file2,configure.pt_foot_file3,configure.px_foot_file3,configure.tar_file3, configure.blast_in_file)
    else:
        shutil.copy2(configure.pt_foot_file2,configure.pt_foot_file3)
        shutil.copy2(configure.px_foot_file2,configure.px_foot_file3)
        shutil.copy2(configure.tar_file2,configure.tar_file3)
    print "include check over"
    configure.progress_simbol=4
    check_file_size(4)
    bind_hairpin(hairpin_file,extra_base_file, configure.pt_foot_hair_file1, configure.px_foot_hair_file1, configure.tar_hair_file1,hair_tar_max)
    print "bind_hairpin over"
    configure.progress_simbol=5
    check_file_size(5)
    self_Tm_check_eachother(configure.pt_foot_hair_file1,configure.px_foot_hair_file1,configure.tar_hair_file1,configure.whole_PTX_tm,configure.pt_foot_hair_file2, configure.px_foot_hair_file2,configure.tar_hair_file2, configure.ptx_melt_file)
    print "self Tm check over"
    configure.progress_simbol=6
    check_file_size(6)
    combine_PT(templates_file,configure.pt_foot_hair_file2,configure.tar_hair_file1,configure.whole_PT1,configure.overhang_pt_tar_max) # bind well-performing template from templates_file and pt foot with haipin from pt_foot_hair_file3
    print "combine_PT is over"
    print "*******",configure.whole_PT1
    configure.progress_simbol=6.5
    check_file_size(6.5)
    self_bind(configure.whole_PT1,configure.whole_PT_self_max,'-pt', configure.whole_PT2)
    print "self_bind PT1 is over"
    configure.progress_simbol=7
    check_file_size(7)
    self_bind(configure.px_foot_hair_file2,configure.whole_PX_self_max,'-px',configure.whole_PX1)
    configure.progress_simbol=8
    check_file_size(8)
    if exclude_line!='':
        combine_check_exclude(configure.whole_PT2, configure.whole_PX1,exclude_line,pt_longer_num,px_longer_num,configure.pt_blast,configure.px_blast,configure.finalptx)
    else:
        combine_pt_px_result(configure.whole_PT2, configure.whole_PX1,configure.finalptx)
    #configure.progress_simbol=9
    #check_file_size(9)
    #print "combine_PT is over"
    #print >>f_p, 'Tm_pt_foot_tar'
    #print >>f_p, Tm_pt_foot_tar
    #print >>f_p, 'tar_include'
    #print >>f_p, tar_include
    #print >>f_p, "Tm_p_foot_hair_tar"
    #print >>f_p, Tm_p_foot_hair_tar
    #print >>f_p, "overhang_bind_pt"
    #print >>f_p, overhang_bind_pt
    #print >>f_p, "whole_self_bind"
    #print >>f_p, whole_self_bind
    f_p.close()
    configure.progress_simbol=10
    check_file_size(10)

def check_file_size(num):
    if num==2:
        try:
           if os.path.getsize(configure.pt_foot_file1) and os.path.getsize(configure.px_foot_file1) and os.path.getsize(configure.tar_file1):
                return
           else:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()
        except:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()
    if num==3:
        try:
           if os.path.getsize(configure.pt_foot_file1) and os.path.getsize(configure.px_foot_file1) and os.path.getsize(configure.tar_file1):
                return
           else:
                QtGui.QMessageBox.critical(self,"Error", "No pt/px combination left after binding pt foot to the target piece, please change concentration or Temperature limit!")
                exit()
        except:
                QtGui.QMessageBox.critical(self,"Error", "No pt/px combination left after binding pt foot to the target piece, please change concentration or Temperature limit!")
                exit()
    if num==4:
        try:
           if os.path.getsize(configure.pt_foot_file1) and os.path.getsize(configure.px_foot_file1) and os.path.getsize(configure.tar_file1):
                return
           else:
                QtGui.QMessageBox.critical(self,"Error", "No pt/px combination left after include list check, please change your include list or blast parameters!")
                exit()
        except:
                QtGui.QMessageBox.critical(self,"Error", "No pt/px combination left after include list check, please change your include list or blast parameters!")
                exit()
    if num==5:
        try:
           if os.path.getsize(configure.pt_foot_file1) and os.path.getsize(configure.px_foot_file1) and os.path.getsize(configure.tar_file1):
                return
           else:
                QtGui.QMessageBox.critical(self,"Error", "No pt/px combination left after binding hairpins to the feet sequence,, please change your haipin sequence!")
                exit()
        except:
                QtGui.QMessageBox.critical(self,"Error", "No pt/px combination left after binding hairpins to the feet sequence,, please change your haipin sequence!")
                exit()
    if num==6:
        try:
           if os.path.getsize(configure.pt_foot_file1) and os.path.getsize(configure.px_foot_file1) and os.path.getsize(configure.tar_file1):
                return
           else:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()
        except:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()
    if num==6.5:
        try:
           if os.path.getsize(configure.whole_PT1):
                return
           else:
                #QtGui.QMessageBox.critical(None,"Error", "Didn't generate complete PT sequence, please try to change parameters!")
                configure.progress_simbol=11
                #QtGui.QMessageBox.critical(None,"Error", "Didn't generate complete PT sequence, please try to change parameters!")
                #exit()
                #cleanup_stop_thread()

        except:
                configure.progress_simbol=11
                #QtGui.QMessageBox.critical(None,"Error", "Didn't generate complete PT sequence, please try to change parameters!")
                #exit()
                #self._Thread__stop()

    if num==7:
        try:
           if os.path.getsize(configure.pt_foot_file1) and os.path.getsize(configure.px_foot_file1) and os.path.getsize(configure.tar_file1):
                return
           else:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()
        except:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()
    if num==8:
        try:
           if os.path.getsize(configure.pt_foot_file1) and os.path.getsize(configure.px_foot_file1) and os.path.getsize(configure.tar_file1):
                return
           else:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()
        except:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()
    if num==9:
        try:
           if os.path.getsize(configure.pt_foot_file1) and os.path.getsize(configure.px_foot_file1) and os.path.getsize(configure.tar_file1):
                return
           else:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()
        except:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()
    if num==10:
        try:
           if os.path.getsize(configure.pt_foot_file1) and os.path.getsize(configure.px_foot_file1) and os.path.getsize(configure.tar_file1):
                return
           else:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()
        except:
                QtGui.QMessageBox.critical(self,"Error", "Didn't generate any pt/px foot combinination, please try to change parameters!")
                exit()


def conbine_final(whole_PT_blast_result,name_seq_pt,whole_PX_blast_result,name_seq_px,whole_PT2,whole_PX1):
    pt_seq={}
    for record in SeqIO.parse(whole_PT2, "fasta") :
	pt_seq[record.id]=record.seq
    px_seq={}
    for record in SeqIO.parse(whole_PX1, "fasta") :
	px_seq[record.id]=record.seq

def self_bind(seq_file, bind_limit,symbol, output_filename): # check if the sequence in the "seq_file" less self bind than "bind_limit", "symbol" ='-pt' or '-px'
    command1='UNAFold.pl '+seq_file+' '+' -n DNA -C 0.000000005 -t 55'+'\n'
    os.system(command1)
    print command1
    command2='ct2rnaml '+seq_file+'.ct'+'\n'
    os.system(command2)
    print command2
    seq={}
    for record in SeqIO.parse(seq_file, "fasta") :
        #print "record.id",record.id
        seq[str(record.id)]=record.seq
    helix_pt=cal_helix(seq_file+'.rnaml')
    f=open(output_filename,'w')
    for i in range(len(helix_pt)):
        temp=helix_pt[i]
        #print temp[2],int(bind_limit)
        if sum(temp[2])<int(bind_limit):
            temp2=temp[0].split(symbol)
            #print "temp2[0]",temp2[0]
            print >>f, '>'+temp2[0]
            print >>f, seq[temp2[0]]
            configure.whole_self_bind[temp2[0]]=temp[2]
    f.close()

def cal_helix(file):
    input_file=open(file, 'r')
    site=-1
    leng=[]
    count=0
    all_cal_helix=[]
    mo_id=""
    line=input_file.readline()
    while (line!=''):
        line=input_file.readline()
	if (line.find('seq-data')!=-1):
	    count=0
	    leng=[]
	    #print "new"
        if (line.find('<molecule id')!=-1):
            mo_id_start=line.find('molecule id=')
	    mo_id_end=line.find('type')
	    mo_id=line[16:mo_id_end-2]
        elif (line.find('End of folding 1 for')!=-1):
            #print "end"
            count=count/2
            #print leng
            all_cal_helix.append([mo_id,count,leng])
            count=0
            leng=[]
            continue
        elif (line.find('helix')!=-1): #count the number of helix
            count=count+1
        elif (line.find('length')!=-1): #count the length of bonds in this helix
            t1=line.find('>')
            t2=line.find('</')
            leng.append(int(line[t1+1:t2]))
    input_file.close()
    return all_cal_helix


def self_Tm_check_eachother(pt_in,px_in,tar_in, lowest_tm,pt_outputfile,px_outputfile,tar_outputfile, p_melt_file):
    ptxfile=open(p_melt_file,'w') # for internal results
    ptx_command='melt.pl -n DNA -t 55 -N 1 -M 0 -C 0.000000005 '+pt_in+' '+px_in
    print ptx_command
    p = os.popen(ptx_command)
    print >>ptxfile, p.read()
    p.close()
    ptxfile.close()
    ptxfile = open(p_melt_file) # this is the .melt file
    ptf= open(pt_in) #pt or px sequence file waiting to be Tm checked
    pxf = open(px_in) #tar sequence file waiting to be Tm checked
    tarf = open(tar_in)
    ptx_bind_line = ptxfile.readline()
    output_pt = open(pt_outputfile,'w')
    output_px = open(px_outputfile,'w')
    output_tar = open(tar_outputfile,'w')
    i=0
    while ptx_bind_line.find('\t')==-1:
	ptx_bind_line = ptxfile.readline()
    print "lastline:", ptx_bind_line
    while 1:
	if ptx_bind_line == "":
		break
	i=i+1
	ptx_bind_line = ptxfile.readline().rstrip()
	ptid = ptf.readline()
	pxid = pxf.readline()
	ptseq = ptf.readline()
	pxseq = pxf.readline()
        tarid =tarf.readline()
        tarseq = tarf.readline()
	if ptx_bind_line !='':
	    ener_pt=ptx_bind_line.split("\t")
	    TM=float(ener_pt[3])
	    #print TM, lowest_tm
	    if (TM<lowest_tm):
		  print >>output_pt,ptid,
		  print >>output_pt,ptseq,
		  print >>output_px,pxid,
		  print >>output_px,pxseq,
                  print >>output_tar,tarid,
                  print >>output_tar,tarseq,
		  configure.Tm_p_foot_hair_tar[ptid]=TM
    ptxfile.close()
    output_pt.close()
    output_px.close()



def self_Tm_check(p_foot_hair_file,tar_foot_hair_file,tm_below,output_pname,p_melt_file): #check the tm when p foot with hairpin binds to the target sequence , p_foot_hair_file is the pt or px sequence file waiting to be Tm checked
    f=open(p_melt_file,'w')
    ptx_command='melt.pl -n DNA -t 55 -N 1 -M 0 -C 0.000000005 '+p_foot_hair_file+' '+tar_foot_hair_file
    print ptx_command
    p = os.popen(ptx_command)
    print >>f, p.read()
    p.close()
    f.close()
    #os.system(command)
    ptxfile = open(p_melt_file) # this is the .melt file
    pf= open(p_foot_hair_file) #pt or px sequence file waiting to be Tm checked
    tarf = open(tar_foot_hair_file) #tar sequence file waiting to be Tm checked
    ptx_bind_line = ptxfile.readline()
    output_p = open(output_pname,'w')
    i=0
    while ptx_bind_line.find('\t')==-1:
	ptx_bind_line = ptxfile.readline()
	#print ptx_bind_line
    ptx_bind_line = ptxfile.readline()
    #print ptx_bind_line
    while 1:
	if ptx_bind_line == "":
		break
	i=i+1
	ptx_bind_line = ptxfile.readline().rstrip()
	pid = pf.readline()
	tarid = tarf.readline()
	pseq = pf.readline()
	tarseq = tarf.readline()
	print "enter tm calculation"
	if ptx_bind_line !='':
            #print ptx_bind_line
	    ptx_bind_line=ptx_bind_line.replace('-','')
	    ener_pt=ptx_bind_line.split("\t")
            #Tm=float(ener_pt[3])
	    TM=float(ener_pt[1])/((float(ener_pt[2])-(1.987*math.log(configure.pt_foot_conc)))/1000)-273.15
	    print TM, tm_below
	    if (TM<tm_below):
		  print >>output_p,pid,
		  print >>output_p,pseq,
		  configure.Tm_p_foot_hair_tar[pid]=TM
    ptxfile.close()
    output_p.close()
    return (output_pname)




def check_binds(seq1,seq2):
    t=alignment.needle(seq1,seq2)
    pieces=t.split()
    leng=[]
    for i in pieces:
	leng.append(len(i))
    #match_num=len(t)-t.count(' ')
    match_num=max(leng)
    return match_num



def bind_hairpin(hairpin_file,extra_base_file, pt_hairpin, px_hairpin, tar_hairpin,hair_tar_max): #bind hairpin to pt foot and px foot

    extra_base={}
    for record in SeqIO.parse(extra_base_file, "fasta") :
	extra_base[record.id]=record.seq

    hairpin={}
    for record in SeqIO.parse(hairpin_file, "fasta") :
	hairpin[record.id]=record.seq
	id2=record.id+'rev'
	hairpin[id2]=revcomp(str(record.seq))

    pt_seq={}
    for record in SeqIO.parse(configure.pt_foot_file3, "fasta") :
	pt_seq[record.id]=record.seq

    px_seq={}
    for record in SeqIO.parse(configure.px_foot_file3, "fasta") :
	px_seq[record.id]=record.seq

    tar_seq={}
    for record in SeqIO.parse(configure.tar_file3, "fasta") :
	tar_seq[record.id]=record.seq



    pt_long_file=open(pt_hairpin,'w')
    px_long_file=open(px_hairpin,'w')
    tar_long_file=open(tar_hairpin,'w')
    iterator=1
    for e in extra_base.keys(): #PT site extra base
	for e2 in extra_base.keys(): #PX site extra base
	    for tar in tar_seq.keys():
		if tar.split('-')[3]==tar.split('-')[4]:
		    gap='0'
		else:
                    start=int(tar.split('-')[3])-int(tar.split('-')[2])
                    end=int(tar.split('-')[4])-int(tar.split('-')[2])
                    #print start, end
		    #gap=tar_seq[start:end]
                    gap='G'
		    #print whether_match(extra_base[e],extra_base[e2],gap)
		if (whether_match(extra_base[e],extra_base[e2],gap)==0):
                    #print "enter1"
		    for h in hairpin.keys():
                            #print check_binds(hairpin[h],tar_seq[tar]),check_binds(revcomp(str(hairpin[h])),tar_seq[tar]), hair_tar_max
			    if check_binds(hairpin[h],tar_seq[tar])<hair_tar_max and check_binds(revcomp(str(hairpin[h])),tar_seq[tar])<hair_tar_max: #not too much hybrid between hairpin and target
                                #print "enter2"
                                str_iterator=str(iterator)
                                partial_name='-'.join(tar.split('-')[1:])
				pt_long_old_id='pt-'+partial_name
				pt_long_new_id=pt_long_old_id+'-'+'e'+'-'+e+'-'+'h'+'-'+h+'-'+str_iterator
				pt_long_seq=hairpin[h]+extra_base[e]+pt_seq[pt_long_old_id]
				#pt_long_seq=str(pt_seq[pt_long_old_id])
				#pt_long_seq='ATATCCATC'
				px_long_old_id='px-'+partial_name
				px_long_new_id=px_long_old_id+'-'+'e'+'-'+e+'-'+'h'+'-'+h+'-'+str_iterator
				px_long_seq=px_seq[px_long_old_id]+extra_base[e]+revcomp(str(hairpin[h]))
                                #px_long_seq='ATATCCATC'
                                #print "partial:", partial_name,pt_long_new_id,pt_long_seq
				print >>pt_long_file, '>'+pt_long_new_id
				print >>pt_long_file, pt_long_seq
				print >>px_long_file, '>'+px_long_new_id
				print >>px_long_file, px_long_seq
				print >>tar_long_file, '>'+tar+'-'+str_iterator
				print >>tar_long_file, tar_seq[tar]
				iterator=iterator+1
    px_long_file.close()
    pt_long_file.close()
    tar_long_file.close()
    return(configure.pt_foot_hair_file1, configure.px_foot_hair_file1,configure.tar_hair_file1)


def whether_match(a,b,c):
    comp=['AT','TA','GC','CG']
    for i in range(len(a)):
	for j in range(len(b)):
	    for h in range(len(c)):
		if ((a[i]+b[j] in comp) | (a[i]+c[h] in comp) | (b[j]+c[h] in comp)):
		    return 1 #there is a match
		else :
		    return 0 # there is no match

def Tm_pfoot(pt_file, px_file, tar_file, pt_tm,px_tm,pt_bind_ok_filename, px_bind_ok_filename, tar_bind_ok_filename, pt_melt, px_melt, pt_foot_conc,px_foot_conc): # to bind pt foot to the target piece, check tm, bind px foot to the target piece, check tm, only output those pairs with a larger tm than input tm criteria
    f=open(pt_melt,'w')
    pt_command='melt.pl -n DNA -t 55 -N 1 -M 0 -C '+str(pt_foot_conc)+' '+tar_file+' '+pt_file
    print "pt_command", pt_command
    p = os.popen(pt_command)
    print >>f, p.read()
    #print "p.read", p.read()
    p.close()
    f.close()
    f=open(px_melt,'w')
    px_command='melt.pl -n DNA -t 55 -N 1 -M 0 -C '+str(px_foot_conc)+' '+tar_file+' '+px_file
    print "pt_command", pt_command
    p = os.popen(px_command)
    print >>f, p.read()
    p.close()
    f.close()
    ###########################
    Tm_calculation(pt_file, px_file, tar_file,float(pt_tm),float(px_tm),pt_bind_ok_filename, px_bind_ok_filename, tar_bind_ok_filename, pt_melt,px_melt, pt_foot_conc,px_foot_conc)


def Tm_calculation(pt_file, px_file, tar_file,ptm,pxm,pt_bind_ok_filename,px_bind_ok_filename, tar_bind_ok_filename, pt_melt,px_melt, pt_foot_conc,px_foot_conc):
#sequence file,ptf and pxf; dH,dH file ptfilename and pxfielname; any pt tm larger than ptm and px tm larger than px tm are chosen
    ptfile = open(pt_melt)
    pxfile = open(px_melt)
    ptf= open(pt_file)
    pxf = open(px_file)
    tarf = open(tar_file)
    temp = open('temp_test.txt','w')
    pt_bind_line = ptfile.readline()
    #print "line", pt_bind_line,"file", pt_melt
    px_bind_line = pxfile.readline()
    #print "line", px_bind_line,"file", px_melt
    output_pt = open(pt_bind_ok_filename,'w')
    output_px = open(px_bind_ok_filename,'w')
    output_tar = open(tar_bind_ok_filename,'w')
    i=0
    while pt_bind_line.find('\t')==-1:
	pt_bind_line = ptfile.readline()
	px_bind_line = pxfile.readline()
    print "last:", pt_bind_line, px_bind_line,
    while 1:
	if pt_bind_line == "":
		break
	i=i+1
	pt_bind_line = ptfile.readline().rstrip()
	px_bind_line = pxfile.readline().rstrip()
	#print "", pxfile, px_bind_line
	ptid = ptf.readline()
	pxid = pxf.readline()
	tarid = tarf.readline()
	ptseq = ptf.readline()
	pxseq = pxf.readline()
	tarseq = tarf.readline()
	#print "enter tm calculation"
	if pt_bind_line !='' and px_bind_line !='':
	    pt_bind_line=pt_bind_line.replace('-','')
	    ener_pt=pt_bind_line.split("\t")
	    px_bind_line=px_bind_line.replace('-','')
	    ener_px=px_bind_line.split("\t")
	    TM_pt=float(ener_pt[1])/((float(ener_pt[2])-(1.987*math.log(pt_foot_conc)))/1000)-273.15
	    TM_px=float(ener_px[1])/((float(ener_px[2])-(1.987*math.log(px_foot_conc)))/1000)-273.15
        print >>temp, ptid,
        print >>temp, pxid,
        print >>temp, TM_pt,
        print >>temp, TM_px,
        if (TM_pt>ptm and TM_px>pxm):
                print "ok", TM_pt, ptm, TM_px, pxm
		print >>output_pt,ptid,
		print >>output_pt,ptseq,
		print >>output_px,pxid,
		print >>output_px,pxseq,
		print >>output_tar,tarid,
		print >>output_tar,tarseq,
        configure.Tm_pt_foot_tar[ptid]=TM_pt
        configure.Tm_px_foot_tar[ptid]=TM_px

    ptfile.close()
    pxfile.close()
    output_pt.close()
    output_px.close()
    output_tar.close()


def ptpx_foot(seq,min_pt,max_pt,min_px,max_px,e,direction): # to get list of the pt and px foot based on the input sequence, pt,px length limit and gap length limit
	gactc_start=[i - 1 for i in range(len(seq)) if seq.startswith('GACTC', i - 1) or seq.startswith('GAGTC', i - 1)]
	s=set([])
	seq_len=len(seq)
	for i in range(seq_len): #start
		for j in range(int(min_pt),int(max_pt)): #pt foot length
			for h in range(0,int(e)): #gap
				for p in range(int(min_px),int(max_px)): #px length
					num=0 #to calculate this range covered how many gactc
					if gactc_start !=[]:
						for each_gactc_start in gactc_start:
							if (i<each_gactc_start) and (i+j+h+p>each_gactc_start):
								num=num+1;
						if num==0: # if none of the gactc is covered in this range
							if i+j+h+p<=seq_len:
								s.add(str(direction)+'-'+str(i)+'-'+str(i+j)+'-'+str(i+j+h)+'-'+str(i+j+h+p))
					else:
						if i+j+h+p<=seq_len:
							s.add(str(direction)+'-'+str(i)+'-'+str(i+j)+'-'+str(i+j+h)+'-'+str(i+j+h+p))
	return(s)

def write_ptpx_foot_file_check_foot(ptpx_foot_list,input_seq,maximum_num, pt_filename,px_filename, tar_filename): # to write the pt px foot into a fast format file, sequece named in a specific format
    ptfile=open(pt_filename,'w')
    pxfile=open(px_filename,'w')
    tarfile=open(tar_filename,'w')
    n=0 # this is the id for pt/px/tar pair
    for tag in ptpx_foot_list:
        eles=tag.split('-')
        i=eles[1]
        j=eles[2]
        h=eles[3]
        p=eles[4]
        num_i=int(i)
        num_j=int(j)
        num_h=int(h)
        num_p=int(p)
        if eles[0]=='p':
            pt=revcomp(input_seq[num_i:num_j])
            px=revcomp(input_seq[num_h:num_p])
            if num_h!=num_j:
                gap=input_seq[num_j:num_h]
            else:
                gap='0'
            #print "pt/px foot binding each other", check_binds(pt,px), maximum_num
            if check_binds(pt,revcomp(px))<maximum_num+1 and (gap.find('G')+gap.find('C')<0):#check if pt and px foot have too much hybridization
                n=n+1
                num=str(n)
                print >>ptfile, '>'+'pt-'+'p'+'-'+i+'-'+j+'-'+h+'-'+p+'-'+num+'\n'+pt
                print >>pxfile, '>'+'px-'+'p'+'-'+i+'-'+j+'-'+h+'-'+p+'-'+num+'\n'+px
                print >>tarfile, '>'+'tar-'+'p'+'-'+i+'-'+j+'-'+h+'-'+p+'-'+num+'\n'+input_seq[num_i:num_p]
        if eles[0]=='n':
            #print "pt/px foot binding each other", check_binds(pt,px), maximum_num
            pt=revcomp(revcomp(input_seq)[num_i:num_j])
            px=revcomp(revcomp(input_seq)[num_h:num_p])
            if num_h!=num_j:
                gap=input_seq[num_j:num_h]
            else:
                gap='0'
            if check_binds(pt,revcomp(px))<maximum_num and (gap.find('G')+gap.find('C')<0):#check if pt and px foot have too much hybridization
                n=n+1
                num=str(n)
                print >>ptfile, '>'+'pt-'+'n'+'-'+i+'-'+j+'-'+h+'-'+p+'-'+num+'\n'+pt
                print >>pxfile, '>'+'px-'+'n'+'-'+i+'-'+j+'-'+h+'-'+p+'-'+num+'\n'+px
                print >>tarfile, '>'+'tar-'+'n'+'-'+i+'-'+j+'-'+h+'-'+p+'-'+num+'\n'+revcomp(input_seq)[num_i:num_p]
    print "tarfile",tarfile
    ptfile.close()
    pxfile.close()
    tarfile.close()
    return(configure.pt_foot_file1,configure.px_foot_file1,configure.tar_file1)

def submit_to_query(line): #from the input blank to generate the entrez query for NCBIWWW.qblast
    strlist=str(line).split(' OR ')
    txid_line=''
    for valist in strlist:
        txid_num=valist[valist.find('(taxid:')+7:valist.find(')')]
        txid='txid'+txid_num+' [ORGN]'
        if txid_line != '':
            txid_line=txid_line+' OR '+txid
        else:
            txid_line=txid
    return txid_line


def combine_check_include(include_line,mismatch_num,tar_file,output_pt_file, output_px_file, output_tar_file, blastfile): #check if the target piece has enough similarity (mismatch)to all organism in the the include list
    (results)=check_blast(tar_file,include_line,mismatch_num, 0,'in',blastfile)
    pt_foot_file3=open(output_pt_file,'w')
    px_foot_file3=open(output_px_file,'w')
    tar_file3=open(output_tar_file,'w')
    i=0
    pt_foot={}
    px_foot={}
    tar={}
    num_organism=len(include_line.split('OR'))
    for record in SeqIO.parse(configure.pt_foot_file2, "fasta") :
        num=record.id.split('-')[len(record.id.split('-'))-1]
        pt_foot[num]='>'+record.id+'\n'+record.seq
    for record in SeqIO.parse(configure.px_foot_file2, "fasta") :
        num=record.id.split('-')[len(record.id.split('-'))-1]
        px_foot[num]='>'+record.id+'\n'+record.seq
    for record in SeqIO.parse(configure.tar_file2, "fasta") :
        num=record.id.split('-')[len(record.id.split('-'))-1]
        tar[num]='>'+record.id+'\n'+record.seq
    if (len(px_foot)!=len(pt_foot) or len(tar)!=len(px_foot)):
        print "error, the amount of px_foot, pt_foot and target sequence are not equal!"
    #print len(results),len(pt_foot),len(px_foot),len(tar),num_organism
    tar_blast={}
    for pair in results.keys():
        if tar_blast.has_key(pair[0]):
            tar_blast[pair[0]]=tar_blast[pair[0]]+1
        else:
            tar_blast[pair[0]]=1
    for j in tar_blast.keys():
        #if tar_blast[j]==num_organism-1: # this needs to be change
        if 1==1:
            tagid=j.split('-')[len(j.split('-'))-1]
            print >>pt_foot_file3,pt_foot[tagid]
            print >>px_foot_file3,px_foot[tagid]
            print >>tar_file3,tar[tagid]
            configure.tar_include[tagid]=tar_blast[j]
    pt_foot_file3.close()
    px_foot_file3.close()
    tar_file3.close()

def combine_pt_px_result(whole_PT2, whole_PX1,final_ptx_filename):
    final_PTX=open(final_ptx_filename,'w')
    pt_final={}
    px_final={}
    for record in SeqIO.parse(whole_PT2, "fasta") :
        pt_final[record.id]='>'+record.id+'\n'+record.seq
    for record in SeqIO.parse(whole_PX1, "fasta") :
        num=record.id.split('-')[len(record.id.split('-'))-1]
        px_final[num]='>'+record.id+'\n'+record.seq
    for pt_fullname in pt_final.keys():
        list_pt_name=pt_fullname.split('-')
        len_pt_name=len(list_pt_name)
        ptx_num=list_pt_name[len_pt_name-4]
        print >>final_PTX,pt_final[pt_fullname]
        print >>final_PTX,px_final[ptx_num]
    print "all done!"
    final_PTX.close()


def combine_check_exclude(whole_PT2, whole_PX1,exclude_line,pt_longer_num,px_longer_num,pt_blast,px_blast,final_ptx_filename): # check if any sequence cross react with "exclude list" organism
    GACTC_YES=0
    (whole_PX_blast_result)=check_blast(whole_PX1,exclude_line,px_longer_num,GACTC_YES,'ex',px_blast)
    print "check_blast PX is over"
    GACTC_YES=1
    (whole_PT_blast_result)=check_blast(whole_PT2,exclude_line,pt_longer_num,GACTC_YES,'ex',pt_blast)
    print "check_blast PT is over"
    configure.progress_simbol=5
    keys_px={}
    keys_px_num={}
    keys_pt={}
    for i in whole_PT_blast_result.keys():
        keys_pt[i[0]]=1 # these are bad pt
    for i in whole_PX_blast_result.keys():
        keys_px[i[0]]=1 # these are bad px
        list=i[0].split('-')
        num=list[len(list)-1]
        keys_px_num[num]=i  # these are bad px
    final_PTX=open(final_ptx_filename,'w')
    pt_final={}
    px_final={}
    print "set final"
    for record in SeqIO.parse(whole_PT2, "fasta") :
        pt_final[record.id]='>'+record.id+'\n'+record.seq
    for record in SeqIO.parse(whole_PX1, "fasta") :
        num=record.id.split('-')[len(record.id.split('-'))-1]
        px_final[num]='>'+record.id+'\n'+record.seq
    print "set final parially done!"
    for pt_fullname in pt_final.keys():
        list_pt_name=pt_fullname.split('-')
        len_pt_name=len(list_pt_name)
        ptx_num=list_pt_name[len_pt_name-4]
        if not keys_pt.has_key(pt_fullname): # if it's not in the "probalby bad" pt list
            if px_final.has_key(ptx_num): # if there is a px linked to this pt
                if not keys_px_num.has_key(ptx_num): # if this px is not in the "probaby bad" px list
                    print >>final_PTX,pt_final[pt_fullname]
                    print >>final_PTX,px_final[ptx_num]
                else:
                    px_fullname=keys_px_num[ptx_num]
                    match_sum_px=0
                    for pair in whole_PX_blast_result.keys():
                        if pair[0]==px_fullname:
                            match_sum_px=match_sum_px+whole_PX_blast_result[pair]
                    if match_sum_px==0: # this px is in the "probably bad" list but not really bad
                        print >>final_PTX,pt_final[pt_fullname]
                        print >>final_PTX,px_final[ptx_num]
        else:
            match_sum_pt=0
            for pair in whole_PT_blast_result.keys():
                if pair[0]==pt_fullname:
                    match_sum_pt=match_sum_pt+whole_PT_blast_result[pair]
            if match_sum_pt==0: # this pt is in the "probably bad" list but not really bad
                if px_final.has_key(ptx_num): # if there is a px linked to this pt
                    if not keys_px_num.has_key(ptx_num): # if this px is not in the "probaby bad" px list
                        print >>final_PTX,pt_final[pt_fullname]
                        print >>final_PTX,px_final[ptx_num]
                    else:
                        px_fullname=keys_px_num[ptx_num]
                        match_sum_px=0
                        for pair_px in whole_PX_blast_result.keys():
                            if pair_px[0]==px_fullname:
                                match_sum_px=match_sum_px+whole_PX_blast_result[pair_px]
                        if match_sum_px==0:
                            print >>final_PTX,pt_final[pt_fullname]
                            print >>final_PTX,px_final[ptx_num]
    print "all done!"
    final_PTX.close()

def include_check(blast_result_file,include_line,seq_name_list,num):
    #print seq_name_list
    include_line2=re.sub('\(taxid:[0-9]+\)','',str(include_line))
    include_list=include_line2.split(' OR ')
    #results = [[1] * len(include_list) for row in range(len(seq_name_list))]
    #results=[[0]* len(include_list) for i in range(len(seq_name_list))]
    #print "len(include_list):",len(include_list)
    results={}
    i=0
    name_seq={}
    for record in NCBIXML.parse(open(blast_result_file)) :
        #print "haha",i, record.query
	if record.alignments :
            name=record.query
	    #print name
            name_seq[i]=name
            #num=split(record.query,"_")[1]
	    #num=int(num,10)
            #leng=split(record.query,"_")[2]
            query_len=record.query_letters # the length of your submitted sequence
            for align in record.alignments :
		for hsp in align.hsps :
		    for index in range(len(include_list)):
			if align.hit_def.lower().find(include_list[index].lower()) >-1:
			    if hsp.identities >= query_len-int(num): # equal or less mistmatch
                                results[(name,index)]=1
        i=i+1
    return (results)

def exclude_check(blast_result_file,exclude_line,seq_name_list,num, GACTC_YES):
    print "exclude_line:", exclude_line
    exclude_line2=re.sub('\(taxid:[0-9]+\)','',str(exclude_line))
    print "exclude_line2:", exclude_line2
    exclude_list=exclude_line2.split(' OR ')
    results = {}
    blast_result_file_handle = open(blast_result_file)
    for record in NCBIXML.parse(blast_result_file_handle) :
	name=record.query
	if record.alignments :
            for align in record.alignments :
		for hsp in align.hsps :
		    for index in range(len(exclude_list)):
			#if align.hit_def.lower().find(exclude_list[index].lower()) >-1:
                        if 1==1:
			    #results[num-1][index]=0
			    if hsp.identities >int(num):
                                if GACTC_YES: # this is for PT
                                    GACTC_start=int(name[name.find('GACTC')+5:])
                                    if (hsp.query_start<GACTC_start,GACTC_start+5<hsp.query_end):
                                        #results[name+'@'+exclude_list[index]]=1
                                        results[(name,index)]=1
                                else: # this is for PX
                                    if hsp.query_start<4:
                                    #results[name+'@'+exclude_list[index]]=1
                                        results[(name,index)]=1
    #print "type::::::::::::::::", type(results)
    return (results)



def check_blast(input_file,taxid_line,num, GACTC_YES,ex_in,blast_filename):
    query_line=submit_to_query(taxid_line)
    input_file = open(input_file,'r')
    typ=''
    seq_name_list=[]
    blast_result_file= open(blast_filename,"w")
    print "input_file", input_file
    for seq_record in SeqIO.parse(input_file, "fasta"):
        if (len(typ)>5000):
            result_handle = NCBIWWW.qblast("blastn", "nr",typ,word_size=13,hitlist_size=100,entrez_query=query_line,expect=10)
            #print typ
            t=result_handle.read()
            blast_result_file.write(t)
            typ=''
            print "5k done!"
        #print "seq", seq_record.format('fasta')
        typ=typ+seq_record.format('fasta')
        seq_name_list.append(seq_record.id)
    blast_result_file.close()
    #print typ
    if typ!='':
        blast_result_file= open(blast_filename,"a")
        result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=query_line,expect=10)
        t2=result_handle.read()
        blast_result_file.write(t2)
        blast_result_file.close()
    input_file.close()
    print "blast job done!"
    if ex_in=='ex':
            (check_blast_results)=exclude_check(blast_filename,taxid_line,seq_name_list,num,GACTC_YES)
    elif ex_in=='in':
            (check_blast_results)=include_check(blast_filename,taxid_line,seq_name_list,num)
    return (check_blast_results)


def combine_PT(ET_file,PT_file,target_file, PT_output,overhang_tar_max): # to combine the PT_file(both foot and hairpin) with good performance template file(ET_file) to get the whole_PT_file
    tar_seq={}
    f=open(PT_output,'w')
    print "start"
    for record in SeqIO.parse(target_file, "fasta") :
        num=record.id.split('-')[len(record.id.split('-'))-1]
	tar_seq[num]=record.seq
    overhang_bind_pt={}
    for ET_record in SeqIO.parse(ET_file, "fasta"):
	for PT_record in SeqIO.parse(PT_file, "fasta"):
		overhang=str(ET_record.seq) #this needs to be fixed, this is just the template,not the whole overhang yet
		num=PT_record.id.split('-')[len(PT_record.id.split('-'))-1]
		if tar_seq.has_key(num):
                    target=tar_seq[num]
                    t=alignment.needle(revcomp(overhang),target)
                    #print "checking"
                    if check_binds(revcomp(overhang),target) <overhang_tar_max and t.find('GAGTC')==-1: # not too much binds with target
                        #print "smaller?",  check_binds(revcomp(overhang),target),overhang_tar_max
                        newseq=ET_record.seq+PT_record.seq
                        newid=PT_record.id+'-'+ET_record.id+'-'+'GACTC'+str(newseq.upper().find('GACTC'))
                        print >>f, '>'+newid
                        print >>f, newseq
                        overhang_bind_pt[newid]=check_binds(revcomp(overhang),target)
                    #else:
                        #print len(revcomp(overhang)),
    #print "end"
    f.close()


def delete_file(filename):
    for file in os.listdir("."):
        if os.path.isfile(file) and file.startswith(filename):
             try:
                  os.remove(file)
             except Exception,e:
                  print e



if __name__ == "__main__":
    target_file='J:\PROXAR\3-2-2012\short_example.txt'
    include_line='Mycobacterium avium (taxid:1764) OR Mycobacterium bovis (taxid:1765)'
    exclude_list='homo sapiens (taxid:9606) OR homo (taxid:9605)'
    hairpin_file='hairpin.txt'
    extra_base_file='extra_base.txt'
    templates_file='good_template.txt'
    exclude_GAGTC=1
    two_direction=1
    PT_foot_min=23
    PT_foot_max=29
    PT_Tm_above='-5'
    PX_foot_min=21
    PX_foot_max=23
    PX_Tm_above='-5'
    max_gap=2
    main_page(target_file,include_line, exclude_list, hairpin_file,extra_base_file, templates_file,exclude_GAGTC,two_direction,PT_foot_min,PT_foot_max,PT_Tm_above,PX_foot_min,PX_foot_max,PX_Tm_above,max_gap)
