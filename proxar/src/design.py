'''
This is for primer desin functions
Created on 2013-1-26

@author: jifeng
'''
#!/usr/bin/python -d
import string
import configure
from Bio import SeqIO
from Bio.Blast import NCBIWWW,NCBIXML
import TmDeltaG
import alignment
class Feet():
    def __init__(self, start_pos1, end_pos1, start_pos2, end_pos2, direction):
        self. pt_start= start_pos1
        self. pt_end = end_pos1
        self. px_start = start_pos2
        self. px_start = end_pos2
        self. direction = direction
        self. pt_tm=0
        self. px_tm=0
    def get_pt(self):
        if self.direction == 'p':
            return revcomp(configure.input_seq[self. pt_start:self. pt_end])
        if self. direction == 'n':
            target=revcomp(configure.input_seq)
            self.seq = revcomp(target[self. pt_start:self. pt_end])
        return self.seq
    def get_px(self):
        if self. direction == 'p':
            return revcomp(configure.input_seq[self. px_start:self. px_end])
        if self. direction == 'n':
            target=revcomp(configure.input_seq)
            self.seq = revcomp(target[self. px_start:self. px_end])  
        return self.seq          

class PTX_long(Feet):
    def __init(self,feet,extrabaset,extrabasex,hairpin):
        self.pt_long=feet
        self.px_long=feet

class primer_set(Feet):
    def __init__(self, feet, extra_base, hairpin, template):
        self.feet=feet
        self.extra_base=extra_base
        self.hairpin=hairpin
        self.template=template

    def bonds(seq1,seq2):
        return MeltingTemp_jf.bonds(qseq, qseq, mono_conc=mono_conc, Mg_conc=diva_conc, dntp_conc=dntp_conc,deltaH=deltaH, deltaS=deltaS)
        
def bind_hairpin_template(blast_ok_feet_list,hairpin_file,extra_base_file):
        pass

def get_foot(seq,start,end,direction):
        if direction == 'p':
            return revcomp(seq[start:end])
        if direction == 'n':
            target=revcomp(seq)
            foot = revcomp(target[start:end])  
            return foot            
def ptpx_foot(seq,min_pt,max_pt,min_px,max_px,max_e,direction,pt_min_tm,pt_max_tm, px_min_tm, px_max_tm,pt_foot_conc,px_foot_conc): # to get list of the pt and px foot based on the input sequence, pt,px length limit and gap length limit
    monovalent_conc=1000
    divalent_conc=0
    dNTP_conc=0.00008
    pt_foot_conc=pt_foot_conc*1000000000
    px_foot_conc=px_foot_conc*1000000000 
    gactc_start=[i - 1 for i in range(len(seq)) if seq.startswith('GACTC', i - 1) or seq.startswith('GAGTC', i - 1)]
    s=set([])
    seq_len=len(seq)
    for i in range(seq_len): #start
        for j in range(int(min_pt),int(max_pt)): #pt foot length
            pt_tm=TmDeltaG.calTm(get_foot(seq,i,i+j,direction), revcomp(get_foot(seq,i,i+j,direction)), monovalent_conc, divalent_conc, pt_foot_conc, dNTP_conc) # mono, div, and dNTP using mM, oligo_conc, oligo_conc nM
            if pt_tm>pt_min_tm and pt_tm<pt_max_tm:
                for h in range(0,max_e): #gap
                    for p in range(int(min_px),int(max_px)): #px length
                        px_tm=TmDeltaG.calTm(get_foot(seq,i,i+j,direction), revcomp(get_foot(seq,i,i+j,direction)), monovalent_conc, divalent_conc, px_foot_conc, dNTP_conc) # mono, div, and dNTP using mM, oligo_conc, oligo_conc nM
                        if px_tm>px_min_tm and px_tm<px_max_tm:
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


def write_fasta_feet(feet_list,filename):
    f=open(filename, 'w')
    p_seq={}
    for feet in feet_list:
        p_seq[feet.pt_start,feet.pt_end, feet.direction]=feet.get_pt
        p_seq[feet.px_start,feet.px_end, feet.direction]=feet.get_px        
    for i in p_seq.keys():
        print >>f, '>', i[0],'-', i[1],'-', i[2]
        print >>f, p_seq[i]
    f.close()

def check_blast_ex(input_filename,taxid_line,num, GACTC_YES,blast_filename,perct):
    query_line=submit_to_query(taxid_line)
    #print "ex", query_line
    input_file = open(input_filename,'r')
    typ=''
    blast_result_file= open(blast_filename,"w")
    for seq_record in SeqIO.parse(input_file, "fasta"):
        if (len(typ)>300):
            result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=query_line,expect=10)
            t=result_handle.read()
            blast_result_file.write(t)
            typ=''
            #print "300 done!"
        typ=typ+seq_record.format('fasta')
    blast_result_file.close()
    blast_result_file= open(blast_filename,"a")
    result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=query_line,expect=10)
    t2=result_handle.read()
    blast_result_file.write(t2)
    blast_result_file.close()
    input_file.close()
    print "blast job done!"
 
def check_blast_in(input_filename,taxid_line,num, GACTC_YES,blast_filename,perct):
    strlist=str(taxid_line).split(' OR ')
    for valist in strlist:
        txid_num=valist[valist.find('(taxid:')+7:valist.find(')')]
        blast_result_file= open(blast_filename+txid_num,"w")
        txid='txid'+txid_num+' [ORGN]'
        typ='' #this is for the input sequence
        input_file = open(input_filename,"r")
        for seq_record in SeqIO.parse(input_file, "fasta"):
            if (len(typ)>200):
                result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=txid,expect=10)
                t=result_handle.read()
                blast_result_file.write(t)
                typ=''
                #print "200 done!"
            typ=typ+seq_record.format('fasta')
            #print "wating", typ,"finished waiting","\n\n"
        if (len(typ)>0):
            #print "working on the leftover"
            result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=txid,expect=10)
            t2=result_handle.read()
            #print typ
            blast_result_file.write(t2)
        input_file.close()
        blast_result_file.close()
    print "blast job done!"


def combine_check_exinclude(input_file,include_line,exclude_line,num,perct,seq_name_list,conserve_mismatch):
    blast_filename='blast_result.xml'
    if exclude_line=='':
        check_blast_results_ex={}
    else:
        check_blast_ex(input_file,exclude_line,num, blast_filename)
        (check_blast_results_ex)=exclude_check(blast_filename,exclude_line,seq_name_list) #if len(results[name])!=0, then good
    if include_line=='':
        check_blast_results_in={}
    else:
        check_blast_in(input_file,include_line,num,blast_filename)
        (check_blast_results_in)=include_check(blast_filename,include_line,seq_name_list,conserve_mismatch) #if len(results[name])==len(number of conserved organism), then good
    #print "ex:,", check_blast_results_ex
    #print "in:",check_blast_results_in
    return (check_blast_results_in,check_blast_results_ex)
        
def exclude_check(blast_result_file,exclude_line,seq_name_list):
    results = {}
    blast_result_file_handle = open(blast_result_file)
    for record in NCBIXML.parse(blast_result_file_handle) :
        name=record.query # name is the submitted sequence name
        results[name]='' #begin with not found anything yet
        if record.alignments :
            for align in record.alignments :
                for hsp in align.hsps :
                    #print "hsp:",hsp.identities, query_len, perct
                    if hsp.identities ==hsp.align_len:
                        if results.has_key(name):
                            if hsp.identities>results[name]:
                                results[name]=hsp.identities
                        else:
                            results[name]=hsp.identities
    return (results) 
    #return(found)
def include_check(blast_result_filename,include_line,seq_name_list,max_mismatch):
    strlist=str(include_line).split(' OR ')
    results={} # intermediate result to show if what organism are conserved in that seq
    results_final={} # final result to show if an seq is conserved or not in all the rquested organism
    for valist in strlist:
        txid_num=valist[valist.find('(taxid:')+7:valist.find(')')]
        blast_result_file= open(blast_result_filename+txid_num,"r")      
        found={}
        for record in NCBIXML.parse(blast_result_file):
            name=record.query
            min_len= record.query_letters-max_mismatch
            if not found.has_key(name): 
                if record.alignments :
                    for align in record.alignments :
                        for hsp in align.hsps :
                            #print "blast: ",hsp.identities,name,query_len,int(num)
                            if hsp.identities == hsp.align_len and hsp.identities>=min_len: # 100% match and has more identities than requirement length of matches
                                found[name]=1 # this valst is conserved in current 
                                if results.has_key(name):
                                    temp=results[name]
                                    temp.append(txid_num)
                                    results[name]=temp
                                else:
                                    temp=[txid_num]
                                    results[name]=temp
                        #print name,results[name]
        blast_result_file.close()
    len_organ=len(strlist)
    for i in results.keys():
        if len(results[i])==len_organ:
            results_final[i]=1
    return (results_final)

def blast_feet_chek(in_blast_result,ex_blast_result, feet_list,pt_minlength, px_minlength):
    feet_ok_list=[]
    for feet in feet_list:
        name1=str(feet. start1)+'-'+str(feet. end1)+'-'+str(feet.direction)
        name2=str(feet. start2)+'-'+str(feet. end2)+'-'+str(feet.direction)
        if name1 in in_blast_result and name2 in in_blast_result: #good conserved sequences
            if ex_blast_result[name1]<pt_minlength and ex_blast_result[name2]<px_minlength: # 
                feet_ok_list.append(feet)
        return feet_ok_list

def revcomp(dna):
    comp = dna.translate(string.maketrans("AGCTagct", "TCGAtcga"))
    lcomp = list(comp)
    lcomp.reverse()
    return  string.join(lcomp,  "")
def bind_hairpin(hairpin_file,extra_base_file,template_file, feet_list, lowest_ptx_long_tm,hair_tar_max): #bind hairpin to pt foot and px foot
    extra_base={}
    template={}
    for record in SeqIO.parse(extra_base_file, "fasta") :
        extra_base[record.id]=record.seq
    hairpin={}
    for record in SeqIO.parse(hairpin_file, "fasta") :
        hairpin[record.id]=record.seq
        id2=record.id+'rev'
        hairpin[id2]=revcomp(str(record.seq))
    for record in SeqIO.parse(template_file,"fasta"):
        template[record.id]=record.seq
    primer_set_list=[]
    target_seq=''
    for e in extra_base.keys(): #PT site extra base
        for e2 in extra_base.keys(): #PX site extra base
            for feet in feet_list:
                if feet.direction=='p':
                    gap=revcomp(configure.input_seq)[feet.end_pos1:feet.start_pos2]
                    target_seq=configure.input_seq[feet.start_pos1:feet.end_pos2]
                else:
                    gap=configure.input_seq[feet.end_pos1:feet.start_pos2]
                    target_seq=revcomp(configure.input_seq[feet.start_pos1:feet.end_pos2])
                if feet.start_pos2+1==feet.end_pos1+1:
                    gap='0'
                if (whether_match(extra_base[e],extra_base[e2],gap)==0):
                    for h in hairpin.keys():
                        if check_binds(hairpin[h],target_seq)<hair_tar_max and check_binds(revcomp(str(hairpin[h])),target_seq)<hair_tar_max: #not too much hybrid between hairpin and target
                            pt_long=hairpin[h]+extra_base[e]+feet.get_pt
                            px_long=feet.get_px+extra_base[e]+revcomp(str(hairpin[h]))
                            if TmDeltaG.calTm(pt_long,px_long)<lowest_ptx_long_tm: #ptx_command='melt.pl -n DNA -t 55 -N 1 -M 0 -C 0.000000005 '+pt_in+' '+px_in
                                for t in template.keys():
                                    align_tar_template=alignment.needle(revcomp(template[t]), target_seq)
                                    if check_binds(revcomp(t),target) <overhang_tar_max and align_tar_template.find('GAGTC')==-1:      feet_long= PTX_long(feet,e,e2,h)
                                        primer_set= primer_set(feet, e,h,t)
                                        primer_set_list.append(primer_set)
    return primer_set_list
def whether_match(a,b,c):
    comp=['AT','TA','GC','CG']
    for i in range(len(a)):
    for j in range(len(b)):
        for h in range(len(c)):
        if ((a[i]+b[j] in comp) | (a[i]+c[h] in comp) | (b[j]+c[h] in comp)):
            return 1 #there is a match
        else :
            return 0 # there is no match
def check_binds(seq1,seq2):
    t=alignment.needle(seq1,seq2)
    pieces=t.split()
    leng=[]
    for i in pieces:
        leng.append(len(i))
    #match_num=len(t)-t.count(' ')
    match_num=max(leng)
    return match_num
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
