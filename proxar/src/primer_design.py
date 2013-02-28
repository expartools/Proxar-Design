'''
#for proxar primer design
Created on Jan 24, 2013

@author: jifeng
'''
#!/usr/bin/python -d
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import TmDeltaG
import re
from  string  import *
import File_jf
import design
import configure



def main_page():
    configure.input_seq=File_jf.check_file(target_file, hairpin_file, extra_base_file,templates_file)
    feet_list=design.ptpx_foot(configure.input_seq,min_pt,max_pt,min_px,max_px,max_e,direction,pt_min_tm,pt_max_tm, px_min_tm, px_max_tm,configure.pt_foot_conc,configure.px_foot_conc)
    design.write_fasta_feet(feet_list,'~feet.fasta')
    (in_blast_result,ex_blast_result)=design.combine_check_exinclude('~feet.fasta',str(include_line),str(exclude_line),int(num),float(perct)/100,seq_name_list)
    blast_ok_feet_list=design.blast_feet_chek(in_blast_result,ex_blast_result, feet_list,pt_min_length, px_min_length)
    primer_set=design.bind_hairpin_template(blast_ok_feet_list,hairpin_file,extra_base_file)



if __name__ == '__main__':
    main_page()
        
