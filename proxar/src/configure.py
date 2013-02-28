import time
progress_simbol=0
input_seq=''


#f1 = NamedTemporaryFile(delete=False,dir='.')
#os.unlink(f1.name)

tempname="~"+str(int(time.time()))
tempname='~1335813823'

pt_foot_file1=tempname+'pt_foot.txt'
px_foot_file1=tempname+'px_foot.txt'
tar_file1=tempname+'tar.txt'
#output of write_ptpx_foot_file_check_foot



pt_bind_filename=tempname+'pt_foot.melt'
px_bind_filename=tempname+'px_foot.melt'
Tm_pt_foot_tar={}
Tm_px_foot_tar={}
#inter result of Tm_pfoot

pt_foot_file2=tempname+'ptok_bind.txt'
px_foot_file2=tempname+'pxok_bind.txt'
tar_file2=tempname+'tarok_bind.txt'
#output of Tm_pfoot


blast_in_file=tempname+'blast_include.xml'
pt_foot_file3=tempname+'ptok_bind_include.txt'
px_foot_file3=tempname+'pxok_bind_include.txt'
tar_file3=tempname+'tarok_bind_include.txt'
tar_include={}
#output of combine_check_include

pt_foot_hair_file1=tempname+'pt_long.txt'
px_foot_hair_file1=tempname+'px_long.txt'
tar_hair_file1=tempname+'tar_long.txt'
#output of bind_hairpin

ptx_melt_file=tempname+'temp_ptx_bind_file.melt'
pt_foot_hair_file2=tempname+'pt_long_thermo_ok1.txt'
px_foot_hair_file2=tempname+'px_long_thermo_ok1.txt'
tar_hair_file2=tempname+'tar_long_thermo_ok1.txt'
Tm_p_foot_hair_tar={}
#output of 2nd self_Tm_check_eachother


whole_PT1=tempname+'wholePT1.txt'
overhang_bind_pt={}
#output of combine_PT

whole_PT2=tempname+'wholePT2.txt'
whole_self_bind={}
#output of 1st self_bind

whole_PX1=tempname+'wholePX1.txt'
#output of 2nd self_bind

pt_blast=tempname+'pt_exlude_blast_result.xml'
px_blast=tempname+'px_exlude_blast_result.xml'
#input of combine_check_exclude

finalptx=tempname+'final_ptx.txt'
#output of combine_check_exclude

#below are advanced settings:
pt_foot_conc=0.0000000005
px_foot_conc=0.000000005
pt_foot_tm=64.5
px_foot_tm=64.5
pt_px_foot_match_max=6
hair_tar_max=5
whole_PTX_tm=13

overhang_pt_tar_max=4
overhang_px_tar_max=4
whole_PT_self_max=10 # need to be add in the advanced page
whole_PX_self_max=10 # need to be add in the advanced page

temper_print='temper_print.txt'
