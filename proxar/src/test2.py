'''
Created on 2013-2-8

@author: jifeng
'''
from Bio.Blast import NCBIWWW, NCBIXML 
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import Entrez, SeqIO
def blast_short(input_filename):
    print "enter"
    #(input_id,input_seq,input_len)=read_fasta(input_filename)
    input_seq=open(input_filename).read()
    if not input_filename: return
    print '\nStarting a BLAST search. This may take awhile...'
    try:
        blast_results = NCBIWWW.qblast('blastn', 'nr',input_seq)
        #save results to a file
        results_filename = '1'+'-blast.xml'
        results_file = open(results_filename, 'w')
        results_file.write(blast_results.read())
        results_file.close()
        blast_results.close()
        print '\nBlast output was written to:\n   '+results_filename
        results_file  = open(results_filename, 'r')
        blast_results = list(NCBIXML.parse(results_file))
        results_file.close()
    except Exception, e:
        print '\nFailed to obtain BLAST search results from NCBI.'
        print_exception(e)
        return
def print_exception(e):
    print '\nException occurred: %s\n%s\n' % (str(type(e)), e.message)

def read_fasta(filename):
    for seq_record in SeqIO.parse(filename, "fasta"):
        return (str(seq_record.id), str(seq_record.seq),len(str(seq_record.seq)))


if __name__ == '__main__':
    blast_short('test_seq.txt')