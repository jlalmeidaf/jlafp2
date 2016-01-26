#antigo mold
import os, sys 
from Modeller_Caller import modeller_caller
class FindTemplates:
    def __init__(self,file_format_pir):
        self.file_format_pir = file_format_pir

    def make_build_profilepy(self):
        script = """\
from modeller import *

log.verbose()
env = environ()

#-- Prepare the input files

#-- Read in the sequence database
sdb = sequence_db(env)
sdb.read(seq_database_file='""" + self.__installfolder__() + """/pdb/pdb_95.pir', seq_database_format='PIR', chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)

#-- Write the sequence database in binary form
sdb.write(seq_database_file='""" + self.__installfolder__() + """/pdb/pdb_95.bin', seq_database_format='BINARY', chains_list='ALL')

#-- Now, read in the binary database
sdb.read(seq_database_file='""" + self.__installfolder__() + """/pdb/pdb_95.bin', seq_database_format='BINARY', chains_list='ALL')

#-- Read in the target sequence/alignment
aln = alignment(env)
aln.append(file='""" + self.file_format_pir + """', alignment_format='PIR', align_codes='ALL')

#-- Convert the input sequence/alignment into
#   profile format
prf = aln.to_profile()

#-- Scan sequence database to pick up homologous sequences
prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, -50), n_prof_iterations=1, check_profile=False, max_aln_evalue=0.01)

#-- Write out the profile in text format
prf.write(file='""" + os.path.dirname(self.file_format_pir) + '/build_profile.prf'"""', profile_format='TEXT')

#-- Convert the profile back to alignment format
aln = prf.to_alignment()

#-- Write out the alignment fileo
aln.write(file='""" + os.path.dirname(self.file_format_pir) + '/build_profilePIR.ali' + """', alignment_format='PIR')
aln.write(file='""" + os.path.dirname(self.file_format_pir) + '/build_profilePAP.ali' + """', alignment_format='PAP')
"""
        arq = file(os.path.dirname(self.file_format_pir) + '/build_profile.py', 'w')
        arq.write(script)
        arq.close()

    def __installfolder__(self):
        '''Returns where the AutoModel was instaled'''
        return os.path.dirname(os.path.realpath(sys.argv[0]))

    def __folder_of_model__(self):
        return(os.path.dirname(self.file_format_pir))

    def run(self):
        self.make_build_profilepy()
        processo = modeller_caller()
        processo.run(self.__folder_of_model__() + '/build_profile.py')