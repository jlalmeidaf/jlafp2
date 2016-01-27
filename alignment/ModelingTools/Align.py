import pdb, os
from Modeller_Caller import modeller_caller
from RunCommand import RunCommand
from Bio import SearchIO
from AlignmentBuilder import AlignmentBuilder
class Align:
    def __init__(self,pdb_folder,modeldir,pdb_file,seq_ali):
        self.seq_ali = modeldir + seq_ali
        self.pdb_folder = pdb_folder
        self.modeldir = modeldir
        self.pdb_file = pdb.pdb(modeldir + pdb_file)
        self.template = modeldir + pdb_file
        self.strname = self.pdb_file.NomedaEstrutura()

        # self.hmmresults = None
        # self.betterid = None
        # self.hmm = None
        self.fasta_aligned_file = None
        self.template_pir = None


    def convert_seqali_pir_to_fasta_formar(self, file_name):
        script = """\
from modeller import *
import os
os.chdir('""" + self.__folder_of_model__() + """')
env = environ()
env.io.hetatm = env.io.water = False
env.io.atom_files_directory = ['""" + self.__folder_of_model__() + """']
code = '""" + self.strname + """'   #   estrutura a ser lida#
mdl = model(env, file=code)

aln = alignment(env)
aln.append_model(mdl, align_codes=code)
aln.write(file='""" + self.__folder_of_model__() + """/template.pir', alignment_format = 'PIR')
aln.write(file='""" + self.__folder_of_model__() + """/template.fasta', alignment_format = 'FASTA')

aln2 = alignment(env)
aln2.append(file='""" + self.__folder_of_model__() + file_name + """', alignment_format = 'PIR')
aln2.write(file='""" + self.__folder_of_model__() + """/seq.fasta', alignment_format='FASTA')

unalinedfasta = file('""" + self.__folder_of_model__() + """/unaligned.fasta', 'w')
templatefasta = file('""" + self.__folder_of_model__() + """/template.fasta', 'r')
seqfasta = file('""" + self.__folder_of_model__() + """/seq.fasta', 'r')

unalinedlist = templatefasta.readlines()
unalinedlist = unalinedlist + seqfasta.readlines()

unalinedfasta.writelines("".join(unalinedlist))
unalinedfasta.close()
"""
        arq = file(self.modeldir + 'get_fasta.py', 'w')
        arq.write(script)
        arq.close()

        self.unaligned_seq = self.__folder_of_model__() + '/unaligned.fasta'
        # self.template_pir = self.__folder_of_model__() + '/template.pir'
        processo = modeller_caller()
        processo.run(self.__folder_of_model__() + 'get_fasta.py')

        self.__get_template_sequence_in_pir_format__()


    def __get_template_sequence_in_pir_format__(self):
        script = """\
from modeller import *
import os
os.chdir('""" + self.__folder_of_model__() + """')
env = environ()
env.io.hetatm = env.io.water = True
env.io.atom_files_directory = './:../atom_files'

code = '""" + self.strname + """'   #   estrutura a ser lida#
mdl = model(env, file=code)

aln = alignment(env)
aln.append_model(mdl, align_codes=code)
aln.write(file='""" + self.__folder_of_model__() + """/str.seq', alignment_format = 'PIR')
"""
        arq = file(self.modeldir + 'get_template.py', 'w')
        arq.write(script)
        arq.close()

        # self.unaligned_seq = self.__folder_of_model__() + '/unaligned.fasta'
        self.template_pir = self.__folder_of_model__() + '/str.seq'
        processo = modeller_caller()
        processo.run(self.__folder_of_model__() + 'get_template.py')

    # def search_an_hmmm(self):
    #       seqfasta = self.modeldir + "/seq.fasta"
    #       result =  self.modeldir + "/hmmsearch.txt"
    #       hmmer = RunCommand(["hmmsearch",self.pfam_database, seqfasta, ">", result])
    #       hmmer.run()
    #       # print "foi"
    #       # seqfasta = modeldir + "/seq.fasta"
    #       # result =  modeldir + "/hmmsearch.txt"
    #       # os.popen("hmmsearch " + self.pfam_database + " " + seqfasta + " > " + result)
    #       self.hmmresults = result

    # def find_better_motif(self):
    #   better_result = None
    #   for qresults in SearchIO.parse(self.hmmresults,'hmmer3-text'):
    #       if(len(qresults.hits)>0 and better_result != None):
    #         print qresults.hits
    #         print qresults.seq_len
    #         print better_result.seq_len
    #         if qresults.seq_len > better_result.seq_len:
    #           print qresults.seq_len
    #           better_result = qresults
    #       else:
    #         if len(qresults.hits)>0:
    #           better_result = qresults
    #   print better_result.id
    #   self.betterid =  better_result.id



    def align_with_muscle(self):
      result = self.modeldir + "seq.ali"
      muscle = RunCommand(["muscle -in ", self.unaligned_seq ," -out ", result])
      muscle.run()
      self.fasta_aligned_file = result

    def convert_fasta_to_pir(self):
      # print self.template
      # print self.fasta_aligned_file
      # print self.template_pir
      try:
        x = AlignmentBuilder(self.modeldir, os.path.basename(self.template), self.fasta_aligned_file,  self.template_pir)
        print "---"
        aliali = self.__folder_of_model__() + 'ali.ali'
        alingFile = file(aliali,'w')
        alingFile.write(x.__get_template_sequence__() + "\n\n") 
        alingFile.write(x.__get_target_sequence__())
        alingFile.close()
        self.aliali = aliali
      except Exception, e:
        print e

    def __folder_of_model__(self):
        return self.modeldir