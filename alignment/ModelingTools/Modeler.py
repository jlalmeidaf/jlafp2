import pdb, os
from Modeller_Caller import modeller_caller
import sequence
class Modeler:
    def __init__(self,pdb_folder,modeldir,pdb_file,ali_ali):
        self.ali_ali = modeldir + ali_ali
        self.pdb_folder = pdb_folder
        self.modeldir = modeldir
        self.seq = sequence.Sequence(self.ali_ali)
        print self.seq.structures_names()
        self.structures = str(self.seq.structures_names()[0:-1])[1:-1]

    def make_get_model_py(self):
        script = """\
from modeller import *
from modeller.automodel import *    # Load the automodel class
import os
os.chdir('""" + self.__folder_of_model__() + """')

log.verbose()

env = environ(rand_seed=-12312)  # To get different models from another script
# directories for input atom files
env.io.hetatm = env.io.water = True


env.io.atom_files_directory = ['""" + self.modeldir + """']
a = automodel(env,
              alnfile='""" + self.ali_ali + """',      # alignment filename (ali.ali)
              knowns=(""" + self.structures + """),             # codes of the templates (variavel da est)
              sequence='seq.ali',           # code of the target (seq.ali)
              assess_methods=assess.GA341)  # request GA341 assessment
a.starting_model= 1                 # index of the first model
a.ending_model  = 5                 # index of the last model
                                    # (determines how many models to calculate)
a.deviation = 4.0                   # has to >0 if more than 1 model

a.make()                            # do homology modelling"""
        arq = file(self.modeldir + 'get_model.py', 'w')
        arq.write(script)
        arq.close()

    def __folder_of_model__(self):
        return self.modeldir

    def model_sequence(self):
        processo = modeller_caller()
        processo.run(self.__folder_of_model__() + 'get_model.py')

    def get_results(self):
        folder = self.modeldir
        mof = 0
        filemof = ""
        for files in os.listdir(folder):
            if (files.startswith('seq.ali.B')):
                    arqr = file(folder + '/' + files, 'r')
                    arqr.readline()
                    linha = arqr.readline()
                    if (mof == 0) or (linha.partition(':      ')[2] < mof):
                            mof = linha.partition(':      ')[2]
                            filemof = folder + '/' + files
                            arqr.close()
        # if not (os.path.exists(folder + '/../results')):
        #         os.mkdir(folder + '/../results')
        # os.popen('cp '+ filemof +  ' ' + folder + '/../results/')
        return filemof