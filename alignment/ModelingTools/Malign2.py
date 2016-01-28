from Modeller_Caller import modeller_caller
import os
class Malign2(object):
    """docstring for Malign"""
    def __init__(self, path):
        self.path = path
        self.myscript = ""

    def create_script_in_folder(self, template, best_sequence):
        script = """from modeller import *
import os

os.chdir('""" + self.path + """')
log.verbose()
env = environ()

env.io.atom_files_directory = './:../atom_files/'

aln = alignment(env)
for (code, chain) in (('""" + template + """', 'A'), ('""" + best_sequence + """', '')):
    mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(mdl, atom_files=code, align_codes=code+chain)

for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50),
               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
               dendrogram_file='ali.tree',
               alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
               feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
               improve_alignment=True, fit=True, write_fit=write_fit,
               write_whole_pdb=whole, output='ALIGNMENT QUALITY')

#aln.write(file='ali.pap', alignment_format='PAP')
aln.write(file='ali.ali', alignment_format='PIR')

aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
           rr_file='$(LIB)/as1.sim.mat', overhang=30,
           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
           gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
           alignment_type='progressive', feature_weights=[0]*6,
           improve_alignment=False, fit=False, write_fit=True,
           write_whole_pdb=False, output='QUALITY')
"""
        script_path = self.path + os.sep + "salign.py"
        loop_script = file(script_path, "w")        
        loop_script.write(script)
        loop_script.close()
        self.myscript = script_path

    def get_model(self):
        processo = modeller_caller()
        processo.run(self.myscript)