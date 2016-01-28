from Modeller_Caller import modeller_caller
import os
class MakeProfile(object):
    """docstring for MakeProfiles"""
    def __init__(self):
        pass
    
    def create_script_in_folder(self,pdb_file):
        self.path = os.path.dirname(pdb_file)
        script = """from modeller import *
from modeller.scripts import complete_pdb
import os

os.chdir('""" + self.path + """')
log.verbose()    # request verbose output
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# directories for input atom files
env.io.atom_files_directory = './:../atom_files'

# read model file
mdl = complete_pdb(env, '""" + pdb_file + """')

s = selection(mdl)
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='""" + os.path.basename(pdb_file)[0:-4] + """.profile',
              normalize_profile=True, smoothing_window=15)
"""
        script_path = self.path + os.sep + "get_profile.py"
        loop_script = file(script_path, "w")        
        loop_script.write(script)
        loop_script.close()
        self.myscript = script_path

    def get_model(self):
        processo = modeller_caller()
        processo.run(self.myscript)