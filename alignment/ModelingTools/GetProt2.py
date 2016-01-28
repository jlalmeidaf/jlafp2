from Modeller_Caller import modeller_caller
import os
class GetProt2(object):
    """docstring for GetProt"""
    def __init__(self, loop_folder):
        self.myscript = ""
        self.path = loop_folder

    def create_script_in_folder(self, template,model):
        template = os.path.basename(template)[0:-4]
        model = os.path.basename(model)[0:-4]
        # myloop = os.path.basename(myloop)[0:-4]
        # self.path = os.path.dirname(model)
        script = """import matplotlib.pyplot as plt
import numpy as np
import modeller
import os

os.chdir('""" + self.path + """')
def get_profile(profile_file, seq):
    '''Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`.'''
    # Read all non-comment and non-blank lines from the file:
    f = file(profile_file)
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            vals.append(float(spl[-1]))
    # Insert gaps into the profile corresponding to those in seq:
    for n, res in enumerate(seq.residues):
        for gap in range(res.get_leading_gaps()):
            vals.insert(n, None)
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals

e = modeller.environ()
a = modeller.alignment(e, file='ali.ali')

try:
	template = get_profile('""" + template + """.profile', a['""" + template + """A'])
	model = get_profile('""" + model + """.profile', a['""" + model + """'])
except Exception, e:
		
	template = get_profile('""" + template + """.profile', a['""" + template + """'])
	model = get_profile('""" + model + """.profile', a['seq.ali'])

# Plot the template and model profiles in the same plot for comparison:
fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.xlabel('Alignment position')
plt.ylabel('DOPE per-residue score')
ax.plot(model, color='red', linewidth=2, label='Model')
ax.plot(template, color='green', linewidth=2, label='Template')
handles, labels = ax.get_legend_handles_labels()
lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1))
ax.grid('on')
fig.savefig('dope_profile_loop', bbox_extra_artists=(lgd,), bbox_inches='tight')
"""
        script_path = self.path + os.sep + "plot_profile.py"
        loop_script = file(script_path, "w")        
        loop_script.write(script)
        loop_script.close()
        self.myscript = script_path

    def get_model(self):
        processo = modeller_caller()
        processo.run(self.myscript)
