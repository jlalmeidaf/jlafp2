from django.shortcuts import render, HttpResponse
from ModelingTools.FindTemplates import FindTemplates
from ModelingTools.TemplateProfile import TemplateProfile
from ModelingTools.GetDataFromPDB import GetDataFromPDB
from ModelingTools.Align import Align
from ModelingTools.Modeler import Modeler
from ModelingTools.Malign2 import Malign2
from ModelingTools.MakeProfile import MakeProfile
from ModelingTools.GetProt2 import GetProt2
import tempfile, os
# Create your views here.

def index(request):
    # context = {'latest_question_list': "oie"}
    return render(request, 'alignment/index.html')

def output(request):
	#cria um arquivo da sequencia
	workdir = tempfile.mkdtemp()
	sequence_file_temp = tempfile.mkstemp(dir = workdir)
	sequence_file_name = sequence_file_temp[1]
	sequence_file = file(sequence_file_name, "w")
	sequence_file.write(request.POST["your_name"])
	sequence_file.close()
	#end#

	#encotra o melhor template
	stepOne = FindTemplates(sequence_file_name)
	stepOne.run()
	profile_of_templates = TemplateProfile(stepOne.profilePRF)
	better_profile = profile_of_templates.getBetterProfile()
	#end#

	#pega o template no site do pdb
	template_manager = GetDataFromPDB(workdir, better_profile.name())
	template_pdb_filename = template_manager.getPDB_File()
	#end#

	#alinhamento inicio
	alignment_manager = Align(workdir + os.sep,workdir + os.sep,  os.path.basename(template_pdb_filename), sequence_file_name)
	alignment_manager.convert_seqali_pir_to_fasta_formar(os.path.basename(sequence_file_name))
	alignment_manager.align_with_muscle()
	alignment_manager.convert_fasta_to_pir()
	#end#

	#modelar inicio
	modeling_manager = Modeler(workdir + os.sep , workdir + os.sep, os.path.basename(template_pdb_filename), os.path.basename(alignment_manager.aliali))
	modeling_manager.make_get_model_py()
	modeling_manager.model_sequence()
	best_model = modeling_manager.get_results()

	#avaliar inicio
	multiple_align = Malign2(workdir + os.sep)
	multiple_align.create_script_in_folder(os.path.basename(template_pdb_filename)[0:-4],os.path.basename(best_model)[0:-4])
	multiple_align.get_model()
	profile_best_model = MakeProfile()
	profile_best_model.create_script_in_folder(workdir + os.sep + os.path.basename(best_model))
	profile_best_model.get_model()
	profile_my_template = MakeProfile()
	profile_my_template.create_script_in_folder(workdir + os.sep + os.path.basename(template_pdb_filename))
	profile_my_template.get_model()
	plot_profiles = GetProt2(workdir + os.sep)
	plot_profiles.create_script_in_folder(template_pdb_filename,best_model)
	plot_profiles.get_model()













	return HttpResponse(best_model)