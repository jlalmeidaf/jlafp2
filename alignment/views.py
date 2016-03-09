from django.shortcuts import render, HttpResponse, redirect
from ModelingTools.FindTemplates import FindTemplates
from ModelingTools.TemplateProfile import TemplateProfile
from ModelingTools.GetDataFromPDB import GetDataFromPDB
# from ModelingTools.Align import Align
from ModelingTools.GPUAlign import GPUAlign
from ModelingTools.Modeler import Modeler
from ModelingTools.Malign2 import Malign2
from ModelingTools.MakeProfile import MakeProfile
from ModelingTools.GetProt2 import GetProt2
from ModelingTools.PDB2 import pdb

from django.core.urlresolvers import reverse
from django.http import JsonResponse

import tempfile, os
import time

workdir = None
template_manager = None
# Create your views here.

def index(request):
    # context = {'latest_question_list': "oie"}
    return render(request, 'alignment/index.html')


def find_templates(request):
	#cria um arquivo da sequencia
	# if not request.POST.has_key("your_name"):
	workdir = tempfile.mkdtemp()
	sequence_file_temp = tempfile.mkstemp(dir = workdir)
	sequence_file_name = sequence_file_temp[1]
	sequence_file = file(sequence_file_name, "w")
	sequence_file.write(">P1;seq\nsequence:seq:::::::0.00: 0.00\n")
	sequence_file.write("\n" +request.POST["your_name"])
	sequence_file.write("*")
	sequence_file.close()
	#end#

	# return render(request, 'alignment/find_template_wait.html')


	#encotra o melhor template
	stepOne = FindTemplates(sequence_file_name)
	stepOne.run()
	profile_of_templates = TemplateProfile(stepOne.profilePRF)
	better_profile = profile_of_templates.getBetterProfile()

	input_sequence = profile_of_templates.list_of_sequences[0]
	#end#

	#pega o template no site do pdb
	template_manager = GetDataFromPDB(workdir, better_profile.name())
	template_pdb_filename = template_manager.getPDB_File()
	# print better_profile.sequence()
	downloaded_pdb = pdb(template_pdb_filename) 
	new_pdb = template_manager.getPDB_File() + "m"
	modified_pdb = file(new_pdb, "w")
	print downloaded_pdb.HetatomsInPDB()
	modified_pdb.write(downloaded_pdb.const({"A":[downloaded_pdb.HetatomsInPDB()]}, ["A"]))
	modified_pdb.close()
	os.rename(new_pdb, template_manager.getPDB_File())
	#end#
	context = {'better_template': better_profile.name(),

	'better_profile_sequence' : (better_profile.sequence()).strip("\n"),
	'better_template_identity' : better_profile.identity(),
	'input_sequence_name' : input_sequence.name(),
	'input_sequence_sequence' : input_sequence.sequence()
	}
	request.session['workdir'] = workdir
	request.session['template_manager'] = template_manager
	request.session['sequence_file_name'] = sequence_file_name
	request.session['template_pdb_filename'] = template_pdb_filename
	reponse = render(request, 'alignment/find_template.html', context)
	return HttpResponse(reponse)


def alignment2(request):
	template_manager = request.session['template_manager']
	workdir = request.session['workdir']
	sequence_file_name = request.session['sequence_file_name']
	# return HttpResponse(x)
	#alinhamento inicio

	template_pdb_filename = request.session['template_pdb_filename']
	alignment_manager = GPUAlign(workdir + os.sep,workdir + os.sep,  os.path.basename(template_pdb_filename), sequence_file_name)
	alignment_manager.convert_seqali_pir_to_fasta_formar(os.path.basename(sequence_file_name))
	alignment_manager.align_with_GPU()
	alignment_manager.convert_fasta_to_pir()
	request.session["alignment_manager"] = alignment_manager
	request.session["template_pdb_filename"] = template_pdb_filename

	file_ = file(alignment_manager.aliali,'r')
	text_ = file_.read()

	context = {'alignment': text_}
	# #end#
	return render(request, 'alignment/alignment.html', context)

def modeling(request):
		#modelar inicio
	workdir = request.session['workdir']
	template_manager = request.session['template_manager']
	template_pdb_filename = request.session['template_pdb_filename']
	alignment_manager = request.session['alignment_manager']
	modeling_manager = Modeler(workdir + os.sep , workdir + os.sep, os.path.basename(template_pdb_filename), os.path.basename(alignment_manager.aliali))
	modeling_manager.make_get_model_py()
	modeling_manager.model_sequence()
	best_model = modeling_manager.get_results()

	file_ = file(best_model,'r')
	text_ = "\n".join(file_.readlines()[1:])
	context = {"modeling" : text_}
	return render(request, 'alignment/model.html', context)

# def output(request):




def conta(request):
    # c.prova(0)
    time.sleep(10)
    redirect = reverse('find_templates')
    return JsonResponse({'redirect': redirect})





# 	#modelar inicio
# 	modeling_manager = Modeler(workdir + os.sep , workdir + os.sep, os.path.basename(template_pdb_filename), os.path.basename(alignment_manager.aliali))
# 	modeling_manager.make_get_model_py()
# 	modeling_manager.model_sequence()
# 	best_model = modeling_manager.get_results()

# 	#avaliar inicio
# 	multiple_align = Malign2(workdir + os.sep)
# 	multiple_align.create_script_in_folder(os.path.basename(template_pdb_filename)[0:-4],os.path.basename(best_model)[0:-4])
# 	multiple_align.get_model()
# 	profile_best_model = MakeProfile()
# 	profile_best_model.create_script_in_folder(workdir + os.sep + os.path.basename(best_model))
# 	profile_best_model.get_model()
# 	profile_my_template = MakeProfile()
# 	profile_my_template.create_script_in_folder(workdir + os.sep + os.path.basename(template_pdb_filename))
# 	profile_my_template.get_model()
# 	plot_profiles = GetProt2(workdir + os.sep)
# 	plot_profiles.create_script_in_folder(template_pdb_filename,best_model)
# 	plot_profiles.get_model()













# 	return HttpResponse(best_model)