from django.shortcuts import render, HttpResponse
from ModelingTools.FindTemplates import FindTemplates
from ModelingTools.TemplateProfile import TemplateProfile
from ModelingTools.GetDataFromPDB import GetDataFromPDB
from ModelingTools.Align import Align
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
	template_sequence_filename = template_manager.getPDB_File()
	#end#

	#alinhamento inicio
	alignment_manager = Align(workdir + os.sep,workdir + os.sep,  os.path.basename(template_sequence_filename), sequence_file_name)
	alignment_manager.convert_seqali_pir_to_fasta_formar(os.path.basename(sequence_file_name))
	alignment_manager.align_with_muscle()
	alignment_manager.convert_fasta_to_pir()







	return HttpResponse(template_sequence_filename)