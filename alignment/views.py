from django.shortcuts import render, HttpResponse
from ModelingTools.FindTemplates import FindTemplates
from ModelingTools.TemplateProfile import TemplateProfile
from ModelingTools.GetDataFromPDB import GetDataFromPDB
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

	template_manager = GetDataFromPDB(os.path.dirname(sequence_file_name), better_profile.name())
	template_sequence_file = template_manager.getPDB_File()


	return HttpResponse(template_sequence_file)