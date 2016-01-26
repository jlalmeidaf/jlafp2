from django.shortcuts import render, HttpResponse
from ModelingTools.FindTemplates import FindTemplates
import tempfile
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

	stepOne = FindTemplates(sequence_file_name)
	stepOne.run()


	return HttpResponse(sequence_file_temp[1])