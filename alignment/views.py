from django.shortcuts import render
# Create your views here.

def index(request):
    context = {'latest_question_list': "oie"}
    return render(request, 'alignment/index.html', context)

def joao(request):
	return HttpResponse("Oi")