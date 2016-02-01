from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    # url(r'output/', views.output, name='output'),
    url(r'^find_templates/$', views.find_templates, name='find_templates'),
    url(r'^alignment2/$', views.alignment2, name='alignment2'),
]