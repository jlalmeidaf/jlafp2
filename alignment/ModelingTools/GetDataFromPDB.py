#data 26 de janeiro de 2016
import urllib2, os

class GetDataFromPDB(object):
	"""docstring for GetDataFromPDB"""
	def __init__(self, path, template_name):
		super(GetDataFromPDB, self).__init__()
		self.path = path
		self.template_name = template_name

	def getPDB_File(self):
		print self.__onlyTemplateName__()
		pdb_file_path = self.path + os.sep + self.__onlyTemplateName__() + ".pdb"
		response = urllib2.urlopen('http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=' + self.__onlyTemplateName__())
		pdb_buffer = response.read()
		pdb_file = file(pdb_file_path, "w")
		pdb_file.write(pdb_buffer)
		pdb_file.close()
		return pdb_file_path
	# http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=1ijk

	def __onlyTemplateName__(self):
		if len(self.template_name) == 4:
			return self.template_name
		else:
			return self.template_name[0:4]