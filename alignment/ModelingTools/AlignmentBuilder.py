import os

from Bio import SearchIO
from Bio import AlignIO
from Bio.PDB import PDBParser
from PirFormat import PirFormat
from modeller import *
from AlignFile import *


class AlignmentBuilder():
	def __init__(self, modeldir, templatePDBFile, alignmentfile, template_pir):
		# self.hmmr3file = hmmr3file
		self.workdir = modeldir
		self.templatePDBFile = templatePDBFile
		self.alignmentfile = alignmentfile
		self.fastaAlignment = self.__get_fasta_alignment_str__()
		self.templatePIR = AlignFile(template_pir)

	def __get_fasta_alignment_str__(self):
		fastaAlignmentFile = AlignIO.read(self.alignmentfile, "fasta") 
		fastaAlignmentList = []
		for record in fastaAlignmentFile:
			fastaAlignmentList.append(record.seq)
		return fastaAlignmentList

	def __get_target_sequence__(self):
		structureName = "seq.ali"
		structureType = "sequence"
		description = ""
		start = ""
		chain_start = ""
		stop = ""
		chain_stop = ""
		sequence = str(self.fastaAlignment[1]).upper() + self.templatePIR.get_all_heteroatoms(0)
		return PirFormat(structureName, structureType, description, start, chain_start, stop, chain_stop, sequence).get()

	def __get_template_sequence__(self):
		#env = environ()
		#mdl = model(env, file=self.templatePDBFile)
		structureName = os.path.basename(self.templatePDBFile).partition(".")[0] #melhorar isso
		parser = PDBParser()
		structure = parser.get_structure(structureName, self.workdir + "/" + self.templatePDBFile)
		model = structure[0] #get model from structure
		listOfChains = model.get_list()
		# structureType = mdl.prottyp
		# description = mdl.name
		# start = mdl.range[0].partition(":")[0]
		#chain_start = mdl.range[0].partition(":")[2]
		# stop = mdl.range[0].partition(":")[0]
		#chain_stop = mdl.range[1].partition(":")[2]
		# structureType = mdl.prottyp
		structureType = "StructureX"
		description = ""
		start = "FIRST"
		chain_start = listOfChains[0].id
		stop = "END"
		chain_stop = listOfChains[-1].id
		sequence = str(self.fastaAlignment[0]).upper() + self.templatePIR.get_all_heteroatoms(0)
		return PirFormat(structureName, structureType, description, start, chain_start, stop, chain_stop, sequence).get()

	def __get_all_heteroatoms__(self, template_pir):
		template_PirFormat = AlignFile(template_pir)
		return template_PirFormat.get_all_heteroatoms()
		

# if __name__ == '__main__':
# 	x = AlignmentBuilder("1bdm.pdb", "aligned.m", 'template.pir')
# 	alingFile = file('ali.ali','w')
# 	alingFile.write(x.__get_template_sequence__() + "\n\n") 
# 	alingFile.write(x.__get_target_sequence__())
# 	alingFile.close()
	# print x.templatePIR.get_all_heteroatoms(0)