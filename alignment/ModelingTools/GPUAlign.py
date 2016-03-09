from Align import Align
from RunCommand import RunCommand
# from PDB2 import pdb
import os
class GPUAlign(Align):
	"""docstring for GPUAlign"""
	def __init__(self,pdb_folder,modeldir,pdb_file,seq_ali):
		self.gpas_aligned_file = ""
		super(GPUAlign, self).__init__(pdb_folder,modeldir,pdb_file,seq_ali)
		# self.aliali = super.aliali
		

	def align_with_GPU(self):
		result = self.modeldir + "out.ali"
		print(os.getcwd()) 
		print self.unaligned_seq
		print result
		gpas = RunCommand(["./gpas/gpas -sm gpas/data/substitution_matrices/BLOSUM62.txt -a nwg -i ", self.unaligned_seq ," -o ", result])
		gpas.run()

		self.gpas_aligned_file = result
		self.convert_gpas_to_fasta()

	def convert_gpas_to_fasta(self):
		gpas_result = file(self.gpas_aligned_file, "r")
		gpas_result_readline = gpas_result.readlines()
		first_seq = gpas_result_readline[10]
		first_seq_name = gpas_result_readline[2].split(" ")[3]
		second_seq = gpas_result_readline[11]
		second_seq_name = gpas_result_readline[3].split(" ")[3]

		fasta_name = self.gpas_aligned_file + ".fasta"
		fasta_file = file(fasta_name, "w")
		fasta_file.write(">" + first_seq_name + "\n")
		fasta_file.write(first_seq)
		fasta_file.write(">" + second_seq_name + "\n")
		fasta_file.write(second_seq)
		self.fasta_aligned_file = fasta_name
		# print fasta_name
		self.seq_ali = fasta_name
		# self.aliali = fasta_name
		return fasta_name


