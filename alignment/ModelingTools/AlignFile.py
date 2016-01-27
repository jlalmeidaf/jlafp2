from modellingfile.PirSequence import *
import os

class AlignFile():
	"""docstring for AlignFile"""
	def __init__(self, align_file):
		# super(AlignFile, self).__init__(align_file, mode)
		self.aliali_file = file(align_file, "r")
		self.tmp_align_file = file(align_file + ".tmp", "w")
		self.selected_sequence = 1
		self.my_ali_file = []
		self.__load_file__()

	def __load_file__(self):
		align_file = self.aliali_file.read()[1:] #fix the error for two return in alignment gerated by modeller
		my_sequence = align_file.split("\n>")
		my_sequence = self.__fix_sequence__(my_sequence)
		self.__create_sequence_list__(my_sequence)

	def __create_sequence_list__(self, list_of_sequence):
		for selected_sequence in list_of_sequence:
			if selected_sequence.endswith("\n"):
				selected_sequence = selected_sequence[:-1]
			pirsequence = PirSequence(selected_sequence)
			self.my_ali_file.append(pirsequence)

	def __fix_sequence__(self, list_of_sequence):
		for i in range(0,len(list_of_sequence)):
			if list_of_sequence[i].startswith("P1"): #find sequences without >
				list_of_sequence[i] = ">" + list_of_sequence[i]
		if list_of_sequence[-1][-1] == "\n":
			list_of_sequence[-1] = list_of_sequence[-1][:-1] #remove \n of last sequence
		return list_of_sequence

	def number_of_sequences(self):
		return len(self.my_ali_file)

	def select_sequence(self, number_of_sequence):
		self.selected_sequence = number_of_sequence - 1

	def show_heteroatoms_of_chain(self, chain):
		return self.my_ali_file[self.selected_sequence].list_heteroatoms(chain)

	def change_heteroatom_of_chain(self, chain, old_heteroatom, new_heteroatom):
		self.my_ali_file[self.selected_sequence].change_heteroatom(chain, old_heteroatom, new_heteroatom)

	def copy_heteroatoms(self,sequence_source,sequence_destination):
		number_of_returns = 1
		sequence_source_content = self.my_ali_file[sequence_source].sequence()
		for posChar in range(0,len(sequence_source_content)):
		
			if not sequence_source_content[posChar].isupper() and not sequence_source_content[posChar] == "-":
				self.my_ali_file[sequence_destination].change_het_in_position(sequence_source_content[posChar], posChar + number_of_returns)

		# for posChar in range(0,len(sequence_source_content)):
		# 	if sequence_source_content[posChar] == "\n":
		# 		number_of_returns = number_of_returns + 1
		# 	elif sequence_source_content[posChar].islower():
		# 		position = posChar + number_of_returns 
		# 		self.my_ali_file[sequence_destination].change_het_in_position(sequence_source_content[posChar], position)


	def write_changes(self):
		for selectedPirsequence in self.my_ali_file:
			# print selectedPirsequence.user_sequence
			self.tmp_align_file.write(selectedPirsequence.user_sequence)
			if selectedPirsequence != self.my_ali_file[-1]: #don't wrap last line
				self.tmp_align_file.write("\n")
		self.__replace_files__()

	def __replace_files__(self):
		name_of_aliali_file = self.aliali_file.name
		name_of_tmp_file = self.tmp_align_file.name

		self.close()

		os.remove(name_of_aliali_file)
		os.rename(name_of_tmp_file, name_of_aliali_file)


		self.aliali_file = file(name_of_aliali_file, "r")
		self.tmp_align_file = file(name_of_tmp_file, "w")

	def close(self):
		self.aliali_file.close()
		self.tmp_align_file.close()

	def get_all_heteroatoms(self, sequence_source):
		my_hets = ""
		sequence_source_content = self.my_ali_file[sequence_source].sequence()
		for posChar in range(0,len(sequence_source_content)):
		
			if not sequence_source_content[posChar].isupper() and not sequence_source_content[posChar] == "-" and not sequence_source_content[posChar] == "\n" and not sequence_source_content[posChar] == "*":
				my_hets = my_hets + sequence_source_content[posChar]

		return my_hets
