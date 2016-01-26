from ProfileSequence import ProfileSequence

class TemplateProfile(object):
	"""docstring for TemplateProfile"""
	def __init__(self, profile_filename):
		super(TemplateProfile, self).__init__()
		self.profile_filename = profile_filename
		self.list_of_sequences = self.__get_sequences__()

	def __get_sequences__(self):
		profile_file = file(self.profile_filename, "r")
		lines = profile_file.readlines()[6:]
		list_of_sequences = []
		for eachLine in lines:
			print eachLine
			list_of_sequences.append(ProfileSequence(eachLine))
		return list_of_sequences

	def getBetterProfile(self):
		better_profile = self.list_of_sequences[0]
		for eachProfile in self.list_of_sequences:
			if eachProfile.identity() > better_profile.identity():
				better_profile = eachProfile
		return better_profile