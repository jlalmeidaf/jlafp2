class PirFormat(object):
	"""docstring for PirFormat"""

	def __init__(self, name, type, description, start, chain_start, stop, chain_stop, sequence):
		self.name = name
		self.type = type
		self.description = description
		self.start = start
		self.chain_start = chain_start
		self.stop = stop
		self.chain_stop = chain_stop
		self.sequence = sequence

	def get(self):
		return ">P1;" + self.name + "\n" + self.type + ":" + self.name + ":" + self.start + ":" \
		 + self.chain_start + ":" + self.stop + ":" + self.chain_stop + ":" + ":" \
		 + self.description + ":" + ":" + "\n" + self.sequence + "*"
		