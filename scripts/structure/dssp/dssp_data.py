from copy import deepcopy

class DsspPdbSummaryChain():
	"""
	A simple class to manage DSSP sumary data.
	"""
	def __init__(self, sequence, secondary_structure, surface_accessibility, 
		surface_accessibility_normalised, mean_surface_accessibility, 
		buried_residues, accessible_residues,offsets):
		"""
		Creates a new DssPdbSummaryChain.
			:param self:
			:param sequence:
			:param secondary_structure:
			:param surface_accessibility:
			:param surface_accessibility_normalised:
			:param mean_surface_accessibility:
			:param buried_residues:
			:param accessible_residues:
		"""
		self.sequence = sequence
		self.construct = sequence.replace("-","")
		self.secondary_structure = secondary_structure
		self.surface_accessibility = surface_accessibility
		self.surface_accessibility_normalised = surface_accessibility_normalised
		self.mean_surface_accessibility = mean_surface_accessibility
		self.buried_residues = buried_residues
		self.accessible_residues = accessible_residues
		self.offsets = offsets
		self.surface_accessibility_difference = []

	def to_dict(self):
		"""
		Creates a dict representation of this object
			:param self: 
		"""
		return {
			"sequence": self.sequence,
			"secondary_structure": self.secondary_structure,
			"surface_accessibility": self.surface_accessibility,
			"surface_accessibility_normalised": self.surface_accessibility_normalised,
			"mean_surface_accessibility": self.mean_surface_accessibility,
			"buried_residues": self.buried_residues,
			"accessible_residues": self.accessible_residues,
			"surface_accessibility_difference": self.surface_accessibility_difference,
			"offsets":self.offsets,
			"construct":self.construct
		}


class DsspPdbChain():
	"""
	A simple class to manage DSSP data for a single chain
	"""
	def __init__(self, pdb_chain):
		self.positions = {}
		self.pdb_chain = pdb_chain
		self.summary = None

	def add_residue(self, offset, data):
		self.positions[offset] = data
	
	def get(self, offset):
		return self.positions[offset]
		
	def get_by_name(self, offset, attribute_name):
		return self.positions[offset][attribute_name]

	def add_summary(self, chain_summary):
		self.summary = chain_summary

	def get_summary(self):
		return self.summary

	def __len__(self):
		return len(self.positions)

	def __contains__(self, item):
		return item in self.positions

	def update_summary_for_sa_diff(self):
		summary = self.get_summary()
		offsets = self.positions.keys()
		offsets_numeric = [int(x)  for x in offsets if x.isnumeric()]
	
		offset_range = range(min(offsets_numeric), max(offsets_numeric)+1)
		sa_diff_summ = "".join([
			self.positions.get(o,{}).get("aa","-") 
			if self.positions.get(o,{}).get("saDifference",-1)<0 else "x" 
			for o in offset_range])
		summary.surface_accessibility_difference = sa_diff_summ
		
		self.add_summary(summary)

	def copy(self):
		"""
		Creates a deep copy of this object.
		Data from positions are deep copied, data from summary is a reference to the original.
			:param self: 
		"""
		copied = DsspPdbChain(self.pdb_chain)
		copied.positions = deepcopy(self.positions)
		copied.summary = self.summary
		return copied

	def to_dict(self):
		"""
		Creates a dictionary representation of this object.
		Intendend to be used to create a json object.
			:param self: 
		"""
		result = {}
		result["positions"] = self.positions
		if self.summary:
			result["summary"] = self.summary.to_dict()
		return result


class DsspDataPdb():
	"""
	A simple class to manage DSSP data for a protein.
	This class is an iterator.
	"""
	def __init__(self):
		"""
		Creates a new DsspDataPdb object.
			:param self: 
		"""
		self._chains = {}
		self._chain_keys = []
		self._count = 0

	def add_chain(self, chain, aDsspPdbChain):
		"""
		Add data of a chain.
			:param self: 
			:param chain: A string used as identifier.
			:param aDsspPdbChain: 
		"""
		self._chains[chain] = aDsspPdbChain

	def add_chain_summary(self, chain, aDsspPdbSummaryChain):
		"""
		Adds the summary to a existing chain.
			:param self: 
			:param chain: 
			:param aDsspPdbSummaryChain: 
		"""
		if chain in self._chains:
			self._chains[chain].add_summary(aDsspPdbSummaryChain)
	
	def chain(self, chain):
		"""
		Retrieves the DsspPdbChain of a chain.
			:param self: 
			:param chain: 
		"""
		if chain in self._chains:
			return self._chains[chain]
		else:
			return None
	
	def summary(self, chain):
		"""
		Retrives the summary of a chain.
			:param self: 
			:param chain: 
		"""
		if chain in self._chains:
			return self._chains[chain].get_summary()
		else:
			return None

	def solvent_accesibility_difference(self, other):
		"""
		Computes de solvent accessibility difference with another DsspDataPdb.
			:param self: 
			:param other: a DsspDataPdb.
		"""
		result = DsspDataPdb()
		for c in self.chains():
			dssp_pdb_chain = self.chain(c).copy()

			for pos in dssp_pdb_chain.positions:
				self_sa = int(dssp_pdb_chain.get_by_name(pos, "solventAccessibility"))
				other_sa = int(other.chain(c).get_by_name(pos, "solventAccessibility")) if c in other and pos in other.chain(c).positions else self_sa
				# Is some weird cases, the positions of every chain in splitted and joined are not
				# equal, it's assumed that there is no change in solvent accessibility in those 
				# cases.
				# Other weird case, is when a pdb has upper and lower chain ids, and the OS is
				# windows that do not differenciate upper case and lower case in file names
				diff_value = self_sa - other_sa if other_sa else 0
				dssp_pdb_chain.positions[pos].update({"saDifference": diff_value})
			dssp_pdb_chain.update_summary_for_sa_diff()
			result.add_chain(c, dssp_pdb_chain)
		return result

	def chains(self):
		"""
		Retrieves the id of all chains.
			:param self: 
		"""
		return self._chains.keys()

	def __iter__(self):
		self._count = 0
		self.chain_keys = self._chains.keys()
		return self

	def next(self):
		if self._count < len(self.chain_keys):
			self._count += 1
			return self._chains[self.chain_keys[self._count-1]]
		else:
			raise StopIteration

	def items(self):
		"""
		Retrieves an iterable for each chain.
			:param self: 
		"""
		return self._chains.items()

	def __contains__(self, item):
		return item in self._chains

	def to_dict(self):
		"""
		Creates a dictionary representation of this object.
		Intendend to be used to create a json object.
			:param self: 
		"""
		return {c: self.chain(c).to_dict() for c in self.chains()}
