import string
from structure.dssp.dssp_data import DsspDataPdb
from structure.dssp.dssp_data import DsspPdbChain
from structure.dssp.dssp_data import DsspPdbSummaryChain


#-----
import logging
logger = logging.getLogger(__name__)
#-----

class Parser(object):
	def __init__(self):
		self.setup_parsing_parameters()

	def setup_parsing_parameters(self):
		self.dssp_headers = {"count": [0, 4],
							"offset": [5, 10],
							"chain": [11, 11],
							"aa": [12, 14],
							"SecondaryStructure": [16],
							"3-turns/helix": [18],
							"4-turns/helix": [19],
							"5-turns/helix": [20],
							"geometricalBend": [21],
							"chirality": [22],
							"betaBridgeLabel1": [23],
							"betaBridgeLabel2": [24],
							"betaBridgePartnerResnum": [25, 28],
							"betaSheetLabel": [29, 32],
							"solventAccessibility": [34, 37],
							"N-H-->O_1": [42, 49],
							"O-->H-N_1": [52, 60],
							"N-H-->O_2": [63, 71],
							"O-->H-N_2": [74, 82],
							"TCO": [84, 90],
							"KAPPA": [91, 96],
							"ALPHA": [97, 102],
							"PHI": [103, 108],
							"PSI": [109, 114],
							"X-CA": [115, 121],
							"Y-CA": [122, 128],
							"Z-CA": [129, 135]}

		# TODO: Locate the source paper of these values
		self.dssp_maximum_accessibility = {
			"A": 124, "B": 157.5, "C": 94, "D": 154, "E": 187, "F": 221,
			"G": 89, "H": 201, "I": 193, "K": 214, "L": 199, "M": 216,
			"N": 161, "P": 130, "Q": 192, "R": 244, "S": 113, "T": 151,
			"V": 169, "W": 264, "Y": 237, "X": 179, "Z": 189.5}

		self.default_sa_normalised = 0.5
		self.cutoff_sa_normalised = 0.25

	def read_string_data(self, dssp_input, chains, pdb_id, use_all_chains=False):
		"""
		Reads the output of a dssp and retrieve the plain data without any post-
		proccess. Input is assumed to be a string.
			:param self:
			:param dssp_input: the result of dssp as string
			:param chains: The chains which information is required to be extracted
			:param pdb_id: The pdb_id that corresponds to the dssp file.
		"""
		dssp = DsspDataPdb()
		atData = False
		
		try:
			dssp_input = dssp_input.strip().decode('utf-8')
			
			for line in dssp_input.strip().split("\n"):
			
				if atData:
					dssp_fields = {}
					for header in self.dssp_headers:
						start = self.dssp_headers[header][0]
						stop = self.dssp_headers[header][-1] + 1
						dssp_fields[header] = line[start:stop].strip()

					if dssp_fields["chain"]:
						pdb_chain = "{}.{}".format(pdb_id, dssp_fields["chain"])
						if use_all_chains or pdb_chain in chains:
							if pdb_chain not in dssp:
								dssp.add_chain(pdb_chain, DsspPdbChain(pdb_chain))
							c_offset = dssp_fields["offset"]
							#if c_offset.isdigit():
							dssp.chain(pdb_chain).add_residue(c_offset, dssp_fields)

				atData = line.strip().startswith("#") or atData
		except:
			logger.error("Error reading DSSP:" + pdb_id)

		return dssp    

	# TODO: make a variant of the function that assume good values for
	#       chains and pdb_id from dssp file content if possible.
	def read_plain_data(self, dssp_path, chains, pdb_id, use_all_chains=False):
		"""
		Reads the output of a dssp and retrieve the plain data without any post-
		proccess.
			:param self:
			:param dssp_path: the path of the input file
			:param chains: The chains which information is required to be extracted
			:param pdb_id: The pdb_id that corresponds to the dssp file.
		"""
		return self.read_string_data(open(dssp_path).read(),chains, pdb_id, use_all_chains)

	def _get_statistics(self, dssp):
		"""
		Get 'statistics' from a DsspDataPdb object.
		the statistic object is only internally to calculate some values for the
		summary.
			:param self:
			:param dssp: a DsspDataPdb object .
		"""
		statistics={pdb_chain:{} for pdb_chain in dssp.chains()}
		for pdb_chain in dssp.chains():
			for header in self.dssp_headers:
				dataDict = {}
				for offset in dssp.chain(pdb_chain).positions:
					dataDict[offset] = dssp.chain(pdb_chain).get(offset)[header]
				statistics[pdb_chain][header] = dataDict
		return statistics

	def _normalize_sa(self, sa, aminoacid):
		try:
			return min(1, sa /
				self.dssp_maximum_accessibility[aminoacid.upper()])
		except:
			return self.default_sa_normalised

	def _append_summary(self, dssp):
		"""
		Append 'summary' dict data to a DsspDataPdb object.
			:param self:
			:param dssp: a dict DsspDataPdb object.
		"""
		def is_buried_from_sa(sa_normalised):
			def f(offset):
				return min(1, sa_normalised[offset]) < self.cutoff_sa_normalised
			return f

		for pdb_chain in dssp.chains():
			sequence = ""
			buried_residues = ""
			accessible_residues = ""
			ss = ""
			sa = {}
			sa_normalised = {}
			statistics = self._get_statistics(dssp)

			if pdb_chain in statistics:
				stats = statistics[pdb_chain]
				offsets = stats["aa"]

				offsets_numeric = [int(x)  for x in offsets.keys() if x.isnumeric()]
				offsets_range = [str(o) for o in range(min(offsets_numeric), max(offsets_numeric) + 1)]

				sequence = "".join([ offsets[str(o)] if str(o) in offsets else "-" for o in offsets_range])
			
				sec_structures = stats.get("SecondaryStructure", {})
				ss = "".join([sec_structures[str(o)]
					if str(o) in sec_structures and sec_structures[o] != "" else "-"
					for o in offsets_range])

				solventAccessibilities = stats.get("solventAccessibility", {})

				sa = {o: float(solventAccessibilities[o])
					for o in offsets if o in solventAccessibilities}

				sa_normalised = {o: self._normalize_sa(sa[o], offsets[o])
					for o in offsets if o in solventAccessibilities}

				is_buried = is_buried_from_sa(sa_normalised)

				buried_residues = [(offsets[str(o)] if is_buried(str(o)) else "x")
					if str(o) in solventAccessibilities else "-"
					for o in offsets_range]

				accessible_residues =  [(offsets[str(o)] if not is_buried(str(o)) else "x")
					if str(o) in solventAccessibilities else "-"
					for o in offsets_range]
				
				dssp.chain(pdb_chain).add_summary(
					DsspPdbSummaryChain(
						sequence, ss, sa, sa_normalised,
						sum(sa.values())/len(sa), buried_residues,
						accessible_residues,offsets))
		return dssp

	def dssp_parse(self, dssp_input, chains, pdb_id, from_string=False, use_all_chains=False):
		"""
		Reads a dssp file and return a dictionary with the data.
			:param self:
			:param dssp_input: The path for the dssp input file, or dssp string if from_string is 
				True.
			:param chains: A collection with the chains which data is required
				to retrieve. Chain names should have a format like "PDB1.A"
				(pdb id followed by a dot and one letter chain id)
			:param pdb_id: The name of the pdb_id that corresponds to the dssp
				file
			: 
		"""
		if from_string:
			dssp = self.read_string_data(dssp_input, chains, pdb_id, use_all_chains)
		else:
			dssp = self.read_plain_data(dssp_input, chains, pdb_id, use_all_chains)
		dssp = self._append_summary(dssp)

		return dssp
