import os, inspect, sys, pprint, re

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBIO import Select
from os.path import splitext

file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../"))
import config_reader

sys.path.append(os.path.join(file_path,"../"))
import config_reader

sys.path.append(os.path.join(file_path,"../data_management/"))
import queryRunner

sys.path.append(os.path.join(file_path,"../utilities/"))
import utilities_error

#-----
import logging
logger = logging.getLogger(__name__)
#-----

class structureUtilities():
	def __init__(self):
		self.options = {}
		self.options['use_only_biological_assembly'] = False
		
		self.data = {}

		config_options = config_reader.load_configeration_options(sections=["general","accessibility"])

		for config_option in config_options:
			self.options[config_option] = config_options[config_option]

	def get_pdb_id(self):
		if isinstance(self.options["pdb_id"],(list)):
			return self.options["pdb_id"][0]
		else:
			return self.options["pdb_id"]	

	
	def get_biological_assembly_chain_comembers(self):
		chain_biological_assembly_members = {}
		
		self.data['biological_assembly'] = self.get_biological_assembly()
		for biological_assembly in self.data['biological_assembly']:
			for chain_a in self.data['biological_assembly'][biological_assembly]["chains"]:
				if chain_a not in chain_biological_assembly_members:
					chain_biological_assembly_members[chain_a] = []
				else:
					continue

				for chain_b in self.data['biological_assembly'][biological_assembly]["chains"]:
					if chain_a != chain_b:
						chain_biological_assembly_members[chain_a].append(chain_b)

		return chain_biological_assembly_members

	def get_biological_assembly(self):
		"""
		Gets the primary biological assemblies for the structure.
			:param self:
		"""

		pdb_id = self.get_pdb_id()
		pdb_structure_data = queryRunner.queryRunner("pdb","download_pdb_structure",{"pdb_id":pdb_id}).run()
		self.options['pdb_file'] = pdb_structure_data['data']

		biological_assemblies = {}

		remark_pattern = re.compile("REMARK 350.+")
		biomolecule_pattern = re.compile("BIOMOLECULE: .+")
		author_determined_pattern = re.compile("AUTHOR DETERMINED BIOLOGICAL UNIT: .+")
		software_determined_pattern = re.compile("SOFTWARE DETERMINED QUATERNARY STRUCTURE: .+")
		software_pattern = re.compile("SOFTWARE USED: .+")
		chains_pattern = re.compile("APPLY THE FOLLOWING TO CHAINS: .+")
		chains_extended_pattern = re.compile("AND CHAINS: .+")
		
		if os.path.exists(self.options['pdb_file'] ):

			pdb_str = open(self.options['pdb_file'] ).read()
		
			remarks_350 = "\n".join(remark_pattern.findall(pdb_str))
			
			for remark_350 in remarks_350.split("REMARK 350 BIOMOLECULE: ")[1:]:
			
				biomolecule = remark_350.split("\n")[0].strip()

				biological_assemblies[biomolecule] = {
					"author_determined":"",
					"software_determined":"",
					"software":"",
					"chains":[]
				}

				author_determined_matches = author_determined_pattern.findall(remark_350)
				software_determined_matches = software_determined_pattern.findall(remark_350)
				software_matches = software_pattern.findall(remark_350)
				chains_matches = chains_pattern.findall(remark_350)
				chains_extended_matches = chains_extended_pattern.findall(remark_350)
				
				if len(author_determined_matches) > 0:
					biological_assemblies[biomolecule]['author_determined'] = author_determined_matches[0].split(":")[1].strip()
				if len(software_determined_matches) > 0:
					biological_assemblies[biomolecule]['software_determined'] = software_determined_matches[0].split(":")[1].strip()
				if len(software_matches) > 0:
					biological_assemblies[biomolecule]['software'] = software_matches[0].split(":")[1].strip()
				if len(chains_matches) > 0:
					for chains_match in chains_matches:
						biological_assemblies[biomolecule]['chains'] += chains_match.split(":")[1].strip().strip(",").replace(" ","").split(",")
				if len(chains_extended_matches) > 0:
					biological_assemblies[biomolecule]['chains'] += chains_extended_matches[0].split(":")[1].strip().strip(",").replace(" ","").split(",")

		if len(biological_assemblies) == 0:
			structure_chain_data = queryRunner.queryRunner("pdb","get_structure_chain_data",{"pdb_id":pdb_id}).run()
			
			return {"1":{
					"author_determined":"Not set. Using all chains",
					"software_determined":"",
					"software":"",
					"chains":[x.split(".")[1] for x in list(structure_chain_data['data'].keys())]
					}
				}
		else:
			return biological_assemblies

	def get_use_chains(self):

		if self.options['pdb_id'] == "alphafold":
			return ["A"]
			
		use_chains = []
		if self.options['use_only_biological_assembly']:
			self.data['biological_assembly'] = self.get_biological_assembly()
			use_chains = self.data['biological_assembly']['1']['chains']
		else:
			self.data['chains_structure'] = queryRunner.queryRunner("pdb","get_chains_structure",{"pdb_id":self.get_pdb_id()}).run()
			
			if 'data' in self.data['chains_structure']:
				use_chains = self.data['chains_structure']['data']
		
		return use_chains

if __name__ == "__main__":
	structureUtilitiesObj = structureUtilities()

	for pdb_id in ["2AST"]:
		structureUtilitiesObj.options['pdb_id'] = pdb_id
		biological_assembly_data = structureUtilitiesObj.get_biological_assembly()
		pprint.pprint(biological_assembly_data)

		#offsets = dssp_data[pdb_id + ".A"]['summary']['surface_accessibility_normalised']
		#accessible_score_sequence = ""
		#for offset in range(min(offsets),max(offsets)):
		#	accessible_score_sequence += str(int(dssp_data[pdb_id + ".A"]['summary']['surface_accessibility_normalised'][offset]*10))

		#print(dssp_data[pdb_id + ".A"]['summary']['sequence'])
		#print(accessible_score_sequence)
		#print("".join(dssp_data[pdb_id + ".A"]['summary']['accessible_residues']))
