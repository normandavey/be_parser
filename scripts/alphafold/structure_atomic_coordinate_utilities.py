import os, inspect, sys, pprint, math, json
from re import L, S
from typing import ValuesView

import Bio

file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../"))
import config_reader

sys.path.append(os.path.join(file_path,"../structure"))
from structure.dssp.dssp import Dssp

sys.path.append(os.path.join(file_path,"../"))
import config_reader

sys.path.append(os.path.join(file_path,"../data_management/"))
import queryRunner

sys.path.append(os.path.join(file_path,"../utilities/"))
import utilities_error

from structure_utilities import structureUtilities 

#-----

import logging
logger = logging.getLogger(__name__)

#-----

# aa - longest distance to the CA from pdb_structures rounded up
# used to define atom search space from residue search space
aa_length = {
	'A': 3, #2.4106069542187374,
	'C': 3.5, #2.811611999622088,
	'D': 4.5, #3.7048607231067634,
	'E': 5.5, #5.01415830778013,
	'F': 5.5, #5.151735722525218,
	'G': 3, #2.396801890982823,
	'H': 5, #4.537050730424373,
	'I': 4.5, #3.9231673159405442,
	'K': 7, #6.340017113045783,
	'L': 4.5, #3.993406255313095,
	'M': 6, #5.238990764086978,
	'N': 4.5, #3.732835597758808,
	'P': 3, #2.4365258715661375,
	'Q': 5.5, #4.987556004508565,
	'R': 8, #7.299905086579639,
	'S': 3, #2.4348306099762866,
	'T': 3, #2.573213847111044,
	'V': 3, #2.5574035671843793,
	'W': 6.5, #6.084074021562211,
	'Y': 6.5, #6.483364873435227
	'X': 10, #6.483364873435227
}

DNA_mapping = {
	" DT":"t",
	" DC":"c",
	" DG":"g",
	" DA":"a",
	"  U":"u",
	"  C":"c",
	"  G":"g",
	"  A":"a"
}

RNA_mapping = {
	" RU":"u",
	" RC":"c",
	" RG":"g",
	" RA":"a",
	"  U":"u",
	"  C":"c",
	"  G":"g",
	"  A":"a"
}

class structureAtomicCoordinateUtilities():
	def __init__(self):
		self.options = {}
		
		self.options['dna_distance_cutoff'] = 20
		self.options['primary_atom_type_cutoff'] = 10
		self.options['heavy_atom_distance_cutoff'] = 6
		self.options['primary_atom_type'] = "CA" # WAS N
		self.options['primary_dna_atom_type'] = {  
			# Details of the label of nucleotides
			# https://bmrb.io/referenc/nomenclature/
			# N labels of nitrogens that are adjacent to the sugar
			"a":"N9",
			"g":"N9",
			"c":"N1",
			"t":"N1",
			"u":"N1"
			}
		self.options['output_format'] = "default"
		self.options['all_atom_search'] = True
		self.options['adjacent_residue_filter_distance'] = 3

		self.data = {}
		self.data["contacts"] = {}

		config_options = config_reader.load_configeration_options(sections=["general","accessibility"])

		for config_option in config_options:
			self.options[config_option] = config_options[config_option]

		self.options['intermolecular_path'] = os.path.join(self.options['data_path'],"pdb","intermolecular")

		if not os.path.exists(self.options['intermolecular_path']):
			os.mkdir(self.options['intermolecular_path'])

	def get_pdb_id(self):
		if isinstance(self.options["pdb_id"],(list)):
			return self.options["pdb_id"][0]
		else:
			return self.options["pdb_id"]	

	def is_DNA(self,residue_code):
		if residue_code in DNA_mapping:
			return True
		else:
			return False

	def is_RNA(self,residue_code):
		if residue_code in RNA_mapping:
			return True
		else:
			return False

	def translateResidueCode(self,residue_code):
		
		if residue_code in DNA_mapping:
			residue = DNA_mapping[residue_code]
		elif residue_code in RNA_mapping:
			residue = RNA_mapping[residue_code]
		else:
			try:
				residue = Bio.PDB.Polypeptide.three_to_one(residue_code)
			except:
				residue = "X"

		return residue
		

	def load_structure(self):
		if self.options['pdb_id'] == "alphafold":
			self.structure = Bio.PDB.PDBParser().get_structure(self.options["pdb_id"], self.options['pdb_file'] )
		else:
			pdb_id = self.get_pdb_id()
			
			logger.debug("Load structure " + pdb_id)
			
			pdb_structure_data = queryRunner.queryRunner("pdb","download_pdb_structure",{"pdb_id":pdb_id}).run()
			self.options['pdb_file'] = pdb_structure_data['data']
			
			self.structure = Bio.PDB.PDBParser().get_structure(self.options["pdb_id"], self.options['pdb_file'] )

	###------------###

	def calculate_aa_lengths(self,protein_chain_a):
		logger.debug("Calculate AA length " + self.get_pdb_id() + " " + protein_chain_a.id )
		
		aa_distances = {}
		for query_res_index in range(0,len(protein_chain_a.child_list)):
			
			query_res = protein_chain_a.child_list[query_res_index]
			query_res_name = self.translateResidueCode(query_res.resname)
			
			if self.options['primary_atom_type'] not in query_res: continue

			for query_res_atom_a in [query_res[self.options['primary_atom_type']]]:
				try:
					for query_res_atom_b in query_res:
						diff_vector_atom  = query_res_atom_a.coord - query_res_atom_b.coord 
						atom_distance = math.sqrt(sum(diff_vector_atom*diff_vector_atom))
						
						if query_res_name not in aa_distances:
							aa_distances[query_res_name] = 0
					
						if atom_distance > aa_distances[query_res_name]:
							aa_distances[query_res_name] = atom_distance	
						
				except:
					utilities_error.printError()

	####----------––----------####
	
	def calculate_chain_molecule_type(self,protein_chain):
	
		try:
			if self.options['pdb_id'] == "alphafold":
				return "Protein"
				
			pdb_structure_data = queryRunner.queryRunner("pdb","get_chains_polymer_type",{"pdb_id":self.get_pdb_id(),"query_chains":[protein_chain]}).run()
			return pdb_structure_data['data'][protein_chain]
		except:
			return "Other"

	###------------######------------######------------######------------######------------######------------###
	###------------######------------######------------######------------######------------######------------###
	###------------######------------######------------######------------######------------######------------###
	###------------######------------######------------######------------######------------######------------###

	def calculate_intramolecular_contacts(self):
		
		pdb_id = self.get_pdb_id()
		self.load_structure()

		logger.debug("Calculating intermolecular interaction statistics for " + pdb_id)

		if pdb_id not in self.data["contacts"]:
			self.data["contacts"][pdb_id] = {}

		structureUtilitiesObj = structureUtilities()
		structureUtilitiesObj.options['pdb_id'] = pdb_id

		use_chains = structureUtilitiesObj.get_use_chains()
		
		
		for protein_chain_a in self.structure[0].child_list:
			if protein_chain_a.id not in use_chains: 
				continue
			
			if protein_chain_a.id not in self.data["contacts"][pdb_id]:
				self.data["contacts"][pdb_id][protein_chain_a.id] = {}
			
			self.calculate_chain_residue_contacts(protein_chain_a,protein_chain_a)
			
		return self.data["contacts"][pdb_id]
		
	###------------###

	def calculate_intermolecular_contacts(self):
		
		pdb_id = self.get_pdb_id()

		tag = "_".join([str(self.options['adjacent_residue_filter_distance']), str(self.options['primary_atom_type']), str(self.options['heavy_atom_distance_cutoff']), str(self.options['primary_atom_type_cutoff'])])
		
		if self.options['all_atom_search']:
			tag += "_all_atom_search"

		out_file = os.path.join(self.options['intermolecular_path'],pdb_id + "." + tag + ".intermolecular.v2.json")

		if not os.path.exists(out_file) or self.options['remake']:
			self.load_structure()

			logger.debug("Calculating intermolecular interaction statistics for " + pdb_id)

			if pdb_id not in self.data["contacts"]:
				self.data["contacts"][pdb_id] = {}

			structureUtilitiesObj = structureUtilities()
			structureUtilitiesObj.options['pdb_id'] = pdb_id
			structureUtilitiesObj.options['use_only_biological_assembly'] = self.options['use_only_biological_assembly']

			use_chains = structureUtilitiesObj.get_use_chains()\

			self.data['biological_assembly'] = structureUtilitiesObj.get_biological_assembly()

			for protein_chain_a in self.structure[0].child_list:
				if protein_chain_a.id not in use_chains: 
					continue
				
				if protein_chain_a.id not in self.data["contacts"][pdb_id]:
					self.data["contacts"][pdb_id][protein_chain_a.id] = {}
				
				for protein_chain_b in self.structure[0].child_list:
					if protein_chain_b.id not in use_chains: continue
					if protein_chain_a != protein_chain_b:
						chain_id_b =  pdb_id + "." + protein_chain_b.id

						if chain_id_b in self.data["contacts"][pdb_id][protein_chain_a.id]: continue

						self.calculate_chain_residue_contacts(protein_chain_a,protein_chain_b)
			
			response = self.data["contacts"][pdb_id]
			logger.debug("Writing " + out_file)

			with open(out_file, 'w') as out_file_handle:
				json.dump(response, out_file_handle)

		else:
			logger.debug("Exists: intermolecular interaction statistics for " + pdb_id + "-> "+ out_file)

		with open(out_file) as outfile:
			response = json.load(outfile)

		return response

	###------------###

	def calculate_interface_chain_contacts_dna(self,query_res,contact_res,primary_distance):

		min_distance = [100,None,None] 

		query_res_offset = query_res.id[1]
		contact_res_offset = contact_res.id[1]
		query_res_name = self.translateResidueCode(query_res.resname)
		contact_res_name = self.translateResidueCode(contact_res.resname)
				
		query_interactor_residue = {}
		contact_interactor_residue = {}

		for query_res_atom in query_res:
			if query_res_atom.element == "H": continue
			for contact_res_atom in contact_res:
				if contact_res_atom.element == "H": continue
				
				diff_vector_atom  = query_res_atom.coord - contact_res_atom.coord 
				atom_distance = math.sqrt(sum(diff_vector_atom*diff_vector_atom))

				if atom_distance < min_distance[0]:
					min_distance = [atom_distance,query_res_atom,contact_res_atom]
		
		if min_distance[0] < self.options['dna_distance_cutoff']:
			
			query_interactor_residue = {
				"query_res_pos":query_res_offset,
				"contact_res_pos":contact_res_offset,
				"query_res_name":query_res_name,
				"contact_res_name":contact_res_name,
				"query_res_atom":min_distance[1].name,
				"contact_res_atom":min_distance[2].name,
				"atom_distance":min_distance[0],
				"primary_atom_distance":primary_distance,
				"primary_atom":self.options['primary_atom_type']
			}

			contact_interactor_residue = {
				"query_res_pos":contact_res_offset,
				"query_res_name":contact_res_name,
				"query_res_atom":min_distance[2].name,
				"contact_res_pos":query_res_offset,
				"contact_res_name":query_res_name,
				"contact_res_atom":min_distance[1].name,
				"atom_distance":min_distance[0],
				"primary_atom_distance":primary_distance,
				"primary_atom":self.options['primary_atom_type']
			}

		return {"query_interactor_residue":query_interactor_residue,"contact_interactor_residue":contact_interactor_residue}

	###------------###

	def calculate_interface_chain_contacts_all_atom(self,query_res,contact_res,primary_distance):

		min_distance = [100,None,None] 

		query_res_offset = query_res.id[1]
		contact_res_offset = contact_res.id[1]
		query_res_name = self.translateResidueCode(query_res.resname)
		contact_res_name = self.translateResidueCode(contact_res.resname)
				
		query_interactor_residue = {}
		contact_interactor_residue = {}

		for query_res_atom in query_res:
			if query_res_atom.element == "H": continue
			for contact_res_atom in contact_res:
				if contact_res_atom.element == "H": continue

				diff_vector_atom  = query_res_atom.coord - contact_res_atom.coord 
				atom_distance = math.sqrt(sum(diff_vector_atom*diff_vector_atom))

				if atom_distance < min_distance[0]:
					min_distance = [atom_distance,query_res_atom,contact_res_atom]
		
		if min_distance[0] < self.options['heavy_atom_distance_cutoff']:
			
			query_interactor_residue = {
				"query_res_pos":query_res_offset,
				"contact_res_pos":contact_res_offset,
				"query_res_name":query_res_name,
				"contact_res_name":contact_res_name,
				"query_res_atom":min_distance[1].name,
				"contact_res_atom":min_distance[2].name,
				"atom_distance":min_distance[0],
				"primary_atom_type":primary_distance,
				"primary_atom":self.options['primary_atom_type']
			}

			contact_interactor_residue = {
				"query_res_pos":contact_res_offset,
				"query_res_name":contact_res_name,
				"query_res_atom":min_distance[2].name,
				"contact_res_pos":query_res_offset,
				"contact_res_name":query_res_name,
				"contact_res_atom":min_distance[1].name,
				"atom_distance":min_distance[0],
				"primary_atom_type":primary_distance,
				"primary_atom":self.options['primary_atom_type']
			}

		return {"query_interactor_residue":query_interactor_residue,"contact_interactor_residue":contact_interactor_residue}

	###------------###

	def calculate_interface_chain_contacts_primary_atom(self,query_res,contact_res,primary_distance):

		query_res_offset = query_res.id[1]
		contact_res_offset = contact_res.id[1]
		query_res_name = self.translateResidueCode(query_res.resname)
		contact_res_name = self.translateResidueCode(contact_res.resname)
		
		query_interactor_residue = {
			"query_res_pos":query_res_offset,
			"query_res_name":query_res_name,
			"query_res_atom":None,
			"contact_res_pos":contact_res_offset,
			"contact_res_name":contact_res_name,
			"contact_res_atom":None,
			"atom_distance":None,
			"primary_atom_type":primary_distance,
			"primary_atom":self.options['primary_atom_type']
		}

		contact_interactor_residue = {
			"query_res_pos":contact_res_offset,
			"query_res_name":contact_res_name,
			"query_res_atom":None,
			"contact_res_pos":query_res_offset,
			"contact_res_name":query_res_name,
			"contact_res_atom":None,
			"atom_distance":None,
			"primary_atom_type":primary_distance,
			"primary_atom":self.options['primary_atom_type']
		}

		return {"query_interactor_residue":query_interactor_residue,"contact_interactor_residue":contact_interactor_residue}

	####----------––----------####
	
	def calculate_chain_residue_contacts(self,protein_chain_a,protein_chain_b):
		pdb_id = self.get_pdb_id() 

		chain_id_a = protein_chain_a.id
		chain_id_b = protein_chain_b.id

		logger.debug("Calculate interface contacts for " + pdb_id + " " + chain_id_a  + " -> " + chain_id_b)
		
		if pdb_id not in self.data["contacts"]:
			self.data["contacts"][pdb_id] = {}

		query_interactor_residues = {}
		contact_interactor_residues = {}

		for query_res_index in range(0,len(protein_chain_a.child_list)):
			query_res = protein_chain_a.child_list[query_res_index]
			query_res_name = self.translateResidueCode(query_res.resname)

			#if query_res.get_full_id()[3][0] != " ": continue
			if not (Bio.PDB.Polypeptide.is_aa(query_res) or query_res in ["a","g","c","t","u"]): continue
			if Bio.PDB.Polypeptide.is_aa(query_res):
				if self.options['primary_atom_type'] not in query_res: 
					logger.error(self.options['primary_atom_type'] + " not in " + str(query_res))
					continue
			
			query_res_offset = query_res.id[1]
			query_res_offset_flanks = [query_res_offset - i for i in range(1,self.options['adjacent_residue_filter_distance'])] + [query_res_offset + i for i in range(1,self.options['adjacent_residue_filter_distance'])]
			
			for contact_res_index in range(0,len(protein_chain_b.child_list)):
			
				try:
					contact_res = protein_chain_b.child_list[contact_res_index]
					contact_res_name = self.translateResidueCode(contact_res.resname)
				
					if not (Bio.PDB.Polypeptide.is_aa(contact_res) or contact_res_name in ["a","g","c","t","u"]): continue
					if Bio.PDB.Polypeptide.is_aa(contact_res):
						if self.options['primary_atom_type'] not in contact_res: 
							logger.error(self.options['primary_atom_type'] + " not in " + str(contact_res))
							continue

					contact_res_offset = contact_res.id[1]
					contact_res_name = self.translateResidueCode(contact_res.resname)
				
					if query_res_offset not in query_interactor_residues: query_interactor_residues[query_res_offset] = {}
					if contact_res_offset not in contact_interactor_residues: contact_interactor_residues[contact_res_offset] = {}
					
					if not Bio.PDB.Polypeptide.is_aa(contact_res) and not Bio.PDB.Polypeptide.is_aa(query_res) : # skip non amino acids
						pass
					elif contact_res == query_res: # skip self contacts
						pass
					elif contact_res_offset in query_res_offset_flanks: # Skip the adjacent amino acids on the linear peptide chain
						pass
					elif self.options['primary_atom_type'] in contact_res and self.options['primary_atom_type'] in query_res: 
						# Check if the right atom is in the amino acid

						diff_vector  = contact_res[self.options['primary_atom_type']].coord - query_res[self.options['primary_atom_type']].coord 
						primary_distance = math.sqrt(sum(diff_vector*diff_vector))

						# Maximum distance between too C-alpha that fit the heavy_atom_distance_cutoff criteria

						if self.options['all_atom_search']:
							primary_distance_check = primary_distance < self.options['heavy_atom_distance_cutoff'] + aa_length[query_res_name] + aa_length[contact_res_name]
						else:
							if self.options['primary_atom_type_cutoff'] == "aa_dependent":
								primary_distance_check = primary_distance < self.options['heavy_atom_distance_cutoff'] + aa_length[query_res_name] + aa_length[contact_res_name]
							else:
								primary_distance_check = primary_distance < self.options['primary_atom_type_cutoff']

						if primary_distance_check: 
							if self.options['all_atom_search']:
								contacts = self.calculate_interface_chain_contacts_all_atom(query_res,contact_res,primary_distance)
								
								if len(contacts["query_interactor_residue"]) > 0:
									if query_res_offset not in query_interactor_residues: query_interactor_residues[query_res_offset] = {}
									query_interactor_residues[query_res_offset][contact_res_offset] = contacts["query_interactor_residue"]

								if len(contacts["contact_interactor_residue"]) > 0:
									if contact_res_offset not in contact_interactor_residues: contact_interactor_residues[contact_res_offset] = {}
									contact_interactor_residues[contact_res_offset][query_res_offset] = contacts["contact_interactor_residue"]
							else:
								contacts = self.calculate_interface_chain_contacts_primary_atom(query_res,contact_res,primary_distance)
								if len(contacts["query_interactor_residue"]) > 0:
									query_interactor_residues[query_res_offset][contact_res_offset] = contacts["query_interactor_residue"]

								if len(contacts["contact_interactor_residue"]) > 0:
									contact_interactor_residues[contact_res_offset][query_res_offset] = contacts["contact_interactor_residue"]
					elif not Bio.PDB.Polypeptide.is_aa(contact_res) and not Bio.PDB.Polypeptide.is_aa(query_res):
						pass
					elif (query_res_name in ["a","g","c","t","u"] and Bio.PDB.Polypeptide.is_aa(contact_res)) or (contact_res_name in ["a","g","c","t","u"] and Bio.PDB.Polypeptide.is_aa(query_res)):
						# Details of the label of nucleotides
						# https://bmrb.io/referenc/nomenclature/
						try:
							primary_distance = None
							if query_res_name in ["a","g","c","t","u"] and Bio.PDB.Polypeptide.is_aa(contact_res):
								diff_vector = contact_res[self.options['primary_atom_type']].coord - query_res[self.options['primary_dna_atom_type'][query_res_name]].coord 
								primary_distance = math.sqrt(sum(diff_vector*diff_vector))

							if contact_res_name in ["a","g","c","t","u"] and Bio.PDB.Polypeptide.is_aa(query_res):
								if self.options['primary_atom_type'] not in query_res: 
									logger.error(self.options['primary_atom_type'] + " not in " + str(query_res))
									continue
								try:
									diff_vector = contact_res[self.options['primary_dna_atom_type'][contact_res_name]].coord - query_res[self.options['primary_atom_type']].coord 
									primary_distance = math.sqrt(sum(diff_vector*diff_vector))
								except:
									logger.error("Distance calculation error: " + str(contact_res)  + " " + str(query_res))

							if primary_distance != None:
								if primary_distance < self.options['dna_distance_cutoff']:
									contacts = self.calculate_interface_chain_contacts_dna(query_res,contact_res,None)
									if len(contacts["query_interactor_residue"]) > 0:
										if query_res_offset not in query_interactor_residues: query_interactor_residues[query_res_offset] = {}
										query_interactor_residues[query_res_offset][contact_res_offset] = contacts["query_interactor_residue"]

									if len(contacts["contact_interactor_residue"]) > 0:
										if contact_res_offset not in contact_interactor_residues: contact_interactor_residues[contact_res_offset] = {}
										contact_interactor_residues[contact_res_offset][query_res_offset] = contacts["contact_interactor_residue"]
						except:
							for atom in query_res:
								logger.error(str(query_res) + " " + str(atom))
							
							for atom in contact_res:
								logger.error(str(contact_res)  + " " + str(atom))

							utilities_error.printError(description=self.options)
					else:
						pass
				except:
					logger.debug("Error @" + str(contact_res_index))
					utilities_error.printError(description=self.options)
					
		
		###################################################

		if chain_id_a not in self.data["contacts"][pdb_id]:
			self.data["contacts"][pdb_id][chain_id_a] = {}

		if chain_id_b not in self.data["contacts"][pdb_id]:
			self.data["contacts"][pdb_id][chain_id_b] = {}
		
		###################################################

		self.data['contacts'][pdb_id][chain_id_a][chain_id_b] = query_interactor_residues
		self.data['contacts'][pdb_id][chain_id_b][chain_id_a] = contact_interactor_residues

		###################################################
		return {"query_interactor_residues":query_interactor_residues,"contact_interactor_residues":contact_interactor_residues}
	
	###------------######------------######------------######------------######------------######------------###
	###------------######------------######------------######------------######------------######------------###
	###------------######------------######------------######------------######------------######------------###
	###------------######------------######------------######------------######------------######------------###


	def use_chain_check(self,chain_id,chain_type="query_chain"):
		use_chain = True

		if chain_type in self.options:
			use_chain = False
			if len(self.options[chain_type]) == 0:
				use_chain = True
			else:
				if isinstance(self.options[chain_type],(str)):
					if chain_id == self.options[chain_type]:
						use_chain = True
				else:
					if chain_id in self.options[chain_type]:
						use_chain = True

		return use_chain

	###----------###

	def calculate_intermolecular_contact_statistics(self):
		pdb_id = self.get_pdb_id()
		self.load_structure()

		self.data["statistics"] = {}
		self.options['all_atom_search'] = False

		logger.debug("Calculating intramolecular interaction statistics for " + pdb_id)
		
		self.data["intramolecular"] = {}
		for protein_chain_a in self.structure[0].child_list:
			if not self.use_chain_check(protein_chain_a.id,chain_type="query_chain"): continue
			self.data["statistics"][protein_chain_a.id ] = {'by_chain':{},'summary':{}}

			for protein_chain_b in self.structure[0].child_list:
				if not self.use_chain_check(protein_chain_b.id,chain_type="use_chains"): continue
				if protein_chain_a.id != protein_chain_b.id:
					mean_local_aa_population = self.calculate_chain_molecular_contact_statistics(protein_chain_a,protein_chain_b)
					if 'data' in mean_local_aa_population:
						self.data["statistics"][protein_chain_a.id ]['by_chain'][protein_chain_b.id] = mean_local_aa_population['data']
					else:
						self.data["statistics"][protein_chain_a.id]['by_chain'][protein_chain_b.id] = mean_local_aa_population['error_type']

			mean_chain_residue_local_aa_counts = []
			for chain_id in self.data["statistics"][protein_chain_a.id]['by_chain']:
				try:
					mean_chain_residue_local_aa_counts.append(self.data["statistics"][protein_chain_a.id]['by_chain'][chain_id]['mean_chain_residue_local_aa_count'])
				except:
					logger.error("Skipping " + pdb_id + " " + protein_chain_a.id + " -> " + chain_id + " mean_chain_residue_local_aa_counts")
				 
			if len(mean_chain_residue_local_aa_counts) > 0:
				self.data["statistics"][protein_chain_a.id]['summary'] = max(mean_chain_residue_local_aa_counts)
				
		return self.data["statistics"]

	###----------###

	def calculate_intramolecular_contact_statistics(self):

		pdb_id = self.get_pdb_id()
		self.load_structure()

		self.data["statistics"] = {}
		self.options['output_format'] = "mean_chain_residue_local_aa_count"
		self.options['all_atom_search'] = False

		logger.debug("Calculating intramolecular interaction statistics for " + pdb_id)
		
		self.data["intramolecular"] = {}
		for protein_chain in self.structure[0].child_list:
			if not self.use_chain_check(protein_chain.id,chain_type="query_chain"): continue
			self.data["statistics"][protein_chain.id ] = {'summary':{}}
			mean_local_aa_population = self.calculate_chain_molecular_contact_statistics(protein_chain)
		
			if 'data' in mean_local_aa_population:
				self.data["statistics"][protein_chain.id ]['summary'] = mean_local_aa_population['data']
			else:
				self.data["statistics"][protein_chain.id ]['summary'] = mean_local_aa_population['error_type']
			
		return self.data["statistics"]

	###------------###

	def calculate_chain_molecular_contact_statistics(self,protein_chain_a,protein_chain_b=None):

		logger.debug("Calculate calculate chain intramolecular contact_statistics for " + protein_chain_a.id)

		###
		
		if protein_chain_b == None:
			protein_chain_b = protein_chain_a

		molecule_type_a = self.calculate_chain_molecule_type(protein_chain_a.id)
		molecule_type_b = self.calculate_chain_molecule_type(protein_chain_b.id)
		interface_contacts = self.calculate_chain_residue_contacts(protein_chain_a,protein_chain_b)
		
		###	

		if 'status' in interface_contacts["query_interactor_residues"]:
			return {"status":"Error","error_type":interface_contacts["query_interactor_residues"]['status']}
		else:
			interactor_residues_summer = {}
			for query_res_id in interface_contacts["query_interactor_residues"]:
				interactor_residues_summer[query_res_id] = len(set(interface_contacts["query_interactor_residues"][query_res_id]))
			
			if len(interactor_residues_summer) != 0:

				if self.options['output_format'] == "detailed_distances":
					mean_local_aa_population = {
						"chain_residue_local_aa_count_total":sum(interactor_residues_summer.values()),
						"chain_residue_count":len(interactor_residues_summer),
						"chain_residue_no_local_aa_count":len(interactor_residues_summer) - list(interactor_residues_summer.values()).count(0),
						"mean_chain_residue_local_aa_count":float(sum(interactor_residues_summer.values()))/len(interactor_residues_summer),
						"interactor_residues_query_res":interface_contacts["query_interactor_residues"],
						"interactor_residues_contact_res":interface_contacts["contact_interactor_residues"]
						}
				elif self.options['output_format'] == "mean_chain_residue_local_aa_count":
					if molecule_type_a == "Protein" or molecule_type_b == "Protein":
						mean_local_aa_population = float(sum(interactor_residues_summer.values()))/len(interactor_residues_summer)
					else:
						mean_local_aa_population = molecule_type_a 
				else:
					mean_local_aa_population = {
						"chain_residue_local_aa_count_total":sum(interactor_residues_summer.values()),
						"chain_residue_count":len(interactor_residues_summer),
						"chain_residue_no_local_aa_count":len(interactor_residues_summer) - list(interactor_residues_summer.values()).count(0),
						"mean_chain_residue_local_aa_count":float(sum(interactor_residues_summer.values()))/len(interactor_residues_summer),
					}

				return {"status":"Success","data":mean_local_aa_population}
			else:
				mean_local_aa_population = False
				return {"status":"Error","error_type":"No data produced. Molecule types = " + molecule_type_a + ":" + protein_chain_a.id + " -> " + molecule_type_b + ":" + protein_chain_b.id }
		

if __name__ == "__main__":
	structureAtomicCoordinateUtilitiesObj = structureAtomicCoordinateUtilities()

	for pdb_id in ["3LV3","3LKN","6MBB","3H8A","1GG6","5J3Q","1HTL","4O9W","6JNR"]:
		structureAtomicCoordinateUtilitiesObj.options['pdb_id'] = pdb_id

		calculate_intramolecular_interaction_statistics_data = structureAtomicCoordinateUtilitiesObj.calculate_intramolecular_interaction_statistics()
		pprint.pprint(calculate_intramolecular_interaction_statistics_data)

		calculate_intermolecular_interaction_statistics_data = structureAtomicCoordinateUtilitiesObj.calculate_intermolecular_interaction_statistics()
		pprint.pprint(calculate_intermolecular_interaction_statistics_data)

