import os, inspect, sys, pprint

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
import utilities_alignment

#-----
import logging
logger = logging.getLogger(__name__)

##-----##

class dsspUtilities():
	def __init__(self):
		self.options = {}
		self.data = {}

		config_options = config_reader.load_configeration_options(sections=["general","accessibility"])

		for config_option in config_options:
			self.options[config_option] = config_options[config_option]

		if not os.path.exists(self.options['dssp_path']):
			logger.error(self.options['dssp_path'] + " does not exist")

	def get_pdb_id(self):
		if isinstance(self.options["pdb_id"],(list)):
			return self.options["pdb_id"][0]
		else:
			return self.options["pdb_id"]	

	def get_dssp_data(self,pdb_id=None,pdb_file=None,alphafold=False):
		logger.debug(pdb_id + " " + pdb_file + " ->" + self.options['dssp_path'])

		if alphafold:
			self.options['chain_ids'] = "A"
		else:
			if pdb_id == None:	
				pdb_id = self.get_pdb_id()
		
			pdb_structure_data = queryRunner.queryRunner("pdb","download_pdb_structure",{"pdb_id":pdb_id}).run()
				
			if pdb_file == None:
				pdb_file = os.path.join(self.options["data_path"],"pdb",pdb_id + ".pdb")

		if 'chain_ids' not in self.options:
			self.options['chain_ids'] = "A"

		dssp_runner = Dssp(self.options['dssp_path'])

		try:
			dssp_monomers, _ = dssp_runner.run(pdb_file, pdb_id, splitted=True,chain_ids=self.options['chain_ids'] )
			dssp_joined, _ = dssp_runner.run(pdb_file, pdb_id, splitted=False,chain_ids=self.options['chain_ids'])
		except:
			utilities_error.printError()
			sys.stderr.write("ERROR: There was a problem running dssp for pdb: {}".format(pdb_file))
			self.data["dssp"] = {}
			raise

		dssp_diff = dssp_joined.solvent_accesibility_difference(dssp_monomers)
		dssp_data = dssp_diff.to_dict()

		for chain_id in dssp_monomers.to_dict().keys():
			surface_accessibility_normalised_split = list(dssp_monomers.to_dict()[chain_id]['summary']['surface_accessibility_normalised'].values())
			mean_surface_accessibility_normalised_split = sum(surface_accessibility_normalised_split)/len(surface_accessibility_normalised_split)
			dssp_data[chain_id]['summary']['mean_surface_accessibility_normalised_split']  = mean_surface_accessibility_normalised_split
		
		return dssp_data

	def get_dssp_complex_data(self):
		return self.get_dssp_data()

	def get_dssp_monomeric_data(self):
		pdb_id = self.get_pdb_id()
	
		pdb_structure_data = queryRunner.queryRunner("pdb","download_pdb_structure",{"pdb_id":pdb_id}).run()
		pdb_file = pdb_structure_data['data']
		
		dssp_runner = Dssp(self.options['dssp_path'])

		try:
			dssp_monomers, _ = dssp_runner.run(pdb_file, pdb_id, splitted=True,chain_ids=self.options['chain_ids'])
		except:
			sys.stderr.write("ERROR: There was a problem running dssp for pdb: {} {}".format(pdb_id,pdb_file))
			self.data["dssp"] = {}
			raise

		dssp_data = dssp_monomers.to_dict()

		for chain_id in dssp_monomers.to_dict().keys():
			surface_accessibility_normalised_split = list(dssp_monomers.to_dict()[chain_id]['summary']['surface_accessibility_normalised'].values())
			mean_surface_accessibility_normalised_split = sum(surface_accessibility_normalised_split)/len(surface_accessibility_normalised_split)
			dssp_data[chain_id]['summary']['mean_surface_accessibility_normalised_split']  = mean_surface_accessibility_normalised_split
		
		return dssp_data

	def get_dssp_sequence_from_offsets(self,offsets):
		offset_keys = list(offsets.keys())
		#print(offset_keys)
		#offset_keys.sort()

		dssp_sequence = ''
		
		numeric_filter = filter(str.isdigit, offset_keys[0])
		last = int("".join(numeric_filter))
		
		for offset_key in offset_keys:
			numeric_filter = filter(str.isdigit, offset_key)
			current = int("".join(numeric_filter))
			if last <= current: 
				pass
			else:
				logger.debug("pdb numberiung out of order",last,offset_key)
				pass 

			dssp_sequence += offsets[offset_key]

		return dssp_sequence

	def align_dssp_sequence_to_uniprot(self,dssp_sequence,uniprot_sequence):
		alignment = utilities_alignment.alignPeptides(dssp_sequence,uniprot_sequence)
		
		return {
			'dssp_sequence_aligned':alignment[0][0],
			'uniprot_sequence_aligned':alignment[0][1]
		}

	def map_dssp_offsets_to_uniprot(self,pdb_id,chain_id,offsets,return_alignment=False):
		logger.debug("map_dssp_offsets_to_uniprot " + pdb_id + " " + chain_id)

		pdb_id = self.get_pdb_id()
		uniprot_mapping_structure_data = queryRunner.queryRunner("pdb","get_sequence_structure",{"pdb_id":pdb_id,"chain_id":chain_id.split(".")[-1]}).run()
		
		offset_keys = list(offsets.keys())
		offset_keys.sort()

		dssp_mapping = {}
		dssp_sequence_alignment = {}

		if 'data' not in uniprot_mapping_structure_data:
			return dssp_mapping

		for uniprot_accession in uniprot_mapping_structure_data['data']:
			
			uniprot_sequence = uniprot_mapping_structure_data['data'][uniprot_accession]['uniprot_sequence'] 
			
			uniprot_sequence_start = uniprot_mapping_structure_data['data'][uniprot_accession]['unp_start']
			uniprot_sequence_end = uniprot_mapping_structure_data['data'][uniprot_accession]['unp_end']
			offset_sequence_start = uniprot_mapping_structure_data['data'][uniprot_accession]['start']['author_residue_number']
			offset_sequence_end = uniprot_mapping_structure_data['data'][uniprot_accession]['end']['author_residue_number']
			
			dssp_sequence = self.get_dssp_sequence_from_offsets(offsets)
			dssp_sequence_alignment = self.align_dssp_sequence_to_uniprot(dssp_sequence,uniprot_sequence)
			
			dssp_offset = 0
			uniprot_offset = uniprot_sequence_start
			
			#print(dssp_sequence)
			#print(uniprot_sequence)
			#print(dssp_sequence_alignment['dssp_sequence_aligned'])
			#print(dssp_sequence_alignment['uniprot_sequence_aligned'])

			for i in range(0,len(dssp_sequence_alignment['dssp_sequence_aligned'])):
				try:
					dssp_sequence_aligned_aa = dssp_sequence_alignment['dssp_sequence_aligned'][i]
					uniprot_sequence_aligned_aa = dssp_sequence_alignment['uniprot_sequence_aligned'][i]
					
					if dssp_offset >= len(offset_keys):
						continue
					
					dssp_offset_pdb_numbering = offset_keys[dssp_offset]

					if not dssp_offset_pdb_numbering.isnumeric():
						pass
					elif int(dssp_offset_pdb_numbering) <=  offset_sequence_start or int(dssp_offset_pdb_numbering) >= offset_sequence_end: 
						if dssp_sequence_aligned_aa != "-": dssp_offset += 1
						continue

					aligned = dssp_sequence_aligned_aa == uniprot_sequence_aligned_aa
					gapped = dssp_sequence_aligned_aa == "-" or uniprot_sequence_aligned_aa == ""

					if not gapped:#dssp_sequence_aligned_aa != "-":
						if dssp_offset_pdb_numbering not in dssp_mapping:
							dssp_mapping[dssp_offset_pdb_numbering] = {}

						dssp_mapping[dssp_offset_pdb_numbering][uniprot_accession] = {
							"uniprot_offset":uniprot_offset,
							"dssp_sequence_aligned_aa":dssp_sequence_aligned_aa,
							"uniprot_sequence_aligned_aa":uniprot_sequence_aligned_aa,
							"gapped":gapped,
							"aligned":aligned
						}

					
					logger.debug("\t".join([
						pdb_id,
						chain_id,
						uniprot_accession,
						str(uniprot_offset),
						str(dssp_offset_pdb_numbering),
						dssp_sequence_aligned_aa,
						uniprot_sequence_aligned_aa,
						str(gapped),
						str(aligned),
						str(uniprot_sequence_start),
						str(uniprot_sequence_end),
						str(dssp_offset)
					]))

					if dssp_sequence_aligned_aa != "-": dssp_offset += 1
					if uniprot_sequence_aligned_aa != "-": uniprot_offset += 1
				except:
					logger.debug("Skipping offset " + str(i))
					raise
		
		logger.debug("Complete")

		if return_alignment:
			return {
				"residue_mapping":dssp_mapping,
				"alignment":dssp_sequence_alignment,
				"region_mapping":uniprot_mapping_structure_data['data']
			}
		else:
			return dssp_mapping
							
	def get_specific_surface_accessibility_data(self):
		self.get_surface_accessibility_data()

	def get_dssp_chain_data(self,pdb_id,pdb_chain,dssp_diff,dssp_monomers,dssp_joined):
		offsets = dssp_monomers.chain(pdb_chain).get_summary().offsets
		dssp_uniprot_mapping = self.map_dssp_offsets_to_uniprot(pdb_id,pdb_chain,offsets)
		
		sa_normalised_diff = dssp_diff.chain(pdb_chain).get_summary().surface_accessibility_normalised if pdb_chain in dssp_diff else {}
		sa_normalised_splitted = dssp_monomers.chain(pdb_chain).get_summary().surface_accessibility_normalised if pdb_chain in dssp_diff else {}
		sa_normalised_joined = dssp_joined.chain(pdb_chain).get_summary().surface_accessibility_normalised if pdb_chain in dssp_diff else {}

		sa_normalised_interface = {}
		sa_normalised_interface_proportion = {}

		offsets = dssp_monomers.chain(pdb_chain).get_summary().offsets
		
		for x in sa_normalised_diff:
			try:
				if sa_normalised_splitted[x] - sa_normalised_joined[x]  != 0:
					sa_normalised_interface[x] = (sa_normalised_splitted[x] - sa_normalised_joined[x])
					sa_normalised_interface_proportion[x] = (sa_normalised_splitted[x] - sa_normalised_joined[x] )/sa_normalised_splitted[x]
				else:
					sa_normalised_interface[x] = 0
					sa_normalised_interface_proportion[x] = 0
			except:
				sa_normalised_interface[x] = 0
				sa_normalised_interface_proportion[x] = 0

		sa_normalised = sa_normalised_diff

		mean_normalised_surface_accessibility = sum(sa_normalised.values())/len(sa_normalised) if pdb_chain in dssp_monomers else 0
		
		chain_details = {}
		chain_details['mean_surface_accessibility'] = mean_surface_accessibility
		chain_details['surface_accessibility_normalised'] = sa_normalised_splitted
		chain_details['surface_accessibility_normalised_complex'] = sa_normalised_joined
		chain_details['surface_accessibility_normalised_interface'] = sa_normalised_splitted
		chain_details['surface_accessibility_normalised_difference'] = sa_normalised_interface
		chain_details['surface_accessibility_normalised_interface_proportion'] = sa_normalised_interface_proportion
		chain_details['sequence'] = dssp_monomers.chain(pdb_chain).get_summary().sequence
		chain_details['offsets'] = dssp_monomers.chain(pdb_chain).get_summary().offsets
		chain_details['secondary_structure'] = dssp_monomers.chain(pdb_chain).get_summary().secondary_structure
		chain_details['accessible_residues'] = dssp_monomers.chain(pdb_chain).get_summary().accessible_residues
		chain_details['buried_residues'] = dssp_monomers.chain(pdb_chain).get_summary().buried_residues
		chain_details['mean_surface_accessibility'] = dssp_monomers.chain(pdb_chain).get_summary().mean_surface_accessibility
		chain_details['mean_normalised_surface_accessibility'] = mean_normalised_surface_accessibility
		chain_details['dssp_uniprot_mapping'] = dssp_uniprot_mapping

		return chain_details

	def get_surface_accessibility_data(self):
			
		pdb_id = self.get_pdb_id()
		pdb_structure_data = queryRunner.queryRunner("pdb","download_pdb_structure",{"pdb_id":pdb_id}).run()
		
		pdb_file = pdb_structure_data['data']
		dssp_runner = Dssp(self.options['dssp_path'])

		try:
			dssp_monomers, _ = dssp_runner.run(pdb_file, pdb_id, splitted=True)
			dssp_joined, _ = dssp_runner.run(pdb_file, pdb_id, splitted=False)
			dssp_diff = dssp_joined.solvent_accesibility_difference(dssp_monomers)
		except:
			logger.error("ERROR: There was a problem running DSSP for pdb: {}".format(pdb_file))
			raise

		if 'chain_id' in self.options:
			pdb_chains = [pdb_id + "." + chain_id for chain_id in self.options['chain_id']]
		else:
			pdb_chains = list(dssp_diff.to_dict().keys())
		

		if self.options['residue_numbering'] == "pdb":
			self.data = {
			"chain_details":{}
			}

			for pdb_chain in pdb_chains:
				self.data['chain_details'][pdb_chain] = self.get_dssp_chain_data(pdb_id,pdb_chain,dssp_diff,dssp_monomers,dssp_joined)
			
			return self.data

		elif self.options['residue_numbering'] == "uniprot":
			self.data = {
				"chain_details":{}
			}
			for pdb_chain in pdb_chains:
				self.data['chain_details'][pdb_chain] = self.get_dssp_chain_data_protein_centric(pdb_id,pdb_chain,dssp_diff,dssp_monomers,dssp_joined)

			return self.data
		else:
			return {"status": "Error","error_type":"You shouldn't be able to reach here!"}

	def get_surface_accessibility_data_bulk(self):
		sa_data = {}
		pdb_ids = self.options["pdb_id"]
		for pdb_id in pdb_ids:
			self.options['pdb_id'] = pdb_id
			sa_data[pdb_id] = self.get_surface_accessibility_data()

		return sa_data
			
	def get_dssp_chain_data_protein_centric(self,pdb_id,pdb_chain,dssp_diff,dssp_monomers,dssp_joined):
		chain_details = self.get_dssp_chain_data(pdb_id,pdb_chain,dssp_diff,dssp_monomers,dssp_joined)
		
		protein_details = {}
		name_mapping = {
				'surface_accessibility_normalised':"sa",
				'surface_accessibility_normalised_complex':"sa_complex",
				'surface_accessibility_normalised_interface':"sa_monomer",
				'surface_accessibility_normalised_difference':"sa_interfaces",
				'surface_accessibility_normalised_interface_proportion':"sa_interfaces_prop"
			}
		data_sources = ['surface_accessibility_normalised',
		'surface_accessibility_normalised_complex',
		'surface_accessibility_normalised_interface',
		'surface_accessibility_normalised_difference',
		'surface_accessibility_normalised_interface_proportion']
	
		for pdb_offset in chain_details['dssp_uniprot_mapping']:
			uniprot_accession = list(chain_details['dssp_uniprot_mapping'][pdb_offset].keys())[0]
			uniprot_offset = chain_details['dssp_uniprot_mapping'][pdb_offset][uniprot_accession]['uniprot_offset']
			if uniprot_accession not in protein_details:
				protein_details[uniprot_accession] = {}

			protein_details[uniprot_accession][uniprot_offset] = chain_details['dssp_uniprot_mapping'][pdb_offset][uniprot_accession]
			protein_details[uniprot_accession][uniprot_offset]['pdb_offset'] = pdb_offset

			for data_source in data_sources:
				protein_details[uniprot_accession][uniprot_offset][name_mapping[data_source]] = chain_details[data_source][pdb_offset]

		return protein_details

	##-----------------##
	##-----------------##


	def get_dssp_attribute_by_offset(self,position_data,attribute):
		dssp_attribute_by_offset = {}
		positions = list(position_data.keys())

		for position in positions:
			dssp_attribute_by_offset[position] = position_data[position][attribute]

		return dssp_attribute_by_offset

	def get_pdb_secondary_structure(self):
			
		pdb_id = self.get_pdb_id()
		pdb_structure_data = queryRunner.queryRunner("pdb","download_pdb_structure",{"pdb_id":pdb_id}).run()
		
		pdb_file = pdb_structure_data['data']
		dssp_runner = Dssp(self.options['dssp_path'])

		try:
			dssp_info, _ = dssp_runner.run(pdb_file, pdb_id, splitted=True)
		except:
			logger.error("ERROR: There was a problem running DSSP for pdb: {}".format(pdb_file))
			raise

		if 'chain_id' in self.options:
			pdb_chains = [pdb_id + "." + chain_id for chain_id in self.options['chain_id']]
		else:
			pdb_chains = list(dssp_info.to_dict().keys())
		
		self.data = {'secondary_structure':{}}
		try:
			for pdb_chain in pdb_chains:
				chain_data = dssp_info.to_dict()[pdb_chain]['summary']
				ss_data = self.define_ss(chain_data['secondary_structure'])
				ss_data['secondary_structure'] = self.get_dssp_attribute_by_offset(dssp_info.to_dict()[pdb_chain]['positions'],"SecondaryStructure")
				ss_data["chain_length"] = len(chain_data['surface_accessibility'].keys())
				self.data['secondary_structure'][pdb_chain] = ss_data
			
			return self.data
		except:
			return {"status": "Error","error_type":"You shouldn't be able to reach here!"}

	def define_ss(self,secondary_structure):
		
		helices = {}
	
		last = secondary_structure[0]
		counter = 1

		compressed_ss = ""
		long_compressed_ss = []
		helix = []
		positions = []

		for i in range(0,len(secondary_structure)):
			if secondary_structure[i]  == last:
				counter += 1
				positions.append(i)
			else:
				compressed_ss += last + str(counter)
		
				if counter >= 9 and last == "H":
					helix.append([last , str(counter),min(positions),max(positions)])
					helix.append([last , str(counter),min(positions),max(positions)])
			
				if (counter >= 4 or (counter >= 3 and last == 'T' )) and last != "-" and last != "" :
					long_compressed_ss.append([last , str(counter),min(positions),max(positions)])
					positions = []
		
				counter = 1
	
			last = secondary_structure[i]
		
		if (counter >= 4 or (counter >= 3 and last == 'T' )) and last != "-" and last != "" :
			long_compressed_ss.append([last , str(counter),min(positions),max(positions)])
			
		helices = helix
		
		ss = "undefined"
		ss_length = "undefined"

		if len(long_compressed_ss) > 1:
			ss = "Multi-partite"
	
		elif len(long_compressed_ss) == 1:
			if long_compressed_ss[0][0] in ["H"]:
				ss = "Alpha-helix"
		
			if long_compressed_ss[0][0] in ["E"]:
				ss = "Strand"
	
			if long_compressed_ss[0][0] in ["B"]:
				ss = "Beta-bridge"
		
			if long_compressed_ss[0][0] in ["G"]:
				ss = "3/10-helix"
		
			if long_compressed_ss[0][0] in ["I"]:
				ss = "Helix-5"
		
			if long_compressed_ss[0][0] in ["T"]:
				ss = "Turn"
		
			if long_compressed_ss[0][0] in ["S"]:
				ss = "Bend"

			ss_length = long_compressed_ss[0][1]
	
		elif len(long_compressed_ss) == 0:
			ss = "Extended"
			ss_length = len(secondary_structure)
			if secondary_structure.count("TT") > 0:
				ss = "Turn"
				ss_length = secondary_structure.count("T") 
				
		ss_composition = {}
		for dssp_code in ['B', 'E', 'G', 'H', 'S', 'T']:
			ss_composition[dssp_code] = secondary_structure.count(dssp_code)

		secondary_structure_info = {
		"ss_type":ss,
		"ss_length":ss_length,
		"ss":secondary_structure,
		"ss_detailed":long_compressed_ss,
		"helices":helices
		}

		return secondary_structure_info

	##-----------------##
	##-----------------##
	
	def get_pdb_uniprot_mapping(self):
		pdb_mapping = {}

		pdb_id = self.get_pdb_id()
		pdb_structure_data = queryRunner.queryRunner("pdb","download_pdb_structure",{"pdb_id":pdb_id}).run()
		
		pdb_file = pdb_structure_data['data']

		dssp_runner = Dssp(self.options['dssp_path'])

		try:
			dssp_monomers, _ = dssp_runner.run(pdb_file, pdb_id, splitted=True)
		except:
			logger.error("ERROR: There was a problem running DSSP for pdb: {}".format(pdb_file))
			raise

		if 'chain_id' in self.options:
			pdb_chains = [pdb_id + "." + chain_id for chain_id in self.options['chain_id']]
		else:
			pdb_chains = list(dssp_monomers.to_dict().keys())
		
		if self.options['residue_numbering'] == "pdb":
			self.data = {
			"chain_details":{}
			}

			for pdb_chain in pdb_chains:
				offsets = dssp_monomers.chain(pdb_chain).get_summary().offsets
				pdb_mapping[pdb_chain] = self.map_dssp_offsets_to_uniprot(pdb_id,pdb_chain,offsets,return_alignment=True)
			#	print(pdb_id,pdb_chain,pdb_mapping[pdb_chain])
			return pdb_mapping
		else:
			return {"status": "Error","error_type":"You shouldn't be able to reach here!"}

if __name__ == "__main__":
	accessibilityUtilitiesOj = dsspUtilities()

	pdb_id = "2AST"
	#pdb_file = "/Users/adminndavey/Documents/Work/projects/Cln2/Figure/ScerCln2_1-376delIns mod1.pdb"
	dssp_data = accessibilityUtilitiesOj.get_dssp_data('2AST')

	accessibilityUtilitiesOj.options['pdb_id'] = pdb_id
	dssp_data = accessibilityUtilitiesOj.get_pdb_secondary_structure()
	
	pprint.pprint(dssp_data)

	#offsets = dssp_data[pdb_id + ".A"]['summary']['surface_accessibility_normalised']
	#accessible_score_sequence = ""
	#for offset in range(min(offsets),max(offsets)):
	#	accessible_score_sequence += str(int(dssp_data[pdb_id + ".A"]['summary']['surface_accessibility_normalised'][offset]*10))

	#print(dssp_data[pdb_id + ".A"]['summary']['sequence'])
	#print(accessible_score_sequence)
	#print("".join(dssp_data[pdb_id + ".A"]['summary']['accessible_residues']))
