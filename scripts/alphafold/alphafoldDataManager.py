import sys,os,inspect,pprint,glob,json,shutil,math

import Bio
import numpy as np

file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../"))
sys.path.append(os.path.abspath(os.path.join(file_path,"../utilities")))
sys.path.append(os.path.abspath(os.path.join(file_path,"../accessibility")))
sys.path.append(os.path.abspath(os.path.join(file_path,"../data_management")))

try:
	sys.path.append(os.path.join(file_path,"../structure/structure_utilities/"))
	from structure_atomic_coordinate_utilities import structureAtomicCoordinateUtilities
except:
	print("Skipping structure_atomic_coordinate_utilities import")
	
import queryRunner

import utilities_error
import utilities_downloader
import utilities_basic

import config_reader

try:
	from dssp_utilities import dsspUtilities
except:
	print("Skipping dssp_utilities import")

try:	
	from simple_mapping_and_accessibility import simpleMappingAccessibilityUtilities
except:
	print("Skipping simple_mapping_and_accessibility import")

try:
	from Bio.PDB import PDBParser
	from Bio.PDB.Polypeptide import *
except:
	print("install Biopython")
	sys.exit()

try:
	import networkx
except:
	print("install networkx")
	sys.exit()

try:
	from Bio.PDB.Polypeptide import three_to_one
except:
	from Bio.PDB.Polypeptide import three_to_index
	from Bio.PDB.Polypeptide import index_to_one

from itertools import groupby
import numpy

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

#-----
import logging
logger = logging.getLogger(__name__)
#-----

class alphafoldDataManager():

	#####------------------------------------------------------------------#####

	def __init__(self):
		self.options = config_reader.load_configeration_options(sections=['general'])
		self.options.update({})

		self.options["smooth_termini"] = False
		self.options["show_smoothing"] = False
		self.options['alphafold_coil_cutoff'] = 30
		self.options["alphafold_short_loop_length_cutoff"] = 10
		self.options["accessibility_score_cutoff"] = 0.55
		self.options["intramolecular_contact_score_cutoff"] = 8
		self.options["smooth_order_length_cutoff"] = 15
		self.options["smooth_disorder_length_cutoff"] = 25
		self.options['adjacent_residue_filter_distance'] = 3
		self.options['alphafold_model_version'] = "4"
		self.options['use_pae_filtering_for_intramolecular_contacts_counts'] = False
		self.options['pae_filtering_for_intramolecular_contacts_counts_cut_off'] = 10
		self.options["filter_inaccessibility"] = True
		self.options["add_atom_details"] = False
		
		self.options['side_chain_contacts_only'] = False
		
		self.options['alphafold_data_path'] = os.path.join(self.options['data_path'],"alphafold")
		self.options['structure_modules_data_path'] = os.path.join(self.options['data_path'],"alphafold","structure_modules")
		self.options['accessibility_data_path'] = os.path.join(self.options['data_path'],"alphafold","accessibility")
		self.options['classification_data_path'] = os.path.join(self.options['data_path'],"alphafold","classification")
		self.options['pdb_path'] = os.path.join(self.options['data_path'],"alphafold","pdb")
		self.options['secondary_structure_path'] = os.path.join(self.options['data_path'],"alphafold","secondary_structure")
		self.options['pLDDT_data_path'] = os.path.join(self.options['data_path'],"alphafold","pLDDT")
		self.options['intramolecular_contacts_data_path'] = os.path.join(self.options['data_path'],"alphafold","intramolecular_contacts")
		self.options['predicted_aligned_error_data_path'] = os.path.join(self.options['data_path'],"alphafold","predicted_aligned_error")

		############################################################################

		if not(os.path.exists(self.options['alphafold_data_path'])):
			os.mkdir(self.options['alphafold_data_path'])

		if not(os.path.exists(self.options['accessibility_data_path'])):
			os.mkdir(self.options['accessibility_data_path'])

		if not(os.path.exists(self.options['structure_modules_data_path'])):
			os.mkdir(self.options['structure_modules_data_path'])

		if not(os.path.exists(self.options['classification_data_path'])):
			os.mkdir(self.options['classification_data_path'])

		if not(os.path.exists(self.options['pdb_path'])):
			os.mkdir(self.options['pdb_path'])

		if not(os.path.exists(self.options['secondary_structure_path'])):
			os.mkdir(self.options['secondary_structure_path'])

		for directory in ['raw','windowed','species','tdt','smoothed', 'tessellation']:
			if not os.path.exists( os.path.join(self.options['accessibility_data_path'],directory)):
				os.mkdir(os.path.join(self.options['accessibility_data_path'],directory))

		if not(os.path.exists(self.options['pLDDT_data_path'])):
			os.mkdir(self.options['pLDDT_data_path'])

		for directory in ['raw','windowed','species','tdt','scaled']:
			if not os.path.exists( os.path.join(self.options['pLDDT_data_path'],directory)):
				os.mkdir(os.path.join(self.options['pLDDT_data_path'],directory))

		if not(os.path.exists(self.options['intramolecular_contacts_data_path'])):
			os.mkdir(self.options['intramolecular_contacts_data_path'])

		for directory in ['smoothed',"windowed"]:
			if not os.path.exists( os.path.join(self.options['intramolecular_contacts_data_path'],directory)):
				os.mkdir(os.path.join(self.options['intramolecular_contacts_data_path'],directory))

		if not(os.path.exists(self.options['predicted_aligned_error_data_path'])):
			os.mkdir(self.options['predicted_aligned_error_data_path'])

		for directory in ['raw','plots','degree']:
			if not os.path.exists( os.path.join(self.options['predicted_aligned_error_data_path'],directory)):
				os.mkdir(os.path.join(self.options['predicted_aligned_error_data_path'],directory))


	############################################################################
	############################################################################
	############################################################################

	def package_alphafold_by_species(self):

		self.options['taxon_id'] = "9606"
		self.options["reviewed"] = True
		self.options["structure_required"] = False

		accessions = queryRunner.queryRunner("uniprot","parse_uniprot_accession_taxa",{"taxon_id":self.options['taxon_id'] ,"reviewed":self.options["reviewed"],"structure_required":self.options["structure_required"]}).run()
		self.options['accession'] = accessions['data']
		logger.debug("Processing " + str(len(self.options['accession'])) + " protein")

		"""
		for filename in os.listdir(os.path.join(self.options['data_path'],"alphafold","pLDDT","raw")):
			if filename.count('accessibility') > 0:
				os.remove(os.path.join(self.options['data_path'],"alphafold","pLDDT","raw",filename))

		for filename in os.listdir(os.path.join(self.options['data_path'],"alphafold","accessibility","raw")):
			if filename.count('pLDDT') > 0:
				os.remove(os.path.join(self.options['data_path'],"alphafold","accessibility","raw",filename))
		"""

			
		try:
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","structure_modules")): os.mkdir(os.path.join(self.options['data_path'],"alphafold","structure_modules"))
		except:
			pass

		try:
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","pLDDT","species")): os.mkdir(os.path.join(self.options['data_path'],"alphafold","pLDDT","species"))
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","pLDDT","species",self.options['taxon_id'])): os.mkdir(os.path.join(self.options['data_path'],"alphafold","pLDDT","species",self.options['taxon_id']))
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","pLDDT","species",self.options['taxon_id'],"raw")): os.mkdir(os.path.join(self.options['data_path'],"alphafold","pLDDT","species",self.options['taxon_id'],"raw"))
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","pLDDT","species",self.options['taxon_id'],"windowed")): os.mkdir(os.path.join(self.options['data_path'],"alphafold","pLDDT","species",self.options['taxon_id'],"windowed"))
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","pLDDT","species",self.options['taxon_id'],"scaled")): os.mkdir(os.path.join(self.options['data_path'],"alphafold","pLDDT","species",self.options['taxon_id'],"scaled"))
		except:
			pass

		try:
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","accessibility","species")): os.mkdir(os.path.join(self.options['data_path'],"alphafold","accessibility","species"))
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id'])): os.mkdir(os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id']))
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id'],"raw")): os.mkdir(os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id'],"raw"))
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id'],"windowed")): os.mkdir(os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id'],"windowed"))
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id'],"smoothed")): os.mkdir(os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id'],"smoothed"))
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id'],"intramolecular_contacts")): os.mkdir(os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id'],"intramolecular_contacts"))
		except:
			pass

		try:
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","secondary_structure","species")): os.mkdir(os.path.join(self.options['data_path'],"alphafold","secondary_structure","species"))
			if not os.path.exists(os.path.join(self.options['data_path'],"alphafold","secondary_structure","species",self.options['taxon_id'])): os.mkdir(os.path.join(self.options['data_path'],"alphafold","secondary_structure","species",self.options['taxon_id']))
		except:
			pass

		counter= {
			"total":len(self.options['accession']),
			"accessibility_out_file_json":0,
			"accessibility_windowed_out_file_json":0,
			"pLDDT_file_json":0,
			"pLDDT_windowed_file_json":0,
			"secondary_structure_file_json":0
		}

		for accession in self.options['accession']:
			logger.debug("Moving " + accession)
			
			accessibility_out_file_json = os.path.join(self.options['data_path'],"alphafold","accessibility","raw", accession + ".alphafold.accessibility.json")
			accessibility_windowed_out_file_json = os.path.join(self.options['data_path'],"alphafold","accessibility","windowed", accession + ".alphafold.windowed_" + str(self.options['window_size']) + ".accessibility.json")
			pLDDT_file_json = os.path.join(self.options['data_path'],"alphafold","pLDDT","raw", accession + ".alphafold.pLDDT.json")
			pLDDT_windowed_file_json = os.path.join(self.options['data_path'],"alphafold","pLDDT","windowed", accession + ".alphafold.windowed_" + str(self.options['window_size']) + ".pLDDT.json")
			pLDDT_file_json = os.path.join(self.options['data_path'],"alphafold","pLDDT","scaled", accession + ".alphafold.scaled.pLDDT.json")
			secondary_structure_file_json = os.path.join(self.options['data_path'],"alphafold","secondary_structure", accession + ".alphafold.ss.json")

			species_accessibility_out_file_json = os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id'],"raw", accession + ".alphafold.accessibility.json")
			species_accessibility_windowed_out_file_json = os.path.join(self.options['data_path'],"alphafold","accessibility","species",self.options['taxon_id'],"windowed", accession + ".alphafold.windowed_" + str(self.options['window_size']) + ".accessibility.json")

			species_pLDDT_file_json = os.path.join(self.options['data_path'],"alphafold","pLDDT","species",self.options['taxon_id'],"raw", accession + ".alphafold.pLDDT.json")
			species_pLDDT_windowed_file_json = os.path.join(self.options['data_path'],"alphafold","pLDDT","species",self.options['taxon_id'],"windowed", accession + ".alphafold.windowed_" + str(self.options['window_size']) + ".pLDDT.json")
			species_pLDDT_file_json = os.path.join(self.options['data_path'],"alphafold","pLDDT","species",self.options['taxon_id'],"scaled", accession + ".alphafold.scaled.pLDDT.json")

			species_secondary_structure_file_json = os.path.join(self.options['data_path'],"alphafold","secondary_structure","species",self.options['taxon_id'],accession + ".alphafold.ss.json")


			if os.path.exists(accessibility_out_file_json) and not os.path.exists(species_accessibility_out_file_json):
				shutil.copy(accessibility_out_file_json,species_accessibility_out_file_json)
				counter["accessibility_out_file_json"] += 1

			if os.path.exists(accessibility_windowed_out_file_json) and not os.path.exists(species_accessibility_windowed_out_file_json):
				shutil.copy(accessibility_windowed_out_file_json,species_accessibility_windowed_out_file_json)
				counter["accessibility_windowed_out_file_json"] += 1

			if os.path.exists(pLDDT_file_json) and not os.path.exists(species_pLDDT_file_json):
				shutil.copy(pLDDT_file_json,species_pLDDT_file_json)
				counter["pLDDT_file_json"] += 1

			if os.path.exists(pLDDT_windowed_file_json) and not os.path.exists(species_pLDDT_windowed_file_json):
				shutil.copy(pLDDT_windowed_file_json,species_pLDDT_windowed_file_json)
				counter["pLDDT_windowed_file_json"] += 1

			if os.path.exists(secondary_structure_file_json) and not os.path.exists(species_secondary_structure_file_json):
				shutil.copy(secondary_structure_file_json,species_secondary_structure_file_json)
				counter["secondary_structure_file_json"] += 1

		return counter

	def make_alphafold_file_by_species(self):
		print("make_alphafold_file_by_species")
		self.options['taxon_id'] = "9606"
		self.options["reviewed"] = True
		self.options["structure_required"] = False

		accessions = queryRunner.queryRunner("uniprot","parse_uniprot_accession_taxa",{"taxon_id":self.options['taxon_id'] ,"reviewed":self.options["reviewed"],"structure_required":self.options["structure_required"]}).run()
		self.options['accession'] = accessions['data']
		#self.options['accession'] = ["P07306","P07311","P07332","P07357","P07359","P07360"]
		"""
		scores = self.calculate_scores()

		for source in ["pLDDT_scaled","pLDDT","accessibility","intramolecular_contacts","secondary_structure"]:
			out_str = ""
			for accession in scores:
				try:
					if isinstance(scores[accession][source],(list)):
						try:
							out_str += accession + "\t" +  ",".join(["%1.3f"%(x) for x in scores[accession][source]]) + "\n"
						except:
							out_str += accession + "\t" +  ",".join(scores[accession][source]) + "\n"
					else:
						out_str += accession + "\t" +  scores[accession] + "\n"
				except:
					utilities_error.printError()
					out_str += accession + "\tmissing\n"

			open(self.options['taxon_id'] + "." + source + '.tdt',"w").write(out_str)
			logger.debug("Writing: " + self.options['taxon_id'] + "." + source + '.tdt')

		"""
		self.options['window_size'] = 20
		for source in ["pLDDT_windowed","accessibility_windowed","pLDDT","accessibility"]:#,"weighted_degree","degree","intramolecular_contacts","pLDDT_scaled","accessibility_smoothed"]:
			
			import matplotlib.pyplot as plt
			plt.figure(figsize= (20,20))

			total_scores = []
			out_str = ""
			counter = 0

			for accession in accessions['data']:
				counter += 1
				logger.debug(counter)
				try:

					sequence_response = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":accession}).run()
					
					###----------------------------###

					if 'data' not in sequence_response: continue
					sequence = sequence_response['data']
					self.options['sequence'] = sequence
					self.options['accession'] = [accession]
					scores = []
					if source == "pLDDT":
						scores = self.get_pLDDT(accession)
					if source == "pLDDT_windowed":
						scores = self.get_pLDDT_windowed_scores(accession)
					if source == "pLDDT_scaled":
						scores =self.get_pLDDT_scaled_scores(accession)
					if source == "accessibility":
						scores =self.get_accessibility(accession)
					if source == "accessibility_windowed":
						scores =self.get_accessibility_windowed_scores(accession)
					if source == "accessibility_smoothed":
						scores =self.get_accessibility_smoothed_scores(accession)
					if source == "secondary_structure":
						scores =self.get_secondary_structure_scores(accession)
					if source == "degree":
						if len(self.options['sequence']) < 2000:
							scores =self.get_degree(accession)
						else:
							scores = []

					###----------------------------###
					
					total_scores += scores

					if len(sequence) != len(scores):
						logger.debug(["sequence length mismatch",counter,accession,len(sequence) != len(scores),len(sequence),len(scores)])
						continue

					if source == "weighted_degree":
						if len(self.options['sequence']) < 2000:
							scores =self.get_weighted_degree(accession)
						else:
							scores = []
					if source == "intramolecular_contacts":
						scores =self.get_intramolecular_contacts_scores(accession)

					if isinstance(scores,(list)):
						out_str += accession + "\t" +  ",".join(["%1.2f"%(x) for x in scores]) + "\n"
					else:
						out_str += accession + "\t" +  scores + "\n"

					###----------------------------###
					
					del scores
					del sequence_response
				except:
					pass
			
			logger.debug("Saving layout")
			plt.hist(total_scores, bins = 50)

			if source in ["accessibility_windowed","pLDDT_windowed"]:
				png_out_path = os.path.join(self.options['data_path'],"alphafold",self.options['taxon_id'] + "." + source + "." + str(self.options['window_size']) + '.png')
			else:
				png_out_path = os.path.join(self.options['data_path'],"alphafold",self.options['taxon_id'] + "." + source + '.png')

			logger.debug("Writing " + png_out_path)
			plt.savefig(png_out_path)
			
			if source in ["accessibility_windowed","pLDDT_windowed"]:
				tdt_out_path = os.path.join(self.options['data_path'],"alphafold",self.options['taxon_id'] + "." + source + "." + str(self.options['window_size']) + '.tdt')
			else:
				tdt_out_path = os.path.join(self.options['data_path'],"alphafold",self.options['taxon_id'] + "." + source + '.tdt')
		
			open(tdt_out_path,"w").write(out_str)

			
			#logger.debug("Writing: " + self.options['taxon_id'] + "." + source + "." + str(self.options['window_size']) + '.tdt')


	def benchmark_alphafold_accessibility(self):

		self.options['taxon_id'] = "9606"
		self.options["reviewed"] = True
		self.options["structure_required"] = False

		accessions = queryRunner.queryRunner("uniprot","parse_uniprot_accession_taxa",{"taxon_id":self.options['taxon_id'] ,"reviewed":self.options["reviewed"],"structure_required":self.options["structure_required"]}).run()
		self.options['accession'] = accessions['data']
		logger.debug("Processing " + str(len(self.options['accession'])) + " protein")

		scores = self.calculate_scores()

		score_distribution = {}
		score_distribution_total = 0

		for i in range(0,11):
			score_distribution[i] = 0

		for accession in scores:
			for i in range(0,len(scores[accession])):
				score_distribution[int(scores[accession][i]*10)] += 1
				score_distribution_total += 1

		for i in range(0,11):
			proportion = float(score_distribution[i])/score_distribution_total
			print(i,"\t",score_distribution[i],"\t","%1.3f"%(proportion),"*"*int(proportion*20))
		return {}

	def benchmark_alphafold_disprot(self):
		out_file_json_path = os.path.join(self.options['data_path'],"alphafold","disprot_benchmarking.json")
		#adjacent_residue_filter_distances = [3,0]
		#intramolecular_contact_score_cutoffs = range(0,15)


		#for adjacent_residue_filter_distance in adjacent_residue_filter_distances:
			#for intramolecular_contact_score_cutoff in intramolecular_contact_score_cutoffs:
				#self.options["intramolecular_contact_score_cutoff"] = intramolecular_contact_score_cutoff
				#self.options['adjacent_residue_filter_distance'] = adjacent_residue_filter_distance

		testing = True

		distributions = {}
		if not os.path.exists(out_file_json_path) or self.options['remake'] or testing:
			disprot_accession_data = queryRunner.queryRunner("pdb","get_disprot_index_accessions",{}).run()
			disprot_data = queryRunner.queryRunner("pdb","get_disprot_by_accession",{'accession':disprot_accession_data['data']}).run()

			regions = {}
			for accession in  disprot_accession_data['data'][0:]:
				try:
					self.options['accession'] = [accession]

					scores = self.calculate_scores()
					#print(accession)
					if accession  not in scores: continue
					#pprint.pprint(disprot_data['data'][accession])

					parse_basic_response = queryRunner.queryRunner("uniprot","parse_basic",{"accession":accession}).run()


					for region in disprot_data['data'][accession]['regions']:
						if region['term_namespace'] == 'Structural state' or region['term_namespace'] == 'Structural transition':
							region_id = accession + "_" + str(region['start']) + "-" + str(region['end'])

							regions[region_id] = {}
							regions[region_id]['accession'] = accession
							regions[region_id]['gene_name'] = parse_basic_response['data']['gene_name']
							regions[region_id]['protein_name'] = parse_basic_response['data']['protein_name']
							regions[region_id]['regions'] = {}

							regions[region_id] = {
								'accession':accession,
								'gene_name':parse_basic_response['data']['gene_name'],
								'protein_name':parse_basic_response['data']['protein_name'],
								"start":str(region['start']),
								"end":str(region['end']),
								"length":str(region['end'] - region['start']),
								"region_id":region['region_id'],
								"term_namespace":region['term_namespace'],
								"term_name":region['term_name'],
								"proviz_link":"http://slim.icr.ac.uk/proviz/proviz.php?tools=alphafold&uniprot_acc=" + accession + "&ali_start=" + str(max(0,region['start'] - 10)) + "&ali_end=" + str(region['end'] + 10) + "&tracks=peptides,Disordered Region,AAA,DisProt:" + region['region_id'] + "," + str(region['start']) + "," + str(region['end'])
							}

							for score_type in ['pLDDT','accessibility',"accessibility_smoothed",'intramolecular_contacts','degree','secondary_structure']:
								if score_type in scores[accession]:
									if score_type in ["secondary_structure"]:
										regions[region_id][score_type] = {
											"-":float(scores[accession][score_type][region['start']:region['end']].count('-'))/len(scores[accession][score_type][region['start']:region['end']]),
											"H":float(scores[accession][score_type][region['start']:region['end']].count('H'))/len(scores[accession][score_type][region['start']:region['end']]),
											"E":float(scores[accession][score_type][region['start']:region['end']].count('E'))/len(scores[accession][score_type][region['start']:region['end']]),
											"G":float(scores[accession][score_type][region['start']:region['end']].count('G'))/len(scores[accession][score_type][region['start']:region['end']]),
											"S":float(scores[accession][score_type][region['start']:region['end']].count('S'))/len(scores[accession][score_type][region['start']:region['end']]),
											"T":float(scores[accession][score_type][region['start']:region['end']].count('T'))/len(scores[accession][score_type][region['start']:region['end']]),
										}
									elif score_type in ["intramolecular_contacts_smoothed","accessibility_smoothed"]:
										if score_type + "_disorder_proportion" not in distributions: distributions[score_type + "_disorder_proportion"] = []
										#print(score_type,"\t",float(scores[accession][score_type][region['start']:region['end']].count(1))/len(scores[accession][score_type][region['start']:region['end']]),"\t",scores[accession][score_type][region['start']:region['end']])
										mean_score = float(scores[accession][score_type][region['start']:region['end']].count(1))/len(scores[accession][score_type][region['start']:region['end']])
										distributions[score_type + "_disorder_proportion"].append("%1.3f"%mean_score)
										regions[region_id][score_type + "_disorder_proportion"] = "%1.3f"%mean_score
									else:
										if score_type + "_mean_score" not in distributions: distributions[score_type + "_mean_score"] = []
										mean_score = float(sum(scores[accession][score_type][region['start']:region['end']]))/len(scores[accession][score_type][region['start']:region['end']])
										distributions[score_type + "_mean_score"].append("%1.3f"%mean_score)
										regions[region_id][score_type + "_mean_score"] = "%1.3f"%mean_score
				except:
					utilities_error.printError()


			logger.debug("Writing " + out_file_json_path)
			utilities_basic.write_to_json(out_file_json_path,regions)
		else:
			logger.debug("Reading " + out_file_json_path)
			regions = json.loads(open(out_file_json_path, 'r').read())


		for distribution in distributions:
			distributions[distribution].sort()

		row = [
				'accession',
				'gene_name',
				'region_id',
				'term_namespace',
				'term_name',
				'start',
				'end',
				'length'
		]

		for distribution in distributions:
			row += [distribution,distribution + " rank"]

		for ss in ["-","H","E","G","S","T"]:
			row.append("ss " + ss + '%')

		print("\t".join(row))

		for region_id in regions:
			try:
				row = [
					regions[region_id]['accession'],
					regions[region_id]['gene_name'],
					regions[region_id]['region_id'],
					regions[region_id]['term_namespace'],
					regions[region_id]['term_name'],
					regions[region_id]['start'],
					regions[region_id]['end'],
					regions[region_id]['length']
				]
				for distribution in distributions:
					row.append(regions[region_id][distribution])
					row.append(str(distributions[distribution].index(regions[region_id][distribution])+1))

				for ss in ["-","H","E","G","S","T"]:
					row.append("%1.1f"%(regions[region_id]["secondary_structure"][ss]*100))

				print("\t".join(row))
			except Exception as e:
				logger.debug(str(e))
		#pprint.pprint(regions_response)

		return regions

	def benchmark_alphafold_disorder(self):
		positive_testing_set = open(os.path.join(self.options['data_path'],"alphafold", "UniProt_accs-IUPred2-positive_testing_set.txt")).read().strip().split("\n")
		negative_testing_set = open(os.path.join(self.options['data_path'],"alphafold", "UniProt_accs-IUPred2-negative_testing_set.txt")).read().strip().split("\n")


		testset_accessions = positive_testing_set + negative_testing_set

		source = "intramolecular_contacts_smoothed"
		sources = ["get_degree"]

		counter = 1

		for source in sources:
			out_str = ""
			for accession in testset_accessions:
				try:
					counter += 1
					logger.debug([counter,accession])
					self.options['accession'] = [accession]

					sequence_response = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":accession}).run()

					if 'data' not in sequence_response: continue
					sequence = sequence_response['data']
					self.options['sequence'] = sequence

					if source == "intramolecular_contacts_smoothed":
						scores = self.get_intramolecular_contacts_smoothed_scores(accession)

					if source == "get_degree":
						scores = self.get_degree(accession)

					if len(scores) > 0:
						try:
							if isinstance(scores,(list)):
								out_str += accession + "\t" +  ",".join(["%1.3f"%x for x in scores]) + "\n"
							else:
								out_str += accession + "\t" +  scores + "\n"
						except:
							print("ERROR1")

					del scores
				except:
					print("ERROR2")

			print(out_str)
			open("testing_set." + source + ".tdt","w").write(out_str)
			logger.debug("Writing: " + "testing_set." + source + ".tdt")

		sys.exit()

		if source in ["intramolecular_contacts_smoothed"]:
			for intramolecular_contact_score_cutoff in range(4,13):
				out_str = ""
				self.options['intramolecular_contact_score_cutoff'] = intramolecular_contact_score_cutoff

				counter = 1

				for accession in testset_accessions:
					try:
						counter += 1
						logger.debug([counter, intramolecular_contact_score_cutoff,len(positive_testing_set + negative_testing_set)])
						self.options['accession'] = [accession]

						sequence_response = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":accession}).run()

						if 'data' not in sequence_response: continue
						sequence = sequence_response['data']
						self.options['sequence'] = sequence


						if source == "intramolecular_contacts_smoothed":
							scores = self.get_intramolecular_contacts_smoothed_scores(accession)


						if source == "get_degree":
							scores = self.get_degree(accession)

						if len(scores) > 0:
							try:
								if isinstance(scores,(list)):
									out_str += accession + "\t" +  ",".join(["%1.3f"%(x) for x in scores]) + "\n"
								else:
									out_str += accession + "\t" +  scores + "\n"
							except:
								print("ERROR")

						#	raise
						del scores
					except:
						print("ERROR")

				print(out_str)
				open("testing_set." + source + ".intramolecular_contact_score_cutoff_" + str(intramolecular_contact_score_cutoff) + '.tdt',"w").write(out_str)
				logger.debug("Writing: " + "testing_set." + source + ".intramolecular_contact_score_cutoff_" + str(intramolecular_contact_score_cutoff) + '.tdt')

		source = "intramolecular_contacts_windowed"
		if source in ["intramolecular_contacts_windowed"]:
			for window_size in [5,10,15,20,25,30,35,40,45,50]:
				out_str = ""
				self.options['window_size'] = window_size
				counter = 1

				for accession in testset_accessions:
					try:
						counter += 1
						logger.debug([counter, intramolecular_contact_score_cutoff,len(positive_testing_set + negative_testing_set)])
						self.options['accession'] = [accession]

						if source == "intramolecular_contacts_windowed":
							scores = self.get_intramolecular_contacts_windowed_scores(accession)

						if len(scores) > 0:
							try:
								if isinstance(scores,(list)):
									out_str += accession + "\t" +  ",".join(["%1.3f"%(x) for x in scores]) + "\n"
								else:
									out_str += accession + "\t" +  scores + "\n"
							except:
								print("ERROR")

						del scores
					except:
						print("ERROR")

				print(out_str)
				open("testing_set." + source + "." + str(window_size) + '.tdt',"w").write(out_str)
				logger.debug("Writing: " + "testing_set." + source + "." + str(window_size) + '.tdt')

		return {}
		sys.exit()

		for window_size in [5,10,15,20,25,30,35,40,45,50]:
			self.options['window_size'] = window_size

			scores = self.calculate_scores()

			for source in ['pLDDT','accessibility',"intramolecular_contacts"]:
				out_str = ""
				for accession in scores:

					if isinstance(scores[accession][source + "_windowed"],(list)):
						out_str += accession + "\t" +  ",".join(["%1.2f"%(x) for x in scores[accession][source + "_windowed"]]) + "\n"
					else:
						out_str += accession + "\t" +  scores[accession] + "\n"

				#print(out_str)

				open("testing_set." + source + "." + str(window_size) + '.tdt',"w").write(out_str)
				logger.debug("Writing: " + "testing_set." + source + "." + str(window_size) + '.tdt')


		for source in ['pLDDT','accessibility',"accessibility_smoothed","intramolecular_contacts","intramolecular_contacts_smoothed"]:
			out_str = ""

			for accession in scores:
				try:
					if isinstance(scores[accession][source],(list)):
						out_str += accession + "\t" +  ",".join(["%1.3f"%(x) for x in scores[accession][source]]) + "\n"
					else:
						out_str += accession + "\t" +  scores[accession][source] + "\n"
				except:
					print("ERROR")

			open("testing_set." + source + '.tdt',"w").write(out_str)
			logger.debug("Writing: " + "testing_set." + source  + '.tdt')


		return {}

	def benchmark_alphafold_motifs(self):

		#motif_instances = queryRunner.queryRunner("database_access","get_instances",{"database_name":"motifs",'columns':'detailed',"additional_params":{"motif_protein_taxon_identifier":"9606"}}).run()
		#motif_instances = queryRunner.queryRunner("database_access","get_instances_by_motif_accession",{"database_name":"motifs",'columns':'detailed',"accession":"P04637"}).run()
		motif_instances = queryRunner.queryRunner("database_access","get_instances",{"database_name":"motifs",'columns':'detailed',"additional_params":{}}).run()
		
		accessions = {
			"motif":[],
			"domain":[]
		}
		
		for motif_instance in motif_instances['data']:
			#if "ELM" in motif_instance['sources']:
				if motif_instance['motif_protein_uniprot_accession'] != None:
					accessions["domain"].append(motif_instance['motif_protein_uniprot_accession'])

				for domain_protein_relation in motif_instance['domain_protein_relation']:
						domain_protein_uniprot_accession = domain_protein_relation['domain_protein_uniprot_accession']

						if domain_protein_uniprot_accession != None:
							domain_start_residue_number = domain_protein_relation['domain_start_residue_number']
							domain_end_residue_number = domain_protein_relation['domain_end_residue_number']
							accessions["motif"].append(domain_protein_uniprot_accession)
				"""
				options = {
					"database_name":"motifs",
					'columns':'detailed',s
					"collapsed":False,
					"include_cols":"other_notes",
					"by_column":"curation_id",
					"by_column_search_term":",".join([str(x) for x in motif_instance['curation_ids']])
				}
				motif_instances_notes = queryRunner.queryRunner("database_access","get_instances",options).run()

				for motif_instances_note in motif_instances_notes['data']:
					domain_protein_uniprot_accession = motif_instances_note['domain_protein_uniprot_accession']
					pprint.pprint(motif_instances_note)
					if domain_protein_uniprot_accession not in accessions and domain_protein_uniprot_accession != None:
						accessions.append(domain_protein_uniprot_accession)
				"""

		self.options['accession'] = accessions["motif"] + accessions["domain"]
		self.options['accession'] = list(set(self.options['accession']))
		self.options['accession'].sort()
		#self.options['accession'].reverse()
		#print("\n".join(self.options['accession']))
		print(len(self.options['accession']))

		data_type_names = {
			'motif':"idr_side",
			'domain':"pocket_side"
		}

		for data_type in ['motif','domain']:
			data_type_name = data_type_names[data_type]
			
			rows = []
			print(len(accessions[data_type]))
			counter = 0
			for accession in list(set(accessions[data_type])):
				counter +=1 
				print(accession,counter,len(list(set(accessions[data_type]))))
				try:
					check_accession = queryRunner.queryRunner("uniprot","check_accession",{"accession":accession}).run()
					if check_accession['status'] == "Error":
						print(accession + " does not exist")
					else:
						sequence_response = queryRunner.queryRunner("uniprot","parse_uniprot_fasta",{"accession":accession}).run()
						#print(sequence_response)
						rows.append(">" + sequence_response['data']['header'] + "\n" + sequence_response['data']['sequence'])
				except:
					utilities_error.printError()

			open(data_type_name + ".list.txt","w").write("\n".join(list(set(accessions[data_type]))))
			open(data_type_name + ".fasta","w").write("\n".join(rows))

		sys.exit()
		logger.debug("Processing " + str(len(self.options['accession'])) + " protein")
		scores = self.calculate_scores()

		#python3 queryManager.py --dataset_type alphafold --verbose True  --task benchmark_alphafold_motifs --debug True

		distribution_score_types = ["pLDDT_windowed'", 'accessibility', 'accessibility_windowed']
		for distribution_score_type in distribution_score_types:

			score_distribution = {
			"motif":{},
			"motif_flanks":{},
			"motif_helices":{},
			"domain":{}
			}

			score_used = {
				"motif":[],
				"domain":[]
			}
			score_distribution_total = {
				"motif":0,
				"domain":0,
				"motif_flanks":0,
				"motif_helices":0,
			}

			for i in range(0,11):
				score_distribution['motif'][i] = 0
				score_distribution['motif_flanks'][i] = 0
				score_distribution['motif_helices'][i] = 0
				score_distribution['domain'][i] = 0

			for motif_instance in motif_instances['data']:
				try:
					options = {
						"database_name":"motifs",
						"collapsed":False,
						"include_cols":"other_notes",
						"by_column":"curation_id",
						"by_column_search_term":",".join([str(x) for x in motif_instance['curation_ids']])
					}

					motif_instances_notes = queryRunner.queryRunner("database_access","get_instances",options).run()

					secondary_structure_type = []
					interaction_type = []
					for motif_instances_note in motif_instances_notes['data']:
						if motif_instances_note['other_notes'] != None:
							if 'secondary structure type' in motif_instances_note['other_notes']:
								secondary_structure_type.append(motif_instances_note['other_notes']['secondary structure type'])

						interaction_type.append(motif_instances_note['interaction_type'])

					#print(motif_instance['motif_protein_uniprot_accession'],motif_instance['motif_protein_uniprot_accession'] in scores,motif_instance['sources'],secondary_structure_type)

					if "ELM" in motif_instance['sources']:
						if motif_instance['motif_protein_uniprot_accession'] in scores:
							#print(scores[motif_instance['motif_protein_uniprot_accession']].keys())
							motif_scores =         scores[motif_instance['motif_protein_uniprot_accession']][distribution_score_type][motif_instance['motif_start_residue_number']:motif_instance['motif_end_residue_number']]
							motif_flank_n_scores = scores[motif_instance['motif_protein_uniprot_accession']][distribution_score_type][motif_instance['motif_start_residue_number']-10:motif_instance['motif_start_residue_number']]
							motif_flank_c_scores = scores[motif_instance['motif_protein_uniprot_accession']][distribution_score_type][motif_instance['motif_end_residue_number']:motif_instance['motif_end_residue_number']+10]

							motif_flank_scores = motif_flank_n_scores  +motif_flank_c_scores

							for i in range(0,len(motif_scores)):
								score_distribution['motif'][int(motif_scores[i]*10)] += 1
								score_distribution_total['motif'] += 1

							for i in range(0,len(motif_flank_scores)):
								score_distribution['motif_flanks'][int(motif_flank_scores[i]*10)] += 1
								score_distribution_total['motif_flanks'] += 1

					if "ELM" in motif_instance['sources']:

							for domain_protein_relation in motif_instance['domain_protein_relation']:
								domain_protein_uniprot_accession = domain_protein_relation['domain_protein_uniprot_accession']

								if motif_instance['motif_protein_uniprot_accession'] in scores:
									if domain_protein_uniprot_accession != None:
										domain_start_residue_number = domain_protein_relation['domain_start_residue_number']
										domain_end_residue_number = domain_protein_relation['domain_end_residue_number']
										domain_scores = []

										if domain_start_residue_number != None and domain_end_residue_number != None and domain_protein_uniprot_accession in scores:
											for i in range(domain_start_residue_number,domain_end_residue_number):
												if domain_protein_uniprot_accession + "_" + str(i) not in score_used:
													score_used['domain'].append(domain_protein_uniprot_accession + "_" + str(i))
													domain_scores.append(scores[domain_protein_uniprot_accession][distribution_score_type][i])

										for i in range(0,len(domain_scores)):
											score_distribution['domain'][int(domain_scores[i]*10)] += 1
											score_distribution_total['domain'] += 1
				except:
					utilities_error.printError()

			#pprint.pprint(score_distribution)
			for source in ['motif','motif_flanks','domain']:
				proportion_sum = 0
				for i in range(0,11):
					proportion = float(score_distribution[source][i])/score_distribution_total[source]
					proportion_sum += proportion
					print(i,"\t","%20s"%distribution_score_type,"\t","%10s"%source,"\t","%10d"%score_distribution[source][i],"\t","%1.3f"%(proportion_sum),"\t","%1.3f"%(proportion),"*"*int(proportion*20))
				print("")

		return {}

	def grab_pdb_file(self,accession):
		url = "https://alphafold.ebi.ac.uk/files/AF-" + accession + "-F1-model_v" + self.options['alphafold_model_version'] + ".pdb"
		out_path = os.path.join(self.options['data_path'],"alphafold", "pdb","AF-" + accession + "-F1-model_v" + self.options['alphafold_model_version'] + ".pdb")
		logger.debug("Downloading: " + url)
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True)
		logger.debug(status)
		return status

	def grab_predicted_aligned_error_file(self,accession):

		pdb_paths = self.get_pdb_paths(accession)
		predicted_aligned_error_paths = []
		for pdb_path in pdb_paths:
			predicted_aligned_error_path = os.path.basename(pdb_path).replace('model_v' + self.options['alphafold_model_version'] + '.pdb','predicted_aligned_error_v' + self.options['alphafold_model_version'] + '.jsonc')
			predicted_aligned_erro_url = os.path.basename(pdb_path).replace('model_v' + self.options['alphafold_model_version'] + '.pdb','predicted_aligned_error_v' + self.options['alphafold_model_version'] + '.json')

			url = "https://alphafold.ebi.ac.uk/files/" + predicted_aligned_erro_url
			out_path = os.path.join(self.options['data_path'],"alphafold","predicted_aligned_error", "raw",predicted_aligned_error_path)

			logger.debug("Grabbing: " + predicted_aligned_erro_url)
			logger.debug("Writing: " + predicted_aligned_error_path)
			predicted_aligned_error_paths.append(out_path)
			logger.debug("Downloading: " + url)
			sessionDownloaderObj = utilities_downloader.sessionDownloader()
			status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True,JSON=True,zipped_JSON=True)
			logger.debug(status)

		return predicted_aligned_error_paths

	def window_scores(self,scores):
		windowed_scores = []
		for i in range(1,len(scores) + 1):
			if self.options['window_size'] == 0:
				windowed_scores.append(scores[str(i)])
			else:
				start_offset = max(0,i - self.options['window_size'])
				end_offset = min(i + self.options['window_size'] + 1,len(scores))

				window_scores = [scores[str(j+1)] for j in range(start_offset,end_offset)]
				windowed_scores.append(sum(window_scores)/len(window_scores))

		return windowed_scores

	def window_scores_list(self,scores):
		windowed_scores = []
		for i in range(0,len(scores)):
			if self.options['window_size'] == 0:
				windowed_scores.append(scores[i])
			else:
				start_offset = max(0,i - self.options['window_size'])
				end_offset = min(i + self.options['window_size'] + 1,len(scores))

				window_scores = [scores[j] for j in range(start_offset,end_offset)]
				windowed_scores.append(sum(window_scores)/len(window_scores))

		return windowed_scores

	def window_scores_asymmetric(self,scores):
		windowed_scores = []
		for i in range(1,len(scores) + 1):
			if self.options['window_size'] == 0:
				windowed_scores.append(scores[str(i)])
			else:
				start_offset = max(0,i - self.options['window_size'])
				end_offset = min(i + self.options['window_size'] + 1,len(scores))

				window_scores_n = [scores[str(j+1)] for j in range(start_offset,i)]
				window_scores_c = [scores[str(j+1)] for j in range(i,end_offset)]
				#print(window_scores_n/len(window_scores_n),window_scores_c/len(window_scores_c))
				#windowed_scores.append(sum(window_scores)/len(window_scores))

		return windowed_scores

	def window_scores_list_asymetric(self,scores):
		windowed_scores = []
		for i in range(0,len(scores)):
			if self.options['window_size'] == 0:
				windowed_scores.append(scores[i])
			else:
				start_offset = max(0,i - self.options['window_size'])
				end_offset = min(i + self.options['window_size'] + 1,len(scores))

				window_scores_n = [scores[j] for j in range(start_offset,i+1)]
				window_scores_c = [scores[j] for j in range(i,end_offset)]

				window_scores_n_mean = sum(window_scores_n)/len(window_scores_n)
				window_scores_c_mean = sum(window_scores_c)/len(window_scores_c)

				window_scores_n_sum_delta = abs(window_scores_n_mean - scores[i])
				window_scores_c_sum_delta = abs(window_scores_c_mean - scores[i])

				if (window_scores_n_sum_delta < window_scores_c_sum_delta):
					score = window_scores_n_mean
					side = "n"
				elif (window_scores_n_sum_delta > window_scores_c_sum_delta):
					score = window_scores_c_mean
					side = "c"
				else:
					score = scores[i]
					side = "-"

				#if self.options['show_smoothing']:
				#	print(i,"\t",side,"%1.3f"%scores[i],"%1.3f"%score,"%20s"%("*"*int(scores[i]*20)),"%20s"%("*"*int(score*20)),"%1.3f"%(window_scores_n_mean-window_scores_c_mean),"%1.3f"%window_scores_n_mean,"%1.3f"%window_scores_c_mean,"%1.3f"%window_scores_n_sum_delta,"%1.3f"%window_scores_c_sum_delta)
				#windowed_scores.append(sum(window_scores)/len(window_scores))

		return windowed_scores

	####-----------------------------------------------------------
	####-----------------------------------------------------------

	def calculate_smoothed_scores(self,scores,smoothing_score_cutoff,cutoff_type="greater_than"):

		self.options["smoothing_score_cutoff"] = smoothing_score_cutoff


		#if self.options['show_smoothing']:
		#	print("seq\t-\t",self.options["sequence"])
		#	print("mask\t-\t",self.output['masked_sequence'])
		self.output = {}
		self.output['masked_binary'] = []
		self.output['scores'] = scores
		for i in range(0,len(scores)):
			if cutoff_type == "greater_than":
				self.output['masked_binary'].append(scores[i] >= self.options["smoothing_score_cutoff"])
			elif cutoff_type == "less_than":
				self.output['masked_binary'].append(scores[i] <= self.options["smoothing_score_cutoff"])

		self.output['masked_binary'][0] = True
		self.output['masked_binary'][-1] = True

		grouped_L = [(k, sum(1 for i in g)) for k,g in groupby(self.output['masked_binary'] )]
		smoothed = False
		replace_structureclass = 0
		while smoothed == False:
			min_score = [self.options["smooth_disorder_length_cutoff"],None,None]

			if grouped_L[0][0] and grouped_L[-1][0] and not self.options["smooth_termini"]:
				iterator = list(range(1,len(grouped_L)-1))
			else:
				iterator = list(range(0,len(grouped_L)))

			if replace_structureclass%2 == 0:iterator.reverse()

			for i in iterator:
				structureclass = grouped_L[i][0]
				length = grouped_L[i][1]

				if length < min_score[0] and structureclass == False and length <= self.options["smooth_disorder_length_cutoff"]:
					min_score = [length,structureclass,i]

				if length < min_score[0] and structureclass == True and length <= self.options["smooth_order_length_cutoff"]:
					min_score = [length,structureclass,i]

				if length == min_score[0] and structureclass == False:
					min_score = [length,structureclass,i]

			if min_score[0] >= self.options["smooth_disorder_length_cutoff"] or min_score[0] == len(self.options["sequence"]):
				break
			

			replace_structureclass = min_score[2]
			range_structureclass =  list(range(max(0,replace_structureclass-1),min(replace_structureclass+2,len(grouped_L))))

			merged_logic = not grouped_L[replace_structureclass][0]
			merged_length = 0

			for i in range_structureclass:
				merged_length += grouped_L[i][1]

			grouped_L[range_structureclass[0]] = (merged_logic,merged_length)

			for i in range_structureclass[1:]:
				del grouped_L[range_structureclass[1]]

			if self.options['show_smoothing']:
				####-----------------------------------------------------------
				####-----------------------------------------------------------
				####-----------------------------------------------------------

				self.output['smoothed_masked_binary'] = []

				for group in grouped_L:
					if group[0] == True:
						self.output['smoothed_masked_binary'] += [True]*group[1]
					if group[0] == False:
						self.output['smoothed_masked_binary'] += [False]*group[1]

				self.output['smoothed_masked_sequence'] = "".join([self.options["sequence"][i] if self.output['smoothed_masked_binary'][i] else "x" for i in range(0,len(self.output['scores']))])

				logger.debug([min_score[1],"\t",min_score[0],"\t",self.output['smoothed_masked_sequence']])

		####-----------------------------------------------------------

		self.output['smoothed_masked_binary'] = []

		for group in grouped_L:
			if group[0] == True:
				self.output['smoothed_masked_binary'] += [True]*group[1]
			if group[0] == False:
				self.output['smoothed_masked_binary'] += [False]*group[1]

		self.output['smoothed_masked_sequence'] = "".join([self.options["sequence"][i] if self.output['smoothed_masked_binary'][i] else "x" for i in range(0,len(self.output['scores']))])

		return [1 if self.output['smoothed_masked_binary'][i] else 0 for i in range(0,len(self.output['scores']))]


	####-----------------------------------------------------------
	####-----------------------------------------------------------

	def calculate_degree(self,accession):
		try:
			matrix = self.read_predicted_aligned_error(accession)

			degrees = []
			if matrix != None:
				offsets = list(matrix.keys())
				offsets.sort()

				g = networkx.Graph()

				g.add_nodes_from(offsets)
				for pos1 in offsets:
					for pos2 in offsets:
						if pos2 in matrix[pos1]:
							if matrix[pos1][pos2] < 5:
								g.add_edge(pos1, pos2)

				for degree in g.degree():
					degrees.append(degree[1])

			return degrees
		except:
			utilities_error.printError()
			return []
	####-----------------------------------------------------------
	####-----------------------------------------------------------

	def calculate_weighted_degree(self,accession):
		try:
			matrix = self.read_predicted_aligned_error(accession)

			weighted_degrees = []
			if matrix != None:
				offsets = list(matrix.keys())
				offsets.sort()

				g = networkx.Graph()

				g.add_nodes_from(offsets)
				for pos1 in offsets:
					for pos2 in offsets:
						if pos2 in matrix[pos1]:
							if abs(pos2 - pos1) < 50:
								weight =  (1/matrix[pos1][pos2])
								g.add_edge(pos1, pos2, weight=weight)

				for degree in g.degree(weight='weight'):
					weighted_degrees.append(degree[1])

			return weighted_degrees
		except:
			utilities_error.printError()
			return []

	####-----------------------------------------------------------
	####-----------------------------------------------------------

	def scale_pLDDT(self,scores):
		scaled_scores = []
		for i in range(0,len(scores)):
			scaled_scores.append(1-(scores[i]/100))
		return scaled_scores

	####-----------------------------------------------------------

	def get_path_accessibility_json(self,accession, tessellation=False,compress=False):
		if tessellation:
			tag = ''
			if self.options['filter_inaccessibility'] == False:
				tag += ".all_residue"
			if self.options['side_chain_contacts_only']:
				tag += ".side_chain_contacts_only"
			if self.options['add_atom_details']:
				tag += ".atom_detail"
			
			format = "json"
			if compress:
				format = "jsonc"

			return os.path.join(self.options['data_path'],"alphafold","accessibility","tessellation", accession + ".alphafold.accessibility" + tag + ".v3." + format)
		
		return os.path.join(self.options['data_path'],"alphafold","accessibility","raw", accession + ".alphafold.accessibility.json")

	def get_path_accessibility_windowed_json(self,accession):
		return os.path.join(self.options['data_path'],"alphafold","accessibility","windowed", accession + ".alphafold.windowed_" + str(self.options['window_size']) + ".accessibility.json")

	def get_path_accessibility_smoothed_json(self,accession):
		tags = [str(self.options['smooth_order_length_cutoff']),str(self.options['smooth_disorder_length_cutoff']), str(self.options["accessibility_score_cutoff"]) ]
		if self.options["smooth_termini"]:
			tags.append("smooth_termini")
			
		return os.path.join(self.options['data_path'],"alphafold","accessibility","smoothed", accession + ".alphafold.smoothed_" + "_".join(tags) + ".accessibility.json")

	def get_path_pLDDT_json(self,accession):
		return os.path.join(self.options['data_path'],"alphafold","pLDDT","raw", accession + ".alphafold.pLDDT.json")

	def get_path_pLDDT_windowed_json(self,accession):
		return os.path.join(self.options['data_path'],"alphafold","pLDDT","windowed", accession + ".alphafold.windowed_" + str(self.options['window_size']) + ".pLDDT.json")

	def get_path_pLDDT_scaled_json(self,accession):
		return os.path.join(self.options['data_path'],"alphafold","pLDDT","scaled", accession + ".alphafold.scaled.pLDDT.json")

	def get_path_secondary_structure_json(self,accession):
		return os.path.join(self.options['data_path'],"alphafold","secondary_structure", accession + ".alphafold.ss.json")

	def get_path_intramolecular_contacts_json(self,accession):
		if self.options['use_pae_filtering_for_intramolecular_contacts_counts']:
			return os.path.join(self.options['data_path'],"alphafold", "intramolecular_contacts","raw",accession + ".alphafold" + ".raw_" + str(self.options['adjacent_residue_filter_distance']) + ".pae_filtered_" + str(self.options['pae_filtering_for_intramolecular_contacts_counts_cut_off'])  +  ".intramolecular_contacts.json")
		else:
			return os.path.join(self.options['data_path'],"alphafold", "intramolecular_contacts","raw",accession + ".alphafold" + ".raw_" + str(self.options['adjacent_residue_filter_distance']) + ".intramolecular_contacts.json")

	def get_path_intramolecular_contacts_smoothed_json(self,accession):
		tags = [str(self.options['smooth_order_length_cutoff']),str(self.options['smooth_disorder_length_cutoff']),str(self.options["intramolecular_contact_score_cutoff"]),str(self.options['adjacent_residue_filter_distance']) ]
		if self.options["smooth_termini"]:
			tags.append("smooth_termini")
		
		if self.options['use_pae_filtering_for_intramolecular_contacts_counts']:
			tags.append(".pae_filtered_" +  str(self.options['pae_filtering_for_intramolecular_contacts_counts_cut_off']) )

		return os.path.join(self.options['data_path'],"alphafold", "intramolecular_contacts","smoothed",accession + ".alphafold.smoothed_" +  "_".join(tags)  + ".intramolecular_contacts.json")
	
	def get_path_intramolecular_contacts_globular_smoothed_json(self,accession):
		tags = [str(self.options['smooth_order_length_cutoff']), str(self.options['smooth_disorder_length_cutoff']), str(self.options["intramolecular_contact_score_cutoff"]), str(self.options['adjacent_residue_filter_distance']) + ".pae_filtered"]
		if self.options["smooth_termini"]:
			tags.append("smooth_termini")
			
		if self.options['use_pae_filtering_for_intramolecular_contacts_counts']:
			tags.append(".pae_filtered_" +  str(self.options['pae_filtering_for_intramolecular_contacts_counts_cut_off'])  )

		return os.path.join(self.options['data_path'],"alphafold", "intramolecular_contacts","smoothed",accession + ".alphafold.smoothed_" +  "_".join(tags)  + ".intramolecular_globular_contacts.json")
			
	def get_path_intramolecular_contacts_windowed_json(self,accession):
		if self.options['use_pae_filtering_for_intramolecular_contacts_counts']:
			return os.path.join(self.options['data_path'],"alphafold", "intramolecular_contacts","windowed",accession + ".alphafold.windowed_" + str(self.options['window_size']) + ".pae_filtered_" +  str(self.options['pae_filtering_for_intramolecular_contacts_counts_cut_off']) + ".intramolecular_contacts.json")
		else:
			return os.path.join(self.options['data_path'],"alphafold", "intramolecular_contacts","windowed",accession + ".alphafold.windowed_" + str(self.options['window_size'])  + ".intramolecular_contacts.json")

	def get_path_predicted_aligned_error_png(self,accession):
		return os.path.join(self.options['data_path'],"alphafold", "predicted_aligned_error","plots",accession + ".predicted_aligned_error.png")

	def get_path_predicted_aligned_error_json(self,accession):
		return os.path.join(self.options['data_path'],"alphafold", "predicted_aligned_error","raw",accession + ".alphafold.predicted_aligned_error.json")

	def get_path_degree_json(self,accession):
		return os.path.join(self.options['data_path'],"alphafold", "predicted_aligned_error","degree",accession + ".alphafold.predicted_aligned_error.degree.json")

	def get_path_weighted_degree_json(self,accession):
		return os.path.join(self.options['data_path'],"alphafold", "predicted_aligned_error","degree",accession + ".alphafold.predicted_aligned_error.weighted_degree.json")

	####-----------------------------------------------------------
	####-----------------------------------------------------------

	def save_accessibility_json(self,accession,scores, tessellation=False):
		out_file_json = self.get_path_accessibility_json(accession, tessellation)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_accessibility_windowed_json(self,accession,scores):
		out_file_json = self.get_path_accessibility_windowed_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_accessibility_smoothed_json(self,accession,scores):
		out_file_json = self.get_path_accessibility_smoothed_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_pLDDT_json(self,accession,scores):
		out_file_json = self.get_path_pLDDT_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_pLDDT_windowed_json(self,accession,scores):
		out_file_json = self.get_path_pLDDT_windowed_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_pLDDT_scaled_json(self,accession,scores):
		out_file_json = self.get_path_pLDDT_scaled_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_secondary_structure_json(self,accession,scores):
		out_file_json = self.get_path_secondary_structure_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_intramolecular_contacts_json(self,accession,scores):
		out_file_json = self.get_path_intramolecular_contacts_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_intramolecular_contacts_smoothed_json(self,accession,scores):
		out_file_json = self.get_path_intramolecular_contacts_smoothed_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_intramolecular_contacts_globular_smoothed_json(self,accession,scores):
		out_file_json = self.get_path_intramolecular_contacts_globular_smoothed_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_intramolecular_contacts_windowed_json(self,accession,scores):
		out_file_json = self.get_path_intramolecular_contacts_windowed_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_degree_json(self,accession,scores):
		out_file_json = self.get_path_degree_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)

	def save_weighted_degree_json(self,accession,scores):
		out_file_json = self.get_path_weighted_degree_json(accession)
		logger.debug("Writing " + out_file_json)
		utilities_basic.write_to_json(out_file_json,scores)
	####-----------------------------------------------------------
	####-----------------------------------------------------------

	def get_accessibility_smoothed_scores(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		accessibility_smoothed_out_file_json = self.get_path_accessibility_smoothed_json(accession)
		if not os.path.exists(accessibility_smoothed_out_file_json) or self.options['remake']:
			accessibility = self.get_accessibility(accession)
			if len(accessibility) == 0:
				return []

			accessibility_smoothed= self.calculate_smoothed_scores(accessibility,self.options["accessibility_score_cutoff"])
			self.save_accessibility_smoothed_json(accession,accessibility_smoothed)
		else:
			if os.path.exists(accessibility_smoothed_out_file_json):
				logger.debug("Reading " + accessibility_smoothed_out_file_json)
				accessibility_smoothed = json.loads(open(accessibility_smoothed_out_file_json, 'r').read())

		return accessibility_smoothed

	def get_accessibility_windowed_scores(self,accession=None):
	
		if accession == None:
			if len(self.options['accession']) == 1:
				accession = self.options['accession'][0]
			else:
				return self.get_accessibility_windowed_scores_list()

		accessibility_windowed_out_file_json = self.get_path_accessibility_windowed_json(accession)
		if not os.path.exists(accessibility_windowed_out_file_json) or self.options['remake']:
			accessibility = self.get_accessibility(accession)

			if 'status' in accessibility:
				return accessibility

			if len(accessibility) == 0:
				return []
			accessibility_windowed = self.window_scores_list(accessibility)
			self.save_accessibility_windowed_json(accession,accessibility_windowed)
		else:
			if os.path.exists(accessibility_windowed_out_file_json):
				logger.debug("Reading " + accessibility_windowed_out_file_json)
				accessibility_windowed = json.loads(open(accessibility_windowed_out_file_json, 'r').read())
		return accessibility_windowed

	def get_accessibility_windowed_scores_list(self):
		data = {}
		for accession in self.options['accession']:
			if accession == None: continue
			data[accession] = self.get_accessibility_windowed_scores(accession)
			#print(accession + "\t" + ",".join(["%1.3f"%x for x in data[accession]]))

		return data
		
	def get_pLDDT_windowed_scores(self,accession=None):
		if accession == None:
			if len(self.options['accession']) == 1:
				accession = self.options['accession'][0]
			else:
				return self.get_pLDDT_windowed_scores_list()

		pLDDT_windowed_file_json = self.get_path_pLDDT_windowed_json(accession)
		if not os.path.exists(pLDDT_windowed_file_json) or self.options['remake']:
			pLDDT = self.get_pLDDT_scaled_scores(accession)

			if len(pLDDT) == 0:
				return []

			pLDDT_windowed = self.window_scores_list(pLDDT)
			self.save_pLDDT_windowed_json(accession,pLDDT_windowed)
		else:
			if os.path.exists(pLDDT_windowed_file_json):
				logger.debug("Reading " + pLDDT_windowed_file_json)
				pLDDT_windowed= json.loads(open(pLDDT_windowed_file_json, 'r').read())
		return pLDDT_windowed

	def get_pLDDT_windowed_scores_list(self):
		data = {}
		for accession in self.options['accession']:
			if accession == None: continue
			data[accession] = self.get_pLDDT_windowed_scores(accession)
			#print(accession + "\t" + ",".join(["%1.3f"%x for x in data[accession]]))

		return data

	def get_pLDDT_scaled_scores(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		pLDDT_scaled_file_json = self.get_path_pLDDT_scaled_json(accession)
		if not os.path.exists(pLDDT_scaled_file_json) or self.options['remake']:
			pLDDT = self.get_pLDDT(accession)
			if len(pLDDT) == 0:
				return []
			pLDDT_scaled = self.scale_pLDDT(pLDDT)
			self.save_pLDDT_scaled_json(accession,pLDDT_scaled)
		else:
			if os.path.exists(pLDDT_scaled_file_json):
				logger.debug("Reading " + pLDDT_scaled_file_json)
				pLDDT_scaled = json.loads(open(pLDDT_scaled_file_json, 'r').read())
		return pLDDT_scaled

	
	def get_secondary_structure_scores(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		secondary_structure_file_json = self.get_path_secondary_structure_json(accession)
		if not os.path.exists(secondary_structure_file_json) or self.options['remake']:
			secondary_structure = self.get_secondary_structure(accession)
			self.save_secondary_structure_json(accession,secondary_structure)
		else:
			if os.path.exists(secondary_structure_file_json):
				logger.debug("Reading " + secondary_structure_file_json)
				secondary_structure = json.loads(open(secondary_structure_file_json, 'r').read())
		return secondary_structure

	def get_intramolecular_contacts_scores(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		intramolecular_contacts_file_json = self.get_path_intramolecular_contacts_json(accession)
		if not os.path.exists(intramolecular_contacts_file_json) or self.options['remake']:
			intramolecular_contacts =self.calculate_intramolecular_contacts(accession)
			self.save_intramolecular_contacts_json(accession,intramolecular_contacts)
		else:
			if os.path.exists(intramolecular_contacts_file_json):
				logger.debug("Reading " + intramolecular_contacts_file_json)
				intramolecular_contacts = json.loads(open(intramolecular_contacts_file_json, 'r').read())
		return intramolecular_contacts

	def get_intramolecular_contacts_windowed_scores(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		intramolecular_contacts_windowed_file_json = self.get_path_intramolecular_contacts_windowed_json(accession)
		if not os.path.exists(intramolecular_contacts_windowed_file_json) or self.options['remake']:
			intramolecular_contacts = self.get_intramolecular_contacts_scores(accession)
			if len(intramolecular_contacts) == 0:
				return []
			intramolecular_contacts_windowed = self.window_scores_list(intramolecular_contacts)
			self.save_intramolecular_contacts_windowed_json(accession,intramolecular_contacts_windowed)
		else:
			if os.path.exists(intramolecular_contacts_windowed_file_json):
				logger.debug("Reading " + intramolecular_contacts_windowed_file_json)
				intramolecular_contacts_windowed =  json.loads(open(intramolecular_contacts_windowed_file_json, 'r').read())
		return intramolecular_contacts_windowed

	def get_intramolecular_contacts_smoothed_scores(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		intramolecular_contacts_smoothed_file_json = self.get_path_intramolecular_contacts_smoothed_json(accession)
		if not os.path.exists(intramolecular_contacts_smoothed_file_json) or self.options['remake']:
			intramolecular_contacts = self.get_intramolecular_contacts_scores(accession)
			if len(intramolecular_contacts) == 0:
				return []
			intramolecular_contacts_smoothed = self.calculate_smoothed_scores(intramolecular_contacts,self.options["intramolecular_contact_score_cutoff"],"less_than")
			self.save_intramolecular_contacts_smoothed_json(accession,intramolecular_contacts_smoothed)
		else:
			if os.path.exists(intramolecular_contacts_smoothed_file_json):
				logger.debug("Reading " + intramolecular_contacts_smoothed_file_json)
				intramolecular_contacts_smoothed = json.loads(open(intramolecular_contacts_smoothed_file_json, 'r').read())

		return intramolecular_contacts_smoothed


	def get_intramolecular_contacts_globular_smoothed_scores(self,accession=None):
		smooth_disorder_length_cutoff = self.options["smooth_disorder_length_cutoff"]
		smooth_order_length_cutoff = self.options["smooth_order_length_cutoff"]
		self.options["smooth_disorder_length_cutoff"] = smooth_order_length_cutoff
		self.options["smooth_order_length_cutoff"] = smooth_disorder_length_cutoff

		if accession == None:
			accession = self.options['accession'][0]

		if 'sequence' not in self.options:
			sequence_response = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":accession}).run()
			self.options['sequence'] = sequence_response['data']

		intramolecular_contacts_smoothed_file_json = self.get_path_intramolecular_contacts_globular_smoothed_json(accession)
		if not os.path.exists(intramolecular_contacts_smoothed_file_json) or self.options['remake']:
			intramolecular_contacts = self.get_intramolecular_contacts_scores(accession)
			
			if len(intramolecular_contacts) == 0:
				return []
			intramolecular_contacts_smoothed = self.calculate_smoothed_scores(intramolecular_contacts,self.options["intramolecular_contact_score_cutoff"],"greater_than")
			self.save_intramolecular_contacts_globular_smoothed_json(accession,intramolecular_contacts_smoothed)
		else:
			if os.path.exists(intramolecular_contacts_smoothed_file_json):
				logger.debug("Reading " + intramolecular_contacts_smoothed_file_json)
				intramolecular_contacts_smoothed = json.loads(open(intramolecular_contacts_smoothed_file_json, 'r').read())

		return intramolecular_contacts_smoothed

	def get_accessibility(self,accession=None, tessellation=False):
		if accession == None:
			accession = self.options['accession'][0]

		accessibility_file_json = self.get_path_accessibility_json(accession, tessellation)
		
		if not os.path.exists(accessibility_file_json) or self.options['remake']:
			
			accessibility = self.calculate_accessibility(accession, tessellation)
			logger.debug("Writing " + accessibility_file_json)
			if 'status' in accessibility:
				return accessibility

			try:
				self.save_accessibility_json(accession,accessibility, tessellation)
			except:
				logger.debug("Error writing file " + accession)
		else:
			if os.path.exists(accessibility_file_json):
				logger.debug("Reading " + accessibility_file_json)
				accessibility = json.loads(open(accessibility_file_json, 'r').read())
		
		return accessibility
		

	def get_pLDDT(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		pLDDT_file_json = self.get_path_pLDDT_json(accession)
		if not os.path.exists(pLDDT_file_json) or self.options['remake']:
			pLDDT =self.calculate_pLDDT(accession)
			self.save_pLDDT_json(accession,pLDDT)
		else:
			if os.path.exists(pLDDT_file_json):
				logger.debug("Reading " + pLDDT_file_json)
				pLDDT = json.loads(open(pLDDT_file_json, 'r').read())

		return pLDDT

	def get_degree(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		degree_file_json = self.get_path_degree_json(accession)
		if not os.path.exists(degree_file_json) or self.options['remake']:
			degree =self.calculate_degree(accession)

			### TODO calculate_degree
			self.save_degree_json(accession,degree)
		else:
			if os.path.exists(degree_file_json):
				logger.debug("Reading " + degree_file_json)
				degree = json.loads(open(degree_file_json, 'r').read())

		return degree

	def get_weighted_degree(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		weighted_degree_file_json = self.get_path_weighted_degree_json(accession)
		if not os.path.exists(weighted_degree_file_json) or self.options['remake']:
			degree =self.calculate_weighted_degree(accession)

			### TODO calculate_degree
			self.save_weighted_degree_json(accession,degree)
		else:
			if os.path.exists(weighted_degree_file_json):
				logger.debug("Reading " + weighted_degree_file_json)
				degree = json.loads(open(weighted_degree_file_json, 'r').read())

		return degree


	####-----------------------------------------------------------
	####-----------------------------------------------------------

	def save_tdt(self,accession,details,scores):
		tdt_str = ""
		row = ["#","Res","SS","SA","SAnorm","SAwin","plot"]
		tdt_str += "\t".join(row) + "\n"
		for i in range(0,len(scores)):
			row = [
				str(i+1),details['sequence'][i],
				details['secondary_structure'][i],
				str(int(details['surface_accessibility'][str(i+1)])),
				"%1.3f"%(details['surface_accessibility_normalised'][str(i+1)]),
				"%1.3f"%scores[i],
				"*"*int(scores[i]*10)
			]
			tdt_str += "\t".join(row) + "\n"

		out_file_tdt = os.path.join(self.options['data_path'],"alphafold","accessibility", "tdt",accession + ".alphafold.accessibility.tdt")
		logger.debug("Writing " + out_file_tdt)
		open(out_file_tdt,"w").write(tdt_str)

	####-----------------------------------------------------------
	####-----------------------------------------------------------

	def calculate_pLDDT(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		try:
			P = PDBParser(QUIET=1)
		except AttributeError:
			P = PDBParser()

		pLDDT_dict = {}
		sequence_response = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":accession}).run()
		sequence = sequence_response['data']

		#

		pdb_paths = self.get_pdb_paths(accession)

		if len(pdb_paths) > 0:
			for pdb_path in pdb_paths:
				construct = ""
				pLDDT = {}
				structure = P.get_structure('structure', pdb_path)
				for atom in structure.get_atoms():
					pLDDT[atom.parent.get_id()[1]] = atom.bfactor

				for res in structure.get_residues():
					resseq = str(res.get_id()[1])

					try:
						one_letter = three_to_one(res.resname)	
					except:
						one_letter = index_to_one(three_to_index(res.resname))
						
					construct +=  one_letter 

				if construct != sequence:
					if sequence.count(construct) == 0 and len(sequence) != len(construct):
						return {
							"status":"Error",
							"error_type":"sequence and model construct do not match",
							"sequence":sequence,
							"construct":construct
						}
					elif sequence.count(construct) == 1:
						protein_offset = len(sequence.split(construct)[0])
					else:
						protein_offset = 0
				else:
					protein_offset = len(sequence.split(construct)[0])

				for i in range(0,len(pLDDT)):
					protein_offset_i = protein_offset + i
					if protein_offset_i not in pLDDT_dict:
						pLDDT_dict[protein_offset_i] = {}

					pLDDT_dict[protein_offset_i][min(i,len(pLDDT)-i)] = pLDDT[i+1]
		else:
			logger.debug(accession + " not found")

		pLDDT_list = []
		
		for i in range(0,len(pLDDT_dict)):
			edge_distances = list(pLDDT_dict[i].keys())
			max_edge_distance = max(edge_distances)
			pLDDT_list.append(pLDDT_dict[i][max_edge_distance])

		return pLDDT_list


	def calculate_accessibility(self, accession=None, tessellation=False):
		if accession == None:
			accession = self.options['accession'][0]

		pdb_paths = self.get_pdb_paths(accession)
		accessibility_dict = {}
		sequence_response = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":accession}).run()
		sequence = sequence_response['data']

		if len(pdb_paths) > 0:
			for pdb_path in pdb_paths:
				#pdb_path = "/home/data/alphafold/pdb/AF-O60673-F10-model_v4.pdb"
				pdb_name = os.path.basename(pdb_path).split(".")[-1]
				
				if tessellation:
					simpleMappingAccessibilityUtilitiesObj = simpleMappingAccessibilityUtilities()

					if self.options['side_chain_contacts_only']:
						tessellation_data = simpleMappingAccessibilityUtilitiesObj.get_tessellation_accessibility(pdb_name, pdb_path, alphafold=True,filter_inaccessibility=self.options['filter_inaccessibility'],skip_atoms = ['CA', 'C', 'N', 'O'])
					else:
						tessellation_data = simpleMappingAccessibilityUtilitiesObj.get_tessellation_accessibility(pdb_name, pdb_path, alphafold=True,filter_inaccessibility=self.options['filter_inaccessibility'],add_atom_details=self.options['add_atom_details'])
					
					if not accessibility_dict:
						accessibility_dict.update(tessellation_data)
						accessibility_dict['chain_details']['A']['data_source'] = {}
						for offset in accessibility_dict['chain_details']['A']['accessible_residues']:
							accessibility_dict['chain_details']['A']['data_source'][offset] = tessellation_data['chain_details']['A']['peptide_offset_shift_centre']
					else:
						try:

							offsets = list(tessellation_data['chain_details']['A']['accessible_residues'].keys())
							offsets.sort()
	
							for offset in offsets:
								if tessellation_data['chain_details']['A']['accessible_residues'][offset]['aa'] != sequence[int(offset)-1]:
									logging.error([
										"amino acid mismatch",
										offset,
										tessellation_data['chain_details']['A']['accessible_residues'][offset]['aa'],
										sequence[int(offset)-1]
									])

								if offset not in accessibility_dict['chain_details']['A']['data_source']:
									accessibility_dict['chain_details']['A']['data_source'][offset] = tessellation_data['chain_details']['A']['peptide_offset_shift_centre']
									for k in accessibility_dict['chain_details']['A']:
										if k in ['data_source','peptide_offset_shift_centre','peptide_offset_shift_start','peptide_offset_shift_end']: continue
										if offset in tessellation_data['chain_details']['A'][k]: 
											accessibility_dict['chain_details']['A'][k][offset] = tessellation_data['chain_details']['A'][k][offset]
										else:
											if tessellation_data['chain_details']['A']['accessible_residues'][offset]['any_atm_score'] > 0:
												logging.error(offset + " not in " + k + " " + str(len(tessellation_data['chain_details']['A'][k])))
								else:
									if abs(int(offset)-tessellation_data['chain_details']['A']['peptide_offset_shift_centre']) < abs(int(offset) - accessibility_dict['chain_details']['A']['data_source'][offset]):
										accessibility_dict['chain_details']['A']['data_source'][offset] = tessellation_data['chain_details']['A']['peptide_offset_shift_centre']
										for k in accessibility_dict['chain_details']['A']:
											if k in ['data_source','peptide_offset_shift_centre','peptide_offset_shift_start','peptide_offset_shift_end']: continue
											if offset in tessellation_data['chain_details']['A'][k]:
												accessibility_dict['chain_details']['A'][k][offset] = tessellation_data['chain_details']['A'][k][offset]
											else:
												if tessellation_data['chain_details']['A']['accessible_residues'][offset]['any_atm_score'] > 0:
													logging.error(offset + " not in " + k + " " + str(len(tessellation_data['chain_details']['A'][k])))
						except:
							utilities_error.printError()
							raise
				else:
					accessibilityUtilitiesOj = dsspUtilities()
					
					accessibility = []
					dssp_data = accessibilityUtilitiesOj.get_dssp_data(pdb_name,pdb_path,alphafold=True)

					if  len(dssp_data['pdb.A']['summary']['construct']) == len(sequence):
						protein_offset = 0
					else:					
						protein_offset = len(sequence.split(dssp_data['pdb.A']['summary']['construct'])[0])

					if sequence.count(dssp_data['pdb.A']['summary']['construct']) == 0 and len(dssp_data['pdb.A']['summary']['construct']) != len(sequence):
						return {"status":"Error","error_type":"UniProt sequence - model sequence mismatch "  + accession + " - " + sequence + " - "  + dssp_data['pdb.A']['summary']['construct'] + " - " + str(len(sequence)) + " - "  + str(len(dssp_data['pdb.A']['summary']['construct']))}

					for i in range(1,len(dssp_data['pdb.A']['summary']['surface_accessibility_normalised']) + 1):
						accessibility.append(dssp_data['pdb.A']['summary']['surface_accessibility_normalised'][str(i)])
					
					for i in range(0,len(accessibility)):
						protein_offset_i = protein_offset + i
						if protein_offset_i not in accessibility_dict:
							accessibility_dict[protein_offset_i] = {}

						accessibility_dict[protein_offset_i][min(i,len(accessibility)-i)] = accessibility[i]

		if tessellation:
			return accessibility_dict
		
		accessibility_list = []
		for i in range(0,len(accessibility_dict)):
			try:
				edge_distances = list(accessibility_dict[i].keys())
				max_edge_distance = max(edge_distances)
				accessibility_list.append(accessibility_dict[i][max_edge_distance])
			except:
				accessibility_list.append(1)
				logger.error(accession)
				utilities_error.printError()
				
		return accessibility_list

	def get_secondary_structure(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		pdb_paths = self.get_pdb_paths(accession)
		secondary_structure_dict = {}
		sequence_response = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":accession}).run()
		sequence = sequence_response['data']

		if len(pdb_paths) > 0:
			for pdb_path in pdb_paths:
				accessibilityUtilitiesOj = dsspUtilities()

				accessibility = []
				pdb_name = os.path.basename(pdb_path).split(".")[-1]
				dssp_data = accessibilityUtilitiesOj.get_dssp_data(pdb_name,pdb_path,alphafold=True)

				construct = dssp_data['pdb.A']['summary']['construct']
				protein_offset = len(construct)

				if construct != sequence:
					if sequence.count(construct) == 0 and len(sequence) != len(construct):
						return {
							"status":"Error",
							"error_type":"sequence and model construct do not match",
							"sequence":sequence,
							"construct":construct
						}
					elif sequence.count(construct) == 1:
						protein_offset = len(sequence.split(construct)[0])
					else:
						protein_offset = 0
				else:
					protein_offset = len(sequence.split(construct)[0])

				secondary_structure = list(dssp_data['pdb.A']['summary']['secondary_structure'])
				#print(pdb_path,protein_offset,sequence.count(construct),len(sequence),len(construct),len(secondary_structure))
				#/home/data/alphafold/pdb/AF-P49454-F4-model_v4.pdb 0 1 3114 600
				for i in range(0,len(secondary_structure)):
					protein_offset_i = protein_offset + i
					if protein_offset_i not in secondary_structure_dict:
						secondary_structure_dict[protein_offset_i] = {}

					secondary_structure_dict[protein_offset_i][min(i,len(secondary_structure)-i)] = secondary_structure[i]
				
		secondary_structure_list = []
		for i in range(0,len(secondary_structure_dict)):
			edge_distances = list(secondary_structure_dict[i].keys())
			max_edge_distance = max(edge_distances)
			secondary_structure_list.append(secondary_structure_dict[i][max_edge_distance])

		return secondary_structure_list


	def get_secondary_structure_collapsed(self,accession=None):
		secondary_structure_list = self.get_secondary_structure()
		secondary_structure_composition = {"residue":{},"region":[]}

		tracker = ["-",0,0]

		for i in range(0,len(secondary_structure_list)):

			if tracker[0] == secondary_structure_list[i]:
				tracker[1] += 1
			else:
				secondary_structure_composition["region"].append(tracker)
				for ii in range(tracker[2],tracker[2] + tracker[1]):
					secondary_structure_composition["residue"][ii] = {tracker[0]:tracker[1]}

				tracker = [secondary_structure_list[i],1,i]

		for ii in range(tracker[2],len(secondary_structure_list)):
			secondary_structure_composition["residue"][ii] = {tracker[0]:tracker[1]}		

		return secondary_structure_composition

	def calculate_intramolecular_contacts(self,accession=None):
		if accession == None:
			accession = self.options['accession'][0]

		pdb_paths = self.get_pdb_paths(accession)
		sequence_response = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":accession}).run()
		sequence = sequence_response['data']

		count = 0

		structureAtomicCoordinateUtilitiesObj = structureAtomicCoordinateUtilities()
		intramolecular_interaction_dict = {}

		if len(pdb_paths) > 1:
			return []

		for pdb_path in pdb_paths:
			structureAtomicCoordinateUtilitiesObj.options['pdb_id'] = "alphafold"
			structureAtomicCoordinateUtilitiesObj.options['pdb_file'] = pdb_path
			structureAtomicCoordinateUtilitiesObj.options['query_chain'] = 'A'
			structureAtomicCoordinateUtilitiesObj.options['output_format'] = "detailed_distances"
			structureAtomicCoordinateUtilitiesObj.options['adjacent_residue_filter_distance'] = self.options['adjacent_residue_filter_distance']

			calculate_intramolecular_interaction_statistics_data = structureAtomicCoordinateUtilitiesObj.calculate_intramolecular_contacts()

			if self.options['use_pae_filtering_for_intramolecular_contacts_counts']:
				predicted_aligned_error = self.read_predicted_aligned_error(accession)
			
			intramolecular_contacts_list = {}
			construct = ""


			offsets = list(calculate_intramolecular_interaction_statistics_data['A']['A'])
			offsets.sort()

			for offset in offsets:
				contact_offset = list(calculate_intramolecular_interaction_statistics_data['A']['A'][offset].keys())

				if len(contact_offset) > 0:
					construct += calculate_intramolecular_interaction_statistics_data['A']['A'][offset][contact_offset[0]]['query_res_name']
				else:
					construct += "."

			protein_offset = len(sequence.split(construct)[0])

			for i in range(1,len(calculate_intramolecular_interaction_statistics_data['A']['A'].keys()) + 1):
				protein_offset_i = protein_offset + i
				if protein_offset_i not in intramolecular_interaction_dict:
					intramolecular_interaction_dict[protein_offset_i] = {}

				intramolecular_contacts_counts = 0
				if self.options['use_pae_filtering_for_intramolecular_contacts_counts']:
					try:
						for ii in calculate_intramolecular_interaction_statistics_data['A']['A'][i].keys():
							if predicted_aligned_error[i-1][ii-1] < self.options['pae_filtering_for_intramolecular_contacts_counts_cut_off']:
								intramolecular_contacts_counts += 1
							else:
								logger.debug("Ignoring low confidence contact: " + str(i)+"-"+str(ii) + " PAE:" + str(predicted_aligned_error[i-1][ii-1]))
					except:
						logger.debug("PAE calculation failed")
						intramolecular_contacts_counts = len(calculate_intramolecular_interaction_statistics_data['A']['A'][i])
				else:
					intramolecular_contacts_counts = len(calculate_intramolecular_interaction_statistics_data['A']['A'][i])

				intramolecular_interaction_dict[protein_offset_i][min(i,len(calculate_intramolecular_interaction_statistics_data['A']['A'].keys())-i)] = intramolecular_contacts_counts 
			"""
			for i in calculate_intramolecular_interaction_statistics_data['A']['A'].keys():
				intramolecular_contacts_list[protein_offset_i].append()
			"""

		intramolecular_contacts_list = []
		offsets = list(intramolecular_interaction_dict.keys())
		offsets.sort()
		for i in offsets:
			edge_distances = list(intramolecular_interaction_dict[i].keys())
			max_edge_distance = max(edge_distances)
			intramolecular_contacts_list.append(intramolecular_interaction_dict[i][max_edge_distance])

		return intramolecular_contacts_list

	####----------––----------####

	def get_pdb_paths(self,accession=None):
		pdb_paths = []

		dirpath =  os.path.join(self.options['data_path'],"alphafold", "pdb")
		if not os.path.exists(dirpath):
			os.mkdir(dirpath)

		error_pdb_path = os.path.join(self.options['data_path'],"alphafold", "pdb", "AF-" + accession + "-F1-model_v" + self.options['alphafold_model_version'] + ".pdb.error")
		gz_glob_path = os.path.join(self.options['data_path'],"alphafold","pdb","AF-" + accession + "-*-model_v" + self.options['alphafold_model_version'] + ".pdb.gz") #-F1-model_v1
		
		for file in glob.glob(gz_glob_path):
			if not os.path.exists(file[:-3]):
				utilities_basic.gunzip_shutil(file,file[:-3])
			pdb_path = file[:-3]
			pdb_paths.append(pdb_path)

		if len(pdb_paths) == 0:
			if not os.path.exists(error_pdb_path):
				status = self.grab_pdb_file(accession)

				if 'status' in status:
					if status['status'] == "Success":
						pdb_paths = [os.path.join(self.options['data_path'],"alphafold", "pdb","AF-" + accession + "-F1-model_v" + self.options['alphafold_model_version'] + ".pdb")]
					else:
						analysis_status = status['status_code']
						open(error_pdb_path,"w").write(str(status))
				else:
					logger.debug(status)
					open(error_pdb_path,"w").write(str(status))
			else:
				logger.error(error_pdb_path)
				analysis_status = "404"
				logger.debug(analysis_status)
				open(error_pdb_path,"w").write(str(analysis_status))

		return pdb_paths



	def calculate_scores(self):

		scores = {}
		counter = 0

		error_log = ""
		for accession in self.options['accession']:
			analysis_status = "Success"

			try:
				counter += 1
				logger.debug("Processing " + str(counter) + " of " + str(len(self.options['accession'])) + " -" + accession)
				error_pdb_path = os.path.join(self.options['data_path'],"alphafold", "pdb", "AF-" + accession + "-F1-model_v" + self.options['alphafold_model_version'] + ".pdb.error")

				error_status_check = True
				if len(self.get_pdb_paths(accession)) == 0:
					error_status_check = False
					analysis_status = "File Missing"
					logger.debug("\t".join(["File Missing",str(len(self.get_pdb_paths(accession))),str(counter),accession,"-","-","-"]))

				elif os.path.exists(error_pdb_path):
					error_status_check = False
					analysis_status = "Status Error"
					logger.debug("\t".join(["Status Error",str(len(self.get_pdb_paths(accession))),str(counter),accession,"-","-","-"]))

				if error_status_check:
					sequence_response = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":accession}).run()
					if 'data' not in sequence_response: continue
					sequence = sequence_response['data']
					self.options['sequence'] = sequence

					scores[accession] = {
						"pLDDT":self.get_pLDDT(accession),
						"pLDDT_windowed":self.get_pLDDT_windowed_scores(accession),
						"pLDDT_scaled": self.get_pLDDT_scaled_scores(accession),
						"accessibility":self.get_accessibility(accession),
						"accessibility_windowed":self.get_accessibility_windowed_scores(accession)	,
						"accessibility_smoothed":self.get_accessibility_smoothed_scores(accession),
						"secondary_structure":self.get_secondary_structure_scores(accession),
						"degree":self.get_degree(accession),
						"weighted_degree":self.get_weighted_degree(accession),
						"intramolecular_contacts":self.get_intramolecular_contacts_scores(accession),
						#"intramolecular_contacts_smoothed":self.get_intramolecular_contacts_smoothed_scores(accession),
						#"intramolecular_contacts_windowed": self.get_intramolecular_contacts_windowed_scores(accession)
					}

					score_lengths = []
					for score_type in scores[accession]:
						score_lengths.append(len(scores[accession][score_type]))
						if len(sequence) != len(scores[accession][score_type]):
							logger.error([accession,score_type,len(sequence), len(scores[accession][score_type])])

					logger.debug(score_lengths)
			except:
				logger.debug(accession)
				logger.debug(utilities_error.printError())

			logger.debug("\t".join([accession,"\t",str(counter),"\t","\t",analysis_status]))

		return scores


	#####------------------------------------------------------------------#####

	def calculate_scores_disprot(self):
		disprot_accession_data = queryRunner.queryRunner("pdb","get_disprot_index_accessions",{}).run()
		disprot_data = queryRunner.queryRunner("pdb","get_disprot_by_accession",{'accession':disprot_accession_data['data']}).run()
		self.options['accession'] = disprot_accession_data['data']
		scores = self.calculate_scores()

	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####


	def get_alphafold_structural_classification(self):
		tag = ".".join([str(t) for t in [self.options['alphafold_short_loop_length_cutoff'],self.options['alphafold_coil_cutoff'],self.options['accessibility_score_cutoff']]])

		classification_out_file_json = os.path.join(self.options['classification_data_path'], self.options['accession'][0] + "." + tag + ".classification.json")
		
		if not os.path.exists(classification_out_file_json) or self.options['remake']:
			alphafold_domain_response = queryRunner.queryRunner("alphafold","get_domains_from_pae_matrix",{"accession":self.options['accession'],"split_accessibility":True,"remake":False}).run()
			secondary_structure_response = queryRunner.queryRunner("alphafold","get_secondary_structure_collapsed",{"accession":self.options['accession'],"remake":False}).run()
			sequence_response = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":self.options['accession']}).run()
			accessibility_response = queryRunner.queryRunner("alphafold","get_accessibility",{"accession":self.options['accession'],"remake":False}).run()

			structural_module_class_data = {}
			for i in range(0,len(sequence_response['data'])):

				try:
					accessibility = accessibility_response['data'][i]
					accessible = accessibility_response['data'][i] > self.options['accessibility_score_cutoff']
				except:
					accessibility = 1
					accessible = True
					logger.error("No accessibility for " + self.options['accession'][0] + ":" + str(i))

				structural_module_class_data[i] = {
					"coil":False,
					"idd":False,
					"domain":False,
					"loop":False,
					"idr_ss":False,
					"unstructured":False,
					"accessible":accessible,
					"island":False,
					"ss":"",
					"ss_context":"",
					"loop_length":0,
					"accessibility":accessibility,
					"module_id":""
				}

			try:
				for i in range(0,len(secondary_structure_response['data']['region'])):
					accessibility_list = []
					for o in range(secondary_structure_response['data']['region'][i][2],secondary_structure_response['data']['region'][i][2] + secondary_structure_response['data']['region'][i][1]):
						accessibility_list.append(accessibility_response['data'][o])
						
					#print(secondary_structure_response['data']['region'][i],self.options['accessibility_score_cutoff'], sum(accessibility_list)/len(accessibility_list))
					
					if sum(accessibility_list)/len(accessibility_list) >= self.options['accessibility_score_cutoff']:
						for o in range(secondary_structure_response['data']['region'][i][2],secondary_structure_response['data']['region'][i][2] + secondary_structure_response['data']['region'][i][1]):
							ss_aa = list(secondary_structure_response['data']['residue'][o].keys())
						
							if "-" in ss_aa:
								structural_module_class_data[o]['unstructured'] = True
							elif 'S' in ss_aa or 'T' in ss_aa or 'B' in ss_aa:
								structural_module_class_data[o]['idr_ss'] = True
							elif 'H' in ss_aa or 'G' in ss_aa or 'E' in ss_aa:
								structural_module_class_data[o]['idr_ss'] = True
			
				for i in range(0,len(sequence_response['data'])):
					if i in secondary_structure_response['data']['residue']:
						if "H" in list(secondary_structure_response['data']['residue'][i].keys()):
							if secondary_structure_response['data']['residue'][i]['H'] >= self.options['alphafold_coil_cutoff']:
								structural_module_class_data[i]['coil'] = True
								structural_module_class_data[i]['idr_ss'] = False
								structural_module_class_data[o]['unstructured'] = False
					
				#logger.error([self.options['accession'],i,len(sequence_response['data']),len(secondary_structure_response['data']['residue'])])
			except:
				logger.error("No SS for " + self.options['accession'][0])
	
			#sys.exit()

			if 'data' in alphafold_domain_response:
				for structural_module in alphafold_domain_response['data']:
					try:
						
						### Coil
						max_coil_length = 0
						for i in range(alphafold_domain_response['data'][structural_module]['start'],alphafold_domain_response['data'][structural_module]['end']):
							try:
								if "H" in list(secondary_structure_response['data']['residue'][i].keys()):
									if secondary_structure_response['data']['residue'][i]['H'] > max_coil_length:
										max_coil_length = secondary_structure_response['data']['residue'][i]['H']
							except:
								logger.error("No SS for " + self.options['accession'][0])
	
			#
						coil_proportion = max_coil_length/(alphafold_domain_response['data'][structural_module]['end'] - alphafold_domain_response['data'][structural_module]['start'])
						
						if coil_proportion > 0.9:
							for i in range(alphafold_domain_response['data'][structural_module]['start'],alphafold_domain_response['data'][structural_module]['end']):
								structural_module_class_data[i]['structural_module_mean_accessibility'] = alphafold_domain_response['data'][structural_module]['mean_accessibility']
								structural_module_class_data[i]['module_id'] = structural_module

							continue
						
						### IDD
						idd_classification = [] 
						for i in range(alphafold_domain_response['data'][structural_module]['start'],alphafold_domain_response['data'][structural_module]['end']):
							idd_classification.append(structural_module_class_data[i]['unstructured'] or structural_module_class_data[i]['idr_ss'] or structural_module_class_data[i]['accessible'])
						
						if idd_classification.count(True)/len(idd_classification) > 0.75:
							for i in range(alphafold_domain_response['data'][structural_module]['start'],alphafold_domain_response['data'][structural_module]['end']):
								structural_module_class_data[i]['structural_module_mean_accessibility'] = alphafold_domain_response['data'][structural_module]['mean_accessibility']
								structural_module_class_data[i]['module_id'] = structural_module
								structural_module_class_data[i]['idd'] = False

							continue

						### Domain
						loop_residues = []
						last_loop_length = 0
						structured_residues = []
						first_region = True

						for i in range(alphafold_domain_response['data'][structural_module]['start'],alphafold_domain_response['data'][structural_module]['end']):

							structural_module_class_data[i]['unstructured'] = False
							structural_module_class_data[i]['idr_ss'] = False

							structural_module_class_data[i]['structural_module_mean_accessibility'] = alphafold_domain_response['data'][structural_module]['mean_accessibility']
							structural_module_class_data[i]['module_id'] = structural_module 

							if structural_module_class_data[i]['structural_module_mean_accessibility'] < self.options['accessibility_score_cutoff']:
								structural_module_class_data[i]['coil'] = False
							
							if i not in structural_module_class_data:
								structural_module_class_data[i] = []

							if i in alphafold_domain_response['data'][structural_module]['residues']:
								loop_length =  len(loop_residues)
								if loop_length > 0:
									for loop_i in loop_residues:
										structural_module_class_data[loop_i]['loop_length'] = len(loop_residues)
									
								if loop_length > 3:
									if first_region == True:
										first_region = False
										if loop_length > 30 and len(structured_residues) < 30:
											for structured_residue_i in structured_residues:
												structural_module_class_data[structured_residue_i]['island'] = True
									
									if last_loop_length > 30 and loop_length > 30 and len(structured_residues) < 30:
										for structured_residue_i in structured_residues:
											structural_module_class_data[structured_residue_i]['island'] = True
										
									last_loop_length =  len(loop_residues)
									loop_residues = []
									structured_residues = []
								
								structural_module_class_data[i]['domain'] = True
								structured_residues.append(i)
							else:
								structural_module_class_data[i]['loop'] = True
								loop_residues.append(i)

						if last_loop_length > 30 and len(structured_residues) < 30:
							for structured_residue_i in structured_residues:
								structural_module_class_data[structured_residue_i]['island'] = True
					except:
						utilities_error.printError()
						raise
			
			structural_module_class_classification = []
			for i in range(0,len(sequence_response['data'])):
				if structural_module_class_data[i]["unstructured"] == True and structural_module_class_data[i]["loop"] == False:
					structural_module_class_classification.append("-")
				elif structural_module_class_data[i]["idr_ss"] == True and structural_module_class_data[i]["loop"] == False:
					structural_module_class_classification.append("s")
				elif structural_module_class_data[i]["coil"] == False and structural_module_class_data[i]["domain"] == False and structural_module_class_data[i]["loop"] == False:
					structural_module_class_classification.append("-")
				elif structural_module_class_data[i]["coil"] == False and structural_module_class_data[i]["domain"] == False and structural_module_class_data[i]["loop"] == True:
					if structural_module_class_data[i]['loop_length'] < self.options['alphafold_short_loop_length_cutoff']:
						structural_module_class_classification.append("l")
					else:
						structural_module_class_classification.append("L")
				elif structural_module_class_data[i]["coil"] == False and structural_module_class_data[i]["domain"] == True:
					if structural_module_class_data[i]['island'] == True:
						structural_module_class_classification.append("d")
					else:
						structural_module_class_classification.append("D")
				elif structural_module_class_data[i]["coil"] == True and structural_module_class_data[i]["domain"] == False:
					structural_module_class_classification.append("C")
				elif structural_module_class_data[i]["coil"] == True:
					structural_module_class_classification.append("c")
				else:
					structural_module_class_classification.append("u")
				
			classification = {
				"structural_module_class_data":structural_module_class_data,
				"structural_module_class_classification":structural_module_class_classification,
				"sequence":sequence_response['data']
			}
		

			logger.debug("Writing " + classification_out_file_json)
			utilities_basic.write_to_json(classification_out_file_json,classification)
		else:
			logger.debug("Reading " + classification_out_file_json)
			classification = json.loads(open(classification_out_file_json, 'r').read())
		
		#print("".join(structural_module_class_classification))
		#sys.exit()
		return classification

	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####

	def set_node_clusters(self,g, clusters):
		for c, v_c in enumerate(clusters):
			for v in v_c:
				# Add 1 to save 0 for external edges
				g.nodes[v]['clusters'] = c + 1

		return g

	def set_edge_clusters(self,g):
		for v, w, in g.edges:
			if g.nodes[v]['clusters'] == g.nodes[w]['clusters']:
				g.edges[v, w]['clusters'] = g.nodes[v]['clusters']
			else:
				g.edges[v, w]['clusters'] = 0

		return g


	def get_color(self,i, r_off=1, g_off=1, b_off=1):
		'''Assign a color to a vertex.'''
		r0, g0, b0 = 0, 0, 0
		n = 16
		low, high = 0.1, 0.9
		span = high - low
		r = low + span * (((i + r_off) * 3) % n) / (n - 1)
		g = low + span * (((i + g_off) * 5) % n) / (n - 1)
		b = low + span * (((i + b_off) * 7) % n) / (n - 1)
		return (r, g, b)

	def read_predicted_aligned_error(self,accession):
		logger.debug("Reading PAE for" + accession)
		predicted_aligned_error_paths = self.grab_predicted_aligned_error_file(accession)

		matrix = {}

		for predicted_aligned_error_path in predicted_aligned_error_paths:
			if os.path.exists(predicted_aligned_error_path):
				for chunk in utilities_basic.read_from_json(predicted_aligned_error_path,zipped=True): #json.loads(open(predicted_aligned_error_path, 'r').read()):
					
					residue1 = chunk['residue1']
					residue2 = chunk['residue2']
					distance = chunk['distance']

					for pos in range(0,len(distance)):
						pos1 = residue1[pos] - 1
						pos2 = residue2[pos] - 1

						if pos1 not in matrix:
							matrix[pos1] = {}

						if pos1 == pos2: continue

						matrix[pos1][pos2] = distance[pos]

		return matrix

	def cluster_predicted_aligned_error(self,accession):
		import networkx
		import pylab as plt
		from networkx.algorithms.community import greedy_modularity_communities

		weight_modifier = 0.5
		matrix = self.read_predicted_aligned_error(accession)

		offsets = list(matrix.keys())
		offsets.sort()

		plt.figure(figsize= (20,20))

		g = networkx.Graph()

		g.add_nodes_from(offsets)
		for pos1 in offsets:
			for pos2 in offsets:
				if pos2 in matrix[pos1]:
					weight =  ((1/matrix[pos1][pos2])**weight_modifier)
					g.add_edge(pos1, pos2, weight=weight)

		clusters = sorted(greedy_modularity_communities(g, weight='weight',resolution=1), key=len, reverse=True)

		g = self.set_node_clusters(g, clusters)
		g = self.set_edge_clusters(g)

		node_color = [self.get_color(g.nodes[v]['clusters']) for v in g.nodes]

		labels = {}
		for offset in offsets:
			labels[offset] = str(offset)

		weights = [g[u][v]['weight'] for u,v in g.edges()]
		pos = networkx.spring_layout(g, k=0.05)

		for degree in g.degree(weight='weight'):
			print(degree[0]+1,"\t","%1.3f"%degree[1],"\t",g.nodes[degree[0]]['clusters'])

		networkx.draw_networkx_edges(g, pos, width=weights, alpha=0.25)
		networkx.draw_networkx_nodes(g, pos, offsets, node_color=node_color)
		networkx.draw_networkx_labels(g, pos, labels, font_size=10, font_color="black")

		out_file_png = self.get_path_predicted_aligned_error_png(accession)
		logger.debug("Writing " + out_file_png)
		plt.savefig(out_file_png)

	#Graph-based community clustering approach to extract protein domains
	#from https://pythonawesome.com/graph-based-community-clustering-approach-to-extract-protein-domains/
	#https://github.com/tristanic/pae_to_domains
	def parse_pae_file(self, accession):

		pae_json_file = self.grab_predicted_aligned_error_file(accession)
		
		data = utilities_basic.read_from_json(pae_json_file[0],zipped=True)
	
		if 'residue1' in data[0]:
			r1, d = data[0]['residue1'],data[0]['distance']
			size = max(r1)
			matrix = numpy.empty((size,size))
			matrix.ravel()[:] = d
		else:
			d = data[0]['predicted_aligned_error']
			size = len(d[0])
			matrix = numpy.array(d)
			

		return matrix

	def domains_from_pae_matrix_networkx(self, pae_power=1, pae_cutoff=5):
		
		accession = self.options['accession'][0]
		structure_modules_out_file_json = os.path.join(self.options['data_path'],"alphafold","structure_modules", accession + ".structure_modules.json")
		
		if not os.path.exists(structure_modules_out_file_json) or self.options['remake']:
			pae_matrix = self.parse_pae_file(accession)

			pLDDT = self.get_pLDDT(accession)
			if self.options['split_accessibility']:
				accessibility = self.get_accessibility(accession)
			
			weights = 1/pae_matrix**pae_power

			g = networkx.Graph()

			edges = numpy.argwhere(pae_matrix < pae_cutoff)
			sel_weights = weights[edges.T[0], edges.T[1]]

			wedges = []
			for (i,j),w in zip(edges,sel_weights):
				if pLDDT[i] >= 70 and pLDDT[j] >= 70:
					wedges.append((i+1,j+1,w))
					
			g.add_weighted_edges_from(wedges)

			split_resolution = self.options['split_resolution']
			clusters = networkx.algorithms.community.greedy_modularity_communities(g, weight='weight', resolution=split_resolution)

			clusters_dict = {}
			for c in clusters:
				if len(c) >= 30:
					cluster_dict = {}
					c_min = int(min(c))
					c_max = int(max(c))
					cluster_dict['start'] = c_min
					cluster_dict['end'] = c_max
					cluster_dict['residues'] = list(c)
					#cluster_dict['residues'] = list(map(int, c))
					cluster_dict['num_residues'] = len(c)
					if self.options['parse_pfam']:
						try:
							pfam_response = queryRunner.queryRunner("uniprot","parse_uniprot_pfam",{"region_start":c_min, "region_end":c_max, "accession":accession}).run()
							if 'data' in pfam_response:
								cluster_dict['pfam'] = pfam_response['data']
						except:
							logger.warning('Failed to parse_uniprot_pfam in domains_from_pae_matrix_networkx: %s', accession)

					total_plddt = 0
					total_pae = 0
					num_pae = 0
					for residue1 in c:
						total_plddt+= pLDDT[residue1-1]
						
						for residue2 in c:
							total_pae+=pae_matrix[residue1-1, residue2-1]
							total_pae+=pae_matrix[residue2-1, residue1-1]
							num_pae+=2
							
					cluster_dict['mean_pLDDT'] = total_plddt / len(c)
					cluster_dict['mean_pae'] = total_pae / num_pae
					
					if self.options['split_accessibility']:
						try:
							total_accessibility = 0
							for residue in c:
								total_accessibility += accessibility[residue-1]
							cluster_dict['mean_accessibility'] = total_accessibility / len(c)
						except:
							logger.warning('Failed to get accessibility in domains_from_pae_matrix_networkx: %s', accession)
					
					clusters_dict[str(c_min)+'_'+str(c_max)] = cluster_dict

			status = utilities_basic.write_to_json(structure_modules_out_file_json,clusters_dict,normalise_json=True)
			logger.debug("Writing " + structure_modules_out_file_json + " " + str(status))
		else:
			clusters_dict = json.loads(open(structure_modules_out_file_json, 'r').read())
			logger.debug("Reading " + structure_modules_out_file_json)
	
		return clusters_dict

if __name__ == "__main__":


	accession = "P20248"

	alphafoldDataManagerObj = alphafoldDataManager()
	alphafoldDataManagerObj.cluster_predicted_aligned_error(accession)

	#alphafoldDataManagerObj.benchmark_alphafold_disprot()



	"""
	from networkx.algorithms import clusters

	clusters = clusters.greedy_modularity_clusters(g, weight='weight', resolution=1.0)
	max_len = max([len(c) for c in clusters])
	clusters = [list(c) + ['']*(max_len-len(c)) for c in clusters]
	print(clusters)
	"""
