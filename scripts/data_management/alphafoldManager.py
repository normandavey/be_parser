import os
import sys

import inspect

file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"../"))
import option_reader

sys.path.append(os.path.join(file_path,"../alphafold"))
from alphafoldDataManager import alphafoldDataManager

from dataManager import dataManager

class alphafoldManager(dataManager):

	##------------------------------------------------------------------##
	## Inherits functions from dataManager
	##------------------------------------------------------------------##

	def setup_data(self):
		self.default_task_options = {
			"task":'help',
			"is_superuser":True,
			"outfile":"",
			"remake":False,
			"verbose":False,
			"debug":False,
			"logfile":False,
			"window_size":15,
			"smooth_termini":False,
			"show_smoothing":True,
			"accessibility_score_cutoff":0.55,
			"intramolecular_contact_score_cutoff":5,
			"smooth_order_length_cutoff":15,
			"smooth_disorder_length_cutoff":25,
			"parse_pfam":False,
			"split_resolution":0.4,
			"split_accessibility":False
		}

		self.allowed_options = [
			"accession",
			"alphafold_coil_cutoff",
			"alphafold_short_loop_length_cutoff",
			"window_size",
			"intramolecular_contact_score_cutoff",
			"smooth_order_length_cutoff",
			"use_pae_filtering_for_intramolecular_contacts_counts",
			"filter_inaccessibility",
			"add_atom_details"
			#"side_chain_contacts_only"
		]

		self.allowed_options_admin = [
		]

		self.task_options = [
			"calculate_scores",
			"get_pLDDT",
			"get_accessibility",
			"get_accessibility_scores",
			"get_degree",
			"get_intramolecular_contacts_scores",
			"get_intramolecular_contacts_smoothed_scores",
			"get_intramolecular_contacts_globular_smoothed_scores",
			"get_intramolecular_contacts_windowed_scores",
			"get_pLDDT_scaled_scores",
			"get_secondary_structure_scores",
			"get_weighted_degree",
			"get_accessibility_smoothed_scores",
			"get_pLDDT_windowed_scores",
			"get_pLDDT_windowed_scores_list",
			"get_pLDDT_scaled_scores_list",
			"get_accessibility_windowed_scores",
			"get_accessibility_windowed_scores_list",
			"get_secondary_structure",
			"get_secondary_structure_collapsed",
			"get_tessellation_accessibility",
			"get_domains_from_pae_matrix",
			"get_alphafold_centrality_scores",
			"get_alphafold_structural_classification",
			"help"
		]

		self.task_options_admin = [
			"calculate_scores_disprot",
			"benchmark_alphafold_disorder",
			"benchmark_alphafold_disprot",
			"benchmark_alphafold_accessibility",
			"benchmark_alphafold_motifs",
			"package_alphafold_by_species",
			"make_alphafold_file_by_species",
		]

		self.required = {
			'all':[],
			"get_alphafold_structural_classification":['accession']
		}

		self.test_options = {
			'accession':"P20248",
			}

		self.required_type = {
			'accession':"string",
			'remake':"bool"
		}

		self.required_format = {
			'accession':"\A([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\-{0,1}[0-9]*\Z"
		}

		# delimiter to split options
		self.list_options = {
			'accession':","
		}

		self.required_valid_options_list = {
			'remake':[True,False],
			'tasks':self.task_options
		}

		self.options = self.default_task_options
		
		for allowed_option in self.allowed_options + self.allowed_options_admin:
			if allowed_option not in self.options:
				self.options[allowed_option] = None

		self.options.update(option_reader.load_commandline_options(self.options,self.options,purge_commandline=True))

	##------------------------------------------------------------------##
	##------------------------------------------------------------------##

	def getData(self,status,task_options):

		###----######----######----###
		#   run_job run_job run_job  #
		###----######----######----###

		data = {}

		dataDownloaderObj = alphafoldDataManager()
		dataDownloaderObj.options.update(task_options)

		if task_options["task"] == "calculate_scores":
			data = dataDownloaderObj.calculate_scores()

		if task_options["task"] == "calculate_scores_disprot":
			data = dataDownloaderObj.calculate_scores_disprot()
		
		if task_options["task"] == "get_accessibility":
			data = dataDownloaderObj.get_accessibility()

		if task_options["task"] == "get_pLDDT":
			data = dataDownloaderObj.get_pLDDT()

		if task_options["task"] == "get_windowed_pLDDT":
			data = dataDownloaderObj.get_windowed_pLDDT()

		if task_options["task"] == "get_accessibility_windowed_scores":
			data = dataDownloaderObj.get_accessibility_windowed_scores()
		
		if task_options["task"] == "get_accessibility_scores":
			data = dataDownloaderObj.get_accessibility_scores()
		
		if task_options["task"] == "get_accessibility_windowed_scores_list":
			data = dataDownloaderObj.get_accessibility_windowed_scores_list()

		if task_options["task"] == "get_degree":
			data = dataDownloaderObj.get_degree()

		if task_options["task"] == "get_intramolecular_contacts_scores":
			data = dataDownloaderObj.get_intramolecular_contacts_scores()

		if task_options["task"] == "get_intramolecular_contacts_smoothed_scores":
			data = dataDownloaderObj.get_intramolecular_contacts_smoothed_scores()
		
		if task_options["task"] == "get_intramolecular_contacts_globular_smoothed_scores":
			data = dataDownloaderObj.get_intramolecular_contacts_globular_smoothed_scores()

		if task_options["task"] == "get_pLDDT_scaled_scores":
			data = dataDownloaderObj.get_pLDDT_scaled_scores()
			
		if task_options["task"] == "get_pLDDT_windowed_scores":
			data = dataDownloaderObj.get_pLDDT_windowed_scores()
			
		if task_options["task"] == "get_pLDDT_windowed_scores_list":
			data = dataDownloaderObj.get_pLDDT_windowed_scores_list()

		if task_options["task"] == "get_secondary_structure_scores":
			data = dataDownloaderObj.get_secondary_structure_scores()
		
		if task_options["task"] == "get_weighted_degree":
			data = dataDownloaderObj.get_weighted_degree()

		if task_options["task"] == "get_accessibility_smoothed_scores":
			data = dataDownloaderObj.get_accessibility_smoothed_scores()

		if task_options["task"] == "get_secondary_structure":
			data = dataDownloaderObj.get_secondary_structure()

		if task_options["task"] == "get_secondary_structure_collapsed":
			data = dataDownloaderObj.get_secondary_structure_collapsed()

		if task_options["task"] == "get_secondary_structure":
			data = dataDownloaderObj.get_secondary_structure()

		if task_options["task"] == "benchmark_alphafold_disorder":
			data = dataDownloaderObj.benchmark_alphafold_disorder()

		if task_options["task"] == "benchmark_alphafold_disprot":
			data = dataDownloaderObj.benchmark_alphafold_disprot()

		if task_options["task"] == "benchmark_alphafold_accessibility":
			data = dataDownloaderObj.benchmark_alphafold_accessibility()

		if task_options["task"] == "benchmark_alphafold_motifs":
			data = dataDownloaderObj.benchmark_alphafold_motifs()

		if task_options["task"] == "package_alphafold_by_species":
			data = dataDownloaderObj.package_alphafold_by_species()

		if task_options["task"] == "make_alphafold_file_by_species":
			data = dataDownloaderObj.make_alphafold_file_by_species()
		
		if task_options["task"] == "get_tessellation_accessibility":
			data = dataDownloaderObj.get_accessibility(tessellation=True)
		
		if task_options["task"] == "get_domains_from_pae_matrix":
			data = dataDownloaderObj.domains_from_pae_matrix_networkx()
		
		if task_options['task'] == 'get_alphafold_centrality_scores':
			alphaFoldPocketScorerObj = AlphaFoldPocketScorer()
			alphaFoldPocketScorerObj.options.update(task_options)
			data = alphaFoldPocketScorerObj.get_pocket_scores()
		
		if task_options['task'] == 'get_alphafold_structural_classification':
			data = dataDownloaderObj.get_alphafold_structural_classification()
			

		return data

##------------------------------------------------------------------##
##------------------------   END    END  ---------------------------##
##------------------------------------------------------------------##

if __name__ == "__main__":

	dataManagerObj = architectureManager()
	dataManagerObj.main()
