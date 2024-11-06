import os
import sys

import inspect

file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"../"))
import option_reader

sys.path.append(os.path.join(file_path,"../accessibility"))
from find_accessibility import accessibilityFinder

from dataManager import dataManager

class accessibiltyManager(dataManager):

	##------------------------------------------------------------------##
	## Inherits functions from dataManager
	##------------------------------------------------------------------##

	def setup_data(self):
		self.default_task_options = {
			"task":'help',
			"is_superuser":True,
			"remake":False,
			"verbose":False,
			'topology_domains':"",
			'split_by_transmembrane_region':True,
			'split_by_topology_regions':True,
			"sa_cutoff":0.5,
			"add_human_readable":False,
			"print_human_readable":False,
			"show_output":False,
			"debug":False,
			"logfile":False
		}

		self.allowed_options = [
			"task",
			"outfile",
			"accession",
			"pdb_id",
			"remake",
			"outfile",
			'topology_domains',
			'split_by_transmembrane_region',
			'split_by_topology_regions',
			"smooth_disorder_length_cutoff",
			"smooth_order_length_cutoff",
			"use_homology_mapping",
			"sa_cutoff",
			"add_human_readable",
			"add_regions",
			"print_human_readable",
			"show_output",
			"output_type",
			"show_smoothing",
			"smooth_termini"
		]

		self.allowed_options_admin = [
		]

		self.task_options = [
			"calculate_accessibility",
			"calculate_disorder",
			"calculate_disorder_proportion",
			"calculate_disorder_by_pdb",
			"calculate_accessibility_consensus",
			"calculate_accessibility_consensus_by_pdb",
			"get_tessellation_accessibility",
			"help"
		]

		self.task_options_admin = []

		self.required = {
			'all':[],
			'calculate_accessibility':['accession'],
			'calculate_accessibility_consensus':['accession'],
			'calculate_accessibility_consensus_by_pdb':["pdb_id"],
			'calculate_disorder':['accession'],
			'calculate_disorder_by_pdb':["pdb_id"]
		}
		
		self.test_options = {
			'accession':"P38634",
			'pdb_id':"4LPA"
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
		
		dataDownloaderObj = accessibilityFinder()
		dataDownloaderObj.options.update(task_options)

		if task_options["task"] == "calculate_accessibility":
			data = dataDownloaderObj.calculate_accessibility()
		
		if task_options["task"] == "calculate_accessibility_consensus":
			data = dataDownloaderObj.calculate_accessibility_consensus()

		if task_options["task"] == "calculate_accessibility_consensus_by_pdb":
			data = dataDownloaderObj.calculate_accessibility_consensus_by_pdb()		

		if task_options["task"] == "calculate_disorder":
			data = dataDownloaderObj.calculate_disorder()

		if task_options["task"] == "calculate_disorder_proportion":
			data = dataDownloaderObj.calculate_disorder_proportion()

		if task_options["task"] == "calculate_disorder_by_pdb":
			data = dataDownloaderObj.calculate_disorder_by_pdb()

		if task_options["task"] == "get_tessellation_accessibility":
			data = dataDownloaderObj.get_tessellation_accessibility()
			

		return data

##------------------------------------------------------------------##
##------------------------   END    END  ---------------------------##
##------------------------------------------------------------------##

if __name__ == "__main__":

	dataManagerObj = accessibiltyManager()
	dataManagerObj.main()
