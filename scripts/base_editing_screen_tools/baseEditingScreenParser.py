import os,json,pprint,sys,inspect,math,random,copy,logging

#-----

file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../"))
import config_reader
import option_reader

sys.path.append(os.path.join(file_path,"../utilities/"))
import utilities_error 
import utilities_log
import utilities_basic

sys.path.append(os.path.join(file_path,"../data_management/"))
import queryRunner

sys.path.append(os.path.join(file_path,"../guide_design/"))
import guideAnnotation


sys.path.append(os.path.join(file_path,"../limma/"))
import limmaHelper

#-----

import base_editing_protein_plots
import base_editing_control_plots

#-----

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.debug("Logging set")

################################################################################
## SECTION I: simbaParser                                                     ##
################################################################################

class baseEditingParser():

	##################################################################################################

	def __init__(self):
		self.set_init_options()
		self.options['debug'] = True
		utilities_log.setupLogger(self.options)

	def set_init_options(self):
		self.options = config_reader.load_configeration_options(sections=['general'])
		
		self.options["job_file"] = False
		self.options["debug"] = False	
		self.options["counts_file"] = ""
		self.options["counts_directory"] = ""
		self.options["region_annotation_file"] = ""
		self.options["region_flanks_length"] = 2
		self.options["remake"] = False
		self.options["testing"] = False
		self.options['normalised_counts_mean_cutoff'] = 1000
		self.options['bootstrap_size'] = 100000
		self.options['plot_controls'] = False
		self.options["r_lib_loc"] = "./"
		self.options['query_time_point'] = 0
		self.options['compare_time_point'] = 18

		self.options['negative_control_labels'] = ["Intergenic control","Non-targeting control"]
		self.options['positive_control_labels'] = ["Essential splice site"]

		self.options["guide_annotation_file_format"] = "generic"
		
		self.options['add_depmap_data'] = True
		self.options["output_include_region_annotation_json"] = True
		self.options["output_include_protein_annotation_json"] = True
		self.options['run_peptools'] = True


		self.options["filter_by_behive"] = False


		self.options['make_region_sliding_window'] = False
		self.options['make_region_tesselation'] = False
		
		##---------------------------------##

		self.options.update(option_reader.load_commandline_options(self.options,{}))

		self.options['initial_options'] = self.options

		self.errors = []

	##---------------------------------##

	def reset_options(self):
		self.options = self.options['initial_options']
		self.options['counts_directory'] = ""
	
	##---------------------------------##
		
	def read_job_file(self):
		if self.options["job_file"]:
			logger.debug("Loading: " + self.options["job_file"])
			self.options.update(json.loads(open(self.options["job_file"]).read()))
		else:
			logger.error("No job_file selected")

	def setupScreen(self):
		
		self.data = {
			'coverage':{},
			'guides_status':{}
		}
		
		##---------------------------------##
		
		if not os.path.exists(self.options['data_directory']):
			print(os.path.abspath(self.options['data_directory']))
			os.mkdir(self.options['data_directory'])
		
		if not os.path.exists(os.path.join(self.options['data_directory'],"regions")):
			os.mkdir(os.path.join(self.options['data_directory'],"regions"))

		if not os.path.exists(os.path.join(self.options['data_directory'],"proteins")):
			os.mkdir(os.path.join(self.options['data_directory'],"proteins"))

		if not os.path.exists(os.path.join(self.options['data_directory'],"tdt")):
			os.mkdir(os.path.join(self.options['data_directory'],"tdt"))

		if not os.path.exists(os.path.join(self.options['data_directory'],"cache")):
			os.mkdir(os.path.join(self.options['data_directory'],"cache"))

		if not os.path.exists(os.path.join(self.options['data_directory'],"coverage")):
			os.mkdir(os.path.join(self.options['data_directory'],"coverage"))

		if not os.path.exists(os.path.join(self.options['data_directory'],"json")):
			os.mkdir(os.path.join(self.options['data_directory'],"json"))

		##---------------------------------##
		
		self.options["design"]["region_annotation"] = {}

		##---------------------------------##

		self.control_types = {
			"Intergenic control":"Intergenic",
			"Essential splice site":"Splice Site",
			"Plasmid backbone sequence":"Backbone",
			"Non-targeting control":"Non-targeting"
		}
		
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
		
	def load_guides_mapping_file(self):
		if 'guide_mapping_file' in self.options:
			logger.debug("load_guides_annotation_file: " + self.options["guide_mapping_file"]) 
			lines = open(self.options["guide_mapping_file"]).read().split("\n")
				
			self.options["design"]['guide_mapping'] = {}
			self.options["design"]['guide_mapping_id'] = {}

			for line in lines[1:]:
				try:
					line_bits = line.split(",")
					self.options["design"]['guide_mapping'][line_bits[0]] = line_bits[1]
					self.options["design"]['guide_mapping_id'][line_bits[1]] = line_bits[0]
				except:
					pass
		
	##################################################################################################
	
	def load_time_points(self):
		lines = open(self.options["counts_file"]).read().strip().split("\n")

		headers = lines[0].split("\t")
		for header_iter in range(0,len(headers)):
			if headers[header_iter] in self.options["design"]["datapoint_labels"]:
				datapoint_label = self.options["design"]["datapoint_labels"][headers[header_iter]]
				if datapoint_label["timepoint"] not in self.data[self.options['screen_name']]['compare_time_points'] and datapoint_label["timepoint"] != self.options['query_time_point']:
					self.data[self.options['screen_name']]['compare_time_points'].append(datapoint_label["timepoint"])

		del lines

	##################################################################################################

	def load_counts_file(self):
		logger.info("Reading " + self.options["counts_file"])

		lines = open(self.options["counts_file"]).read().strip().split("\n")
		count_file_format = "v1"

		if "count_file_format" in self.options:
			count_file_format = self.options["count_file_format"]

		headers = lines[0].split("\t")
		uniprot_check = {}

		grna_statistics = {
			"observed":0,
			"contains_zero_value":0,
			"contains_low_count_value":0
		}

		gene_name_mapping = {}
		essentiality = {}

		for line in lines[1:]:
			try:
				line_bits = line.split("\t")

				##-------------------------##
	
				if 'guide_mapping_file' in self.options:
					if line_bits[0] in self.options["design"]['guide_mapping']:
						guide_sequence = self.options["design"]['guide_mapping'][line_bits[0]] 
					else:
						if count_file_format in ['v1','v2','v3']:
							guide_sequence = line_bits[0]
						if count_file_format in ['v4','v5']:
							guide_sequence = line_bits[1]
				else:
					if count_file_format in ['v1','v2','v3']:
						guide_sequence = line_bits[0]
					if count_file_format in ['v4','v5']:
						guide_sequence = line_bits[1]
				
				if guide_sequence not in self.data[self.options['screen_name']]['ngs_data']:
					self.data[self.options['screen_name']]['ngs_data'][guide_sequence] = {
						"data":{},
						"data_t0_normalised":{},
						"data_log2_cpm_normalised":{}
					}

				contains_zero_value = False
				contains_low_count_value = False

				for header_iter in range(0,len(headers)):
					if headers[header_iter] == "": continue
					if headers[header_iter] in self.options["design"]["datapoint_labels"]:
						datapoint_label = self.options["design"]["datapoint_labels"][headers[header_iter]]
						if datapoint_label["timepoint"] not in [self.options['query_time_point'],self.options['compare_time_point']]: continue
						if datapoint_label["timepoint"] not in self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['data']:
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['data'][datapoint_label["timepoint"]] = {}

						if float(line_bits[header_iter]) == 0:
							contains_zero_value = True

						if float(line_bits[header_iter]) < 30 and datapoint_label["timepoint"] == self.options['query_time_point']:
							contains_low_count_value = True

						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['data'][datapoint_label["timepoint"]][datapoint_label["replicate"]] = float(line_bits[header_iter])
					else:
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence][headers[header_iter]] = line_bits[header_iter]
				
				if contains_zero_value:
					grna_statistics["contains_zero_value"] += 1
					self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['contains_zero_value'] = True
				else:
					if 'contains_zero_value' not in self.data[self.options['screen_name']]['ngs_data'][guide_sequence]:
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['contains_zero_value'] = False
				
				if contains_low_count_value:
					grna_statistics["contains_low_count_value"] += 1
					self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['contains_low_count_value'] = True
				else:
					if 'contains_low_count_value' not in self.data[self.options['screen_name']]['ngs_data'][guide_sequence]:
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['contains_low_count_value'] = False
					#continue

				guide_info = self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['Gene']
				del self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['Gene']

				guide_info_bits = guide_info.split("|")
				
				##-------------------------##
	
				if count_file_format == 'v1':
					if len(guide_info_bits) == 1:
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = True
					
						if guide_info_bits[0].count("Intergenic control"):
							control_type = "Intergenic control"
						elif guide_info_bits[0].count("Non-targeting control"):
							control_type = "Non-targeting control"
						elif guide_info_bits[0].count("Essential splice site"):
							control_type = "Essential splice site"
						elif guide_info_bits[0].count("Plasmid backbone sequence"):
							control_type = "Plasmid backbone sequence"
						elif guide_info_bits[0].count("backbone sequence"):
							control_type = "Plasmid backbone sequence"
						else:
							control_type = "Unknown: " + guide_info_bits[0]

						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['seq'] = guide_sequence
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] = control_type

						if control_type not in self.data[self.options['screen_name']]['controls']:
							self.data[self.options['screen_name']]['controls'][control_type] = []

						self.data[self.options['screen_name']]['controls'][control_type].append(guide_sequence)
					else:
						accession = guide_info_bits[0]
						accession,uniprot_check = self.update_obsolete_accesion(accession,uniprot_check)
						
						
						if "filter_accession" in self.options:
							if accession not in self.options["filter_accession"]:
								del self.data[self.options['screen_name']]['ngs_data'][guide_sequence]
								continue
						
						try:
							pam_type = self.options["design"]['guide_annotation'][guide_sequence]['pam_type']
							annotation = self.options["design"]['guide_annotation'][guide_sequence]['annotation']
						except:
							pam_type = "-"
							annotation = "-"

						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = False
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['seq'] = guide_sequence
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['pam_type'] = pam_type
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation'] = annotation
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = False
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] = False
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession'] = accession 
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['gene'] = guide_info_bits[1]
						
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_start'] = guide_info_bits[2]
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_end'] = guide_info_bits[3]
					
						if len(guide_info_bits[4]) > 0:
							cbe_edit_bits = guide_info_bits[4].split(":")
							if len(cbe_edit_bits) < 4: continue
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['cbe_edit_mutation'] = cbe_edit_bits[1]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['cbe_edit_wildtype_peptide'] = cbe_edit_bits[2]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['cbe_edit_mutant_peptide'] = cbe_edit_bits[3]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['cbe_edit_mutation_details'] = self.process_mutations(cbe_edit_bits[1])
						
						if len(guide_info_bits[5]) > 0:
							abe_edit_bits = guide_info_bits[5].split(":")
							if len(abe_edit_bits) < 4: continue
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['abe_edit_mutation'] = abe_edit_bits[1]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['abe_edit_wildtype_peptide'] = abe_edit_bits[2]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['abe_edit_mutant_peptide'] = abe_edit_bits[3]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['abe_edit_mutation_details'] = self.process_mutations(abe_edit_bits[1])	
				
				elif count_file_format == 'v3':
					if len(guide_info_bits) == 1:
						accession = guide_info_bits[0]
						accession,uniprot_check = self.update_obsolete_accesion(accession,uniprot_check)
						
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = True
					
						if guide_info_bits[0].count("Intergenic control"):
							control_type = "Intergenic control"
						elif guide_info_bits[0].count("Non-targeting control"):
							control_type = "Non-targeting control"
						elif guide_info_bits[0].count("Essential splice site"):
							control_type = "Essential splice site"
						elif guide_info_bits[0].count("Plasmid backbone sequence"):
							control_type = "Plasmid backbone sequence"
						elif guide_info_bits[0].count("backbone sequence"):
							control_type = "Plasmid backbone sequence"
						else:
							control_type = "Unknown: " + guide_info_bits[0]
							
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] = control_type
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['seq'] = guide_sequence

						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['gene'] = control_type
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession'] = accession

						if control_type not in self.data[self.options['screen_name']]['controls']:
							self.data[self.options['screen_name']]['controls'][control_type] = []

						self.data[self.options['screen_name']]['controls'][control_type].append(guide_sequence)
					else:
						accession = guide_info_bits[0]
						accession,uniprot_check = self.update_obsolete_accesion(accession,uniprot_check)

						if "filter_accession" in self.options:
							if accession not in self.options["filter_accession"]:
								del self.data[self.options['screen_name']]['ngs_data'][guide_sequence]
								continue
							
						try:
							pam_type = self.options["design"]['guide_annotation'][guide_sequence]['pam_type']
							annotation = self.options["design"]['guide_annotation'][guide_sequence]['annotation']
						except:
							pam_type = "-"
							annotation = "-"

						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = False
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['seq'] = guide_sequence
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['pam_type'] = pam_type
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation'] = annotation
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = False
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] = False
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession'] = accession 
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['gene'] = guide_info_bits[1]
					
						try:
							try:
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_start'] = guide_info_bits[2].split(":")[2].split("-")[0]
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_end'] = guide_info_bits[2].split(":")[2].split("-")[0]
							except:
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_start'] = guide_info_bits[3].split(":")[2].split("-")[0]
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_end'] = guide_info_bits[3].split(":")[2].split("-")[0]
						except:
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_start'] = guide_info_bits[2]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_end'] = guide_info_bits[2]
							
						
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['abe_edit_mutation'] = guide_info_bits[2]
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['cbe_edit_mutation'] = guide_info_bits[3]
					
					
				elif count_file_format == 'v2':
					if 'annotation' not in self.data[self.options['screen_name']]['ngs_data'][guide_sequence]:
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation'] = ""
				
					if self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation'].count("Control") > 0:
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = True
						
						if self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation'].count("Intergenic control"):
							control_type = "Intergenic control"
						elif self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation'].count("Non-targeting control"):
							control_type = "Non-targeting control"
						elif self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation'].count("Essential splice site"):
							control_type = "Essential splice site"
						elif self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation'].count("Plasmid backbone sequence"):
							control_type = "Plasmid backbone sequence"
						elif self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation'].count("backbone sequence"):
							control_type = "Plasmid backbone sequence"
						else:
							control_type = "Unknown: " +self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation']

						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] = control_type
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['seq'] = guide_sequence

						if control_type not in self.data[self.options['screen_name']]['controls']:
							self.data[self.options['screen_name']]['controls'][control_type] = []

						self.data[self.options['screen_name']]['controls'][control_type].append(guide_sequence)
					else:
						accession = self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation']

						if accession not in uniprot_check:
							try:
								uniprot_check_response = queryRunner.queryRunner("uniprot","check_uniprot_entry",{"accession":accession}).run()

								if 'updated_accession' in uniprot_check_response:
									uniprot_check[accession] = uniprot_check_response['updated_accession']
									logger.error("Updating obsolete UniProt accession: " + accession + " " + uniprot_check_response['updated_accession'])
							except:
								raise
						
						if accession in uniprot_check: accession = uniprot_check[accession]

						if "filter_accession" in self.options:
							if accession not in self.options["filter_accession"]:
								del self.data[self.options['screen_name']]['ngs_data'][guide_sequence]
								continue
							
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = False
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['seq'] = guide_sequence
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] = False

						try:
							if len(self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['Predicted_ABE_edit']) > 2:
								abe_edit_bits = self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['Predicted_ABE_edit'].split(":")
								peptide_start = min(list(self.process_mutations(abe_edit_bits[1]).keys()))
								peptide_end = max(list(self.process_mutations(abe_edit_bits[1]).keys()))

							if len(self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['Predicted_CBE_edit']) > 2:
								cbe_edit_bits = self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['Predicted_CBE_edit'].split(":")
								peptide_start = min(list(self.process_mutations(cbe_edit_bits[1]).keys()))
								peptide_end = max(list(self.process_mutations(cbe_edit_bits[1]).keys()))
						except:
							peptide_start = "1"
							peptide_end = "1"

						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession'] = accession
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['gene'] = self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['Genename']
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_start'] = peptide_start
						self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_end'] = peptide_end
						
						if len(self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['Predicted_CBE_edit']) > 2:
							cbe_edit_bits = self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['Predicted_CBE_edit'].split(":")
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['cbe_edit_mutation'] = cbe_edit_bits[1]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['cbe_edit_wildtype_peptide'] = cbe_edit_bits[2]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['cbe_edit_mutant_peptide'] = cbe_edit_bits[3]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['cbe_edit_mutation_details'] = self.process_mutations(cbe_edit_bits[1])
						
						if len(self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['Predicted_ABE_edit']) > 2:
							abe_edit_bits = self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['Predicted_ABE_edit'].split(":")
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['abe_edit_mutation'] = abe_edit_bits[1]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['abe_edit_wildtype_peptide'] = abe_edit_bits[2]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['abe_edit_mutant_peptide'] = abe_edit_bits[3]
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['abe_edit_mutation_details'] = self.process_mutations(abe_edit_bits[1])	
				
				elif count_file_format in ['v4','v5']:
					if count_file_format == 'v4':
						guide_sequence = guide_info

					if guide_sequence in self.options["design"]['guide_annotation']:

						guide_info_bits = self.options["design"]['guide_annotation'][guide_sequence]['annotation'].split("|")
						
						non_target_type = False
						if 'type' in self.options["design"]['guide_annotation'][guide_sequence]:
							if self.options["design"]['guide_annotation'][guide_sequence]['type'] not in ['targets']:
								non_target_type = True
						
						if non_target_type:
							gene_name = guide_info_bits[0]
							
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = True

							control_type = self.options["design"]['guide_annotation'][guide_sequence]['type']

							if self.options["design"]['guide_annotation'][guide_sequence]['type'] == 'control_splice_sites':
								control_type = "Essential splice site " + guide_info_bits[1] + guide_info_bits[2]
							elif self.options["design"]['guide_annotation'][guide_sequence]['type'] == 'splice_sites':
								control_type = "Control splice site"
							elif self.options["design"]['guide_annotation'][guide_sequence]['type'] == 'non_targetting_control':
								control_type = "Non-targeting control"
							elif self.options["design"]['guide_annotation'][guide_sequence]['type'] == 'intergenic_control':
								control_type = "Intergenic control"
							else:
								control_type = "Unknown: " + guide_info_bits[0]

							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] = control_type
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['seq'] = guide_sequence

							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['gene'] = self.options["design"]['guide_annotation'][guide_sequence]['annotation']

							if gene_name not in gene_name_mapping:
								gene_name_mapping_instance = queryRunner.queryRunner("uniprot","get_uniprot_mapping_details",{"identifier":[gene_name],"mapping_from":"GENENAME","mapping_to":"ACC"}).run()["data"]
								if gene_name in gene_name_mapping_instance:
									gene_name_mapping[gene_name] = gene_name_mapping_instance[gene_name]['Entry']
								else:
									gene_name_mapping[gene_name] = "-"

							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession'] = gene_name_mapping[gene_name] 

							if control_type not in self.data[self.options['screen_name']]['controls']:
								self.data[self.options['screen_name']]['controls'][control_type] = []

							self.data[self.options['screen_name']]['controls'][control_type].append(guide_sequence)

						elif self.options["design"]['guide_annotation'][guide_sequence]['annotation'].count("|") == 0:
							accession = guide_info_bits[0]
							accession,uniprot_check = self.update_obsolete_accesion(accession,uniprot_check)
							
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = True
						
							if guide_info_bits[0].count("Intergenic control"):
								control_type = "Intergenic control"
							elif guide_info_bits[0].count("Non-targeting control"):
								control_type = "Non-targeting control"
							elif guide_info_bits[0].count("Essential splice site"):
								control_type = "Essential splice site"
							elif guide_info_bits[0].count("Plasmid backbone sequence"):
								control_type = "Plasmid backbone sequence"
							elif guide_info_bits[0].count("backbone sequence"):
								control_type = "Plasmid backbone sequence"
							else:
								control_type = "Unknown: " + guide_info_bits[0]
							
							
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] = control_type
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['seq'] = guide_sequence

							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['gene'] = control_type
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession'] = accession

							if control_type not in self.data[self.options['screen_name']]['controls']:
								self.data[self.options['screen_name']]['controls'][control_type] = []

							self.data[self.options['screen_name']]['controls'][control_type].append(guide_sequence)
						else:
							accession = guide_info_bits[0]
							
							accession,uniprot_check = self.update_obsolete_accesion(accession,uniprot_check)

							if "filter_accession" in self.options:
								if accession not in self.options["filter_accession"]:
									del self.data[self.options['screen_name']]['ngs_data'][guide_sequence]
									continue
								
							try:
								pam_type = self.options["design"]['guide_annotation'][guide_sequence]['pam_type']
								annotation = self.options["design"]['guide_annotation'][guide_sequence]['annotation']
							except:
								pam_type = "-"
								annotation = "-"

							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = False
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['seq'] = guide_sequence
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['pam_type'] = pam_type
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['annotation'] = annotation
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control'] = False
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] = False
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession'] = accession 
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['gene'] = guide_info_bits[1]
							try:
								try:
									self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_start'] = guide_info_bits[2].split(":")[2].split("-")[0]
									self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_end'] = guide_info_bits[2].split(":")[2].split("-")[0]
								except:
									self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_start'] = guide_info_bits[3].split(":")[2].split("-")[0]
									self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_end'] = guide_info_bits[3].split(":")[2].split("-")[0]
							except:
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_start'] = guide_info_bits[2]
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['peptide_end'] = guide_info_bits[2]
							
							try:
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['abe_edit_mutation'] = guide_info_bits[2]
							except:
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['abe_edit_mutation'] = ""

							try:
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['cbe_edit_mutation'] = guide_info_bits[3]
							except:
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['cbe_edit_mutation'] = ""

						
						if self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession'] not in essentiality:
							hap1_gene_effect = {}
							mean_gene_effect = ""
							
							try:
								mean_gene_effect = queryRunner.queryRunner("depmap","get_CRISPR_gene_dependency",{"accession":self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession']}).run()['data']['mean_gene_effect']
							except:
								pass

							try:
								hap1_gene_effect = queryRunner.queryRunner("depmap","get_hap1_gene_dependency_data",{"accession":self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession']}).run()['data']
								if len(hap1_gene_effect) == 0:
									hap1_gene_effect = {'gene_name': '', 'p_value': '-', 'q_value': '-', 'ratio': '-', 'selected': ''}
							except:
								pass
							
							essentiality[self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession']] = {
								"depmap":mean_gene_effect,
								"hap1":hap1_gene_effect
							}
						
						try:
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['essentiality_hap1'] = essentiality[self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession']]['hap1']['selected']
							if self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] == "Control splice site" and self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['essentiality_hap1'] == "YES":
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] = "Control splice site essential"
						except:
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['essentiality_hap1'] = "-"

						try:
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['essentiality_depmap'] = essentiality[self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession']]['depmap']
							if self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] == "Control splice site" and float(self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['essentiality_depmap']) < -0.5:
								self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control_type'] = "Control splice site essential"
						except:
							self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['essentiality_depmap'] = "-"

					else:
						del self.data[self.options['screen_name']]['ngs_data'][guide_sequence]
						logger.error(guide_sequence + " not in design")
			except:
				logger.debug("ERROR - : " + guide_sequence + " " + self.options['screen_name'])
				utilities_error.printError()

		del lines

	##################################################################################################
	
	def load_guides_annotation_file(self):
		logger.debug("load_guides_annotation_file: " + self.options["guide_annotation_file"]) 
		lines = open(self.options["guide_annotation_file"]).read().split("\n")
			
		self.options["design"]['guide_annotation'] = {}
		self.data['coverage']["design_guide_total"] = 0

		for line in lines[1:]:
			try:
				if self.options["guide_annotation_file_format"] == "reporter":
					line_bits = line.split("\t")
					self.options["design"]['guide_annotation'][line_bits[2]] = {
						"pam_type":line_bits[3].strip('"'),
						"annotation":line_bits[4].strip('"'),
						"type":line_bits[8]
					}
				elif self.options["guide_annotation_file_format"] == "reporter_v2":
					line_bits = line.split("\t")
					self.options["design"]['guide_annotation'][line_bits[2]] = {
						"pam_type":"NGN",
						"annotation":line_bits[5].strip('"'),
						"type":line_bits[8]
					}
				else:
					line_bits = line.split("\t")
					self.options["design"]['guide_annotation'][line_bits[2]] = {
						"pam_type":line_bits[3].strip('"'),
						"annotation":line_bits[4].strip('"')
					}

				self.data['coverage']["design_guide_total"] += 1
			except:
				pass
		
		logger.debug("Guide count:" + str(self.data['coverage']["design_guide_total"]))
	
	##################################################################################################

	def load_region_annotation_file(self):
		logger.debug("load_region_annotation_file")
		lines = open(self.options["region_annotation_file"]).read().split("\n")

		logger.debug("Region count:" + str(len(lines)))

		headers = lines[0].split("\t")
		headers_to_list = ['dataset','pdb','pmid','source','specficity_class','type','slimprint_p']
		
		region_annotation = {}

		protein_regions = {}
		protein_residues = {}

		for line in lines[1:]:
			try:
				line_bits = line.split("\t")
				region_id = "_".join([line_bits[2],line_bits[4],line_bits[5]])
		
				region_tmp = {}

				
				if line_bits[2] not in protein_regions:
					protein_regions[line_bits[2]] = {}
					protein_residues[line_bits[2]] = {}

				for offset in list(range(int(line_bits[4]),int(line_bits[5])+1)):
					if offset not in protein_residues[line_bits[2]]:
						protein_residues[line_bits[2]][offset] = []

					protein_residues[line_bits[2]][offset].append(region_id)
				
				protein_regions[line_bits[2]][region_id] = list(range(int(line_bits[4]),int(line_bits[5])+1))
			except:
				utilities_error.printError()

		#---##---##---##---##---##---##---##---##---##---##---##---##---##---##---##---#
		### COLLAPSE OVERLAPPING REGIONS 
		#---##---##---##---##---##---##---##---##---##---##---##---##---##---##---##---#
		
		region_mapping_id = {}
		for accession in protein_residues:
			for group in self.find_consecutive_groups(list(protein_residues[accession].keys())):
				update_region_id = "_".join([accession,str(min(group)),str(max(group))])
				for offset in group:
					for region_id in protein_residues[accession][offset]:
						region_mapping_id[region_id] = {"id":update_region_id,"peptide_start":str(min(group)),"peptide_end":str(max(group))}

		for line in lines[1:]:
			try:
				line_bits = line.split("\t")
				region_id = "_".join([line_bits[2],line_bits[4],line_bits[5]])

				if int(line_bits[1]) == 0:
					logger.error("No mutations: " + line)
					continue
			
				region_use_id = region_mapping_id[region_id]['id']

				if region_use_id not in region_annotation:
					region_annotation[region_use_id] = {}
					region_annotation[region_use_id]["regions"] = {}
					for header_iter in range(0,len(headers)):
						if headers[header_iter] in headers_to_list:
							region_annotation[region_use_id][headers[header_iter]] = [line_bits[header_iter]]
						else:
							region_annotation[region_use_id][headers[header_iter]] = line_bits[header_iter]

				elif region_use_id in self.options["design"]['region_annotation']:
					for header_iter in range(0,len(headers)):
						if headers[header_iter] in headers_to_list:
							region_annotation[region_use_id][headers[header_iter]] += [line_bits[header_iter]]


				region_annotation[region_use_id]['peptide_start'] = region_mapping_id[region_id]['peptide_start']
				region_annotation[region_use_id]['peptide_end'] = region_mapping_id[region_id]['peptide_end']

				#---##---##---##---##---##---##---##---##---##---##---##---##---##---##---##---#

				region_tmp = {}
				for header_iter in range(0,len(headers)):
					region_tmp[headers[header_iter]] = line_bits[header_iter]

				region_annotation[region_use_id]["regions"][region_id] = region_tmp
			except:
				utilities_error.printError()
		
		for region_annotation_id in region_annotation:
			self.options["design"]['region_annotation'][region_annotation_id] = region_annotation[region_annotation_id]
			
		logger.debug("Region count unique:" + str(len(self.options["design"]['region_annotation'])))
		
		self.add_guide_id_to_region_annotation()
		

	##################################################################################################

	def plot_count_data(self):
		baseEditingControlPlotterObj = base_editing_control_plots.baseEditingControlPlotter()
		baseEditingControlPlotterObj.data = self.data[self.options['screen_name']]['ngs_data']
		baseEditingControlPlotterObj.options['data_directory'] = self.options['data_directory']
		baseEditingControlPlotterObj.options['screen_name'] = self.options["original_screen_name"] 
		baseEditingControlPlotterObj.options['compare_time_point'] = self.options['compare_time_point']

		baseEditingControlPlotterObj.plot_count_data()

	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
			
	def normalise_t0_files(self):
		reference_timepoint = self.options['query_time_point']
		for guide_id in self.data[self.options['screen_name']]['ngs_data']:
			for timepoint in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data']:
				self.data[self.options['screen_name']]['ngs_data'][guide_id]['data_t0_normalised'][timepoint] = {}
				for replicate in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][timepoint]:
					try:
						t_start =     self.data[self.options['screen_name']]['ngs_data'][guide_id]['data_log2_cpm_normalised'][timepoint][replicate]
						t_reference = self.data[self.options['screen_name']]['ngs_data'][guide_id]['data_log2_cpm_normalised'][reference_timepoint][replicate]
						data_t0_normalised =  t_start - t_reference
						self.data[self.options['screen_name']]['ngs_data'][guide_id]['data_t0_normalised'][timepoint][replicate] = data_t0_normalised
					except:
						self.data[self.options['screen_name']]['ngs_data'][guide_id]['data_t0_normalised'][timepoint][replicate] = 0

	##################################################################################################

	def plot_t0_normalised_data(self):
		baseEditingControlPlotterObj = base_editing_control_plots.baseEditingControlPlotter()
		baseEditingControlPlotterObj.data = self.data[self.options['screen_name']]['ngs_data']
		baseEditingControlPlotterObj.options['data_directory'] = self.options['data_directory']
		baseEditingControlPlotterObj.options['screen_name'] = self.options["original_screen_name"] 
		baseEditingControlPlotterObj.options['compare_time_point'] = self.options['compare_time_point']

		baseEditingControlPlotterObj.plot_control_t0_normalised_data()

	##################################################################################################
				
	def normalise_cpm_log2_files(self):
		sample_total_counts = {}

		for guide_id in self.data[self.options['screen_name']]['ngs_data']:
			for timepoint in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data']:
				if timepoint not in sample_total_counts:sample_total_counts[timepoint] = {}
				for replicate in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][timepoint]:
					if replicate not in sample_total_counts[timepoint]:sample_total_counts[timepoint][replicate] = 0 
					sample_total_counts[timepoint][replicate] += self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][timepoint][replicate]

				
		for guide_id in self.data[self.options['screen_name']]['ngs_data']:
			for timepoint in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data']:
				self.data[self.options['screen_name']]['ngs_data'][guide_id]['data_log2_cpm_normalised'][timepoint] = {}
				for replicate in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][timepoint]:
					try:
						self.data[self.options['screen_name']]['ngs_data'][guide_id]['data_log2_cpm_normalised'][timepoint][replicate] = math.log2((self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][timepoint][replicate]/sample_total_counts[timepoint][replicate])*1000000)
					except:
						self.data[self.options['screen_name']]['ngs_data'][guide_id]['data_log2_cpm_normalised'][timepoint][replicate] = 0

	##################################################################################################
	
	def process_mutations(self,mutation_formatted):
		mutations = {}
		for mutation in mutation_formatted.split("__"):
			mutation_bits = mutation.split("->")
			mutations[int(mutation_bits[0][1:])] = {
				"wt":mutation_bits[0][0],
				"mut":mutation_bits[1][0]
			}

		return mutations

	##################################################################################################
				
	def process_count_statistics(self):
		logger.debug("process_count_statistics initialised")
		for guide_id in self.data[self.options['screen_name']]['ngs_data']:
			try:
				self.data[self.options['screen_name']]['ngs_data'][guide_id]["normalised_counts"] = []
				self.data[self.options['screen_name']]['ngs_data'][guide_id]["counts"] = []				
				
				for timepoint in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data']:
					if timepoint == self.options['compare_time_point']:		
						for replicate in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][timepoint]:
							self.data[self.options['screen_name']]['ngs_data'][guide_id]["normalised_counts"].append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['data_t0_normalised'][timepoint][replicate])
							self.data[self.options['screen_name']]['ngs_data'][guide_id]["counts"].append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][timepoint][replicate])
				
				normalised_counts_mean = sum(self.data[self.options['screen_name']]['ngs_data'][guide_id]["normalised_counts"])/len(self.data[self.options['screen_name']]['ngs_data'][guide_id]["normalised_counts"])
				if normalised_counts_mean == 0: continue

				self.data[self.options['screen_name']]['ngs_data'][guide_id]["normalised_counts_mean"] = math.log2(normalised_counts_mean)

				counts_mean = sum(self.data[self.options['screen_name']]['ngs_data'][guide_id]["counts"])/len(self.data[self.options['screen_name']]['ngs_data'][guide_id]["counts"])
				if counts_mean == 0: continue
				self.data[self.options['screen_name']]['ngs_data'][guide_id]["counts_mean"] = counts_mean
			except:
				logger.error(guide_id + " failed")
				pprint.pprint(self.data[self.options['screen_name']]['ngs_data'][guide_id]["normalised_counts"])
	
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
		
	def check_guide_for_filtering(self,guide_id,return_reason=False):
		if self.data[self.options['screen_name']]['ngs_data'][guide_id]['contains_low_count_value']: 
			if return_reason: return 'contains_low_count_value'
			return True
		if self.data[self.options['screen_name']]['ngs_data'][guide_id]['contains_zero_value']: 
			if return_reason: return 'contains_zero_value'
			return True
		
		if 'guide_off_target' in self.data[self.options['screen_name']]:
			if guide_id in self.data[self.options['screen_name']]['guide_off_target' ]:
				try:
					if self.data[self.options['screen_name']]['guide_off_target'][guide_id]['0'] > 5:
						if return_reason: return 'guide_off_target'
						return True
				except:
					logger.error(guide_id + " not in guide_off_target")
		else:
			logger.error(guide_id + " not in guide_off_target")

		if self.options["filter_by_behive"]:
			if 'guide_efficiency' in self.data[self.options['screen_name']]:
				if guide_id in self.data[self.options['screen_name']]['guide_efficiency']:
					if 'uniprot_offsets' not in self.data[self.options['screen_name']]['guide_efficiency'][guide_id]:
						if return_reason: return 'beehive_no_edits'
						return True
					if len(self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['uniprot_offsets']) == 0:
						if return_reason: return 'beehive_no_edits'
						return True
					if 0 in self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['uniprot_offsets']:
						if return_reason: return 'beehive_pre_cterm_edits'
						return True
					if 1 in self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['uniprot_offsets']:
						if return_reason: return 'beehive_cterm_edits'
						return True
			else:
				logger.error(guide_id + " not in guide_efficiency")
			
		return False
	
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	
	def get_mutation_tag(self):
		if self.options["editor_type"] == "ABE":
			mutation_tag = 'abe_edit_mutation'
		if self.options["editor_type"] == "CBE":
			mutation_tag = 'cbe_edit_mutation'
		if self.options["editor_type"] == "Dual":
			mutation_tag = 'dual_edit_mutation'
		
		return mutation_tag
	
	def annotate_guide_status(self):
		for region_id in self.options["design"]['region_annotation']:
			self.options["design"]['region_annotation'][region_id]['guides_status'] = {}
			if "guide_ids" in self.options["design"]['region_annotation'][region_id]:
				for guide_id in self.options["design"]['region_annotation'][region_id]["guide_ids"]:
					guide_check = self.check_guide_for_filtering(guide_id,return_reason=True)

					if region_id.count("Q9Y5A9") > 0:
						print(region_id,region_id.count("Q9Y5A9"),guide_id,guide_check)
						print(guide_id,self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change_guide'])

					if guide_check != False:
						self.data['guides_status'][guide_id] = guide_check
						continue
					else:
						self.data['guides_status'][guide_id] = ""

					try:
						mutation_tag =self.get_mutation_tag()

						if mutation_tag not in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
							self.data['guides_status'][guide_id] = "no_predicted_edits"
							continue
							
						if self.data[self.options['screen_name']]['ngs_data'][guide_id][mutation_tag] == "":
							self.data['guides_status'][guide_id] = "predicted_edits_empty"
							continue
					
						if "limma_scores" in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
							pass
						else:
							self.data['guides_status'][guide_id] = "no_limma"
							continue
					except:
						utilities_error.printError()

				for guide_id in self.options["design"]['region_annotation'][region_id]["guide_ids"]:
					self.options["design"]['region_annotation'][region_id]['guides_status'][guide_id] = self.data['guides_status'][guide_id]
		
	#################################################################################################	

	def mw_region_statistics(self):
		logger.debug("mw_region_statistics initialised")

		from scipy.stats import mannwhitneyu

		log_fold_changes = []
		p_value_adjs = []

		logger.debug(len(self.data[self.options['screen_name']]['ngs_data']))

		for guide_id in self.data[self.options['screen_name']]['ngs_data']:  
			if self.check_guide_for_filtering(guide_id): continue
			
			if 'control' not in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
				logger.error("Control status not in:" + guide_id)
				continue
			
			if self.data[self.options['screen_name']]['ngs_data'][guide_id]['control'] == False:
				try:
					mutation_tag =self.get_mutation_tag()

					if mutation_tag not in self.data[self.options['screen_name']]['ngs_data'][guide_id]:continue
					if self.data[self.options['screen_name']]['ngs_data'][guide_id][mutation_tag] == "":continue

					if "limma_scores" in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
						log_fold_changes.append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'])
						p_value_adjs.append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['p_value_adj'])
				except:
					logger.error(guide_id + " p_value_adjs")	
		
		for region_id in self.options["design"]['region_annotation']:
			if "guide_ids" in self.options["design"]['region_annotation'][region_id]:
				log_fold_changes_sample = []
				p_value_adjs_sample = []
				unfiltered_guides = []
				
				for guide_id in self.options["design"]['region_annotation'][region_id]["guide_ids"]:
					guide_check = self.check_guide_for_filtering(guide_id,return_reason=True)
					
					if guide_check != False:
						continue

					try:
						mutation_tag = self.get_mutation_tag()

						if mutation_tag not in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
							continue
							
						if self.data[self.options['screen_name']]['ngs_data'][guide_id][mutation_tag] == "":
							continue
						
						if "limma_scores" in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
							log_fold_changes_sample.append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'])
							p_value_adjs_sample.append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['p_value_adj'])
							unfiltered_guides.append(guide_id)
						else:
							continue
					except:
						utilities_error.printError()

				### Check for sufficient data for stats calculation
				if len(unfiltered_guides) <=1:
					self.options["design"]['region_annotation'][region_id].update({
						"mw_p":1,
						"mw_p_log":0,
						"mean_fold_change":0,
						"unfiltered_guides":unfiltered_guides,
						"fold_changes":log_fold_changes_sample,
						"guides":{}
					})

					for guide_id in self.options["design"]['region_annotation'][region_id]["guide_ids"]:
						self.options["design"]['region_annotation'][region_id]["guides"][guide_id] = self.data[self.options['screen_name']]['ngs_data'][guide_id]
						self.data[self.options['screen_name']]['ngs_data'][guide_id]["region_mw_p"] = 1

					continue

				### Stats calculation
				mean_fold_change = sum(log_fold_changes_sample)/len(log_fold_changes_sample)
				
				try:				
					U1, mw_p = mannwhitneyu(log_fold_changes_sample,log_fold_changes)
					mw_p_log = math.log10(mw_p)
				except:
					logger.error([log_fold_changes_sample,log_fold_changes])
		
				self.options["design"]['region_annotation'][region_id].update({
					"mw_p":mw_p,
					"mw_p_log":mw_p_log,
					"mean_fold_change":mean_fold_change,
					"unfiltered_guides":unfiltered_guides,
					"fold_changes":log_fold_changes_sample,
					"guides":{}
				})
				
				for guide_id in self.options["design"]['region_annotation'][region_id]["guide_ids"]:
					self.options["design"]['region_annotation'][region_id]["guides"][guide_id] = self.data[self.options['screen_name']]['ngs_data'][guide_id]
					self.data[self.options['screen_name']]['ngs_data'][guide_id]["region_mw_p"] = mw_p

		logger.debug("mw_region_statistics completed")

	#################################################################################################	

	def bootstrap_guide_statistics(self):
		logger.debug("bootstrap_guide_statistics")
		try:
			positive_controls = []
			negative_controls = []
			background = []

			for guide_id in self.data[self.options['screen_name']]['ngs_data']:
				if self.check_guide_for_filtering(guide_id): continue
				if 'control' not in self.data[self.options['screen_name']]['ngs_data'][guide_id]: continue
				if self.data[self.options['screen_name']]['ngs_data'][guide_id]['control']:
					if self.data[self.options['screen_name']]['ngs_data'][guide_id]['control_type'] in self.options['positive_control_labels']:
						positive_controls.append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'])
					elif self.data[self.options['screen_name']]['ngs_data'][guide_id]['control_type'] in self.options['negative_control_labels']:
						negative_controls.append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'])
					else:
						logger.error([len(positive_controls),len(negative_controls),self.data[self.options['screen_name']]['ngs_data'][guide_id]['control_type'], self.options['positive_control_labels'],self.options['negative_control_labels']])
				else:
					background.append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'])

			if len(positive_controls) == 0 or len(negative_controls) == 0:
				logger.error("Controls not available - check control labels - Try adding the corect control names")
				sys.exit()

			positive_controls.sort()
			negative_controls.sort()
			background.sort()

			from scipy.stats import bootstrap
			import numpy as np
			import bisect
			
			bootstrap_negative_controls = bootstrap((negative_controls,), np.median, confidence_level=0.99, random_state=1, method='percentile')
			ci_l, ci_u = bootstrap_negative_controls.confidence_interval
			bootstrap_positive_controls = bootstrap((positive_controls,), np.median, confidence_level=0.99, random_state=1, method='percentile')
			
			logger.debug(["bootstrap_negative_controls:",bootstrap_negative_controls])
			logger.debug(["bootstrap_negative_controls:",ci_l, ci_u])
			logger.debug(["bootstrap_positive_controls:",bootstrap_positive_controls])

			positive_controls_means_by_len = {}
			negative_controls_means_by_len = {}
			background_means_by_len = {}

			logger.debug("bootstrap_guide_statistics: annotate_regions")
			for region_id in self.options["design"]['region_annotation']:
				try:
					positive_test = []

					self.options["design"]['region_annotation'][region_id]["positive_controls_bootstraps"] = -1
					self.options["design"]['region_annotation'][region_id]["negative_controls_bootstraps"] = -1
					self.options["design"]['region_annotation'][region_id]["background_bootstraps"] = -1
				
					if "guide_ids" in self.options["design"]['region_annotation'][region_id]:
						for guide_id in self.options["design"]['region_annotation'][region_id]["guide_ids"]:
							if self.check_guide_for_filtering(guide_id): continue
							try:
								positive_test.append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'])
							except:
								pass
					else:
						continue
						
					if len(positive_test) == 0:
						continue

					vector_length = len(positive_test)
					positive_test_mean = sum(positive_test)/len(positive_test)

					positive_controls_means = []
					negative_controls_means = []
					background_means = []
					
					if vector_length not in positive_controls_means_by_len:
						for i in range(0,self.options['bootstrap_size']):
							positive_controls_mean = sum(random.choices(positive_controls, k=vector_length))/vector_length
							positive_controls_means.append(positive_controls_mean)
							negative_controls_mean = sum(random.choices(negative_controls, k=vector_length))/vector_length
							negative_controls_means.append(negative_controls_mean)
							backgrounds_mean = sum(random.choices(background, k=vector_length))/vector_length
							background_means.append(backgrounds_mean)

						positive_controls_means_by_len[vector_length] = positive_controls_means
						negative_controls_means_by_len[vector_length] = negative_controls_means
						background_means_by_len[vector_length] = background_means
					else:
						positive_controls_means = positive_controls_means_by_len[vector_length] 
						negative_controls_means = negative_controls_means_by_len[vector_length]
						background_means = background_means_by_len[vector_length]
					
					positive_controls_means.sort()
					negative_controls_means.sort()
					background_means.sort()
				
					try:
						self.options["design"]['region_annotation'][region_id]["positive_controls_bootstraps"] = bisect.bisect_left(positive_controls_means,positive_test_mean)/self.options['bootstrap_size']
						self.options["design"]['region_annotation'][region_id]["negative_controls_bootstraps"] = bisect.bisect_left(negative_controls_means,positive_test_mean)/self.options['bootstrap_size']
						self.options["design"]['region_annotation'][region_id]["background_bootstraps"] = bisect.bisect_left(background_means,positive_test_mean)/self.options['bootstrap_size']
					except:
						utilities_error.printError()
				except:
					utilities_error.printError()

			logger.debug("bootstrap_guide_statistics: ngs_data")
			import statistics

			negative_controls_stdev = statistics.stdev(negative_controls)
			negative_controls_mean = sum(negative_controls)/len(negative_controls)
			
			for guide_id in self.data[self.options['screen_name']]['ngs_data']:
				try:
					if self.check_guide_for_filtering(guide_id): continue
					self.data[self.options['screen_name']]['ngs_data'][guide_id]['negative_controls_rank'] = bisect.bisect_left(negative_controls,self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'])/len(negative_controls)
					self.data[self.options['screen_name']]['ngs_data'][guide_id]['positive_controls_rank'] = bisect.bisect_left(positive_controls,self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'])/len(positive_controls)
					self.data[self.options['screen_name']]['ngs_data'][guide_id]['background_rank'] = bisect.bisect_left(background,self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'])/len(background)
					self.data[self.options['screen_name']]['ngs_data'][guide_id]['z_score'] = (self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'] - negative_controls_mean)/negative_controls_stdev
					self.data[self.options['screen_name']]['ngs_data'][guide_id]['below_ci_l'] = str(self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'] < ci_l)
					self.data[self.options['screen_name']]['ngs_data'][guide_id]['above_ci_u'] = str(self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores']['log_fold_change'] > ci_u)
					self.data[self.options['screen_name']]['ngs_data'][guide_id]['outside_ci'] = str(self.data[self.options['screen_name']]['ngs_data'][guide_id]['below_ci_l'] or self.data[self.options['screen_name']]['ngs_data'][guide_id]['above_ci_u'])
				except Exception as e:
					logger.debug(["Error",guide_id,e])	
				
			logger.debug("bootstrap_guide_statistics completed")
		except:
			logger.error("bootstrap_guide_statistics failed")
			utilities_error.printError()
			raise

	##################################################################################################
				
	def limma_guide_statistics(self):

		logger.debug("limma_guide_statistics initialised")

		limma_data = {}
		wt_timepoint = self.options['compare_time_point']
		mut_timepoint = self.options['query_time_point']
		
		##-------------------------------------------------------------------------##

		for guide_id in self.data[self.options['screen_name']]['ngs_data']:
			try:
				limma_data[guide_id] = {
					'guide_id':guide_id,
					'fold_change':None,
					'mutant':[],
					'wildtype':[],
					'wildtype':[],
					'wildtype_counts':[],
					'mutant_counts':[],
				}

				for replicate in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][wt_timepoint]:
					limma_data[guide_id]['wildtype'].append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['data_log2_cpm_normalised'][wt_timepoint][replicate])
					limma_data[guide_id]['wildtype_counts'].append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][wt_timepoint][replicate])
					
				for replicate in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][mut_timepoint]:
					limma_data[guide_id]['mutant'].append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['data_log2_cpm_normalised'][mut_timepoint][replicate])
					limma_data[guide_id]['mutant_counts'].append(self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][mut_timepoint][replicate])

			except Exception as e:
				logger.debug(e)	

		self.data[self.options['screen_name']]['limma_scores'] = limmaHelper.limma(limma_data,r_lib_loc=self.options["r_lib_loc"],force=False)

		for guide_id in self.data[self.options['screen_name']]['limma_scores']:
			self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma_scores'] = self.data[self.options['screen_name']]['limma_scores'][guide_id]

		logger.debug("limma_guide_statistics completed")

	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
		

	def update_obsolete_accesion(self,accession,uniprot_check):
		return [accession,uniprot_check]
		
	def format_guides(self,missing_oligos=[]):
		protein_guides = {}
		oligo_list = []
		
		for guide_sequence in list(self.data[self.options['screen_name']]['ngs_data'].keys()):
			try:
				if len(missing_oligos) > 0:
					if guide_sequence not in missing_oligos: continue

				if self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['control']: continue
				if self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession'] not in protein_guides:
					protein_guides[self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession']] = []
				
				protein_guides[self.data[self.options['screen_name']]['ngs_data'][guide_sequence]['accession']].append(guide_sequence)
				oligo_list.append(guide_sequence)
			except:
				logger.error(guide_sequence)
		
		logger.debug("Oligo count " + str(len(oligo_list)) + " " + str(len(list(self.data[self.options['screen_name']]['ngs_data'].keys()))))
		
		return {
			"protein_guides":protein_guides,
			"oligo_list":oligo_list
		}

	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
		
	def annotate_guides(self,retry=False):
		missing_oligos = []
		formatted_guides = self.format_guides(missing_oligos=missing_oligos)
		protein_guides = formatted_guides['protein_guides']
		oligo_list = formatted_guides['oligo_list']
		
		logger.debug("Oligo count " + str(len(oligo_list)) + " " + str(len(list(self.data[self.options['screen_name']]['ngs_data'].keys()))))
		
		if len(oligo_list) == 0: return

		self.guideAnnotationObj = guideAnnotation.guideAnnotation()

		if self.options['editor_type'] == "ABE":
			self.guideAnnotationObj.options['behive_base_editor'] = 'ABE8'
			self.guideAnnotationObj.options['behive_celltype'] = 'mES'

		if self.options['editor_type'] == "CBE":
			self.guideAnnotationObj.options['behive_base_editor'] = 'BE4'	
			self.guideAnnotationObj.options['behive_celltype'] = 'mES'

		if self.options['editor_type'] == "Dual":
			self.guideAnnotationObj.options['behive_base_editor'] = 'ABE8'	
			self.guideAnnotationObj.options['behive_celltype'] = 'mES'

		self.guideAnnotationObj.options['editors'] = [self.options['editor_type']]
		self.guideAnnotationObj.options['protein_gRNAs'] = protein_guides
		self.guideAnnotationObj.options['species'] = "9606"
		self.guideAnnotationObj.options['protein_centric'] = False
		
		###

		self.annotate_guide_edit_mapping(oligo_list)
		edit_oligo_list = list(self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]])
		
		missing_guides = list(set(oligo_list).difference(edit_oligo_list))
		logger.error(str(len(missing_guides)) + " Missing Guides")
		logger.error(missing_guides)
		
		self.annotate_guide_off_target(edit_oligo_list)
		self.annotate_guide_off_target_coordinates(edit_oligo_list)
		self.annotate_guide_efficiency(edit_oligo_list)
		self.annotate_guide_accessiblity(edit_oligo_list)
		self.annotate_guide_snps(edit_oligo_list)
		self.annotate_guide_features(edit_oligo_list)
		self.annotate_guide_motifs(edit_oligo_list)
		self.annotate_guide_conservation(edit_oligo_list)
		self.process_guide_annotation(edit_oligo_list)
		
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###

	def annotate_guide_edit_mapping(self,oligo_list,retry=False):
		tags = self.options["editor_type"]
		guide_annotation_json_path = os.path.join(self.options['data_directory'],"cache", self.options['screen_name'] + "." + tags + ".guide_edit_mapping.json")
		
		if not os.path.exists(guide_annotation_json_path) :
			logger.debug("Writing: " + guide_annotation_json_path)
			self.data[self.options['screen_name']]['guide_annotation'] = self.guideAnnotationObj.annotate_protein_guides()
			utilities_basic.write_to_json(guide_annotation_json_path,self.data[self.options['screen_name']]['guide_annotation'])
		else:
			self.data[self.options['screen_name']]['guide_annotation'] = utilities_basic.read_from_json(guide_annotation_json_path)
			logger.debug("Reading: " + guide_annotation_json_path)

		if self.options["editor_type"] in self.data[self.options['screen_name']]['guide_annotation']:
			logger.debug([guide_annotation_json_path,len(self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]]), len(oligo_list)])
			if len(self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]]) != len(oligo_list):
				logger.error("Guides missing - " + str(len(self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]]))  + " != " + str(len(oligo_list)) + " - " + guide_annotation_json_path)
				logger.error("Guides missing - " + guide_annotation_json_path)
				
	##-----------------------##

	def annotate_guide_off_target_coordinates(self,oligo_list,retry=False):
		guide_off_target_json_path = os.path.join(self.options['data_directory'],"cache", self.options['screen_name'] + ".guide_off_target_coordinates.json")
		
		if not os.path.exists(guide_off_target_json_path):
			logger.debug("Writing: " + guide_off_target_json_path)
			self.data[self.options['screen_name']]['guide_off_target_coordinates'] = queryRunner.queryRunner("ensembl","find_sequence_human_genome",{"oligo_sequence":oligo_list,"counts_only":True,"counts_output_type":"coordinates","remake":False}).run()
			if 'data' not in self.data[self.options['screen_name']]['guide_off_target_coordinates']:
				pprint.pprint(self.data[self.options['screen_name']]['guide_off_target_coordinates'])

			self.data[self.options['screen_name']]['guide_off_target_coordinates'] = self.data[self.options['screen_name']]['guide_off_target_coordinates']['data']
			utilities_basic.write_to_json(guide_off_target_json_path,self.data[self.options['screen_name']]['guide_off_target_coordinates'])
			
		else:
			logger.debug("Reading: " + guide_off_target_json_path)
			self.data[self.options['screen_name']]['guide_off_target_coordinates'] = utilities_basic.read_from_json(guide_off_target_json_path)
			
		if len(self.data[self.options['screen_name']]['guide_off_target_coordinates']) != len(oligo_list):
			logger.error([guide_off_target_json_path,len(self.data[self.options['screen_name']]['guide_off_target_coordinates']),len(oligo_list)])
			logger.error("Guides missing - " + guide_off_target_json_path)

			if retry == False: 
				pass
			else:
				logger.error("Guides missing - " + guide_off_target_json_path)
				logger.error("Retry failed")
				
	##-----------------------##
				
	def annotate_guide_off_target(self,oligo_list,retry=False):
		guide_off_target_json_path = os.path.join(self.options['data_directory'],"cache", self.options['screen_name'] + ".guide_off_target.json")
		guide_off_target_no_complexity_filtering_json_path = os.path.join(self.options['data_directory'],"cache", self.options['screen_name'] + ".no_complexity_filtering.guide_off_target.json")
		
		if not os.path.exists(guide_off_target_json_path):
			logger.debug("Writing: " + guide_off_target_json_path)
			self.data[self.options['screen_name']]['guide_off_target'] = queryRunner.queryRunner("ensembl","find_sequence_human_genome",{"oligo_sequence":oligo_list,"counts_only":True,"remake":False}).run()
			if 'data' not in self.data[self.options['screen_name']]['guide_off_target']:
				pprint.pprint(self.data[self.options['screen_name']]['guide_off_target'])

			self.data[self.options['screen_name']]['guide_off_target'] = self.data[self.options['screen_name']]['guide_off_target']['data']
			utilities_basic.write_to_json(guide_off_target_json_path,self.data[self.options['screen_name']]['guide_off_target'])
			
		else:
			logger.debug("Reading: " + guide_off_target_json_path)
			self.data[self.options['screen_name']]['guide_off_target'] = utilities_basic.read_from_json(guide_off_target_json_path)
			
		if len(self.data[self.options['screen_name']]['guide_off_target']) != len(oligo_list):
			logger.error([guide_off_target_json_path,len(self.data[self.options['screen_name']]['guide_off_target']),len(oligo_list)])
			logger.error("Guides missing - " + guide_off_target_json_path)

			if not os.path.exists(guide_off_target_no_complexity_filtering_json_path):
				logger.debug("Writing off target with no complextity filtering: " + guide_off_target_no_complexity_filtering_json_path)
				missing_oligo_list = list(set(oligo_list).difference(self.data[self.options['screen_name']]['guide_off_target'].keys()))
				guide_off_target_no_complexity_filtering = queryRunner.queryRunner("ensembl","find_sequence_human_genome",{"oligo_sequence":missing_oligo_list,"counts_only":True,"use_filtering":False,"chunks_size":10,"remake":True}).run()
				
				self.data[self.options['screen_name']]['guide_off_target'].update(guide_off_target_no_complexity_filtering['data'])
				utilities_basic.write_to_json(guide_off_target_no_complexity_filtering_json_path, guide_off_target_no_complexity_filtering)
			else:
				logger.debug("Reading: " + guide_off_target_no_complexity_filtering_json_path)
				guide_off_target_no_complexity_filtering = utilities_basic.read_from_json(guide_off_target_no_complexity_filtering_json_path)
				self.data[self.options['screen_name']]['guide_off_target'].update(guide_off_target_no_complexity_filtering['data'])
			
			if retry == False: 
				pass
			else:
				logger.error("Guides missing - " + guide_off_target_json_path)
				logger.error("Retry failed")

	##-----------------------##
		
	def annotate_guide_efficiency(self,oligo_list,retry=False):
		tags = self.guideAnnotationObj.options['behive_base_editor'] + "." + self.guideAnnotationObj.options['behive_celltype']
		guide_efficiency_json_path = os.path.join(self.options['data_directory'],"cache", self.options['screen_name'] + "." + tags + ".guide_efficiency.json")
		
		if not os.path.exists(guide_efficiency_json_path):
			logger.debug("Writing: " + guide_efficiency_json_path)
			self.data[self.options['screen_name']]['guide_efficiency'] = self.guideAnnotationObj.annotate_guides_efficiency()
			utilities_basic.write_to_json(guide_efficiency_json_path,self.data[self.options['screen_name']]['guide_efficiency'])
		else:
			logger.debug("Reading: " + guide_efficiency_json_path)
			self.data[self.options['screen_name']]['guide_efficiency'] = utilities_basic.read_from_json(guide_efficiency_json_path)
			
		if len(self.data[self.options['screen_name']]['guide_efficiency']) != len(oligo_list):
			logger.error([guide_efficiency_json_path,len(self.data[self.options['screen_name']]['guide_efficiency']), len(oligo_list)])
			logger.error("Guides missing - " + guide_efficiency_json_path)
			
			if retry == False: 
				pass
			else:
				logger.error("Guides missing - " + guide_efficiency_json_path)
				logger.error("Retry failed")
	
	##-----------------------##
		
	def annotate_guide_accessiblity(self,oligo_list,retry=False):
		tags = self.options["editor_type"]
		guide_accessibility_json_path = os.path.join(self.options['data_directory'],"cache", self.options['screen_name'] + "." + tags + ".guide_accessibility.json")

		if (not os.path.exists(guide_accessibility_json_path)):
			logger.debug("Writing: " + guide_accessibility_json_path)
			self.data[self.options['screen_name']]['guide_accessibility'] = self.guideAnnotationObj.annotate_guide_accessiblity(annotated_gRNAs=self.data[self.options['screen_name']]['guide_annotation'])
			utilities_basic.write_to_json(guide_accessibility_json_path,self.data[self.options['screen_name']]['guide_accessibility'])
		else:
			logger.debug("Reading: " + guide_accessibility_json_path)
			self.data[self.options['screen_name']]['guide_accessibility'] = utilities_basic.read_from_json(guide_accessibility_json_path)
			
		if len(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]]) != len(oligo_list):
			logger.error([guide_accessibility_json_path,len(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]]), len(oligo_list)])
			logger.error("Guides missing - " + guide_accessibility_json_path)
			if retry == False: 
				pass
			else:
				logger.error("Guides missing - " + guide_accessibility_json_path)
				logger.error("Retry failed")

	##-----------------------##
				
	def annotate_guide_conservation(self,oligo_list,retry=False):
		tags = self.options["editor_type"]
		guide_conservation_json_path = os.path.join(self.options['data_directory'],"cache", self.options['screen_name'] + "." + tags + ".guide_conservation.json")
		
		if (not os.path.exists(guide_conservation_json_path)):
			logger.debug("Writing: " + guide_conservation_json_path)
			self.data[self.options['screen_name']]['guide_conservation'] = self.guideAnnotationObj.annotate_guide_conservation(annotated_gRNAs=self.data[self.options['screen_name']]['guide_annotation'])
			utilities_basic.write_to_json(guide_conservation_json_path,self.data[self.options['screen_name']]['guide_conservation'])
		else:
			self.data[self.options['screen_name']]['guide_conservation'] = utilities_basic.read_from_json(guide_conservation_json_path)
			logger.debug("Reading: " + guide_conservation_json_path)

		logger.debug([guide_conservation_json_path,len(self.data[self.options['screen_name']]['guide_conservation'][self.options["editor_type"]]), len(oligo_list)])
			
	##-----------------------##
 
	def annotate_guide_features(self,oligo_list,retry=False):
		tags = self.options["editor_type"]
		guide_features_json_path = os.path.join(self.options['data_directory'],"cache", self.options['screen_name'] + "." + tags + ".guide_features.json")
		
		if (not os.path.exists(guide_features_json_path)):
			logger.debug("Writing: " + guide_features_json_path)
			self.data[self.options['screen_name']]['guide_features'] = self.guideAnnotationObj.annotate_guide_features(annotated_gRNAs=self.data[self.options['screen_name']]['guide_annotation'])
			utilities_basic.write_to_json(guide_features_json_path,self.data[self.options['screen_name']]['guide_features'])
		else:
			self.data[self.options['screen_name']]['guide_features'] = utilities_basic.read_from_json(guide_features_json_path)
			logger.debug("Reading: " + guide_features_json_path)

		logger.debug([guide_features_json_path,len(self.data[self.options['screen_name']]['guide_features']), len(oligo_list)])
		
	##-----------------------##

	def annotate_guide_motifs(self,oligo_list,retry=False):
		tags = self.options["editor_type"]
		guide_motifs_json_path = os.path.join(self.options['data_directory'],"cache", self.options['screen_name'] + "." + tags + ".guide_motifs.json")
		
		if (not os.path.exists(guide_motifs_json_path)):
			logger.debug("Writing: " + guide_motifs_json_path)
			self.data[self.options['screen_name']]['guide_motifs'] = self.guideAnnotationObj.annotate_guide_motifs(annotated_gRNAs=self.data[self.options['screen_name']]['guide_annotation'])
			utilities_basic.write_to_json(guide_motifs_json_path,self.data[self.options['screen_name']]['guide_motifs'])
		else:
			self.data[self.options['screen_name']]['guide_motifs'] = utilities_basic.read_from_json(guide_motifs_json_path)
			logger.debug("Reading: " + guide_motifs_json_path)

		logger.debug([guide_motifs_json_path,len(self.data[self.options['screen_name']]['guide_motifs']), len(oligo_list)])
		
	##-----------------------##
				
	def annotate_guide_snps(self,oligo_list,retry=False):
		tags = self.options["editor_type"]
		guide_snps_json_path = os.path.join(self.options['data_directory'],"cache", self.options['screen_name'] + "." + tags + ".guide_snps.json")
		
		if (not os.path.exists(guide_snps_json_path)):
			logger.debug("Writing: " + guide_snps_json_path)
			self.data[self.options['screen_name']]['guide_snps'] = self.guideAnnotationObj.annotate_guide_snps(annotated_gRNAs=self.data[self.options['screen_name']]['guide_annotation'])
			utilities_basic.write_to_json(guide_snps_json_path,self.data[self.options['screen_name']]['guide_snps'])
		else:
			self.data[self.options['screen_name']]['guide_snps'] = utilities_basic.read_from_json(guide_snps_json_path)
			logger.debug("Reading: " + guide_snps_json_path)

		logger.debug([guide_snps_json_path,len(self.data[self.options['screen_name']]['guide_snps'][self.options["editor_type"]]), len(oligo_list)])
		
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
		
	def process_guide_annotation(self,oligo_list):
		for guide_id in oligo_list:
			if guide_id not in self.data[self.options['screen_name']]['ngs_data']: continue

			#-----#
   
			structural_offsets = []
			structural_module_class = []
			structural_topology = []

			alphafold_accessibility = []
			alphafold_disorder = []
			alphafold_pLDDT = []

			structural_tessellation_accessibility = []
			structural_tessellation_accessibility_binary = []
			structural_tessellation_accessibility_any_atom = []
			structural_tessellation_accessibility_sidechain_atom = []
			structural_tessellation_accessibility_type = []

			if guide_id in self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]]:
				for offset in self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]][guide_id]:
					structural_offsets.append(str(offset))
					alphafold_accessibility.append(str(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]][guide_id][offset]['alphafold_accessibility']))
					alphafold_disorder.append(str(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]][guide_id][offset]['alphafold_disorder']))
					alphafold_pLDDT.append(str(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]][guide_id][offset]['alphafold_pLDDT']))
					structural_module_class.append(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]][guide_id][offset]['structural_module_class'])
					structural_tessellation_accessibility.append(str(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]][guide_id][offset]['alphafold_accessibility']))
					structural_tessellation_accessibility_binary.append(str(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]][guide_id][offset]['tessellation_accessibility']))
					structural_tessellation_accessibility_any_atom.append(str(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]][guide_id][offset]['tessellation_accessibility_any_atom']))
					structural_tessellation_accessibility_sidechain_atom.append(str(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]][guide_id][offset]['tessellation_accessibility_sidechain_atom']))
					structural_tessellation_accessibility_type.append(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]][guide_id][offset]['tessellation_accessibility_type'])
					structural_topology.append(self.data[self.options['screen_name']]['guide_accessibility'][self.options["editor_type"]][guide_id][offset]['topology'])
			else:
				if not self.data[self.options['screen_name']]['ngs_data'][guide_id]['control']:
					logger.error(guide_id + " not in guide_accessibility")			
			#-----#
					
			conservation_raw = []
			conservation_rlc = []
		
			if guide_id in self.data[self.options['screen_name']]['guide_conservation'][self.options["editor_type"]]:
				for offset in self.data[self.options['screen_name']]['guide_conservation'][self.options["editor_type"]][guide_id]:
					structural_offsets.append(str(offset))
					conservation_raw.append("%1.3f"%(self.data[self.options['screen_name']]['guide_conservation'][self.options["editor_type"]][guide_id][offset]['conservation_raw']))
					conservation_rlc.append("%1.3f"%(self.data[self.options['screen_name']]['guide_conservation'][self.options["editor_type"]][guide_id][offset]['conservation_rlc']))
			else:
				if not self.data[self.options['screen_name']]['ngs_data'][guide_id]['control']:
					logger.error(guide_id + " not in guide_accessibility")			
			
			#-----#
					
			strand =  ""
			peptide_start = ""
			peptide_end = ""
			protein_changes =  ""
			nucleotide_edits =  ""
			genome_edited =  ""
			proteome_edited =  ""
			splice_site =  ""
			splice_site_info = ""
			
			if guide_id in self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]]:
				if len(self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]][guide_id]['difference']) > 0:
					offsets = [int(x) for x in self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]][guide_id]['difference'].keys()]

					if len(offsets) > 0:
						peptide_start = min(offsets)
						peptide_end = max(offsets)

				strand = self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]][guide_id]['strand']
				protein_changes = ",".join(self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]][guide_id]['changes'])
				nucleotide_edits = str(len(self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]][guide_id]['edits']))
				genome_edited = str(self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]][guide_id]['genome_edited'])
				proteome_edited = str(self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]][guide_id]['proteome_edited'])
				splice_site = str(self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]][guide_id]['splice_site'])
			
				if self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]][guide_id]['splice_site']:
					splice_site_info = "exon:" + str(self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]][guide_id]['splice_site_info']['boundary_counter']) + ":" + self.data[self.options['screen_name']]['guide_annotation'][self.options["editor_type"]][guide_id]['splice_site_info']['splice_site_type']
				else:
					splice_site_info = ""

			else:
				if not self.data[self.options['screen_name']]['ngs_data'][guide_id]['control']:
					logger.error(guide_id + " not in guide_annotation")			

			#-----#
					
			guide_off_target_exact = ""
			guide_off_target_mismatches = ""
			
			if guide_id in self.data[self.options['screen_name']]['guide_off_target']:
				guide_off_target_exact = str(self.data[self.options['screen_name']]['guide_off_target'][guide_id]["0"] if "0" in self.data[self.options['screen_name']]['guide_off_target'][guide_id] else 0)
				guide_off_target_mismatches = "|".join([str(i) + ":" + str(self.data[self.options['screen_name']]['guide_off_target'][guide_id][str(i)]) if str(i) in self.data[self.options['screen_name']]['guide_off_target'][guide_id] else str(i) + ":" + "0" for i in range(0,4)])
			else:
				if not self.data[self.options['screen_name']]['ngs_data'][guide_id]['control']:
					logger.error(guide_id + " not in guide_off_target")			

			#-----#

			behive_most_likely_edit = ""
			behive_most_likely_edit_efficiency = ""
			behive_most_likely_single_edit = ""
			behive_peptide_start = ""
			behive_peptide_end = ""
			behive_logit_score = ""
			behive_edits = []
				
			if guide_id in self.data[self.options['screen_name']]['guide_efficiency']:
				offsets = []

				behive_logit_score = self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['behive_logit_score']
				for edit in self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['gRNA_edits_pooled_efficiency_by_aa']:
					if len(edit) > 0:
						offsets.append(int(edit[2:].split("-")[0]))
				
				if len(offsets) > 0:
					behive_peptide_start = min(offsets)
					behive_peptide_end = max(offsets)

				for edit in self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['gRNA_edits_pooled_efficiency_by_aa']:
					behive_edit = self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['gRNA_edits_pooled_efficiency_by_aa'][edit]

					if edit == "": edit = ""
					behive_edits.append(edit + ":" + "%1.3f"%behive_edit['behive_predicted_frequency_aggregated'])
				
				edit = self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['maximum_frequency']['all']["edit"] 
				if edit == "": edit = ""
				behive_most_likely_edit = edit
				behive_most_likely_edit_efficiency = "%1.3f"%self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['maximum_frequency']['all']["frequency"]

				edit = self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['maximum_frequency']['single']["edit"] 
				if edit == "": edit = ""
				behive_most_likely_single_edit = edit + ":" + "%1.3f"%self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['maximum_frequency']['single']["frequency"]
			else:
				if not self.data[self.options['screen_name']]['ngs_data'][guide_id]['control']:
					pass#logger.error(guide_id + " not in guide_efficiency")

			#-----#
			
			guide_snps = []
			if len(self.data[self.options['screen_name']]['guide_snps'][self.options["editor_type"]]) > 0:
				if guide_id in self.data[self.options['screen_name']]['guide_snps'][self.options["editor_type"]]:
					for offset in self.data[self.options['screen_name']]['guide_snps'][self.options["editor_type"]][guide_id]:
						for snp in self.data[self.options['screen_name']]['guide_snps'][self.options["editor_type"]][guide_id][offset]:
							row = [
								snp['description']['wildtype'] + "->" + snp['description']['mutated'] + "@" + str(offset),
								snp['description']['type'],
								",".join(snp['description']['clinical_significance']),
								",".join(snp['description']['disease']) if 'disease' in snp['description'] else '-',
								",".join(snp['description']['xrefs']),
								",".join([str(x) for x in snp['description']['pmids']]) if 'pmids' in snp['description'] else '-'
							]
							guide_snps.append("|".join(row))

			#-----#
   
			guide_features = {
				"guide_feature_domain":[],
				"guide_feature_structure":[],
				"guide_feature_isoform":[],
				"guide_feature_modification":[],
				"guide_feature_region":[],
				"guide_feature_topology":[],
				"guide_feature_mutagenesis":[],
				"guide_feature_other":[]
			}
			
			if len(self.data[self.options['screen_name']]['guide_features'][self.options["editor_type"]]) > 0:
				if guide_id in self.data[self.options['screen_name']]['guide_features'][self.options["editor_type"]]:
					for offset in self.data[self.options['screen_name']]['guide_features'][self.options["editor_type"]][guide_id]:
						for feature in self.data[self.options['screen_name']]['guide_features'][self.options["editor_type"]][guide_id][offset]:
							try:	
								if feature['type'] in ['snp','motif','secondary_structure']: continue
								description = ""
								if feature['type'] == 'domain':
									description = [feature['description']['domain_name'],feature['description']['name'],feature['source_id'],feature['external_id']]
								elif feature['type'] == 'modification':
									if 'description' in feature['description']: 
										description = [feature['description']['name'] + " - " + feature['description']['description'],"source:" + feature['source_id']]
									else:
										description = [feature['description']['name'],"source:" + feature['source_id']]
								elif feature['type'] == 'isoform':
									if 'wildtype' in feature:
										description = [feature['description']['description'],feature['wildtype'],feature['mutated']]
									else:
										description = [feature['description']['description'],"-","-"]
								elif feature['type'] == 'structure':
									description = [feature['source_id'] + ":" + feature['external_id']]
								elif feature['type'] == 'switch':
									feature['type']  = "other"
									description = ["switch" + ":" + feature['description']['type'] + ":" + feature['description']['subtype']]
								else:
									description = [feature['description']['description']]
									
								if "guide_feature_" + feature['type'] not in guide_features:
									logger.error("guide_feature_" + feature['type'])
								else:
									try:
										row = [
											str(feature['start']),
											str(feature['stop']),
										] + [str(x) for x in description] + [
											"PMID:" + ",".join([str(pmid) for pmid in feature['description']['pmids']]) if 'pmids' in feature['description'] else "-"
										]
										
										guide_features["guide_feature_" + feature['type']].append("|".join(row))
									except:
										logger.error(row)
							except:
								logger.error(feature)
			#-----#   

			guide_motifs = []
			if len(self.data[self.options['screen_name']]['guide_motifs'][self.options["editor_type"]]) > 0:
				if guide_id in self.data[self.options['screen_name']]['guide_motifs'][self.options["editor_type"]]:
					for motif in self.data[self.options['screen_name']]['guide_motifs'][self.options["editor_type"]][guide_id]:
						row = [
							str(motif['motif_start_residue_number']),
							str(motif['motif_end_residue_number']),
							",".join(motif['interaction_types']),
							motif['motif_sequence'],
							",".join(motif['binding_proteins']),
							",".join(motif['specificity_classes']),
							",".join(motif['pmid'])
						]

						guide_motifs.append("|".join(row))
		
			#-----#
   
			guide_annotation = {}
			guide_annotation['strand'] = strand
			guide_annotation['peptide_start'] = peptide_start
			guide_annotation['peptide_end'] = peptide_end
			guide_annotation['protein_changes'] = protein_changes
			guide_annotation['nucleotide_edits'] = nucleotide_edits
			guide_annotation['genome_edited'] = genome_edited
			guide_annotation['proteome_edited'] = proteome_edited


			guide_annotation['behive_logit_score'] = behive_logit_score
			guide_annotation['behive_most_likely_edit'] = behive_most_likely_edit
			guide_annotation['behive_most_likely_edit_match'] = str(protein_changes==behive_most_likely_edit)
			guide_annotation['behive_most_likely_edit_efficiency'] = behive_most_likely_edit_efficiency
			guide_annotation['behive_most_likely_single_edit'] = behive_most_likely_single_edit
			guide_annotation['behive_edits'] = "|".join(behive_edits)
			guide_annotation['behive_peptide_start'] = behive_peptide_start
			guide_annotation['behive_peptide_end'] = behive_peptide_end

			guide_annotation['splice_site'] = splice_site
			guide_annotation['splice_site_info'] = splice_site_info
			guide_annotation['guide_off_target_exact'] = guide_off_target_exact
			guide_annotation['guide_off_target_mismatches'] = guide_off_target_mismatches

			guide_annotation['conservation_raw'] = "|".join(conservation_raw)
			guide_annotation['conservation_rlc_p'] = "|".join(conservation_rlc)

			guide_annotation['structural_module_class'] = "|".join(structural_module_class)
			guide_annotation['structural_offsets'] = "|".join(structural_offsets)
			guide_annotation['alphafold_accessibility'] = "|".join(alphafold_accessibility)
			guide_annotation['alphafold_disorder'] = "|".join(alphafold_disorder)
			guide_annotation['alphafold_pLDDT'] = "|".join(alphafold_pLDDT)
			guide_annotation['structural_tessellation_accessibility'] = "|".join(structural_tessellation_accessibility)
			guide_annotation['structural_tessellation_accessibility_any_atom'] = "|".join(structural_tessellation_accessibility_any_atom)
			guide_annotation['structural_tessellation_accessibility_sidechain_atom'] = "|".join(structural_tessellation_accessibility_sidechain_atom)
			guide_annotation['structural_tessellation_accessibility_type'] = "|".join(structural_tessellation_accessibility_type)
			guide_annotation['structural_topology'] = "|".join(structural_topology)

			guide_annotation['guide_snps'] = ";".join(guide_snps)
			guide_annotation['guide_motifs'] = ";".join(guide_motifs)
			guide_annotation['guide_feature_domain'] = ";".join(guide_features['guide_feature_domain'])
			guide_annotation['guide_feature_structure'] = ";".join(guide_features['guide_feature_structure'])
			guide_annotation['guide_feature_isoform'] = ";".join(guide_features['guide_feature_isoform'])
			guide_annotation['guide_feature_modification'] = ";".join(guide_features['guide_feature_modification'])
			guide_annotation['guide_feature_region'] = ";".join(guide_features['guide_feature_region'])
			guide_annotation['guide_feature_topology'] = ";".join(guide_features['guide_feature_topology'])
			guide_annotation['guide_feature_mutagenesis'] = ";".join(guide_features['guide_feature_mutagenesis'])
			guide_annotation['guide_feature_other'] = ";".join(guide_features['guide_feature_other'])
			
			self.data[self.options['screen_name']]['ngs_data'][guide_id].update(guide_annotation)
			
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
		
	def find_consecutive_groups(self,lst):
		consecutive_groups = []
		current_group = []
		lst.sort()

		for i, num in enumerate(lst):
			if i == 0 or num != lst[i - 1] + 1:
				if current_group:
					consecutive_groups.append(current_group)
					current_group = []
			current_group.append(num)

		if current_group:
			consecutive_groups.append(current_group)

		return consecutive_groups
		
	##################################################################################################

	def make_region_sliding_window(self):
		logger.debug("make_region_sliding_window")

		guide_mapping = {}
		for guide_id in self.data[self.options['screen_name']]['ngs_data']:
			try:
				if self.check_guide_for_filtering(guide_id): continue
				if self.data[self.options['screen_name']]['ngs_data'][guide_id]['control']: continue
				accession = self.data[self.options['screen_name']]['ngs_data'][guide_id]['accession']
				
				peptide_start = self.data[self.options['screen_name']]['ngs_data'][guide_id]['behive_peptide_start']
				peptide_end =  self.data[self.options['screen_name']]['ngs_data'][guide_id]['behive_peptide_end']

				if accession not in guide_mapping:
					response = queryRunner.queryRunner("uniprot","parse_basic",{"accession":accession}).run()
					guide_mapping[accession] = {
						"guides":{},
						"gene":self.data[self.options['screen_name']]['ngs_data'][guide_id]['gene'],
						"sequence":response['data']['sequence']
					}
				
				for offset in range(peptide_start,peptide_end+1):
					if offset not in guide_mapping[accession]["guides"]:
						guide_mapping[accession]["guides"][offset] = []
					
					guide_mapping[accession]["guides"][offset].append(guide_id)
			except:
				logger.debug("make_region_sliding_window " + guide_id) 

		##------------------------##

		self.options['region_window'] = 5
		for accession in guide_mapping:
			for offset in guide_mapping[accession]["guides"]:
				region_annotation = {}
				region_annotation['accession'] = accession
				region_annotation['gene'] = guide_mapping[accession]['gene']
				region_annotation['peptide_start'] = max(0,offset-self.options['region_window'])
				region_annotation['peptide_end'] = offset+self.options['region_window']
				
				region_annotation['peptide'] = guide_mapping[accession]['sequence'][region_annotation['peptide_start']:region_annotation['peptide_end'] +1]
				flanking_guide_offsets = [x for x in range(region_annotation['peptide_start'],region_annotation['peptide_end']+1) if x in guide_mapping[accession]["guides"]]

				region_annotation['guide_ids'] = []
				for flanking_guide_offset in flanking_guide_offsets:
					for guide_id in guide_mapping[accession]["guides"][flanking_guide_offset]:
						if guide_id not in region_annotation['guide_ids']:
							region_annotation['guide_ids'].append(guide_id)
				
				region_annotation_id = "_".join([accession,str(region_annotation['peptide_start']),str(region_annotation['peptide_end'])])
				self.options["design"]['region_annotation'][region_annotation_id] = region_annotation
	
	##################################################################################################

	def make_region_tesselation(self):
		logger.debug("make_region_tesselation")

		guide_mapping = {}
		counter = 0
		for guide_id in self.data[self.options['screen_name']]['ngs_data']:
			
				try:
					if self.check_guide_for_filtering(guide_id):
						continue
					if self.data[self.options['screen_name']]['ngs_data'][guide_id]['control']: 
						continue
					if guide_id not in self.data[self.options['screen_name']]['guide_efficiency']:
						continue
				
					accession = self.data[self.options['screen_name']]['ngs_data'][guide_id]['accession']
					offsets = self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['uniprot_offsets']
			
					if accession not in guide_mapping:
						response = queryRunner.queryRunner("uniprot","parse_basic",{"accession":accession}).run()
						guide_mapping[accession] = {
							"guides":{},
							"gene":self.data[self.options['screen_name']]['ngs_data'][guide_id]['gene'],
							"sequence":response['data']['sequence']
						}

					for offset in offsets:
						if offset not in guide_mapping[accession]["guides"]:
							guide_mapping[accession]["guides"][offset] = []
						
						counter += 1
						guide_mapping[accession]["guides"][offset].append(guide_id)
				except:
					logger.error("make_region_3d " + guide_id) 
					raise

		##------------------------##

		for accession in guide_mapping:
			response = queryRunner.queryRunner("alphafold","get_tessellation_accessibility",{"accession":accession}).run()

			if response['status'] == 'Error':
				pprint.pprint(response)
		
			for offset in guide_mapping[accession]["guides"]:
				region_annotation = {}
				region_annotation['accession'] = accession
				region_annotation['gene'] = guide_mapping[accession]['gene']
				region_annotation['peptide_start'] = offset
				region_annotation['peptide_end'] = offset
				region_annotation['direct_neighbors'] = []
				region_annotation.update(response['data']['chain_details']['A']['accessible_residues'][str(offset)])

				flanking_guide_offsets = [str(offset)]
				if str(offset) in response['data']['chain_details']['A']['direct_neighbors']:
					flanking_guide_offsets += response['data']['chain_details']['A']['direct_neighbors'][str(offset)]
					region_annotation['direct_neighbors'] = response['data']['chain_details']['A']['direct_neighbors'][str(offset)]

				region_annotation['peptide'] = guide_mapping[accession]['sequence'][region_annotation['peptide_start']:region_annotation['peptide_end'] +1]
				region_annotation['guide_ids'] = copy.deepcopy(guide_mapping[accession]["guides"][offset])
				region_annotation['guide_offsets'] = {}
			
				for flanking_guide_offset in flanking_guide_offsets:
					if int(flanking_guide_offset) in guide_mapping[accession]["guides"]:
						region_annotation['guide_offsets'][flanking_guide_offset] = {}
						for guide_id in guide_mapping[accession]["guides"][int(flanking_guide_offset)]:
							if guide_id not in region_annotation['guide_ids']:
								region_annotation['guide_ids'].append(guide_id)
							
							region_annotation['guide_offsets'][flanking_guide_offset][guide_id] = self.data[self.options['screen_name']]['guide_efficiency'][guide_id]['gRNA_edits_pooled_efficiency_by_aa']

				region_annotation_id = "_".join([accession,str(region_annotation['peptide_start']),str(region_annotation['peptide_end'])])
				self.options["design"]['region_annotation'][region_annotation_id] = region_annotation


	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
		
	def add_guide_id_to_region_annotation(self):
		logger.debug("add_guide_id_to_region_annotation")
	
		for guide_id in self.data[self.options['screen_name']]['ngs_data']:
			for region_id in self.options["design"]['region_annotation']:
				try:
					if 'accession' not in self.data[self.options['screen_name']]['ngs_data'][guide_id]: continue
					if 'accession' not in self.options["design"]['region_annotation'][region_id]: continue
					
					if self.data[self.options['screen_name']]['ngs_data'][guide_id]['accession'] == self.options["design"]['region_annotation'][region_id]['accession']:
						add_guides = False

						guide_id_regions_offsets = []
						if 'cbe_edit_mutation_details' in  self.data[self.options['screen_name']]['ngs_data'][guide_id]:
							guide_id_regions_offsets += list(self.data[self.options['screen_name']]['ngs_data'][guide_id]['cbe_edit_mutation_details'].keys())
						if 'abe_edit_mutation_details' in  self.data[self.options['screen_name']]['ngs_data'][guide_id]:
							guide_id_regions_offsets += list(self.data[self.options['screen_name']]['ngs_data'][guide_id]['abe_edit_mutation_details'].keys())

						if len(guide_id_regions_offsets) == 0: continue

						guide_id_regions_start = min(guide_id_regions_offsets)
						guide_id_regions_end = max(guide_id_regions_offsets)
						
						if guide_id_regions_start >= int(self.options["design"]['region_annotation'][region_id]['peptide_start']) and guide_id_regions_start <= int(self.options["design"]['region_annotation'][region_id]['peptide_end']):					
							add_guides = True
						if guide_id_regions_end >= int(self.options["design"]['region_annotation'][region_id]['peptide_start']) and guide_id_regions_end <= int(self.options["design"]['region_annotation'][region_id]['peptide_end']):					
							add_guides = True

						if add_guides:
							if region_id in self.options["design"]['region_annotation']:
								if "guide_ids" not in self.options["design"]['region_annotation'][region_id]:
									self.options["design"]['region_annotation'][region_id]["guide_ids"] = []
							
							self.options["design"]['region_annotation'][region_id]["guide_ids"].append(guide_id)

						self.data[self.options['screen_name']]['ngs_data'][guide_id]['region'] = region_id
				except:
					logger.error([region_id,guide_id])
					utilities_error.printError()
		
		logger.debug("Completed add_guide_id_to_region_annotation")

	##################################################################################################

	def load_region_annotation(self):
		if "run_peptools" in self.options:
			if self.options["run_peptools"] == False:
				return 
			
		offsets = []

		for region_id in self.options["design"]['region_annotation']:
			try:
				offsets.append(
					[
						self.options["design"]['region_annotation'][region_id]['accession'],
						str(self.options["design"]['region_annotation'][region_id]['peptide_start']),
						str(self.options["design"]['region_annotation'][region_id]['peptide_end']),
					]
				)
			except:
				logger.error("Annotation failed for " + region_id)
				logger.error(self.options["design"]['region_annotation'][region_id])
		
		offsets.sort()
		
		for i in range(0,len(offsets),10000):
			offsets_formatted = []
			for offset in offsets[i:i+10000]:
				offsets_formatted.append({"accession":offset[0],"start":int(offset[1]),"end":int(offset[2])})
				
			options = {
				"debug":True,
				"toolname":"peptools",
				"iupred_cutoff":0,
				"flank": 10,
				"taxon_id": 9606,
				"binding_domains": [],
				"binding_proteins": [],
				"offset":offsets_formatted,
				"shared_pval_cutoff": 0.0001,
				"cached": True,
				"check_shared_domains_interactions": False,
				"check_shared_proteins_interactions": False,
				"check_shared_proteins": False,
				"check_shared_domains": False,
				"motif":"",
				'database_name':"peptools",
				"is_superuser":True,
			}
			
			response = queryRunner.queryRunner("database_access","run_peptools_offset_annotation",options).run()
			
			if 'data' in response:
				for region_peptools_annotation in response['data']['result']['result']:
					region_id = region_peptools_annotation['ProteinAcc'] + "_" + str(region_peptools_annotation['SeqStart']) + "_" + str(region_peptools_annotation['SeqStop'])

					if region_id not in self.options["design"]['region_annotation']:
						self.options["design"]['region_annotation'][region_id] = {"peptools":region_peptools_annotation}
					else:
						self.options["design"]['region_annotation'][region_id]["peptools"] = region_peptools_annotation
					
	##################################################################################################
		
	def annotate_regions(self):
		logger.debug("annotate_regions: intialised")
		mean_gene_effect_dict = {}
		hap1_gene_effect_dict = {}

		green_scale = ["0A6921","1A8828","429B46","94C58C","D4E8D1"]
		red_scale = ["FF3333","FF7070","FF9999","FFD6D6","FFEBEB"]

		green_scale.reverse()
		red_scale.reverse()

		for region_id in self.options["design"]['region_annotation']:
			if "filter_accession" in self.options:
				if self.options["design"]['region_annotation'][region_id]['accession'] not in self.options["filter_accession"]:
					return

			try:
				if "guide_ids" in self.options["design"]['region_annotation'][region_id]:
					sample_name = []
					snps = []
					modifications = []
					motifs = []
					mutation_track = []

					accessibility = 0.0
					significant_guides = 0
					confident_guides = 0
					guides_count = 0

					for guide_id in self.options["design"]['region_annotation'][region_id]["guide_ids"]:
						try:
							
							mutation_tag = self.get_mutation_tag()

							if self.options["design"]['region_annotation'][region_id]['guides_status'][guide_id] == "":	
								if 'accessibility' in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
									accessibility = self.data[self.options['screen_name']]['ngs_data'][guide_id]['accessibility']
								
								if self.data[self.options['screen_name']]['ngs_data'][guide_id][mutation_tag] != "":
									guides_count += 1

									if 'motif' in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
										for motif_instance in self.data[self.options['screen_name']]['ngs_data'][guide_id]['motif']:
											motifs.append(str(motif_instance['name'].replace("SLiM","Unclassified")) + "@" + str(motif_instance['description']['start']) + ":" +  str(motif_instance['description']['stop']) )

									if 'accessibility' in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
										accessibility = self.data[self.options['screen_name']]['ngs_data'][guide_id]['accessibility']
									
									if "outside_ci" in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
										if self.data[self.options['screen_name']]['ngs_data'][guide_id]['outside_ci'] == "True":
											confident_guides += 1
									
									colour = "AAA"
									if guide_id in self.data[self.options['screen_name']]['limma_scores']:
										if self.data[self.options['screen_name']]['limma_scores'][guide_id]['p_value'] < 0.01:
											significant_guides += 1
								
										if self.data[self.options['screen_name']]['limma_scores'][guide_id]['log_fold_change'] < 0:
											colour = red_scale[int(min(len(red_scale)-1,abs(self.data[self.options['screen_name']]['limma_scores'][guide_id]['log_fold_change'])*2-1))]
										else:
											colour = green_scale[int(min(len(green_scale)-1,abs(self.data[self.options['screen_name']]['limma_scores'][guide_id]['log_fold_change'])*2-1))]

										if self.data[self.options['screen_name']]['limma_scores'][guide_id]['p_value'] > 0.01:
											colour = "AAA"
										
										if self.data[self.options['screen_name']]['limma_scores'][guide_id]['p_value'] < 0.01:
											if 'snps' in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
												for snp in self.data[self.options['screen_name']]['ngs_data'][guide_id]['snps']:
													
													if snp['name'].count("frameshift") == 0 and snp['name'].count("stop") == 0:
														if 'Disease' in snp['description']:
															snps += [snp['name'] + ":" +  snp['description']['mutation'] + ":" +  snp['description']['Disease'] + "@" + str(snp['description']['start'])] 
														else:
															snps += [snp['name'] + ":" +  snp['description']['mutation'] + ":-@" + str(snp['description']['start'])] 

											
											if 'modification' in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
												modifications += [x['name'] + "@" + str(x['description']['start']) for x in self.data[self.options['screen_name']]['ngs_data'][guide_id]['modification']]
										

											sample_name.append(self.data[self.options['screen_name']]['ngs_data'][guide_id][mutation_tag] + ':fc=' + "%1.2g"%self.data[self.options['screen_name']]['limma_scores'][guide_id]['log_fold_change'] + ":p=" + "%1.2g"%self.data[self.options['screen_name']]['limma_scores'][guide_id]['p_value'])
									
									else:
										logger.debug(guide_id + " not in limma scores")
										logger.debug(self.data[self.options['screen_name']]['ngs_data'][guide_id])

									try:
										track_mutation_start = min(list(self.data[self.options['screen_name']]['ngs_data'][guide_id][mutation_tag + '_details'].keys()))
										track_mutation_end = max(list(self.data[self.options['screen_name']]['ngs_data'][guide_id][mutation_tag + '_details'].keys()))
										track_mutations = [value['mut'] for value in self.data[self.options['screen_name']]['ngs_data'][guide_id][mutation_tag + '_details'].values()]
		
										mutation_track.append("peptides,Mutations," + colour + "," + "".join(track_mutations) + "," + str(track_mutation_start) + "," + str(track_mutation_end))
									except:
										logger.error(guide_id + " ProViz")
							else:
								try:
									if self.data[self.options['screen_name']]['limma_scores'][guide_id]['p_value'] < 0.01:
										sample_name.append('OutOfWindow:fc=' + "%1.2g"%self.data[self.options['screen_name']]['limma_scores'][guide_id]['log_fold_change'] + ":p=" + "%1.2g"%self.data[self.options['screen_name']]['limma_scores'][guide_id]['p_value'])
								except:
									pass
						except:
							utilities_error.printError()
				
					motifs = list(set(motifs))
					modifications = list(set(modifications))
					snps = list(set(snps))

					###------------------------------------------###

					mean_gene_effect = -1

					if self.options['add_depmap_data']:
						if self.options["design"]['region_annotation'][region_id]['accession'] not in mean_gene_effect_dict:
							mean_gene_effect = 100
							hap1_gene_effect = {'gene_name': '', 'p_value': '-', 'q_value': '-', 'ratio': '-', 'selected': ''}
							try:
								mean_gene_effect_data = queryRunner.queryRunner("depmap","get_CRISPR_gene_dependency",{"accession":self.options["design"]['region_annotation'][region_id]['accession']}).run()
								mean_gene_effect = mean_gene_effect_data['data']['mean_gene_effect']
							except:
								pass

							try:
								hap1_gene_effect = queryRunner.queryRunner("depmap","get_hap1_gene_dependency_data",{"accession":self.options["design"]['region_annotation'][region_id]['accession']}).run()['data']
								if len(hap1_gene_effect) == 0:
									hap1_gene_effect = {'gene_name': '', 'p_value': '-', 'q_value': '-', 'ratio': '-', 'selected': ''}
							except:
								pass
							
							mean_gene_effect_dict[self.options["design"]['region_annotation'][region_id]['accession']] = mean_gene_effect
							hap1_gene_effect_dict[self.options["design"]['region_annotation'][region_id]['accession']] = hap1_gene_effect
						else:
							mean_gene_effect = mean_gene_effect_dict[self.options["design"]['region_annotation'][region_id]['accession']]
							hap1_gene_effect = hap1_gene_effect_dict[self.options["design"]['region_annotation'][region_id]['accession']]
					
					###------------------------------------------###

					self.options["design"]['region_annotation'][region_id]['accession'] = self.options["design"]['region_annotation'][region_id]['accession']
					self.options["design"]['region_annotation'][region_id]['gene'] = self.options["design"]['region_annotation'][region_id]['gene']
					self.options["design"]['region_annotation'][region_id]['peptide_start'] = self.options["design"]['region_annotation'][region_id]['peptide_start']
					self.options["design"]['region_annotation'][region_id]['peptide_end'] = self.options["design"]['region_annotation'][region_id]['peptide_end']
					self.options["design"]['region_annotation'][region_id]['peptide'] = self.options["design"]['region_annotation'][region_id]['peptide']
					self.options["design"]['region_annotation'][region_id]['mean_gene_effect'] = mean_gene_effect
					self.options["design"]['region_annotation'][region_id]['hap1_gene_effect'] = hap1_gene_effect
					self.options["design"]['region_annotation'][region_id]["accessibility"] = accessibility
					self.options["design"]['region_annotation'][region_id]['sample_name'] = sample_name
					self.options["design"]['region_annotation'][region_id]["motifs"] = motifs
					self.options["design"]['region_annotation'][region_id]["snps"] = snps
					self.options["design"]['region_annotation'][region_id]["modifications"] = modifications
					self.options["design"]['region_annotation'][region_id]['significant_guides'] = significant_guides
					self.options["design"]['region_annotation'][region_id]['confident_guides'] = confident_guides
					self.options["design"]['region_annotation'][region_id]['guides_count'] = guides_count
					self.options["design"]['region_annotation'][region_id]["mutation_track"] = mutation_track
					self.options["design"]['region_annotation'][region_id]["mean_gene_effect"] = mean_gene_effect

					if guides_count > 0:
						self.options["design"]['region_annotation'][region_id]['significant_guides_proportion'] = significant_guides/guides_count
						self.options["design"]['region_annotation'][region_id]['confident_guides_proportion'] = confident_guides/guides_count
					else:
						self.options["design"]['region_annotation'][region_id]['significant_guides_proportion'] = 0
						self.options["design"]['region_annotation'][region_id]['confident_guides_proportion'] = 0

			except:
				logger.error("annotate_regions: " + region_id + " failed")
				utilities_error.printError()
		
		logger.debug("annotate_regions: completed")
	
	##################################################################################################
	
	def load_mutation_annotation(self):

		if "run_peptools" in self.options:
			if self.options["run_peptools"] == False:
				return 
			
		offsets = []
		guide_id_mapping = {}
		
		guide_ids = list(self.data[self.options['screen_name']]['ngs_data'])
		guide_ids.sort()

		for guide_id in guide_ids:
			if 'accession' in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
				accession = self.data[self.options['screen_name']]['ngs_data'][guide_id]['accession']
				if 'abe_edit_mutation' in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
					for mutation in self.data[self.options['screen_name']]['ngs_data'][guide_id]['abe_edit_mutation'].split("__"):
						offsets.append([accession,mutation[1:].split("->")[0]])
						if accession + "_" + mutation[1:].split("->")[0] not in guide_id_mapping:
							guide_id_mapping[accession + "_" + mutation[1:].split("->")[0]] = []

						guide_id_mapping[accession + "_" + mutation[1:].split("->")[0]].append(guide_id)

				if 'cbe_edit_mutation' in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
					for mutation in self.data[self.options['screen_name']]['ngs_data'][guide_id]['cbe_edit_mutation'].split("__"):
						offsets.append([accession,mutation[1:].split("->")[0]])

						if accession + "_" + mutation[1:].split("->")[0] not in guide_id_mapping:
							guide_id_mapping[accession + "_" + mutation[1:].split("->")[0]] = []
							
						guide_id_mapping[accession + "_" + mutation[1:].split("->")[0]].append(guide_id)
		
		for i in range(0,len(offsets),10000):
			offsets_formatted = []
			for offset in offsets[i:i+10000]:
				offsets_formatted.append({"accession":offset[0],"start":int(offset[1]),"end":int(offset[1])})
				
			options = {
			"debug":True,
			"toolname":"peptools",
			"iupred_cutoff":0,
			"flank": 10,
			"taxon_id": 9606,
			"binding_domains": [],
			"binding_proteins": [],
			"offset":offsets_formatted,
			"shared_pval_cutoff": 0.0001,
			"cached": True,
			"check_shared_domains_interactions": False,
			"check_shared_proteins_interactions": False,
			"check_shared_proteins": False,
			"check_shared_domains": False,
			"motif":"",
			'database_name':"peptools",
			"is_superuser":True,
			}
			
			response = queryRunner.queryRunner("database_access","run_peptools_offset_annotation",options).run()
			
			if 'data' in response:
				for mutation in response['data']['result']['result']:
					for guide_id in guide_id_mapping[mutation['ProteinAcc'] + "_" + str(mutation['SeqStart'])]:

						self.data[self.options['screen_name']]['ngs_data'][guide_id]['accessibility'] = mutation['Accessibility']
						self.data[self.options['screen_name']]['ngs_data'][guide_id]['accessibility_raw'] = mutation['Accessibility_raw']
						self.data[self.options['screen_name']]['ngs_data'][guide_id]['conservation_metazoa'] = mutation['Conservation metazoa']['score']

						if 'SNP' in mutation:
							snps = []
							for snp in mutation['SNP']:
								if snp['description']['distance'] == 0:
									snps.append(snp)

							self.data[self.options['screen_name']]['ngs_data'][guide_id]['snps'] = snps
						else:
							self.data[self.options['screen_name']]['ngs_data'][guide_id]['snps'] = []
						
						if 'Motif' in mutation:
							motifs = []
							for motif in mutation['Motif']:
								if motif['description']['distance'] == 0:
									motifs.append(motif)

							self.data[self.options['screen_name']]['ngs_data'][guide_id]['motif'] = motifs
						else:
							self.data[self.options['screen_name']]['ngs_data'][guide_id]['motif'] = []

						if 'Modification' in mutation:	
							modifications = []
							for modification in mutation['Modification']:
								if modification['description']['distance'] == 0:
									modifications.append(modification)
							
							self.data[self.options['screen_name']]['ngs_data'][guide_id]['modification'] = modifications
						else:
							self.data[self.options['screen_name']]['ngs_data'][guide_id]['modification'] = []

			else:
				logger.error(response)

	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------######------------------------######------------------------######------------------------###
		
	def make_web_page_data(self):
		##-------------------------------------------------------------------------##
		##
		## Make the data for the web page
		##
		##-------------------------------------------------------------------------##

		limma_mutation_statistics_data = {
			"points":[],
			"range":{
				"fold_change":[1000,-1000],
				"p_value":[1000,-1000],
			}
		}
		
		for guide_id in self.data[self.options['screen_name']]['limma_scores']:
			if 'abe_edit_mutation' in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
				accessibility = -1
				if 'accessibility' in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
					accessibility = self.data[self.options['screen_name']]['ngs_data'][guide_id]['accessibility']

				self.data[self.options['screen_name']]['ngs_data'][guide_id]['limma'] = self.data[self.options['screen_name']]['limma_scores'][guide_id]

				annotation = {
					"fold_change":self.data[self.options['screen_name']]['limma_scores'][guide_id]['log_fold_change'],
					"p_value":self.data[self.options['screen_name']]['limma_scores'][guide_id]['log_p_value'],
					"accessibility": accessibility,
					"accession":self.data[self.options['screen_name']]['ngs_data'][guide_id]['accession'],
					"gene_name":self.data[self.options['screen_name']]['ngs_data'][guide_id]['gene'] ,
					"peptide_start":self.data[self.options['screen_name']]['ngs_data'][guide_id]['peptide_start'],
					"peptide_end":self.data[self.options['screen_name']]['ngs_data'][guide_id]['peptide_end'],
					'mutation':self.data[self.options['screen_name']]['ngs_data'][guide_id]['abe_edit_mutation']
				}

				if self.data[self.options['screen_name']]['limma_scores'][guide_id]['log_p_value'] > 1:
					limma_mutation_statistics_data['points'].append(annotation)

				if annotation['fold_change'] < limma_mutation_statistics_data["range"]["fold_change"][0]: limma_mutation_statistics_data["range"]["fold_change"][0] = annotation['fold_change']
				if annotation['fold_change']> limma_mutation_statistics_data["range"]["fold_change"][1]: limma_mutation_statistics_data["range"]["fold_change"][1] = annotation['fold_change']
				if annotation['p_value'] < limma_mutation_statistics_data["range"]["p_value"][0]: limma_mutation_statistics_data["range"]["p_value"][0] = annotation['p_value']
				if annotation['p_value'] > limma_mutation_statistics_data["range"]["p_value"][1]: limma_mutation_statistics_data["range"]["p_value"][1] = annotation['p_value']	
		
		json_path = "/home/data/base_editing/"  + self.options['screen_name'] + "." + str(self.options['compare_time_point']) + ".base_editing_mutation_scatterplot.json"
		logger.debug("Writing " + json_path)
		
		with open(json_path,'w') as outfile:
			json.dump(limma_mutation_statistics_data, outfile)

		return {"status":"Success"}
	
	##################################################################################################
				
	def write_guide_annotation_file(self):
		
		p_limma_tag =  [
			"log_p_value",
			"log_p_value_adj",
			"log_fold_change"
		]

		tags = [
			'seq',
			'strand',
			'gene',
			'accession',
			'essentiality_hap1',
			'essentiality_depmap',
			'peptide_start',
			'peptide_end',
			'splice_site',
			'splice_site_info',
			"control",
			"control_type",
			'protein_changes',
			'nucleotide_edits',
			'genome_edited',
			'proteome_edited',
			'behive_peptide_start',
			'behive_peptide_end',
			'behive_logit_score',
			'behive_most_likely_edit',
			'behive_most_likely_edit_efficiency',
			'behive_most_likely_edit_match',
			'behive_most_likely_single_edit',
			'behive_edits',
			'guide_off_target_exact',
			'guide_off_target_mismatches',
			'conservation_raw',
			'conservation_rlc_p',
			'alphafold_accessibility',
			'alphafold_disorder',
			'alphafold_pLDDT',
			'structural_module_class',
			'structural_offsets',
			'structural_tessellation_accessibility',
			'structural_tessellation_accessibility_any_atom',
			'structural_tessellation_accessibility_sidechain_atom',
			'structural_tessellation_accessibility_type',
			'structural_topology',
 			#'cbe_edit_wildtype_peptide',
			#'cbe_edit_mutant_peptide',
 			#'cbe_edit_mutation',
 			#'abe_edit_wildtype_peptide',
		 	#'abe_edit_mutant_peptide',
 			#'abe_edit_mutation',
			'normalised_counts_mean',
		]

		limma_tag = [
			'sample_1_TPM',
			'sample_2_TPM',
		]

		count_tags = [
			"z_score",
			"positive_controls_rank",
			"negative_controls_rank",
			"background_rank",
			"sample_1_counts_mean",
			"sample_2_counts_mean",
			"sample_1_counts",
			"sample_2_counts",
		#	"ratio"
		]

		annotation_tags = [
			'guide_snps',
			'guide_motifs',
			'guide_feature_domain',
			'guide_feature_structure',
			'guide_feature_isoform',
			'guide_feature_modification',
			'guide_feature_region',
			'guide_feature_topology',
			'guide_feature_mutagenesis',
			'guide_feature_other'
		]

		region_tags = [
			"source",
			"type",
			"didi_identifier",
			"pmid",
			"dataset",
			"specficity_class",
			"pdb",
			"slimprint_p"
		]

		#--------------------------------------------------------#

		rows = ["\t".join(p_limma_tag + tags+limma_tag+count_tags+annotation_tags+region_tags)]

		for guide_id in self.data[self.options['screen_name']]['ngs_data']:
			if 'output_include_controls' in self.options:
				if self.options['output_include_controls'] == False:
					if self.data[self.options['screen_name']]['ngs_data'][guide_id]['control']:
						continue
			
			try:
				if self.data[self.options['screen_name']]['ngs_data'][guide_id]['contains_low_count_value']: continue
				if self.data[self.options['screen_name']]['ngs_data'][guide_id]['contains_zero_value']: continue

				if "split_by_accesion" in self.options:
					if self.data[self.options['screen_name']]['ngs_data'][guide_id]['control']: continue

				row = []
				region_id = guide_id
				
				for tag in p_limma_tag:
					try:
						row.append("%1.2f"%self.data[self.options['screen_name']]['limma_scores'][guide_id][tag])
					except:
						row.append("")
						utilities_error.printError()

				if self.data[self.options['screen_name']]['ngs_data'][guide_id]['control'] == False:
					try:
						region_id = "_".join([self.data[self.options['screen_name']]['ngs_data'][guide_id]['accession'],str(self.data[self.options['screen_name']]['ngs_data'][guide_id]['peptide_start']),str(self.data[self.options['screen_name']]['ngs_data'][guide_id]['peptide_end'])])
					except:
						region_id = "_".join([self.data[self.options['screen_name']]['ngs_data'][guide_id]['annotation'],str(self.data[self.options['screen_name']]['ngs_data'][guide_id]['peptide_start']),str(self.data[self.options['screen_name']]['ngs_data'][guide_id]['peptide_end'])])
				
				for tag in tags:
					try:
						row.append(str(self.data[self.options['screen_name']]['ngs_data'][guide_id][tag]))
					except:
						row.append("")

						annotation_tags

				for tag in limma_tag:
					if tag == "sample_1_TPM":
						tag = 'mutant'

					if tag == "sample_2_TPM":
						tag = 'wildtype'

					try:
						if tag in ['mutant','wildtype']:
							row.append("|".join(["%1.2f"%x for x in self.data[self.options['screen_name']]['limma_scores'][guide_id][tag]]))
						else:
							row.append("%1.2f"%self.data[self.options['screen_name']]['limma_scores'][guide_id][tag])
					except:
						row.append("")
						utilities_error.printError()
			
				try:
					row.append("%1.5f"%self.data[self.options['screen_name']]['ngs_data'][guide_id]["z_score"])
					row.append("%1.5f"%self.data[self.options['screen_name']]['ngs_data'][guide_id]["positive_controls_rank"])
					row.append("%1.5f"%self.data[self.options['screen_name']]['ngs_data'][guide_id]["negative_controls_rank"])
					row.append("%1.5f"%self.data[self.options['screen_name']]['ngs_data'][guide_id]["background_rank"])
				except:
					row.append("")
					row.append("")
					row.append("")
					row.append("")

				try:
					counts = []
					
					for sample in self.data[self.options['screen_name']]['ngs_data'][guide_id]["data"]:
						counts.append(str(int(sum([x for x in self.data[self.options['screen_name']]['ngs_data'][guide_id]["data"][sample].values()])/len(self.data[self.options['screen_name']]['ngs_data'][guide_id]["data"][sample]))))

					for sample in self.data[self.options['screen_name']]['ngs_data'][guide_id]["data"]:
						counts.append("|".join([str(int(x)) for x in self.data[self.options['screen_name']]['ngs_data'][guide_id]["data"][sample].values()]))

					row += counts
				except:
					row.append("")
					row.append("")
					row.append("")
					row.append("")

				try:
					row.append("|".join(["%1.2f"%x for x in self.data[self.options['screen_name']]['ngs_data'][guide_id]["counts"]]))
				except:
					row.append("")

				try:
					row.append("|".join(["%1.2f"%x for x in self.data[self.options['screen_name']]['ngs_data'][guide_id]["normalised_counts"]]))
				except:
					row.append("")

				for tag in annotation_tags:
					try:
						row.append(str(self.data[self.options['screen_name']]['ngs_data'][guide_id][tag]))
					except:
						row.append("")


				for region_tag in region_tags:
					try:
						row.append(str(self.options["design"]['region_annotation'][region_id][region_tag]))
					except:
						row.append("")
				
				rows.append("\t".join(row))
			except:
				logger.error("Error adding row for " + guide_id)
				utilities_error.printError()
				#raise
				
		tdt_path = os.path.join(self.options['data_directory'],"tdt", self.options['screen_name'] + "." + str(self.options['query_time_point']) + "-" + str(self.options['compare_time_point']) + ".guide_annotation.tdt")
		logger.info("Writing:" + tdt_path)
		open(tdt_path,'w').write("\n".join(rows))
		
	##################################################################################################

	def write_regions_results_json(self):
		if self.options["output_include_region_annotation_json"]:
			for region_id in self.options["design"]['region_annotation']:

				json_path = os.path.join(self.options['data_directory'], "regions",region_id + "." + self.options['screen_name'] + "." + str(self.options['query_time_point']) + "-" + str(self.options['compare_time_point']) + ".region.annotation.json")
				#logger.debug("Writing: " + json_path)

				regions_data = {
					#"screen":self.data[self.options['screen_name']]['regions'][region_id],
					"design":self.options["design"]['region_annotation'][region_id]
				}

				with open(json_path,'w') as outfile:
					json.dump(regions_data, outfile)
		
	def write_protein_results_json(self):
		if self.options["output_include_protein_annotation_json"]:
			guide_ids = list(self.data[self.options['screen_name']]['ngs_data'])
			guide_ids.sort()

			protein_guide_data = {}
			for guide_id in guide_ids:
				if 'accession' in self.data[self.options['screen_name']]['ngs_data'][guide_id]:
					if self.data[self.options['screen_name']]['ngs_data'][guide_id]['control']: continue
					accession = self.data[self.options['screen_name']]['ngs_data'][guide_id]['accession']

					if accession not in protein_guide_data:
						protein_guide_data[accession] = {}

					protein_guide_data[accession][guide_id] = self.data[self.options['screen_name']]['ngs_data'][guide_id]

			for accession in protein_guide_data:
				json_path = os.path.join(self.options['data_directory'], "proteins",accession + "." + self.options['screen_name'] + "." + str(self.options['query_time_point']) + "-" + str(self.options['compare_time_point']) + ".protein.annotation.json")
				concise_json_path = os.path.join(self.options['data_directory'], "proteins",accession + "." + self.options['screen_name'] + "." + str(self.options['query_time_point']) + "-" + str(self.options['compare_time_point']) + ".protein.annotation.concise.json")
				
				with open(json_path,'w') as outfile:
					json.dump(protein_guide_data[accession], outfile)

				######

				concise_tags = ["peptide_start","peptide_end","guide_off_target_exact","splice_site"]
				concise_limma_tags = ["log_fold_change","log_p_value_adj"]

				response = queryRunner.queryRunner("database_access","get_attributes",{"accession":accession,"database_name":"database"}).run()
			
				concise_protein_guide_data = {
					"guides":{},
					"alphafold_accessibility":response['data']['alphafold_accessibility'],
					"alphafold_secondary_structure":response['data']['alphafold_secondary_structure'],
					"conservation":response['data']['Conservation_raw']['metazoa']
				}

				for guide_id in protein_guide_data[accession]:
					concise_protein_guide_data["guides"][guide_id] = {}
					for concise_tag in concise_tags:
						if concise_tag in protein_guide_data[accession][guide_id]:
							concise_protein_guide_data["guides"][guide_id][concise_tag] = protein_guide_data[accession][guide_id][concise_tag]
						else:
							concise_protein_guide_data["guides"][guide_id][concise_tag] = ""
					
					for concise_limma_tag in concise_limma_tags:
						if concise_limma_tag in protein_guide_data[accession][guide_id]["limma_scores"]:
							concise_protein_guide_data["guides"][guide_id][concise_limma_tag] = protein_guide_data[accession][guide_id]["limma_scores"][concise_limma_tag]
						else:
							concise_protein_guide_data["guides"][guide_id][concise_limma_tag] = ""
					
				with open(concise_json_path,'w') as outfile:
					json.dump(concise_protein_guide_data, outfile)

	##################################################################################################
 
	def write_results_json(self):
		json_path = os.path.join(self.options['data_directory'], "json",self.options['screen_name'] + "." + str(self.options['query_time_point']) + "-" + str(self.options['compare_time_point']) + ".region.annotation.json")
		logger.debug("Writing: " + json_path)
		with open(json_path,'w') as outfile:
			json.dump(self.options["design"]['region_annotation'], outfile)

		json_path = os.path.join(self.options['data_directory'], "json",self.options['screen_name'] + "." + str(self.options['query_time_point']) + "-" + str(self.options['compare_time_point']) + ".guide.annotation.json")
		logger.debug("Writing: " + json_path)
		with open(json_path,'w') as outfile:
			json.dump(self.data[self.options['screen_name']]['ngs_data'], outfile)		

		json_path = os.path.join(self.options['data_directory'],"json",self.options['screen_name'] + "." + str(self.options['query_time_point']) + "-" + str(self.options['compare_time_point']) + ".guide.controls.json")
		logger.debug("Writing: " + json_path)
		with open(json_path,'w') as outfile:
			json.dump(self.data[self.options['screen_name']]['controls'], outfile)		

	##################################################################################################

	def write_region_annotation(self):
		if self.options["output_include_region_annotation_json"]:
			header = "\t".join([
				"Region",
				'Accession',
				'Gene Name',
				'Start',
				'End',
				'Peptide',
				'Source',
				'Dataset',
				'Motif Class',
				'Motif Type',
				'PMID',
				'PDB',
				"DepMap Mean Gene Effect",
				"HAP1 Gene Effect Ratio",
				"HAP1 Gene Effect Q-value",
				"HAP1 Gene Effect Required",
				'Alphafold Accessibility',
				'Mean Fold Change (log2)',
				'M-W p-value (-log10)',
				"Positive Controls Bootstraps - Depleted",
				"Negative Controls Bootstraps - Depleted",
				"Background Bootstraps - Depleted",
				"Positive Controls Bootstraps - Enriched",
				"Negative Controls Bootstraps - Enriched",
				"Background Bootstraps - Enriched",
				"Guides Count",
				"Guides Count",
				"Off Target Genome Matches Count Mean",
				"Significant Guides Count",
				"Guides Count Region",
				"Guides Count Region Stats",
				"Significant Guides Count",
				"Confident Guides Count",
				"Significant Guides Proportion",
				"Significant Guides Proportion",
				"Sample Name",
				"Motifs",
				"Mutagenesis",
				"SNPs",
				"Modifications",
				"Extended Peptide",
				"Accessibility Raw",
				"AF SS classification",
				"AF SS"
			])


			rows = []

			try:
				source = "|".join(self.options["design"]['region_annotation'][region_id]['source'])
				dataset = "|".join(self.options["design"]['region_annotation'][region_id]['dataset'])
				specficity_class = "|".join(self.options["design"]['region_annotation'][region_id]['specficity_class'])
				region_type = "|".join(self.options["design"]['region_annotation'][region_id]['type'])
				pmid = "|".join(self.options["design"]['region_annotation'][region_id]['pmid'])
				pdb = "|".join(self.options["design"]['region_annotation'][region_id]['pdb'])
			except:
				source = ""
				dataset = ""
				specficity_class = ""
				region_type = ""
				pmid = ""
				pdb = ""
			
			for region_id in self.options["design"]['region_annotation']:
				row = []

				off_target_count_mean = 0
				guide_count_region_annotation = 0
				guide_count_stats = 0
				significant_guide_count = 0

				if 'guides' in self.options["design"]['region_annotation'][region_id]:
					guide_count_stats = len(self.options["design"]['region_annotation'][region_id]['guides'])
					
				if 'guide_ids' in self.options["design"]['region_annotation'][region_id]:
					guide_count_region_annotation = len(self.options["design"]['region_annotation'][region_id]["guide_ids"])
					for guide_id in self.options["design"]['region_annotation'][region_id]["guide_ids"]:
						try:
							if self.options["design"]['region_annotation'][region_id]['guides'][guide_id]['limma_scores']['log_p_value'] > 2: 
								significant_guide_count += 1
						except:
							pass
				
					off_target_counts = []
					for guide_id in self.options["design"]['region_annotation'][region_id]["guide_ids"]:
						try:
							off_target_counts.append(int(self.data[self.options['screen_name']]['ngs_data'][guide_id]['guide_off_target_exact']))
						except:
							logger.error(guide_id + " no guide_off_target_exact")
				
					off_target_count_mean = -1
					try:
						off_target_count_mean = sum(off_target_counts)/len(off_target_counts)
					except:
						logger.error(region_id + " off_target_count_mean")
				
				try:
					motifs = []
					snps = []
					mutagenesis = []
					modifications = []

					if 'peptools' in self.options["design"]['region_annotation'][region_id]:
						if 'Motif' in self.options["design"]['region_annotation'][region_id]['peptools']:
							for motif_instance in self.options["design"]['region_annotation'][region_id]['peptools']['Motif']:
								motifs.append(str(motif_instance['name'].replace("SLiM","Unclassified")) + "@" + str(motif_instance['description']['start']) + ":" +  str(motif_instance['description']['stop']) )

				
						if 'SNP' in self.options["design"]['region_annotation'][region_id]['peptools']:
							for snp in self.options["design"]['region_annotation'][region_id]['peptools']['SNP']:
								
								if snp['name'].count("frameshift") == 0 and snp['name'].count("stop") == 0:
									if 'Disease' in snp['description']:
										snps += [snp['name'] + ":" +  snp['description']['mutation'] + ":" +  snp['description']['Disease'] + "@" + str(snp['description']['start'])] 
									else:
										snps += [snp['name'] + ":" +  snp['description']['mutation'] + ":-@" + str(snp['description']['start'])] 

						
						if 'Mutagenesis' in self.options["design"]['region_annotation'][region_id]['peptools']:
							try:
								for x in self.options["design"]['region_annotation'][region_id]['peptools']['Mutagenesis']:
									if "mutation" in x['description']:
										mutagenesis += [x['name'] + "@" + str(x['description']['start']) + ":" + str(x['description']['mutation'])  + ":" + str(x['name']) ]
									else:
										mutagenesis += [x['name'] + "@" + str(x['description']['start']) + ":-:" + str(x['name'])]
							except:
								pprint.pprint(self.options["design"]['region_annotation'][region_id]['peptools']['Mutagenesis'])

						if 'Modification' in self.options["design"]['region_annotation'][region_id]['peptools']:
							try:
								modifications += [x['name'] + "@" + str(x['description']['start']) + ":" + str(x['description']['description']) + ":" + str(x['description']['enzymes'])  for x in self.options["design"]['region_annotation'][region_id]['peptools']['Modification']]
							except:
								pprint.pprint(self.options["design"]['region_annotation'][region_id]['peptools']['Modification'])

					row = [
						region_id,
						self.options["design"]['region_annotation'][region_id]['accession'],
						self.options["design"]['region_annotation'][region_id]['gene'] if 'gene' in self.options["design"]['region_annotation'][region_id] else "-",
						str(self.options["design"]['region_annotation'][region_id]['peptide_start']),
						str(self.options["design"]['region_annotation'][region_id]['peptide_end']),
						str(self.options["design"]['region_annotation'][region_id]['peptide']) if 'peptide' in self.options["design"]['region_annotation'][region_id] else "-",
						source,
						dataset,
						specficity_class,
						region_type,
						pmid,
						pdb,
						str(self.options["design"]['region_annotation'][region_id]['mean_gene_effect']) if 'mean_gene_effect' in self.options["design"]['region_annotation'][region_id] else "-",
						str(self.options["design"]['region_annotation'][region_id]['hap1_gene_effect']['ratio']) if 'hap1_gene_effect' in self.options["design"]['region_annotation'][region_id] else "-",
						str(self.options["design"]['region_annotation'][region_id]['hap1_gene_effect']['q_value'])  if 'hap1_gene_effect' in self.options["design"]['region_annotation'][region_id] else "-",
						str(self.options["design"]['region_annotation'][region_id]['hap1_gene_effect']['selected']) if 'hap1_gene_effect' in self.options["design"]['region_annotation'][region_id] else "-",
						"%1.3f"%self.options["design"]['region_annotation'][region_id]["accessibility"]  if 'accessibility' in self.options["design"]['region_annotation'][region_id] else "-",
						"%1.3f"%self.options["design"]['region_annotation'][region_id]['mean_fold_change']  if 'mean_fold_change' in self.options["design"]['region_annotation'][region_id] else "-",
						"%1.3f"%self.options["design"]['region_annotation'][region_id]['mw_p_log']  if 'mw_p_log' in self.options["design"]['region_annotation'][region_id] else "-",
						"%1.3g"%self.options["design"]['region_annotation'][region_id]['positive_controls_bootstraps'] if 'positive_controls_bootstraps' in self.options["design"]['region_annotation'][region_id] else "-",
						"%1.3g"%self.options["design"]['region_annotation'][region_id]['negative_controls_bootstraps'] if 'negative_controls_bootstraps' in self.options["design"]['region_annotation'][region_id] else "-",
						"%1.3g"%self.options["design"]['region_annotation'][region_id]['background_bootstraps'] if 'background_bootstraps' in self.options["design"]['region_annotation'][region_id] else "-",
						"%1.3g"%(1-self.options["design"]['region_annotation'][region_id]['positive_controls_bootstraps']) if 'positive_controls_bootstraps' in self.options["design"]['region_annotation'][region_id] else "-",
						"%1.3g"%(1-self.options["design"]['region_annotation'][region_id]['negative_controls_bootstraps']) if 'negative_controls_bootstraps' in self.options["design"]['region_annotation'][region_id] else "-",
						"%1.3g"%(1-self.options["design"]['region_annotation'][region_id]['background_bootstraps']) if 'background_bootstraps' in self.options["design"]['region_annotation'][region_id] else "-",
						str(guide_count_region_annotation),
						str(guide_count_stats),
						str(off_target_count_mean),
						str(significant_guide_count),
						str(self.options["design"]['region_annotation'][region_id]['guides_count']) if 'guides_count' in self.options["design"]['region_annotation'][region_id] else "-",
						str(self.options["design"]['region_annotation'][region_id]['significant_guides']) if 'significant_guides' in self.options["design"]['region_annotation'][region_id] else "-",
						str(self.options["design"]['region_annotation'][region_id]['confident_guides']) if 'confident_guides' in self.options["design"]['region_annotation'][region_id] else "-",
						"%1.3f"%(self.options["design"]['region_annotation'][region_id]['significant_guides_proportion']) if 'significant_guides_proportion' in self.options["design"]['region_annotation'][region_id] else "-",
						"%1.3f"%(self.options["design"]['region_annotation'][region_id]['confident_guides_proportion']) if 'confident_guides_proportion' in self.options["design"]['region_annotation'][region_id] else "-",
						"|".join(self.options["design"]['region_annotation'][region_id]['sample_name']) if 'sample_name' in self.options["design"]['region_annotation'][region_id] else "-",
						"|".join(motifs),
						"|".join(mutagenesis),
						"|".join(snps),
						"|".join(modifications)
					]

					try:
						row.append(self.options["design"]['region_annotation'][region_id]['peptools']['Hit'])
					except:
						row.append("-")

					try:
						row.append("%1.3f"%self.options["design"]['region_annotation'][region_id]['peptools']['Accessibility_raw'])
					except:
						row.append("-")
					
					try:
						if 'formatted' in self.options["design"]['region_annotation'][region_id]['peptools']['AF_SS_classification']:
							row.append(self.options["design"]['region_annotation'][region_id]['peptools']['AF_SS_classification']['formatted'][0]['name'])
						else:
							row.append("-")
					except:
						row.append("-")
					
					try:
						row.append(self.options["design"]['region_annotation'][region_id]['peptools']['AF_SS_classification']['per_residue'])
					except:
						row.append("-")
						
					try:
						rows.append("\t".join(copy.deepcopy(row)))
					except:
						logger.error("Error formatting row for " + region_id)
						utilities_error.printError()
				except:
					logger.error("Error adding row for " + region_id)
					utilities_error.printError()
					rows.append(region_id + "\t" + str(utilities_error.getError()))

			tdt_path = os.path.join(self.options['data_directory'], "tdt",self.options['screen_name'] + "." + str(self.options['query_time_point']) + "-" + str(self.options['compare_time_point']) + ".region_annotation.tdt")
			logger.info("Writing: " + tdt_path)
			open(tdt_path,'w').write("\n".join([header] + rows))

	##################################################################################################

	def write_region_server_json(self):
		mw_region_statistics_data = {
			"points":[],
			"range":{
				"fold_change":[1000,-1000],
				"p_value":[1000,-1000],
			}
		}

		tags = [
			'sample_name',
			'accession',
			'gene',
			'peptide_start',
			'peptide_end',
			'peptide',
			'mean_fold_change',
			'mw_p_log',
			'mean_gene_effect',
			'accessibility',
			'source',
			'dataset',
			'specficity_class',
			'pmid',
			'pdb',
			'mutations',
			'modifications',
			'snps',
			'motifs',
			'negative_controls_bootstraps',
			'significant_guides',
			'confident_guides',
			'mutation_track'
		]

		for region_id in self.options["design"]['region_annotation']:
			try:
				if 'significant_guides' not in self.options["design"]['region_annotation'][region_id]:
					logger.error("significant_guides not in region_id: " + region_id)
					continue

				if self.options["design"]['region_annotation'][region_id]['significant_guides'] >= 1 or self.options["design"]['region_annotation'][region_id]['mw_p_log'] >= 2:
					try:
						logger.debug("Adding " + region_id)
						
						regions_data = {}
						regions_data['total_guides'] = len(self.options["design"]['region_annotation'][region_id]['guide_ids'])

						most_significant_guide = {"log_p_value":0,"log_fold_change":0}
						for guide_id in self.options["design"]['region_annotation'][region_id]['guides']:
							try:
								if most_significant_guide['log_p_value'] > self.data[self.options['screen_name']]['limma_scores'][guide_id]['log_p_value']:
									
									most_significant_guide = {
										"log_p_value":self.data[self.options['screen_name']]['limma_scores'][guide_id]['log_p_value'],
										"log_fold_change":self.data[self.options['screen_name']]['limma_scores'][guide_id]['log_fold_change']
										}
							except:
								pass
						
						regions_data['most_significant_guide'] = most_significant_guide

						for tag in tags:
							try:
								regions_data[tag] = self.options["design"]['region_annotation'][region_id][tag]
							except:
								pass
						
						mw_region_statistics_data['points'].append(regions_data)

						if self.options["design"]['region_annotation'][region_id]['mean_fold_change'] < mw_region_statistics_data["range"]["fold_change"][0]: 
							mw_region_statistics_data["range"]["fold_change"][0] = self.options["design"]['region_annotation'][region_id]['mean_fold_change']
						if self.options["design"]['region_annotation'][region_id]['mean_fold_change'] > mw_region_statistics_data["range"]["fold_change"][1]: 
							mw_region_statistics_data["range"]["fold_change"][1] = self.options["design"]['region_annotation'][region_id]['mean_fold_change']
						if self.options["design"]['region_annotation'][region_id]['mw_p_log'] < mw_region_statistics_data["range"]["p_value"][0]: 
							mw_region_statistics_data["range"]["p_value"][0] = self.options["design"]['region_annotation'][region_id]['mw_p_log']
						if self.options["design"]['region_annotation'][region_id]['mw_p_log'] > mw_region_statistics_data["range"]["p_value"][1]:
							mw_region_statistics_data["range"]["p_value"][1] = self.options["design"]['region_annotation'][region_id]['mw_p_log']
					except:
						logger.error(utilities_error.getError())
			except:
				logger.error(region_id)
				logger.error(utilities_error.getError())
				
				
		json_path =  os.path.join(self.options['data_directory'],self.options['screen_name'] + "." + str(self.options['query_time_point']) + "-" + str(self.options['compare_time_point']) + ".base_editing_region_scatterplot.json")
		logger.debug("Writing: " + json_path)
		with open(json_path,'w') as outfile:
			json.dump(mw_region_statistics_data, outfile)

	###------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------###

	def write_protein_plots(self):
		base_editing_protein_plotsObj = base_editing_protein_plots.baseEditingProteinPlots()
		base_editing_protein_plotsObj.options['plot_guide_data'] = self.data[self.options['screen_name']]['ngs_data']
		base_editing_protein_plotsObj.options['sample'] = self.options['screen_name'] + "_" + str(self.options['query_time_point']) + "-" + str(self.options['compare_time_point'])
		base_editing_protein_plotsObj.options['protein_png_output_directory_path'] = os.path.join(self.options["data_directory"],"protein_plots")

		if not os.path.exists(base_editing_protein_plotsObj.options['protein_png_output_directory_path']):
			os.mkdir(base_editing_protein_plotsObj.options['protein_png_output_directory_path'])
		
		base_editing_protein_plotsObj.plot_proteins()
	
	###------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------###

	def plot_control_data(self):
		baseEditingControlPlotterObj = base_editing_control_plots.baseEditingControlPlotter()
		baseEditingControlPlotterObj.data = self.data[self.options['screen_name']]['ngs_data']
		baseEditingControlPlotterObj.options['data_directory'] = self.options['data_directory']
		baseEditingControlPlotterObj.options['screen_name'] = self.options["original_screen_name"] 
		baseEditingControlPlotterObj.options['compare_time_point'] = self.options['compare_time_point']

		try:
			if "control_boxplot_keys" in self.options and "control_boxplot_labels" in self.options:
				baseEditingControlPlotterObj.plot_control_boxplot(sorted_keys=self.options['control_boxplot_keys'],sorted_keys_label=self.options['control_boxplot_labels'])
			else:
				baseEditingControlPlotterObj.plot_control_boxplot()
		except:
			logging.error('Skipping plot_control_boxplot')

		try:
			if "control_roc_curve_labels" in self.options:
				baseEditingControlPlotterObj.plot_control_roc_curves(positive_control_labels=self.options["control_roc_curve_labels"])
			else:
				baseEditingControlPlotterObj.plot_control_roc_curves()
		except:
			logging.error('Skipping plot_control_roc_curves')		
			
	###------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------###
	###------------------------######------------------------######------------------------###

	def calculate_coverage(self):
		self.data['coverage']["sequencing_guide_total"] = len(self.data[self.options['screen_name']]['ngs_data'])
		self.data['coverage']["timepoints"] = {}
		
		for guide_id in self.data[self.options['screen_name']]['ngs_data']:
			for timepoint in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data']:
				if timepoint not in self.data['coverage']["timepoints"]: self.data['coverage']["timepoints"][timepoint] = {}
				for replicate in self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][timepoint]:
					if replicate not in self.data['coverage']["timepoints"][timepoint]: self.data['coverage']["timepoints"][timepoint][replicate] = 0
					if self.data[self.options['screen_name']]['ngs_data'][guide_id]['data'][timepoint][replicate] > 0:
						self.data['coverage']["timepoints"][timepoint][replicate] += 1

		logger.debug("COVERAGE")
		logger.debug(self.data['coverage'])

		for timepoint in self.data['coverage']["timepoints"]:
			for replicate in self.data['coverage']["timepoints"][timepoint]:
				row = [
					str(timepoint),
					str(replicate),
					str(self.data['coverage']["timepoints"][timepoint][replicate]),
					str(self.data['coverage']["sequencing_guide_total"]),
					str(self.data['coverage']["design_guide_total"]),
					"%1.4f"%(self.data['coverage']["timepoints"][timepoint][replicate]/self.data['coverage']["sequencing_guide_total"]),
					"%1.4f"%(self.data['coverage']["timepoints"][timepoint][replicate]/self.data['coverage']["design_guide_total"]),
				]

				logger.debug("\t".join(row))

		json_path = os.path.join(self.options['data_directory'], "coverage",self.options['screen_name'] + ".guide.coverage.json")
		logger.debug("Writing: " + json_path)
		with open(json_path,'w') as outfile:
			json.dump(self.data['coverage'], outfile)		


	##################################################################################################
	##################################################################################################
	##################################################################################################
	
	def process_screen(self):
		self.read_job_file()

		errors = []
		self.options["original_screen_name"] = self.options["screen_name"]

		if "split_by_accesion" in self.options:
			screen_name = self.options["screen_name"]

			logger.debug(screen_name)
			for filter_accession in self.options["split_by_accesion"]:
				try:
					self.reset_options()
					response = queryRunner.queryRunner("uniprot","parse_basic",{"accession":filter_accession}).run()
				
					self.options["screen_name"] = response['data']['gene_name'] + "." + screen_name 
					self.options["filter_accession"] = filter_accession
					self.process_data()
				except:
					utilities_error.printError()
					logger.error(filter_accession + " processing failed")
					errors.append(filter_accession)
					raise
					
		else:
			self.process_data()

		if len(errors) > 0:
			logger.error("Failed accessions:")
			logger.error(errors)

	def process_data(self):
		self.setupScreen()

		if os.path.exists(self.options['counts_directory']):
			counts_files = os.listdir(self.options['counts_directory'])
		else:
			if isinstance( self.options['counts_file'], list):
				counts_files = self.options['counts_file']
				self.options['counts_directory'] = os.path.abspath(os.path.dirname(self.options['counts_file'][0]))
			else:
				counts_files = [self.options['counts_file']]
				self.options['counts_directory'] = os.path.abspath(os.path.dirname(self.options['counts_file']))

		self.data[self.options['screen_name']] = {
			"ngs_data":{},
			"controls":{},
			"metrics":{},
			"compare_time_points":[]
		}

		self.load_guides_mapping_file()
		self.load_guides_annotation_file()
		
		for counts_file	in counts_files:
			self.options['counts_file'] = os.path.join(self.options['counts_directory'],counts_file)

			if self.options['counts_file'].split(".")[-1] == "txt":
				logger.debug("Processing:" + self.options['counts_file'])
				self.load_time_points()
				self.load_counts_file()
				self.normalise_cpm_log2_files()
				self.normalise_t0_files() # Wrong order?
				self.plot_count_data()
				self.plot_t0_normalised_data()
		
		###----------------------###
  
		self.calculate_coverage()
		self.annotate_guides()
		
		if self.options['make_region_sliding_window']: 
			self.make_region_sliding_window()

		if self.options['make_region_tesselation']: 
			self.make_region_tesselation()
		 
		###----------------------###

		if os.path.exists(self.options['region_annotation_file']):
			self.load_region_annotation_file()
			self.load_region_annotation()

		###----------------------###

		self.limma_guide_statistics()
		self.mw_region_statistics()
		self.bootstrap_guide_statistics()
		self.annotate_guide_status()

		###----------------------###

		self.plot_control_data()
		self.annotate_regions()
		
		if 'output_include_protein_plots' in self.options:
			if self.options['output_include_protein_plots']:
				self.write_protein_plots()

		###----------------------###

		self.write_guide_annotation_file()
		self.write_region_annotation()
		self.write_results_json()
		self.write_protein_results_json()
		self.write_regions_results_json()

################################################################################
## SECTION II: MAIN PROGRAM                                                  ##
################################################################################

if __name__ == "__main__":

	baseEditingParserObj = baseEditingParser()
	baseEditingParserObj.process_screen()

################################################################################
## END OF SECTION II                                                          ##
################################################################################

