import sys,os,inspect,pprint,json,time

file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../"))
sys.path.append(os.path.abspath(os.path.join(file_path,"../utilities")))

sys.path.append(os.path.join(file_path,"../data_management/"))
import queryRunner

import utilities_error
import utilities_downloader
import utilities_rest

import config_reader
import option_reader

#-----
import logging
logger = logging.getLogger(__name__)
#-----

class motifQueryRestManager():
	
	#####------------------------------------------------------------------#####
	
	def __init__(self):

		self.options = config_reader.load_configeration_options(sections=['general'])
		self.options.update({
			"peptide_padding":3,
			"collapsed":True,
			"by_column":None,	
			"by_column_search_term":None,
			"data_sources":None,
			"group_by":None,
			"columns":'detailed',
			"include_cols":None,
			"is_prototype_structure":None,
			"additional_params":None,
			"is_active_version":None,
			"curator_instance_assessment":"TP",
			"api_endpoint":None,
			"motif_database_source":"motifs"
			})
				
	############################################################################
	############################################################################
	############################################################################

	def set_columns_required(self):
		if self.options['columns'] == 'detailed':
			self.options['include_cols'] = [
				'completeness',
				'curation_ids',
				'motif_proteins',
				'domain_proteins',
				'domain_protein_relation',
				'curators',
				'sources',
				'references',
				'pdb_ids',
				'motif_classes',
				'specificity_classes',
				'pockets',
				'curator_assessment',
				'details',
				'other_notes',
				'interaction_types'
			]
			self.options['columns'] = [
				'specificity_classes',
				'motif_classes'
			]

		if self.options['columns'] == 'compact':
				self.options['columns'] = []

	def set_api_path(self):
		if self.options["api_endpoint"] != None:
			self.options['api_path'] = '/restapi/' + self.options['motif_database_source'] + '/' + self.options["api_endpoint"]
		else:
			if self.options['collapsed']:
				self.options['api_path'] = '/restapi/' + self.options['motif_database_source'] + '/collapsed_instances'
			else:
				self.options['api_path'] = '/restapi/' + self.options['motif_database_source'] + '/instances'
				
	def make_query(self):
		logger.debug("Making query from params")
		params = []
		
		if self.options['data_sources'] != None:
			params.append("sources=" + self.options['data_sources'])

		if self.options['by_column'] != None and self.options['by_column_search_term'] != None:
			params.append(self.options['by_column'] + "=" + str(self.options['by_column_search_term']))
		
		if self.options['group_by'] != None and self.options['collapsed']:
			params.append("group_by=" + self.options['group_by'])

		if self.options['columns'] != None and self.options['collapsed']:
			params.append("columns=" + ",".join(self.options['columns']))

		if self.options['include_cols'] != None and self.options['collapsed']:
			params.append("include_cols=" + ",".join(self.options['include_cols']))

		if self.options['is_prototype_structure'] != None and self.options['collapsed']:
			params.append("is_prototype_structure=" + str(self.options['is_prototype_structure']))

		if self.options['is_active_version'] != None and self.options['collapsed']:
			params.append("is_active_version=" + str(self.options['is_active_version']))
		
		if self.options['curator_instance_assessment'] != None and self.options['collapsed']:
			params.append("curator_instance_assessment=" + str(self.options['curator_instance_assessment']))

		if self.options['additional_params'] != None:
			for additional_param in self.options['additional_params']:
				params.append(additional_param + "=" + str(self.options['additional_params'][additional_param]))

		self.options['params'] = "&".join(params)

	############################################################################

	def run_motifs_query(self):

		self.set_columns_required()
		self.set_api_path()
		self.make_query()

		return self.run_query()

	def run_query(self):
		url = self.options['api_path'] + "?" + self.options['params']
		logger.debug("Submitting:" + url)

		sessionDownloaderObj = utilities_rest.RestSubmission()
		sessionDownloaderObj.server_url = self.options['server_url']
		response = sessionDownloaderObj.run_query(url,submission="GET")
		logger.debug("Received response from:" + url)
		return response

	def run_update(self):
		logger.debug("Submitting:" + self.options['api_path'] )

		sessionDownloaderObj = utilities_rest.RestSubmission()
		
		if "username" not in self.options:return {"status":"Error","error_type":"Username/Password issue"}
		if "password" not in self.options:return {"status":"Error","error_type":"Username/Password issue"}

		if self.options["username"] != "" and self.options["password"] != "":
			logger.debug("Loggin in with usename and password: " + self.options["username"])
			sessionDownloaderObj.username = self.options["username"] #""
			sessionDownloaderObj.password = self.options["password"] #""

			login_token_response = sessionDownloaderObj.loginToken()
			sessionDownloaderObj.token = login_token_response['access_token']

			sessionDownloaderObj.server_url = self.options['server_url']
			response = sessionDownloaderObj.run_query(self.options['api_path'] ,self.options['params'],submission="POST")
			logger.debug("Received response from:" + self.options['api_path'] )
			return response
		else:
			logger.error("Username/Password issue")
			return {"status":"Error","error_type":"Username/Password issue"}

	#####------------------------------------------------------------------#####
		
	def check_query(self):
		return {"status":"Success"}

	#####------------------------------------------------------------------#####	
	
	def get_instances(self):
		return self.run_motifs_query()

	def get_instances_without_specificity_class(self):
		#http://slim.icr.ac.uk/restapi/motifs/collapsed_instances?group_by=pdb&specificity_acc=None
		self.options['group_by'] = "pdb"
		self.options['columns'] = "compact"
		
		self.options['additional_params'] = {
			"specificity_acc":"None",
		}
		
		return self.run_motifs_query()

	def get_instances_with_structural_information(self):
		#http://slim.icr.ac.uk/restapi/motifs/collapsed_instances?group_by=pdb&pocket_acc=None&include_cols=pockets
		self.options['group_by'] = "pdb"
		self.options['columns'] = "compact"
		self.options['collapsed'] = False
		self.options['include_cols'] = ["pockets",'other_notes']
		self.options['additional_params'] = {
			"source":"PDB"
		}
		
		return self.run_motifs_query()
		
	def get_instances_without_pocket_class(self):
		#http://slim.icr.ac.uk/restapi/motifs/collapsed_instances?group_by=pdb&pocket_acc=None&include_cols=pockets
		self.options['group_by'] = "pdb"
		self.options['columns'] = "compact"
		self.options['collapsed'] = False
		self.options['include_cols'] = ["pockets",'other_notes']
		self.options['additional_params'] = {
			"pocket_acc":"None",
			"source":"PDB"
		}
		
		return self.run_motifs_query()

	def get_prototype_structures(self):
		#http://slim.icr.ac.uk/restapi/motifs/pdb_motifs?is_prototype_structure=true&columns=specificity_classes&limit=2
		self.options["api_endpoint"] = "pdb_motifs"
		self.options['columns'] = ['specificity_classes']

		self.options['additional_params'] = {
			"is_prototype_structure":True,
		}
		return self.run_motifs_query()


	def get_instances_with_structures(self):
		#http://slim.icr.ac.uk/restapi/motifs/collapsed_instances?group_by=pdb&include_cols=pockets,specificity_classes&limit=10
		self.options['group_by'] = "pdb"
		self.options['columns'] = "compact"
		self.options['include_cols'] = ["pockets","specificity_classes"]
		
		return self.run_motifs_query()

	def get_instances_by_binding_protein_region(self):
		self.options['additional_params'] = {
			"domain_protein_uniprot_accession":self.options['accession'],
			"domain_start":self.options['domain_start'],
			"domain_end":self.options['domain_end']
		}

		logger.debug(self.options['additional_params'])
		return self.run_motifs_query()

	def get_instances_by_motif_accession(self):
		self.options['additional_params'] = {
			"motif_protein_uniprot_accession":self.options['accession'],
		}
		logger.debug(self.options['additional_params'] )
		return self.run_motifs_query()

	def get_instances_by_domain_accession(self):
		self.options['additional_params'] = {
			"domain_protein_uniprot_accession":self.options['accession'],
		}
		logger.debug(self.options['additional_params'] )
		return self.run_motifs_query()

	def get_instances_by_specificity_class(self):
		self.options['additional_params'] = {
			"specificity_class":self.options['specificity_class'],
		}
		logger.debug(self.options['additional_params'] )
		return self.run_motifs_query()

	def get_instances_by_specificity_acc(self):
		self.options['additional_params'] = {
			"specificity_acc":self.options['specificity_acc'],
		}
		logger.debug(self.options['additional_params'] )
		return self.run_motifs_query()

	def get_instances_by_elm_class(self):
		self.options['additional_params'] = {
			"motif_class":self.options['elm_class'],
		}
		logger.debug(self.options['additional_params'] )
		return self.run_motifs_query()
	
	def get_instances_by_pfam_id(self):
		self.options['additional_params'] = {
			"domain_pfam_accession":self.options['pfam_id'],
		}
		logger.debug(self.options['additional_params'] )
		return self.run_motifs_query()

	def get_instances_by_curation_id(self):
		self.options['group_by'] = "motif"
		self.options['include_cols'] = ["curation_ids","motif_proteins","domain_proteins","references","pdb_ids","motif_classes","specificity_classes","pocket"]

		self.options['additional_params'] = {
			"curation_id":self.options['curation_id']
		}

		logger.debug(self.options['additional_params'] )
		return self.run_motifs_query()

	def get_instances_pdb_pocket_fasta(self):
		fasta = ""

		self.options['group_by'] = "pdb"
		self.options['columns'] = "compact"
		self.options['collapsed'] = False
		self.options['include_cols'] = ["pockets",'other_notes']
		self.options['additional_params'] = {
			"source":"PDB"
		}
		
		counter = 0
		instances = self.run_motifs_query()
		for instance in instances:
			try:
				counter += 1
				print(instance['interaction_pdb'][0],counter,len(instances))
				
				sequences = queryRunner.queryRunner("pdb","get_sequence_structure",{"pdb_id":instance['interaction_pdb'][0],"chain_id":instance['other_notes']['domain chain']}).run()
				for protein_accession in sequences['data']:
					
					fasta += ">" + " ".join([
						instance['interaction_pdb'][0] + "." + instance['other_notes']['domain chain'],
						protein_accession, 
						str(sequences['data'][protein_accession][0]['unp_start']), 
						str(sequences['data'][protein_accession][0]['unp_end'])
					])+ "\n"

					fasta +=sequences['data'][protein_accession][0]['construct'] + "\n"
			except:
				pass
		print(fasta)
		open("motif_pocket_pdb.fasta","w").write(fasta)
		return {}

	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####

	def get_pocket_class(self):

		if self.options["database_source"] == "ELM":
			self.options['api_path'] = 'restapi/elm/pockets'
		else:
			self.options['api_path'] = 'restapi/motifs/pockets'

		self.options['params'] = {}
		self.options['additional_params'] = {
			"pocket_id":self.options['pocket_id']
		}

		params = []
		for additional_param in self.options['additional_params']:
			params.append(additional_param + "=" + str(self.options['additional_params'][additional_param]))

		self.options['params'] = "&".join(params)
		pocket = self.run_query()
		return pocket

	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####

	def get_specificity_class_pssm_data(self):
		self.options['api_path'] = 'restapi/motifs/pssms'
		self.options['params'] = {}
		self.options['additional_params'] = {
			"specificity_acc":self.options['specificity_id'],
			"is_main":True,
			"padding":self.options['peptide_padding']
		}

		if "scoring_method" in self.options:
			self.options['additional_params']["scoring_method"] = self.options["scoring_method"].replace(" ","%20")

		params = []
		for additional_param in self.options['additional_params']:
			params.append(additional_param + "=" + str(self.options['additional_params'][additional_param]))

		self.options['params'] = "&".join(params)
		pssm = self.run_query()
		return pssm
	
	def get_specificity_class_pssm_data(self):
		self.options['api_path'] = 'restapi/motifs/pssms'
		self.options['params'] = {}
		self.options['additional_params'] = {
			"ELM_identifier":self.options['ELM_identifier'],
			"is_main":True,
			"padding":self.options['peptide_padding']
		}

		if "scoring_method" in self.options:
			self.options['additional_params']["scoring_method"] = self.options["scoring_method"].replace(" ","%20")

		params = []
		for additional_param in self.options['additional_params']:
			params.append(additional_param + "=" + str(self.options['additional_params'][additional_param]))

		self.options['params'] = "&".join(params)
		pssm = self.run_query()
		return pssm

	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####

	def get_specificity_class_by_region(self):

		self.options['api_path'] = 'restapi/elm/instances'
		
		self.options['pfam_mapping']  = {}
		if 'accession' in self.options and 'region_start' in self.options and 'region_end' in self.options and 'pfam_id' not in self.options:
			self.options['pfam_mapping']  = queryRunner.queryRunner("uniprot","parse_region_annotation",{"accession":self.options['accession'],"region_start":self.options['region_start'],'region_end':self.options['region_end']}).run()
			self.options['pfam_id'] = list(self.options['pfam_mapping']['data']['Pfam'].keys())
		else:
			if 'pfam_id' in self.options:
				self.options['pfam_id'] = [self.options['pfam_id']]

		specificity_classes = {
			"protein":{}
		}

		self.options['additional_params'] = {
			"domain_protein_uniprot_accession":self.options['accession'],
			"group_by":"specificity_class",
			"include_cols":"domain_pfams,pockets,motif_classes"
		}
		
		logger.debug(self.options['additional_params'])
		specificity_classes_responses = self.run_motifs_query()
		for specificity_class in specificity_classes_responses:
			specificity_classes['protein'][specificity_class['specificity_acc']] = specificity_class

		if 'pfam_id' in self.options:
			#specificity_classes["pfam_id"] = self.options['pfam_id']
			specificity_classes["protein_domain"] = {}
			specificity_classes["domain"] = {}

			for pfam_id in self.options['pfam_id']:
				self.options['additional_params'] = {
					"domain_protein_uniprot_accession":self.options['accession'],
					"domain_pfam_accession":pfam_id,
					"group_by":"specificity_class",
					"include_cols":"domain_pfams,pockets,motif_classes"
				}
				
				logger.debug(self.options['additional_params'])
				specificity_classes_responses = self.run_motifs_query()

				for specificity_class in specificity_classes_responses:
					specificity_classes['protein_domain'][specificity_class['specificity_acc']] = specificity_class


				self.options['additional_params'] = {
					"domain_pfam_accession":pfam_id,
					"group_by":"specificity_class",
					"include_cols":"domain_pfams,pockets,motif_classes"
				}
				
				logger.debug(self.options['additional_params'])
				specificity_classes_responses = self.run_motifs_query()

				for specificity_class in specificity_classes_responses:
					specificity_classes['domain'][specificity_class['specificity_acc']] = specificity_class

		return specificity_classes
	
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####

	def get_specificity_class_xref_data(self):
		self.options['api_path'] = 'restapi/elm/specificity'
		self.options['params'] = {}
		self.options['additional_params'] = {
			"columns":"xrefs,specificity_id"
		}

		params = []
		for additional_param in self.options['additional_params']:
			params.append(additional_param + "=" + str(self.options['additional_params'][additional_param]))

		self.options['params'] = "&".join(params)
		specificity_class_xref = self.run_query()
		return specificity_class_xref


	def get_specificity_class_xref(self):
		xrefs = {}
		specificity_class_xref = self.get_specificity_class_xref_data()
		for specificity_class in specificity_class_xref:
			if specificity_class['specificity_acc'] not in list(xrefs.keys()):
				xrefs[specificity_class['specificity_acc']] = {
					"specificity_id":specificity_class['specificity_id'],
					"ELM_identifier":[],
					"ELM_accession":[]
				}

			for xref in specificity_class['xrefs']:
				if xref['source'] == 'ELM':
					xrefs[specificity_class['specificity_acc']]["ELM_identifier"].append(xref['external_id'])
					xrefs[specificity_class['specificity_acc']]["ELM_accession"] += xref['other_names']

		return xrefs

	def get_specificity_class_reverse_xref(self):
		xrefs = {}
		specificity_class_xref = self.get_specificity_class_xref_data()
		for specificity_class in specificity_class_xref:
			for xref in specificity_class['xrefs']:
				try:
					for other_name in xref['other_names']:
						xrefs[other_name] = {
							"specificity_id":xref['external_id'],
							"DIDI_identifier":[specificity_class['specificity_id']],
							"DIDI_accession":[specificity_class['specificity_acc']]
						}
				except:
					logger.error(xref)

		return xrefs
		
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####

	def get_specificity_class(self):
		self.options['additional_params'] = {
			"specificity_id":self.options['specificity_id']
		}
		
		pocket = self.get_specificity_class_by_additional_params()
		return pocket
	
	def get_specificity_class_by_specificity_acc(self):
		self.options['additional_params'] = {
			"specificity_acc":self.options['specificity_acc']
		}
		
		pocket = self.get_specificity_class_by_additional_params()
		return pocket
	
	
	def get_specificity_class_by_additional_params(self):
		self.options['api_path'] = 'restapi/motifs/specificity'
		self.options['params'] = {}

		params = []
		for additional_param in self.options['additional_params']:
			params.append(additional_param + "=" + str(self.options['additional_params'][additional_param]))

		self.options['params'] = "&".join(params)
		pocket = self.run_query()
		return pocket

	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####

	def get_alignment_by_specificity_acc(self):
		import tempfile
		url = os.path.join(self.options['server_url'],"restapi/motifs/alignment?is_main=True&is_active=True&padding=" + str(self.options['peptide_padding']) + "&specificity_acc=" + self.options['specificity_acc'])
		out_path = os.path.join(tempfile.gettempdir(),self.options['specificity_acc'] + "." + str(self.options['peptide_padding']) + ".json")

		logger.debug(url)	
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True)

		alignment = {}
		with open(out_path) as data_file:
			alignment = json.load(data_file)

		return alignment

	def get_alignment_by_elm_class(self):
		import tempfile
		url = os.path.join(self.options['server_url'],"restapi/elm/alignment?is_main=True&padding=" + str(self.options['peptide_padding']) + "&elm_identifier=" + self.options['elm_class'])
		out_path = os.path.join(tempfile.gettempdir(),self.options['elm_class'] + "." + str(self.options['peptide_padding']) + ".json")
					
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True)

		alignment = {}
		with open(out_path) as data_file:
			alignment = json.load(data_file)

		return alignment

	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####


	def add_tag_by_accession(self):
		response = {}
		
		
		for accession in self.options['accessions']:
			self.options['accession'] = accession
			motif_proteins = self.get_instances_by_motif_accession()
			domain_proteins = self.get_instances_by_domain_accession()


			for instance in motif_proteins:
				self.add_tag_to_instance(instance)
				
			for instance in domain_proteins:
				self.add_tag_to_instance(instance)
				
			response[accession] = {"status":"Success"} 

		return response

	def add_tag_to_instance(self,instance):
		for curation_id in instance['curation_ids']:
			self.options['curation_id'] = curation_id

			self.options['collapsed'] = False
			curation_instance_data = self.get_instances_by_curation_id()

			curation_instance = curation_instance_data[-1]

			self.options['collapsed'] = True
			self.options['api_path'] = '/restapi/curation/add/curation_version/' + str(curation_id)
			
			try:
				curation_instance["tags"] = list(set(curation_instance["tags"] + self.options["tags"]))
			except:
				logger.error("Error merging tags " + str(curation_instance["tags"]))
				curation_instance["tags"] = self.options["tags"]

			curation_instance["update_description"] = "added new tag"
			curation_instance["interaction_kd"] = {}

			self.options['params'] = json.dumps(curation_instance)
			self.run_update()
		
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	
	def add_pocket_class(self):
		
		self.options['api_path'] = '/restapi/motifdb/pocket/add'

		pocket_options ={
			"pocket_id": "STAG2_RAD21_pocket",
			"prototype_structure": "6QNX",
			"name": "STAG2-RAD21 composite pocket",
			"name_synonyms": [
				"YxF pocket"
			],
			"domain_ids": [
				""
			],
			"taxonomic_range": [
				"Metazoa",
				"Fungi"
			]
		}

		self.options['params'] = json.dumps(pocket_options)

		pprint.pprint(pocket_options)
		response = self.run_update()
		pprint.pprint(response)

	def delete_pocket_class(self):		
		self.options['api_path'] = '/restapi/motifdb/pocket/delete/' + self.option['pocket_acc']
		response = self.run_update()
		pprint.pprint(response)

	def add_specificity_class(self):
	
		class_options = {
			"review_description":"YxF cohesin binding pocket abstract review_description",
			"molecular_description":"*Cohesin* catalyses the folding of the genome into loops that are anchored by CTCF1. The molecular mechanism of how cohesin and CTCF structure the 3D genome has remained unclear. Here we show that a segment within the CTCF N terminus interacts with the SA2-SCC1 subunits of human cohesin. We report a crystal structure of *SA2-SCC1* in complex with CTCF at a resolution of 2.7 Å, which reveals the molecular basis of the interaction. We demonstrate that this interaction is specifically required for CTCF-anchored loops and contributes to the positioning of cohesin at CTCF binding sites. A similar motif is present in a number of established and newly identified cohesin ligands, including the cohesin release factor WAPL2,3. Our data suggest that CTCF enables the formation of chromatin loops by protecting cohesin against loop release. These results provide fundamental insights into the molecular mechanism that enables the dynamic regulation of chromatin folding by cohesin and CTCF.",
			"short_description":"YxF cohesin binding pocket **molecular** details short_description",
			"pocket_acc": "POCKET0000000163",
			"specificity_id": "STAG2_RAD21_YxF_pocket",
			"prototype_structure": "6QNX",
			"name": "YxF cohesin binding pocket",
			"name_synonyms": [
				"YxF motif"
			],
			"functional_class": [
				"Binding"
				#{Targetting}
				#{Modification}
				#{Trafficking}
				#{Binding,Targetting}
				#{Degradation}
				#{Docking}
				#{Cleavage}
				#{Binding}
			],
			#"tags": "string",
			"termini_binding": "False",
			"binding_functional_class": "Ligand",
			#"trafficking_from_location": "",
			#"trafficking_to_location": "string",
			"targetting_location": "Cohesin complex",
			#"docking_enzyme_type": "string",
			"consensus": "[YF].F",
			"taxonomic_range": [
				"Metazoa",
				"Fungi"
			],
			"taxonomic_range_missing": [
			],
			"core_regex": "[YF].F",
			"xrefs": [
				{}
			],
			"pmids": [
				"31905366"
			],
			"modifications": [
				{}
			]
		}
		self.options['specificity_id'] = 'STAG2_RAD21_YxF_pocket'
		response = self.get_specificity_class()
		
		if len(response) > 0:
			self.options['api_path'] = '/restapi/motifdb/specificity/update/' + response[0]['specificity_acc']
			
			for class_option in class_options:
				update = False
				class_options_update = {}

				if (class_option not in response[0]) and ("specificity_" + class_option not in response[0]): 
					class_options_update = {
						"column": class_option,
						"value": str(class_options[class_option])
					}
					update = True
				else:
					if class_option  in response[0]:
						response_value = response[0][class_option]
					if "specificity_" + class_option  in response[0]:
						response_value = response[0]["specificity_" + class_option]

					if response_value!= class_options[class_option]:
						class_options_update = {
							"column": class_option,
							"value": str(class_options[class_option])
						}
						update = True

				if update:
					logger.debug("Updating entry: " + response[0]['specificity_acc'])
					logger.debug(class_options_update)
					self.options['params'] = json.dumps(class_options_update)
					update_response = self.run_update()
					pprint.pprint(update_response)
		else:	
			self.options['api_path'] = '/restapi/motifdb/specificity/add'
			self.options['params'] = json.dumps(class_options)
			response = self.run_update()

	def delete_specificity_class(self):		
		self.options['api_path'] = '/restapi/motifdb/specificity/delete/' + self.options['specificity_acc']
		self.options['params'] = {}
		response = self.run_update()
		pprint.pprint(response)
		
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####

	def add_instance(self):
		"""
		data =  {
			"motif_protein_uniprot_accession": motif_protein_acc,
			"motif_start_residue_number": 399,
			"motif_end_residue_number": 405,
			"motif_sequence": "DSGLSMS",
			"source": "delete_test",
			"pmid": "20048001",
			"domain_start_residue_number": 1,
			"domain_end_residue_number": 210,
			"domain_protein_uniprot_accession": "Q9Y297",
			"domain_pfam_accession": "PF00400",
			"interaction_pdb": None,
			"interaction_type": "SLiM",
			"evidence_classes": None,
			"evidence_methods":  ['MI_0096', 'MI_0943', 'MI_0107', 'MI_0007', 'MI_0663', 'MI_0051', 'MI_0114'],
			"evidence_logic": None,
			"curator_instance_assessment": None,
			"conditional_motif": False,
			"conditional_motif_data": None,
			"other_notes":  None,
			"confidence":  None,
			"reliabilities": None,
			"elm_accession": None,
			"specificity_acc": "IDIC000019",
			"motif_instance_id": None,
			"interaction_kd": {"value":7.14,"min":None,"max":None,"sd":0.39},
			"throughput": "high", ## NEW
			"cell_cycle_data": {
				"cell_state_start":"G1",
				"cell_state_end": "S",
				"cell_state_range_start": 5,
				"cell_state_range_end": 10
				},
			"modification_data": [
				{
					"is_synthetic_motif":False,
					"name": "Phospho peptide",
					"protein_id": motif_protein_acc,
					"protein_start":None,
					"protein_end":None,
					"wildtype_sequence": None,
					"description":None,
					"interaction_kd_value":None,
					"interaction_kd_sd":None,
					"pdb":None,
					"modified_residues": [
						{
						"modification_type": "Phosphorylation",
						"modification_effect":"Activating",
						"modified_protein_position": 400,
						"modified_residue": None,
						"wildtype_residue": "S",
						"modyfing_enzyme": "P48730",
						"notes": "If not phosphorylated, interaction was abolish",
						"is_required": True,
						"ptm_id":None
						}
					]
				}
			]
		}
		"""

		self.options['api_path'] = '/restapi/curation/add/curation' 
		self.options['params'] = json.dumps(self.options['dataset_instance'])
		is_running = True
		
		while is_running:
			resp = self.run_update()
			
			if resp['status'] in ['Submitted']:
				print(resp)
			
			if resp['status'] in ['Finished','Not completed']:
				is_running = False
				print(resp)
				print(resp['curation_id'])
			elif resp['status'] in ['Error']:
				print("ERROR")
				print(resp)
				is_running = False
			else:
				print(resp)
				time.sleep(5)

		return resp

	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####
	#####------------------------------------------------------------------#####

