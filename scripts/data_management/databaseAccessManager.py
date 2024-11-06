import os
import sys

import inspect

file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"../"))
import option_reader

sys.path.append(os.path.join(file_path,"../database_access/"))

from dataManager import dataManager

### python motifManager.py --task test --data_sources viruses --filter_by "motif taxonomy" --filter_term Viruses --filter_logic contains --verbose True --task get_statistics  --data_type "motif species"  --statistics_detailed True

class databaseAccessManager(dataManager):

	##------------------------------------------------------------------##
	## Inherits functions from dataManager
	##------------------------------------------------------------------##

	def setup_data(self):
		self.default_rest_options = {
			"task":'help',
			"is_superuser":True,
			"collapsed":None,
			"remake":False,
			"verbose":False,
			"debug":False,
			"logfile":False,
			"server_url":'http://slim.icr.ac.uk'
		}

		self.allowed_options = [
			"task",
			"outfile",
			"username",
			"password",
			"database_name",
			"columns",
			'accession',
			'region_start',
			'region_end',
			'accessions',
			"additional_params",
			"server_url",
			#Database
   			"residue_centric",
			#Motifs
			"motif_database_source",
			"by_column",	
			"peptide_padding",
			"by_column_search_term",
			"collapsed",
			"data_sources",
			"group_by",
			"include_cols",
			"is_prototype_structure",
			"is_active_version",
			"additional_params",
			"data_sources",
			'elm_class',
			"specificity_class",
			"specificity_id",
			"specificity_acc",
			"pocket_id",
			"pocket_acc",
			'pfam_id',
			'domain_start',
			'domain_end',
			"group_by",
			"curation_id",
			'tags',		
			'dataset_instance',	
			"skip_feature_annotations",
			#Consensus
			'peptide',
			#Interactions
			'detection_methods',
			'use_STRING',
			"table"
			#GENERAL PEPTOOLS
			'toolname',
			'run_enrichment',
			"taxon_id",
			"search_type",
			"iupred_cutoff",
			"flank",
			"motif",
			"binding_domains",
			"binding_proteins",
			"masked_aa",
			"check_shared_domains_interactions", ## check onnly if PPI exists
			"check_shared_proteins_interactions", ## check onnly if PPI exists
			"check_shared_domains",## check shared functional annotataions
			"check_shared_proteins", ## check shared functional annotataions
			"shared_pval_cutoff",
			"cached", ## unless you want to rerun job
			"query_database",
			#PEPTOOLS
			"peptides",
			"offset",
			#PSSMSEARCH
			"pssm",
			"pssm_conservation_cutoff",
			"taxonomic_range",
			"scoring_method",
			"mask_cutoff",
			"pssm_pval_method",
			"pssm_pval_cutoff",
			"pssm_motif",
			"filter_by_regexp",
			"motif",
			"pssm_peptides",
			"hits_cutoff",
			#Datasets
			"dataset_id",
			"dataset_name",
			"convert_dataset_id",
			'group_id',
			"dataset_user",
			"dataset_instance_ids",
			"dataset_instance_data",
			"dataset_instance_file",
			"dataset_instance",
			#USERS
			"user_username",
			"user_password",
			#PSSM
			"pssm_id",
			"aln_id",
			"params",
			"username",
			"password"
		]

		self.allowed_options_admin = [
			"motif_instance"
		]

		self.task_options_admin = [
			]

		self.task_options = [
			#Database
			"get_features",
			"get_attributes",
			"get_architecture",
			#Motifs
			"get_instances_by_elm_class",
			#'get_instances_by_specificity_class',
			'get_instances_by_specificity_acc',
			'get_instances_by_motif_accession',
			'get_instances_by_domain_accession',
			'get_instances_by_pfam_id',
			'get_instances_by_binding_protein_region',
			'get_instances_by_curation_id',
			"get_prototype_structures",
			'get_alignment_by_elm_class',
			'get_alignment_by_specificity_acc',
			'get_pocket_class',
			'get_specificity_class',
			'get_specificity_class_by_specificity_acc',
			'get_specificity_class_pssm_data',
			'get_specificity_class_by_region',
			'get_specificity_class_xref',
			'get_specificity_class_reverse_xref',
			#Consensus
			'get_consensus_matches',
			'get_consensus_from_pssm',
			#Interactions
			'get_interactions',
			'get_interactions_by_accession',
			"run_peptools_offset_annotation",
			"run_peptools_offset_enrichment",
			#PepTools
			"run_protein_re_annotation",
			"run_protein_pssm_annotation",
		 	"run_protein_pssm_pssm_annotation",
			'run_peptools_enrichment',
			'run_peptools_annotation',
			"run_pssmsearch_annotation",
			"run_pssmsearch_enrichment",
			"run_slimsearch_annotation",
			"run_slimsearch_enrichment",
			"run_go_enrichment",
			#PSSM
			'browse_pssm_alignments',
			'browse_pssms',
		]

		self.task_options_admin = [
			#PSSM
			'submit_alignment',
			'remove_alignment',
			'make_pssm',
			'add_pssm',
			'remove_pssm',
			#Motifs
			'get_instances',
			"get_instances_without_specificity_class",
			"get_instances_without_pocket_class",
			"get_instances_with_structural_information",
			"get_instances_with_structures",
			"get_instances_pdb_pocket_fasta",
			'add_tag_by_accession',
			'add_specificity_class',
			'delete_specificity_class',
			'add_pocket_class',
			'delete_pocket_class',
			'add_instance',
			#USERS
			"get_authentication_token",
			"activate_user",
			"inactivate_user",
			"update_users_password",
			"get_user_info",
			"get_users_info",
			"delete_user",
			#Datasets
			'list_user_groups',
			'list_my_user_groups',
			'add_dataset',
			'add_dataset_instances',
			'get_dataset_instances',
			'delete_dataset',
			'get_dataset_details',
			'get_dataset_details_by_name',
			'add_dataset_from_file',
			'create_user_group',
			'delete_user_groups',  			
			'list_user_groups',
			'list_my_user_groups',
			'list_user_groups_by_group_id',
			'add_username_to_user_groups',
			'remove_username_from_user_groups'
		]

		self.task_skip_queue = [
			"help"
		]

		self.required = {
			'all':['database_name'],
			"get_features":['accession'],
			"get_attributes":['accession'],
			"get_architecture":['accession'],
			"add_tag_by_accession":['tags','accessions'],
			"get_instances":[],
			"add_instance":['dataset_instance'],
			'get_instances_by_elm_class':['elm_class'],
			'get_alignment_by_elm_class':['elm_class'],
			'get_alignment_by_specificity_acc':['specificity_acc'],
			'get_instances_by_specificity_class':['specificity_class'],
			'get_instances_by_specificity_acc':['specificity_acc'],
			'get_instances_by_motif_accession':['accession'],
			'get_instances_by_domain_accession':['accession'],
			'get_instances_by_binding_protein_region':['accession','domain_start','domain_end'],
			'get_instances_by_pfam_id':['pfam_id'],
			'get_instances_by_curation_id':['curation_id'],
			'get_consensus_matches':['peptide'],
			'get_interactions':['accession'],
			'get_interactions_by_accession':['accession'],
			'run_peptools_enrichment':['peptides','taxon_id'],
			'run_peptools_annotation':['peptides','taxon_id'],
			'run_peptools_offset_annotation':['offset','taxon_id'],
			'run_peptools_offset_enrichment':['offset','taxon_id'],
			'run_pssm_enrichment':['pssm','taxon_id'],
			'run_pssm_annotation':['pssm','taxon_id'],
			"run_protein_re_annotation":["query_database"],
			"run_protein_pssm_annotation":["query_database"],
		 	"run_protein_pssm_pssm_annotation":["query_database"],
			'delete_user_groups':['group_id'],	
			'list_user_groups_by_group_id':['group_id'],
			'remove_username_from_user_groups':['group_id','dataset_user'],
			"activate_user":["user_username"],
			"inactivate_user":["user_username"],
			"update_users_password":["user_username","user_password"],
			"get_authentication_token":['username','password'],
			"get_user_info":["user_username"],
			"get_users_info":[],
			"delete_user":["user_username"],
			'add_specificity_class':['specificity_class'],
			'get_specificity_class_pssm_data':['specificity_id'],
			'delete_specificity_class':['specificity_id'],
			'get_specificity_class':['specificity_id'],
			'get_specificity_class_by_specificity_acc':['specificity_acc'],
			'add_pocket_class':['pocket_class'],
			'delete_pocket_class':['pocket_id'],
			'get_pocket_class':['pocket_id'],
			'get_instances_by_specificity_acc':['specificity_acc'],
			'get_consensus_from_pssm':['pssm']
		}

		self.test_options = {
			'accession':"P20248",
			'pfam_id':"PF00400",
			'elm_class':'LIG_APCC_ABBA_1',
			'specificity_class':'IDIC000192',
			'specificity_acc':'IDIC000192',
			'specificity_id':'14-3-3_phosphopocket',
			'pocket_id':'IDIP000001',
			"database_name":"data_sources",
			"database_source":"motifs",
			"group_by":"motif",
			"domain_start":250,
			"domain_end":500,
			"database_name":"motifs",
			"curation_id":"25214,25540",
			"taxon_id":"9606",
			"peptide":'METGSLFTSGIKRHLKDKRISKT',
			"peptides":['METGSLFTSGIKRHLKDKRISKT'],
			"offset":[{"start":10,"end":20,"accession":"P20248"}]
		}

		self.required_type = {
			'data_sources':"string",
			'group_by':"string"
		}

		self.convert_type = {
		}

		self.required_format = {
			'accession':"\A([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\-{0,1}[0-9]*\Z",
			'accessions':"\A([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\-{0,1}[0-9]*\Z",
			'detection_methods':"\AMI_[0-9]{4}\Z"
		}

		# delimiter to split options
		self.list_options = {
			"data_sources":",",
			"accessions":",",
			"filter_by":",",
			"data_type":",",
			'detection_methods':",",
			'peptides':",",
			'tags':",",
			'dataset_instance_ids':","
		}

		self.required_valid_options_list = {
			'remake':[True,False],
			'group_by':['motif','domain_proteins',"pdb"],
			'database_name':['users','motifs','interactions','peptools','datasets','consensus',"pssm","database"],
			'by_column':[
				"specificity_acc",
				"elm_accession",
				"motif_class",
				"motif_taxonomic_range",
				"motif_protein_uniprot_accession" ,
				"motif_protein_taxon_identifier",
				"motif_protein_family",
				"domain_protein_uniprot_accession" ,
				"domain_protein_family", 
				"domain_pfam_accession",
				"clan_id",
				"curation_id",
				"username"
			],
			'tasks':self.task_options
		}

		self.options = self.default_rest_options

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

		if task_options["database_name"] == "users":
			import utilities_users

			dataDownloaderObj = utilities_users.UserRestSubmission()
			dataDownloaderObj.options.update(task_options)
	
			if task_options["task"] == "activate_user":
				data = dataDownloaderObj.activateUser()
			
			if task_options["task"] == "inactivate_user":
				data = dataDownloaderObj.inactivateUser()

			if task_options["task"] == "update_users_password":
				data = dataDownloaderObj.updateUsersPassword()

			if task_options["task"] == "get_user_info":
				data = dataDownloaderObj.getUserInfo()

			if task_options["task"] == "get_users_info":
				data = dataDownloaderObj.getUsers()

			if task_options["task"] == "delete_user":
				data = dataDownloaderObj.deleteUser()

			if task_options["task"] == "get_authentication_token":
				data = dataDownloaderObj.get_authentication_token()

		if task_options["database_name"] == "database":
			import database_query_rest
			dataDownloaderObj = database_query_rest.databaseQueryRestManager()
			dataDownloaderObj.options.update(task_options)

			if task_options["task"] == "get_features":
				data = dataDownloaderObj.get_features()
			
			if task_options["task"] == "get_attributes":
				data = dataDownloaderObj.get_attributes()

		if task_options["database_name"] == "pssm":
			import pssm_query_rest
			dataDownloaderObj = pssm_query_rest.pssmQueryRestManager()
			dataDownloaderObj.options.update(task_options)
			
			if task_options["task"] == "browse_pssm_alignments":
				data = dataDownloaderObj.browse_pssm_alignments()
			if task_options["task"] == "browse_pssms":
				data = dataDownloaderObj.browse_pssms()
			if task_options["task"] == "make_pssm":
				data = dataDownloaderObj.make_pssm()
			if task_options["task"] == "add_pssm":
				data = dataDownloaderObj.add_pssm()
			if task_options["task"] == "remove_pssm":
				data = dataDownloaderObj.remove_pssm()
			if task_options["task"] == "submit_alignment":
				data = dataDownloaderObj.submit_alignment()
			if task_options["task"] == "remove_alignment":
				data = dataDownloaderObj.remove_alignment()
			
		if task_options["database_name"] == "motifs":
			import motif_query_rest

			dataDownloaderObj = motif_query_rest.motifQueryRestManager()
			dataDownloaderObj.options.update(task_options)
			
			if task_options["task"] == "get_instances":
				data = dataDownloaderObj.get_instances()
			elif task_options["task"] == "get_instances_by_elm_class":
				data = dataDownloaderObj.get_instances_by_elm_class()
			elif task_options["task"] == "get_instances_by_specificity_class":
				data = dataDownloaderObj.get_instances_by_specificity_class()
			elif task_options["task"] == "get_instances_by_specificity_acc":
				data = dataDownloaderObj.get_instances_by_specificity_acc()
			elif task_options["task"] == "get_instances_by_motif_accession":
				data = dataDownloaderObj.get_instances_by_motif_accession()
			elif task_options["task"] == "get_instances_by_domain_accession":
				data = dataDownloaderObj.get_instances_by_domain_accession()
			elif task_options["task"] == "get_instances_by_binding_protein_region":
				data = dataDownloaderObj.get_instances_by_binding_protein_region()
			elif task_options["task"] == "get_instances_by_pfam_id":
				data = dataDownloaderObj.get_instances_by_pfam_id()
			elif task_options["task"] == "get_instances_by_curation_id":
				data = dataDownloaderObj.get_instances_by_curation_id()
			elif task_options["task"] == "get_instances_without_specificity_class":
				data = dataDownloaderObj.get_instances_without_specificity_class()
			elif task_options["task"] == "get_instances_without_pocket_class":
				data = dataDownloaderObj.get_instances_without_pocket_class()
			elif task_options["task"] == "get_instances_with_structural_information":
				data = dataDownloaderObj.get_instances_with_structural_information()
			elif task_options["task"] == "get_prototype_structures":
				data = dataDownloaderObj.get_prototype_structures()
			elif task_options["task"] == "get_instances_with_structures":
				data = dataDownloaderObj.get_instances_with_structures()
			elif task_options["task"] == "get_instances_pdb_pocket_fasta":
				data = dataDownloaderObj.get_instances_pdb_pocket_fasta()
			elif task_options["task"] == "add_tag_by_accession":
				data = dataDownloaderObj.add_tag_by_accession()		
			elif task_options["task"] == "get_alignment_by_elm_class":
				data = dataDownloaderObj.get_alignment_by_elm_class()	
			elif task_options["task"] == "get_alignment_by_specificity_acc":
				data = dataDownloaderObj.get_alignment_by_specificity_acc()	
			elif task_options["task"] == "add_pocket_class":
				data = dataDownloaderObj.add_pocket_class()
			elif task_options["task"] == "delete_pocket_class":
				data = dataDownloaderObj.delete_pocket_class()
			elif task_options["task"] == "add_specificity_class":
				data = dataDownloaderObj.add_specificity_class()	
			elif task_options["task"] == "delete_specificity_class":
				data = dataDownloaderObj.delete_specificity_class()		
			elif task_options["task"] == "get_specificity_class":
				data = dataDownloaderObj.get_specificity_class()
			elif task_options["task"] == "get_specificity_class_by_specificity_acc":
				data = dataDownloaderObj.get_specificity_class_by_specificity_acc()		
			elif task_options["task"] == "get_pocket_class":
				data = dataDownloaderObj.get_pocket_class()		
			elif task_options["task"] == "get_specificity_class_xref":
				data = dataDownloaderObj.get_specificity_class_xref()	
			elif task_options["task"] == "get_specificity_class_reverse_xref":
				data = dataDownloaderObj.get_specificity_class_reverse_xref()		
			elif task_options["task"] == "get_specificity_class_by_region":
				data = dataDownloaderObj.get_specificity_class_by_region()		
			elif task_options["task"] == "get_specificity_class_pssm_data":
				data = dataDownloaderObj.get_specificity_class_pssm_data()	
			elif task_options["task"] == "add_instance":
				data = dataDownloaderObj.add_instance()	
				
			del dataDownloaderObj
			
		if task_options["database_name"] == "consensus":
			import consensus_query_rest

			dataDownloaderObj = consensus_query_rest.consensusQueryRestManager()
			dataDownloaderObj.options.update(task_options)

			if task_options["task"] == "get_consensus_matches":
				data = dataDownloaderObj.get_consensus_matches()

			if task_options["task"] == "get_consensus_from_pssm":
				data = dataDownloaderObj.get_consensus_from_pssm()

		if task_options["database_name"] == "interactions":
			import interaction_query_rest

			dataDownloaderObj = interaction_query_rest.interactionQueryRestManager()
			dataDownloaderObj.options.update(task_options)

			if task_options["task"] == "get_interactions":
				data = dataDownloaderObj.get_interactions()
			
			if task_options["task"] == "get_interactions_by_accession":
				data = dataDownloaderObj.get_interactions_by_accession()

		if task_options["database_name"] == "peptools":
			import peptools_query_rest

			dataDownloaderObj = peptools_query_rest.peptoolsQueryRestManager()
			dataDownloaderObj.options.update(task_options)

			
			if task_options["task"] == "run_protein_re_annotation":
				data = dataDownloaderObj.run_protein_re_annotation()

			if task_options["task"] == "run_protein_pssm_annotation":
				data = dataDownloaderObj.run_protein_pssm_annotation()

			if task_options["task"] == "run_protein_pssm_pssm_annotation":
				data = dataDownloaderObj.run_protein_pssm_pssm_annotation()

			if task_options["task"] == "run_peptools_annotation":
				data = dataDownloaderObj.run_peptools_annotation()
	
			if task_options["task"] == "run_peptools_enrichment":
				data = dataDownloaderObj.run_peptools_enrichment()

			if task_options["task"] == "run_pssmsearch_annotation":
				data = dataDownloaderObj.run_pssmsearch_annotation()

			if task_options["task"] == "run_pssmsearch_enrichment":
				data = dataDownloaderObj.run_pssmsearch_enrichment()
			
			if task_options["task"] == "run_slimsearch_annotation":
				data = dataDownloaderObj.run_slimsearch_annotation()

			if task_options["task"] == "run_slimsearch_enrichment":
				data = dataDownloaderObj.run_slimsearch_enrichment()

			if task_options["task"] == "run_peptools_offset_annotation":
				data = dataDownloaderObj.run_peptools_offset_annotation()

			if task_options["task"] == "run_peptools_offset_enrichment":
				data = dataDownloaderObj.run_peptools_offset_enrichment()

			if task_options["task"] == "run_go_enrichment":
				data = dataDownloaderObj.run_go_enrichment()
		
		if task_options["database_name"] == "datasets":
			import dataset_query_rest

			dataDownloaderObj = dataset_query_rest.datasetQueryRestManager()
			dataDownloaderObj.options.update(task_options)

			if task_options["task"] == "list_user_groups":
				data = dataDownloaderObj.list_user_groups()
			
			if task_options["task"] == "list_my_user_groups":
				data = dataDownloaderObj.list_my_user_groups()

			if task_options["task"] == "add_dataset":
				data = dataDownloaderObj.add_dataset()

			if task_options["task"] == "delete_dataset":
				data = dataDownloaderObj.delete_dataset()
			
			if task_options["task"] == "get_dataset_details":
				data = dataDownloaderObj.get_dataset_details()

			if task_options["task"] == "get_dataset_details_by_name":
				data = dataDownloaderObj.get_dataset_details_by_name()
				
			if task_options["task"] == "add_dataset_from_file":
				data = dataDownloaderObj.add_dataset_from_file()
			
			if task_options["task"] == "add_dataset_instances":
				data = dataDownloaderObj.add_dataset_instances()

			if task_options["task"] == "get_dataset_instances":
				data = dataDownloaderObj.get_dataset_instances()
				
			if task_options["task"] == "create_user_group":
				data = dataDownloaderObj.create_user_group()

			if task_options["task"] == "delete_user_groups":
				data = dataDownloaderObj.delete_user_groups()

			if task_options["task"] == "list_user_groups_by_group_id":
				data = dataDownloaderObj.list_user_groups_by_group_id()

			if task_options["task"] == "add_username_to_user_groups":
				data = dataDownloaderObj.add_username_to_user_groups()

			if task_options["task"] == "remove_username_from_user_groups":
				data = dataDownloaderObj.remove_username_from_user_groups()

		return data
		
##------------------------------------------------------------------##
##------------------------   END    END  ---------------------------##
##------------------------------------------------------------------##

if __name__ == "__main__":
	dataManagerObj = databaseAccessManager()
	dataManagerObj.main()
