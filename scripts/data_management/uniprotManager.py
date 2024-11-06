import os
import sys

import inspect

file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"../"))
import option_reader

sys.path.append(os.path.join(file_path,"../downloaders"))
import uniprotDownloader

from dataManager import dataManager

class uniprotManager(dataManager):

	##------------------------------------------------------------------##
	## Inherits functions from dataManager
	##------------------------------------------------------------------##

	def setup_data(self):
		self.default_task_options = {
			"task":'help',
			"is_superuser":True,
			"remake":False,
			"update":False,
			"verbose":False,
			"reviewed":False,
			"mapping_from":'UniProtKB_AC-ID',
			"mapping_to":"UniProtKB",
			"mapping_format":"json",
			"format":"list",
			"flank_length":0,
			"uniref_identity":0.5,
			"debug":False,
			"logfile":False,
			"include_version": False,
			"structure_required":False,
			"protein_existence": False
		}

		self.allowed_options = [
			"task",
			"update",
			'verbose',
			"outfile",
			"accession",
			"region_start",
			"region_end",
			"pdb_id",
			"taxon_id",
			"taxonomy",
			"identifier",
			"reviewed",
			"regions",
			"flank_length",
			"mapping_from",
			"mapping_to",
			"mapping_format",
			"format",
			"pmid",
			'parse_features',
			'parse_protein_region',
			'feature_type',
			'use_features',
			'parse_keywords',
			'parse_snps',
			'parse_attributes',
			'parse_disorder',
			'parse_isoforms',
			'parse_sequence',
			'parse_generic',
			'parse_db_xref',
			"remake",
			'parseCommandline',
			'uniref_identity',
			'structure_required',
			'keyword_id',
			"remake_age",
			"include_version",
			"protein_existence"
		]

		self.allowed_options_admin = [
		]

		self.task_options = [
			"parse_uniprot",
			"parse_uniprot_bulk",
			"parse_uniprot_fasta",
			"parse_uniprot_by_pdb",
			"parse_mobidb",
			"parse_sequence",
			"parse_uniprot_pfam",
			"parse_attributes",
			"parse_features",
			"parse_basic",
			"parse_name",
			"parse_domains",
			"parse_go",
			"parse_go_detailed",
			"parse_isoforms",
			"parse_organism_host",
			"parse_keywords",
			"parse_all_keywords",
			"parse_localisation",
			"parse_pdb",
			"parse_secondary_accessions",
			"parse_snps",
			"parse_region_annotation",
			"parse_region_sequence",
			"parse_region_features",
			"parse_secondary_structure",
			"parse_mutagenesis_site",
			"parse_diseases",
			"parse_db_xref",
			"parse_pathways",
			"parse_drugbank",
			"parse_ptms",
			"parse_regions_of_interest",
			"get_uniprot_mapping",
			"get_uniprot_mapping_details",
			"parse_taxon_identifier_information",
			"parse_taxon_identifier_information_from_uniprot",
			"parse_uniprot_accession_taxa",
			"parse_taxonomy_information",
			"parse_uniprot_uniref_members",
			"parse_uniprot_uniref_clustername",
			"parse_closest_reviewed_protein_uniref",
			"parse_ptm_vocabulary",
			"check_accession",
			"check_uniprot_entry",
			"check_uniprot_versions_taxa",
			"get_protein_topology",
			"get_protein_topology_string",
			'get_keyword_proteins',
			"parse_history",
			"help"
		]

		self.task_options_admin = [
			'check_protein_accession_files'
		]

		self.required = {
			'all':[],
			"parse_uniprot":["accession"],
			"parse_uniprot_bulk":["accession"],
			"parse_uniprot_by_pdb":["pdb_id"],
			"parse_mobidb":["accession"],
			"parse_sequence":["accession"],
			"parse_uniprot_pfam":["accession"],
			"parse_attributes":["accession"],
			"parse_features":["accession"],
			"parse_basic":["accession"],
			"parse_name":["accession"],
			"parse_domains":["accession"],
			"parse_go":["accession"],
			"parse_isoforms":["accession"],
			"parse_organism_host":["accession"],
			"parse_keywords":["accession"],
			"parse_all_keyword":[],
			"parse_localisation":["accession"],
			"parse_pdb":["accession"],
			"parse_go":["accession"],
			"parse_go_detailed":["accession"],
			"parse_secondary_accessions":["accession"],
			"parse_snps":["accession"],
			"parse_region_annotation":["accession","region_start","region_end"],
			"parse_region_sequence":["regions"],
			"parse_region_features":["regions"],
			"get_uniprot_mapping":["identifier"],
			"get_uniprot_mapping_details":["identifier"],
			"parse_taxon_identifier_information":["taxon_id"],
			"parse_taxon_identifier_information_from_uniprot":["taxon_id"],
			"parse_uniprot_accession_taxa":["taxon_id"],
			"parse_taxonomy_information":["taxonomy"],
			"parse_uniprot_uniref_members":["accession"],
			"parse_uniprot_uniref_clustername":["accession"],
			"parse_closest_reviewed_protein_uniref":["accession"],
			"parse_ptm_vocabulary":[],
			"check_accession":["accession"],
			"get_protein_topology":["accession"],
			"get_protein_topology_string":["accession"],
			'get_keyword_proteins':["keyword_id"],
			'check_protein_accession_files':["accession"],
			"parse_history": ["accession"],
			"check_uniprot_versions_taxa":["taxon_id"],
		}

		self.test_options = {
			'accession':"P06400",
			'regions':"P04637;10;34,P04600;100;200,P20248;1;156",
			'pdb_id':"2AST",
			'taxon_id':"7227",
			'taxonomy':"Insecta",
			'keyword_id':"KW-0195",
			'identifier':"CCNA2_HUMAN",
			"structure_required":False,
			"region_start":1,
			"region_end":100
		}

		self.required_type = {
			'accession':"string",
			'identifier':"string",
			'remake':"bool",
			'flank_length':"int"
		}

		self.required_format = {
			'accession':"\A([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\-{0,1}[0-9]*\Z"
		}

		# delimiter to split options
		self.list_options = {
			"identifier":",",
			"accession":",",
			"regions":",",
			"pdb_id":",",
			"use_features":","

		}

		self.required_valid_options_list = {
			'remake':[True,False],
			'tasks':self.task_options
		}

		self.options.update(self.default_task_options)

		for allowed_option in self.allowed_options + self.allowed_options_admin:
			if allowed_option not in self.options:
				self.options[allowed_option] = None

		if 'parseCommandline' not in self.options: self.options['parseCommandline'] = True
		if 'parseCommandline' in self.options:
			if self.options['parseCommandline'] == None:
				self.options['parseCommandline'] = True

		if self.options['parseCommandline']:
			self.options.update(option_reader.load_commandline_options(self.options,self.options))

		for k, v in self.default_task_options.items():
			if k not in self.options:continue
			if isinstance(v, bool):
				if not isinstance(self.options[k], bool):
					if str(self.options[k]).lower() in ['true','t']:
						self.options[k] = True
					else:
						self.options[k] = False



	##------------------------------------------------------------------##

	def getData(self,status,task_options):

		###----######----######----###
		#   run_job run_job run_job  #
		###----######----######----###

		data = {}

		dataDownloaderObj = uniprotDownloader.uniprotDownloader()
		dataDownloaderObj.options.update(task_options)

		if task_options["task"] == "parse_uniprot":
			data = dataDownloaderObj.parseUniProt(task_options['accession'][0])
		if task_options["task"] == "parse_uniprot_bulk":
			data = dataDownloaderObj.parseUniProtBulk(task_options['accession'])
		if task_options["task"] == "parse_uniprot_fasta":
			data = dataDownloaderObj.parseUniProtFasta(task_options['accession'][0])
		if task_options["task"] == "parse_uniprot_by_pdb":
			data = dataDownloaderObj.parseUniProtByPDB(task_options['pdb_id'][0])
		if task_options["task"] == "parse_mobidb":
			data = dataDownloaderObj.parseMobiDB(task_options['accession'][0])
		if task_options["task"] == "parse_sequence":
			data = dataDownloaderObj.parseSequence(task_options['accession'][0])
		if task_options["task"] == "parse_uniprot_pfam":
			data = dataDownloaderObj.parseUniProtPfam(task_options['accession'][0])
		if task_options["task"] == "parse_attributes":
			data = dataDownloaderObj.parseAttributes(task_options['accession'][0])
		if task_options["task"] == "parse_features":
			data = dataDownloaderObj.parseFeatures(task_options['accession'][0])
		if task_options["task"] == "parse_basic":
			data = dataDownloaderObj.parseBasic(task_options['accession'][0])
		if task_options["task"] == "parse_name":
			data = dataDownloaderObj.parseName(task_options['accession'][0])
		if task_options["task"] == "parse_domains":
			data = dataDownloaderObj.parseDomains(task_options['accession'][0])
		if task_options["task"] == "parse_go":
			data = dataDownloaderObj.parseGO(task_options['accession'][0])
		if task_options["task"] == "parse_go_detailed":
			data = dataDownloaderObj.parseGODetailed(task_options['accession'][0])
		if task_options["task"] == "parse_organism_host":
			data = dataDownloaderObj.parseOrganismHost(task_options['accession'][0])
		if task_options["task"] == "parse_isoforms":
			data = dataDownloaderObj.parseIsoforms(task_options['accession'][0])
		if task_options["task"] == "parse_keywords":
			data = dataDownloaderObj.parseKeywords(task_options['accession'][0])
		if task_options['task'] == 'parse_all_keywords':
			data = dataDownloaderObj.parseKeywordsDetails()
		if task_options["task"] == "parse_localisation":
			data = dataDownloaderObj.parseLocalisation(task_options['accession'][0])
		if task_options["task"] == "parse_pdb":
			data = dataDownloaderObj.parsePDBs(task_options['accession'][0])
		if task_options["task"] == "parse_secondary_accessions":
			data = dataDownloaderObj.parseSecondaryAccessions(task_options['accession'][0])
		if task_options["task"] == "parse_snps":
			data = dataDownloaderObj.parseSNPs(task_options['accession'][0])
		if task_options["task"] == "parse_region_sequence":
			data = dataDownloaderObj.parse_region_sequence(task_options['regions'],flank_length=task_options['flank_length'])
		if task_options["task"] == "parse_region_features":
			data = dataDownloaderObj.parse_region_features(task_options['regions'],flank_length=task_options['flank_length'])
		if task_options["task"] == "parse_region_annotation":
			data = dataDownloaderObj.parse_region_annotation()
		if task_options["task"] == "get_uniprot_mapping":
			data = dataDownloaderObj.get_uniprot_mapping(task_options['identifier'], mapping_from = task_options['mapping_from'], mapping_to = task_options['mapping_to'])
		if task_options["task"] == "get_uniprot_mapping_details":
			data = dataDownloaderObj.get_uniprot_mapping_details(task_options['identifier'],mapping_from=task_options["mapping_from"],mapping_to=task_options["mapping_to"],format=task_options["mapping_format"])
		if task_options["task"] == "parse_taxon_identifier_information":
			data = dataDownloaderObj.parse_taxon_identifier_information()
		if task_options['task'] == 'parse_taxon_identifier_information_from_uniprot':
			data = dataDownloaderObj.parse_taxon_identifier_information_from_uniprot()
		if task_options["task"] == "parse_taxonomy_information":
			data = dataDownloaderObj.parse_taxonomy_information()
		if task_options["task"] == "parse_secondary_structure":
			data = dataDownloaderObj.parseSecondaryStructure(task_options['accession'][0])
		if task_options["task"] == "parse_mutagenesis_site":
			data = dataDownloaderObj.parseMutagenesis(task_options['accession'][0])
		if task_options["task"] == "parse_ptms":
			data = dataDownloaderObj.parsePTMs(task_options['accession'][0])
		if task_options["task"] == "parse_regions_of_interest":
			data = dataDownloaderObj.parseRegionsOfInterest(task_options['accession'][0])
		if task_options["task"] == "parse_diseases":
			data = dataDownloaderObj.parseDiseases(task_options['accession'][0])
		if task_options["task"] == "parse_uniprot_accession_taxa":
			data = dataDownloaderObj.parse_uniprot_accession_taxa(task_options['taxon_id'],reviewed=task_options['reviewed'],structures=task_options['structure_required'], include_version = task_options['include_version'], protein_existence = task_options["protein_existence"],format=task_options['format'] )
		if task_options["task"] == "parse_uniprot_uniref_members":
			data = dataDownloaderObj.getUniProtUnirefMembers(task_options['accession'][0],task_options['uniref_identity'])
		if task_options["task"] == "parse_uniprot_uniref_clustername":
			data = dataDownloaderObj.getUniProtUnirefClusterName(task_options['accession'][0],task_options['uniref_identity'])
		if task_options["task"] == "parse_closest_reviewed_protein_uniref":
			data = dataDownloaderObj.getClosestReviewedProteinUniRef(task_options['accession'][0])
		if task_options["task"] == "check_accession":
			data = dataDownloaderObj.check_accession(task_options['accession'][0])
		if task_options["task"] == "check_uniprot_entry":
			data = dataDownloaderObj.checkUniProtEntry(task_options['accession'][0])
		if task_options["task"] == "parse_db_xref":
			data = dataDownloaderObj.parseDbXref(task_options['accession'][0])
		if task_options["task"] == "parse_pathways":
			data = dataDownloaderObj.parsePathways(task_options['accession'][0])
		if task_options["task"] == "parse_drugbank":
			data = dataDownloaderObj.parseDrugbank(task_options['accession'][0])
		if task_options['task'] == 'parse_ptm_vocabulary':
			data = dataDownloaderObj.parsePTMvocabulary()
		if task_options["task"] == "get_protein_topology":
			data = dataDownloaderObj.get_protein_topology(task_options['accession'][0])
		if task_options["task"] == "get_protein_topology_string":
			data = dataDownloaderObj.get_protein_topology_string(task_options['accession'][0])
		if task_options["task"] == "get_keyword_proteins":
			data = dataDownloaderObj.get_keyword_proteins()
		if task_options["task"] == "check_protein_accession_files":
			data = dataDownloaderObj.check_protein_accession_files()
		if task_options["task"] == "parse_history":
			data = dataDownloaderObj.parseProteinHistory(accession= task_options['accession'][0])
		if task_options["task"] == "check_uniprot_versions_taxa":
			data = dataDownloaderObj.check_uniprot_versions_taxa(task_options['taxon_id'])



		return data

	##------------------------------------------------------------------##
	##------------------------   END    END  ---------------------------##
	##------------------------------------------------------------------##

if __name__ == "__main__":

	dataManagerObj = uniprotManager()
	dataManagerObj.main()
