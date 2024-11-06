import os
import sys

import inspect

file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"../"))
import option_reader

sys.path.append(os.path.join(file_path,"../downloaders"))
import ensemblDownloader

from dataManager import dataManager


class ensemblManager(dataManager):

	##------------------------------------------------------------------##
	## Inherits functions from dataManager
	##------------------------------------------------------------------##

	def setup_data(self):
		self.default_task_options = {
			"task":'help',
			"sequence_information_type":"protein",
			"is_superuser":True,
			"accession":None,
			"verbose":False,
			"debug":False,
			"remake":False,
			"annotate_all_exact_matches":False,
			"logfile":False,
		}

		self.allowed_options = [
			"task",
			"outfile",
			"accession",
			"ensembl_gene_id",
			"ensembl_transcript_id",
			"gene_name",
			"chromosome",
			"region_start",
			"region_end",
			"feature_type",
			"flanks",
			"peptide_start",
			"peptide_end",
			"peptide_sequence",
			"oligo_sequence",
			"annotate_all_exact_matches",
			"exact_matches_only",
			"counts_only",
			"counts_output_type",
			"sequence_information_type",
			"sequence_length",
			"sequence_samples",
			"sequence_requirements",
			"sequence_flanks_3'",
			"sequence_flanks_5'",
			"dna_region_coordinates",
			"evalue",
			"remake",
			"detailed",
			"chunks_size",
			"use_filtering"
		]

		self.allowed_options_admin = [
		]

		self.task_options = [
			"get_ensembl_data_for_region",
			"get_ensembl_data_for_gene",
			"grab_ensembl_gene_id_by_genename",
			"grab_ensembl_gene_id_by_uniprot_xref",
			"grab_ensembl_transcript_id_by_genename_peptide",
			"grab_ensembl_transcript_id_by_genename",
			"grab_ensembl_transcript_id_by_uniprot_xref",
			"grab_exon_data",
			"grab_gene_genename",
			"grab_gene_transcripts",
			"grab_protein_exons",
			"grab_gene_transcript_ids",
			"grab_gene_transcript_information",
			"grab_transcript_gene_exons",
			"grab_transcript_protein_exons",
			"grab_transcript_protein",
			"grab_transcript_sequence_information",
			"grab_transcript_swissprot_xref",
			"grab_transcript_uniprot_xref",
			"grab_genomic_region_features",
			"find_sequence_human_genome",
			"grab_random_sequence",
			"grab_paralogues",
			"grab_ortholgues",
			"map_oligo_to_gene",
			"grab_dna_regions",
			"help"
		]

		self.task_options_admin = []

		self.required = {
			'all':[],
			'grab_genomic_region_features':["feature_type","chromosome","region_start","region_end"]
		}

		self.test_options = {
			'accession':"P20248",
			"peptide_start":10,
			"peptide_end":20,
			"peptide_sequence":"ATREAGSALLA",
			"flanks":1,
			"gene_name":"CCNA2",
			"ensembl_gene_id":"ENSG00000145386",
			"ensembl_protein_id":"ENSP00000274026",
			"ensembl_transcript_id":"ENST00000274026",
			"oligo_sequence":"CCCTAGCTCATACTCAG"
			}

		self.required_type = {
			'accession':"string",
			'remake':"bool"
		}

		self.convert_type = {
		}

		self.required_format = {
			'accession':"\A([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\-{0,1}[0-9]*\Z"
		}

		# delimiter to split options
		self.list_options = {
			'accession':",",
			'oligo_sequence':","
		}

		self.required_valid_options_list = {
			'remake':[True,False],
			'tasks':self.task_options,
			"find_sequence_human_genome":['oligo_sequence'],
			"map_oligo_to_gene":['oligo_sequence']
		}

		self.options = self.default_task_options

		for allowed_option in self.allowed_options + self.allowed_options_admin:
			if allowed_option not in self.options:
				self.options[allowed_option] = None

		self.options.update(option_reader.load_commandline_options(self.options,self.options))

	##------------------------------------------------------------------##
	##------------------------------------------------------------------##

	def getData(self,status,task_options):

		###----######----######----###
		#   run_job run_job run_job  #
		###----######----######----###

		data = {}

		dataDownloaderObj = ensemblDownloader.ensemblDownloader()
		dataDownloaderObj.options.update(task_options)


		if task_options["task"] == "get_ensembl_data_for_region":
			data = dataDownloaderObj.get_ensembl_data_for_region()
		if task_options["task"] == "get_ensembl_data_for_gene":
			data = dataDownloaderObj.get_ensembl_data_for_gene()
		if task_options["task"] == "grab_ensembl_gene_id_by_genename":
			data = dataDownloaderObj.grab_ensembl_gene_id_by_genename()
		if task_options["task"] == "grab_ensembl_gene_id_by_uniprot_xref":
			data = dataDownloaderObj.grab_ensembl_gene_id_by_uniprot_xref()
		if task_options["task"] == "grab_ensembl_transcript_id_by_genename_peptide":
			data = dataDownloaderObj.grab_ensembl_transcript_id_by_genename_peptide()
		if task_options["task"] == "grab_ensembl_transcript_id_by_genename":
			data = dataDownloaderObj.grab_ensembl_transcript_id_by_genename()
		if task_options["task"] == "grab_ensembl_transcript_id_by_uniprot_xref":
			data = dataDownloaderObj.grab_ensembl_transcript_id_by_uniprot_xref()
		if task_options["task"] == "grab_exon_data":
			data = dataDownloaderObj.grab_exon_data()
		if task_options["task"] == "grab_gene_genename":
			data = dataDownloaderObj.grab_gene_genename()
		if task_options["task"] == "grab_gene_transcripts":
			data = dataDownloaderObj.grab_gene_transcripts()
		if task_options["task"] == "grab_gene_transcript_information":
			data = dataDownloaderObj.grab_gene_transcript_information()
		if task_options["task"] == "grab_gene_transcript_ids":
			data = dataDownloaderObj.grab_gene_transcript_ids()
		if task_options["task"] == "grab_transcript_gene_exons":
			data = dataDownloaderObj.grab_transcript_gene_exons()
		if task_options["task"] == "grab_transcript_protein_exons":
			data = dataDownloaderObj.grab_transcript_protein_exons()
		if task_options["task"] == "grab_protein_exons":
			data = dataDownloaderObj.get_exon_protein_data()
		if task_options["task"] == "grab_transcript_protein":
			data = dataDownloaderObj.grab_transcript_protein()
		if task_options["task"] == "grab_transcript_sequence_information":
			data = dataDownloaderObj.grab_transcript_sequence_information()
		if task_options["task"] == "grab_transcript_swissprot_xref":
			data = dataDownloaderObj.grab_transcript_swissprot_xref()
		if task_options["task"] == "grab_transcript_uniprot_xref":
			data = dataDownloaderObj.grab_transcript_uniprot_xref()
		if task_options["task"] == "grab_genomic_region_features":
			data = dataDownloaderObj.grab_genomic_region_features()
		if task_options["task"] == "find_sequence_human_genome":
			data = dataDownloaderObj.find_sequence_human_genome()
		if task_options["task"] == "map_oligo_to_gene":
			data = dataDownloaderObj.map_oligo_to_gene()
		if task_options["task"] == "grab_paralogues":
			data = dataDownloaderObj.grab_paralogues()
		if task_options["task"] == "grab_ortholgues":
			data = dataDownloaderObj.grab_ortholgues()
		if task_options["task"] == "grab_random_sequence":
			data = dataDownloaderObj.grab_random_sequence()
		if task_options["task"] == "grab_dna_regions":
			data = dataDownloaderObj.grab_dna_regions()
			
			
		
		return data

##------------------------------------------------------------------##
##------------------------   END    END  ---------------------------##
##------------------------------------------------------------------##

if __name__ == "__main__":

	dataManagerObj = ensemblManager()
	dataManagerObj.main()
