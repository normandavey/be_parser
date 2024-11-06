import os, sys, inspect, json, pprint, hashlib, string, subprocess, time, copy,random,shutil

from Bio import pairwise2
				
file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../"))
import config_reader
import option_reader

sys.path.append(os.path.join(file_path,"../utilities"))
import utilities_downloader
import utilities_error
import utilities_codons
import utilities_basic
import utilities_alignment

sys.path.append(os.path.join(file_path,"../data_management/"))
import queryRunner

#-----

import zipfile

#-----
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
#-----

class ensemblDownloader():
	##------------------------------------------------------------------##

	def __init__(self):
		self.options = {
			"remake":False,
			"species":"homo_sapiens"
		}

		self.options.update(config_reader.load_configeration_options(sections=["general"]))
		self.options.update(option_reader.load_commandline_options(self.options,{}))
		
		self.options['annotate_isoforms'] = True
		self.options['all_isoforms'] = True
		self.options['annotate_isoforms'] = True
		self.options['get_protein_data'] = True
		self.options['get_exon_data'] = False
		self.options['exact_matches_only'] = False
		self.options['transcript_centric'] = False
		self.options['use_filtering'] = True
		self.options['remove_scaffolds'] = False
		self.options['counts_only'] = False
		self.options['counts_output_type'] = "numeric"
		self.options['evalue'] = "1000"
		self.options['remake_age'] = 1800
		self.options['chunks_size'] = 100
		self.options['sequence_length'] = 20
		self.options['sequence_samples'] = 1
		self.options['sequence_requirements'] = {}
		
		self.paths = {}
		self.paths['ensembl_path'] = os.path.join(self.options['data_path'],"ensembl")
		self.paths['exons_path'] = os.path.join(self.paths['ensembl_path'],"exons")
		self.paths['homology_path'] = os.path.join(self.paths['ensembl_path'],"homology")
		self.paths['cds_genomic_alignment_path'] = os.path.join(self.paths['ensembl_path'],"cds_genomic_alignment")
		self.paths['dna_region_coordinates_path'] = os.path.join(self.paths['ensembl_path'],"dna_region_coordinates")
		self.paths['gene_id_path'] = os.path.join(self.paths['ensembl_path'],"gene_id")
		self.paths['ensembl_protein_path'] = os.path.join(self.paths['ensembl_path'],"ensembl_protein")
		self.paths['ensembl_protein_exons_path'] = os.path.join(self.paths['ensembl_path'],"ensembl_protein_exons")
		self.paths['ensembl_gene_exons_path'] = os.path.join(self.paths['ensembl_path'],"ensembl_gene_exons")
		self.paths['ensembl_gene_path'] = os.path.join(self.paths['ensembl_path'],"ensembl_gene")
		self.paths['gene_genename_path'] = os.path.join(self.paths['ensembl_path'],"gene_name")
		self.paths['uniprot_xref_path'] = os.path.join(self.paths['ensembl_path'],"uniprot_xref")
		self.paths['sequence_information_path'] = os.path.join(self.paths['ensembl_path'],"sequence_information")
		
		self.paths['genomes_dir_path'] = os.path.join(self.options['data_path'],"genomes")
		self.paths['genomes_region_path'] = os.path.join(self.options['data_path'],"genomes","regions")
		self.paths['genomes_guides_path'] = os.path.join(self.options['data_path'],"genomes","gRNA")
		self.paths['genomes_guides_mapping_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","gRNA_mapping")
		self.paths['genomes_guides_fasta_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","gRNA_fasta")
		self.paths['genomes_guides_blast_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","blast")
		self.paths['genomes_guides_blast_counts_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","blast","counts")
		self.paths['genomes_guides_blast_raw_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","blast","raw")
		self.paths['genomes_guides_blast_coordinates_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","blast","counts_coordinates")
		self.paths['genome_path'] = os.path.join(self.options['data_path'],"genomes","Homo_sapiens.GRCh38.dna.primary_assembly.collapsed.fa")

		self.options.update(self.paths)

		for path in self.paths:
			if path in ['genome_path']: continue
			if not os.path.exists(self.paths[path]):
				os.mkdir(self.paths[path])
 
	######### paralogues orthologues homology ###############
	######### paralogues orthologues homology ###############
	######### paralogues orthologues homology ###############
	######### paralogues orthologues homology ###############

	def grab_homologues(self,homology_type="paralogues"):
		sequence_information_type = "condensed"
		url = "https://rest.ensembl.org/homology/id/" + self.options['ensembl_gene_id'].split(".")[0] + "?type=" + homology_type + ";format=condensed;content-type=application/json" 
		out_path = os.path.join(self.options['ensembl_path'],'homology',self.options['accession'][0] + "." + homology_type + ".json")

		if 'sequence_information_type' in self.options:
			if self.options['sequence_information_type'] == "full":
				sequence_information_type = "full"
				url = "https://rest.ensembl.org/homology/id/" + self.options['ensembl_gene_id'].split(".")[0] + "?type=" + homology_type + ";format=full;content-type=application/json" 
				out_path = os.path.join(self.options['ensembl_path'],'homology',self.options['accession'][0] + "." + homology_type + ".full.json")
			
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'],headers={ "Accept":"application/json"})
		
		with open(out_path) as outfile:
			json_data = json.load(outfile)

		if 'homologies' in json_data[0]:
			homologue_data = {}
			similarity = {}
			for homologue in json_data[0]['homologies']:
				if sequence_information_type == "condensed":
					self.options['ensembl_gene_id'] = homologue['id']
				else:
					self.options['ensembl_gene_id'] = homologue['target']['id']
					import utilities_alignment
					similarity = utilities_alignment.sequenceSimilarityScores(homologue['source']['align_seq'],homologue['target']['align_seq'])
					
				gene_uniprot_xrefs = self.grab_gene_uniprot_xref()
				
				for gene_uniprot_xref in gene_uniprot_xrefs:
					response_db_xref = queryRunner.queryRunner("uniprot","parse_basic",{"accession":gene_uniprot_xref['primary_id']}).run()
					if response_db_xref['data']['reviewed']:
						homologue_data[gene_uniprot_xref['primary_id']] = gene_uniprot_xref

						if sequence_information_type == "full":
							homologue_data[gene_uniprot_xref['primary_id']].update(similarity)

			return homologue_data
		else:
			return {}

	def grab_paralogues(self):
		ensembl_gene_ids = self.grab_ensembl_gene_id_by_uniprot_xref(as_list=True)
		self.options['ensembl_gene_id']  = ensembl_gene_ids[0].split(".")[0]
		homologue_data = self.grab_homologues(homology_type="paralogues")
		return homologue_data

	def grab_orthologues(self):
		ensembl_gene_ids = self.grab_ensembl_gene_id_by_uniprot_xref(as_list=True)
		self.options['ensembl_gene_id']  = ensembl_gene_ids[0]
		homologue_data = self.grab_homologues(homology_type="orthologues")
		return homologue_data
 
	######### START EXON FUNCTIONS ###############
	######### START EXON FUNCTIONS ###############
	######### START EXON FUNCTIONS ###############
	######### START EXON FUNCTIONS ###############

	def grab_exon_data(self):
		url = "https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&accession=" + self.options['accession'][0]
		out_path = os.path.join(self.options['ensembl_path'],'exons',self.options['accession'][0] + ".exons.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'],headers={ "Accept" : "application/json"})
		
		if os.path.exists(out_path):
			with open(out_path) as outfile:
				json_data = json.load(outfile)

			return json_data
		else:
			logging.error(status)
			return {}
	

	def grab_protein_gene_id(self):
		url = "https://www.ebi.ac.uk/proteins/api/genecentric/"  + self.options['accession'][0]
		out_path = os.path.join(self.options['ensembl_path'],'gene_id',self.options['accession'][0] + ".gene_id.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'],headers={ "Accept" : "application/json"})
		
		if os.path.exists(out_path):
			with open(out_path) as outfile:
				json_data = json.load(outfile)
		else:
			logging.error(status)
			return {}
	
		gene_id = []

		if 'relatedGene' in json_data:
			for relatedGene in json_data['relatedGene']:
				if 'geneNameType' in relatedGene:
					if relatedGene['geneNameType'] == 'Ensembl':
						if relatedGene['geneName'] not in gene_id:
							gene_id.append(relatedGene['geneName'])

		return gene_id


	def grab_exon_details(self):
		self.options['accession'] = self.options['accession'][0]
		response_db_xref = queryRunner.queryRunner("uniprot","parse_db_xref",{"accession":self.options['accession'].split("-")[0]}).run()
		ensembl_gene_ids = []
	
		transcript_ensembl_exon_data = {}

		if 'Ensembl' in response_db_xref['data']:
			for transcript_id in response_db_xref['data']['Ensembl']:
				if response_db_xref['data']['Ensembl'][transcript_id]['gene ID'] not in ensembl_gene_ids:
					ensembl_gene_ids.append(response_db_xref['data']['Ensembl'][transcript_id]['gene ID'])

		if len(ensembl_gene_ids) > 0:
			for ensembl_gene_id in ensembl_gene_ids:
				
				logger.debug(ensembl_gene_id)
				self.options['ensembl_gene_id'] = ensembl_gene_id

				transcript_ensembl_protein_ids = self.grab_gene_transcript_ids()
				transcript_gene_exons = self.grab_transcript_gene_exons()

				for transcript_gene_exon in transcript_gene_exons:
					transcript_ensembl_exon_data[transcript_gene_exon['exon_id']] = {
						"in_transcript":False
					}
					
					for tag in ['start','end','ensembl_end_phase','ensembl_phase','constitutive']:
						if tag in ['start','end']:
							transcript_ensembl_exon_data[transcript_gene_exon['exon_id']]["genomic_" + tag] = transcript_gene_exon[tag] 
						else:
							transcript_ensembl_exon_data[transcript_gene_exon['exon_id']][tag] = transcript_gene_exon[tag] 

		return transcript_ensembl_exon_data


	def process_isoform_exon_data(self):

		uniprot_isoform_exon_match_found = False
		
		try:
			response_db_xref = queryRunner.queryRunner("uniprot","parse_db_xref",{"accession":self.options['accession'][0].split("-")[0]}).run()
			query_sequence = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":self.options['accession']}).run()['data']

			ensembl_gene_ids = []
			protein_list = []
			isoforms_data = {}

			transcript_ensembl_protein_data = {}
			transcript_ensembl_exon_data = {}
			transcript_uniprot_mapping = {}

			ensembl_gene_ids = self.grab_ensembl_gene_id_by_genename()

			if len(ensembl_gene_ids) > 0:
				for ensembl_gene_id in ensembl_gene_ids:
					logger.info("> "*10 + ensembl_gene_id)
					self.options['ensembl_gene_id'] = ensembl_gene_id
					transcript_ensembl_protein_ids = self.grab_gene_transcript_ids()
					transcript_gene_exons = self.grab_transcript_gene_exons()

					#-----------------------------------------------------------------------------------------#
					
					for transcript_gene_exon in transcript_gene_exons:
						try:
							self.options['dna_region_coordinates'] =  ':'.join([transcript_gene_exon['seq_region_name'],str(transcript_gene_exon['start']),str(transcript_gene_exon['end']),str(transcript_gene_exon['strand'])]) 
							
							transcript_ensembl_exon_data[transcript_gene_exon['exon_id']] = {
								"in_transcript":False
							}
							grab_dna_regions = self.grab_dna_regions()

							transcript_ensembl_exon_data[transcript_gene_exon['exon_id']]["dna_sequence"] = grab_dna_regions[0]['seq'] 
							
							for tag in ['start','end','ensembl_end_phase','ensembl_phase']:
								if tag in ['start','end']:
									try:	
										transcript_ensembl_exon_data[transcript_gene_exon['exon_id']]["genomic_" + tag] = transcript_gene_exon[tag] 
									except:
										transcript_ensembl_exon_data[transcript_gene_exon['exon_id']]["genomic_" + tag] = "-"
								else:
									try:	
										transcript_ensembl_exon_data[transcript_gene_exon['exon_id']][tag] = transcript_gene_exon[tag] 
									except:
										transcript_ensembl_exon_data[transcript_gene_exon['exon_id']][tag] = "-"
						except:
							logging.error("Exon not added:" + transcript_gene_exon['exon_id'] + " " + str(transcript_gene_exon))
							logging.error(utilities_error.getError())

					#-----------------------------------------------------------------------------------------#

					options = {
						"identifier":transcript_ensembl_protein_ids,
						"mapping_from":"Ensembl_Transcript",
						"mapping_to":"ACC"
					}

					response_mapping = queryRunner.queryRunner("uniprot","get_uniprot_mapping_details",options).run()
					
					#-----------------------------------------------------------------------------------------#
					
					for transcript_ensembl_protein_id in transcript_ensembl_protein_ids:
						self.options['ensembl_transcript_id'] = transcript_ensembl_protein_id
						transcript_protein = self.grab_transcript_protein()

						if 'id' not in transcript_protein:
							logging.error("No identifier in " + transcript_ensembl_protein_id)
							continue
						
						self.options['ensembl_protein_id'] = transcript_protein['id']
						transcript_protein_exons = self.grab_transcript_protein_exons()
						
						uniprot_accession = ""
						
						if 'Ensembl' in response_db_xref['data']:
							if transcript_ensembl_protein_id in response_db_xref['data']['Ensembl']:
								if 'molecule' in response_db_xref['data']['Ensembl'][transcript_ensembl_protein_id]:
									uniprot_accession = response_db_xref['data']['Ensembl'][transcript_ensembl_protein_id]['molecule']
								
						if uniprot_accession == "":
							transcript_protein_uniprot_xref = self.grab_transcript_swissprot_xref()
						
							if len(transcript_protein_uniprot_xref) != 0:
								uniprot_accession = transcript_protein_uniprot_xref[0]['primary_id']

						if uniprot_accession == "":
							transcript_protein_uniprot_xref = self.grab_transcript_uniprot_xref()
							if len(transcript_protein_uniprot_xref) != 0:
								uniprot_accession = transcript_protein_uniprot_xref[0]['primary_id']
	
						if 'data' in response_mapping:
							if transcript_ensembl_protein_id in response_mapping['data']:
								if 'Entry' in response_mapping['data'][transcript_ensembl_protein_id]:
									uniprot_accession = response_mapping['data'][transcript_ensembl_protein_id]['Entry'] 
						
						protein_list.append(uniprot_accession)

						if uniprot_isoform_exon_match_found == False:
							uniprot_isoform_exon_match_found = transcript_protein['seq'] == query_sequence

						if uniprot_accession == "":
							logging.error("Cannot map accession for transcript - setting as query " + ensembl_gene_id + " " + transcript_ensembl_protein_id + " " + self.options['accession'][0])
							uniprot_accession = self.options['accession'][0]

						if uniprot_accession == self.options['accession'][0] or self.options['all_isoforms'] or (uniprot_accession.split("-")[-1] == "1" and self.options['accession'] == uniprot_accession.split("-")[0]):
							response_sequence = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":uniprot_accession}).run()

							if 'updated_accession' in response_sequence:
								response_sequence = queryRunner.queryRunner("uniprot","parse_sequence",{"accession":response_sequence['updated_accession']}).run()

							try:
								if self.options['debug'] == True:
									if response_sequence['data'] != transcript_protein['seq']: 
										
										logger.debug(uniprot_accession + " - " + transcript_ensembl_protein_id + " does not match uniprot sequence")
										
										if 1 in [response_sequence['data'].count(transcript_protein['seq']),transcript_protein['seq'].count(response_sequence['data'])]:
											logger.debug(uniprot_accession + " ->" + transcript_ensembl_protein_id + " - Uniprot extension")
											alignment = utilities_alignment.alignPeptides(response_sequence['data'],transcript_protein['seq'])
											"""
											print(transcript_ensembl_protein_id)
											print("".join(alignment[0][0]))
											print("".join(alignment[0][1]))
											print(query_sequence)
											#print([alignment[0][0][i] == alignment[0][1][i] for i in range(0,len(alignment[0][0]))].count("True")/(len(alignment[0][0]) - alignment[0][0].count("-")))
											#pass
											"""
										else:	
											logger.debug(uniprot_accession + " ->" + transcript_ensembl_protein_id + " - Uniprot mismatch")
											alignment = utilities_alignment.alignPeptides(response_sequence['data'],transcript_protein['seq'])
											"""
											print(transcript_ensembl_protein_id)
											print("".join(alignment[0][0]))
											print("".join(alignment[0][1]))
											print(query_sequence)
											"""
											#print([alignment[0][0][i] == alignment[0][1][i] for i in range(0,len(alignment[0][0]))].count("True"))
											#print([alignment[0][0][i] == alignment[0][1][i] for i in range(0,len(alignment[0][0]))].count("True")/(len(alignment[0][0]) - alignment[0][0].count("-")))
										
									else:
										logger.debug(uniprot_accession + " ->" + transcript_ensembl_protein_id + " - Exact " + uniprot_accession)
										"""
										print(transcript_ensembl_protein_id)
										print(response_sequence['data'])
										print(transcript_protein['seq'])
										print(query_sequence)
										"""
							except:
								logger.debug(uniprot_accession + " - " + transcript_ensembl_protein_id + " does not return data")
								#pprint.pprint(response_sequence)
								#continue

							transcript_ensembl_protein_id_info = {
								"ensembl_gene_id":ensembl_gene_id,
								"transcript_ensembl_protein_id":transcript_ensembl_protein_id,
								"uniprot_accession":uniprot_accession,
								"transcript_protein_exons":transcript_protein_exons,
								"transcript_sequence":transcript_protein['seq']
							}
							
							if self.options['ensembl_transcript_id'] in transcript_uniprot_mapping:
								try:
									if transcript_uniprot_mapping[self.options['ensembl_transcript_id']]['molecule'][-2:] == "-1":
										transcript_ensembl_protein_id_info["uniprot_accession_isoform"] = transcript_uniprot_mapping[self.options['ensembl_transcript_id']]['molecule'][:-2]
									else:
										transcript_ensembl_protein_id_info["uniprot_accession_isoform"] = transcript_uniprot_mapping[self.options['ensembl_transcript_id']]['molecule']

									uniprot_accession = transcript_ensembl_protein_id_info["uniprot_accession_isoform"]
								except:
									transcript_ensembl_protein_id_info["uniprot_accession_isoform"] = uniprot_accession
							else:
								transcript_ensembl_protein_id_info["uniprot_accession_isoform"] = uniprot_accession
							
							if self.options['annotate_isoforms']:
								uniprot_accession_primary_isoform = uniprot_accession.split("-")[0]

								if uniprot_accession_primary_isoform not in isoforms_data:
									response_isoforms = queryRunner.queryRunner("uniprot","parse_isoforms",{"accession":uniprot_accession_primary_isoform}).run()
								
									if 'data' in response_isoforms:	
										if 'isoforms' in response_isoforms['data']:	
											isoforms_data[uniprot_accession_primary_isoform] = response_isoforms['data']['isoforms']
							
							for transcript_protein_exon_iter in range(0,len(transcript_protein_exons)):
								try:
									transcript_protein_exon = copy.deepcopy(transcript_protein_exons[transcript_protein_exon_iter])
									transcript_protein_exon['sequence'] = transcript_protein['seq'][transcript_protein_exon['start']-1:transcript_protein_exon['end']]
								
									if "proteins" not in transcript_ensembl_exon_data[transcript_protein_exon['id']]:
										transcript_ensembl_exon_data[transcript_protein_exon['id']]["proteins"] = {}
								
									transcript_ensembl_exon_data[transcript_protein_exon['id']]["proteins"][uniprot_accession] = transcript_protein_exon['rank']
									transcript_ensembl_exon_data[transcript_protein_exon['id']]["sequence"] = transcript_protein_exon['sequence']
									transcript_ensembl_exon_data[transcript_protein_exon['id']]["in_transcript"] = True
									transcript_ensembl_exon_data[transcript_protein_exon['id']]["full_protein_coding_exon"] = True
									transcript_ensembl_exon_data[transcript_protein_exon['id']]["protein_coding_start"] = 0
									transcript_ensembl_exon_data[transcript_protein_exon['id']]["protein_coding_end"] = len(transcript_ensembl_exon_data[transcript_protein_exon['id']]['dna_sequence'])
									transcript_ensembl_exon_data[transcript_protein_exon['id']]["exon_length"] = len(transcript_ensembl_exon_data[transcript_protein_exon['id']]['dna_sequence'])

									if transcript_ensembl_exon_data[transcript_protein_exon['id']]['ensembl_phase'] == 0:
										transcript_ensembl_exon_data[transcript_protein_exon['id']]["translated_sequence"]  = utilities_codons.translate_oligo(transcript_ensembl_exon_data[transcript_protein_exon['id']]['dna_sequence'][transcript_ensembl_exon_data[transcript_protein_exon['id']]['ensembl_phase']:])
										transcript_ensembl_exon_data[transcript_protein_exon['id']]["cds_sequence"] = transcript_ensembl_exon_data[transcript_protein_exon['id']]['dna_sequence']
									elif transcript_ensembl_exon_data[transcript_protein_exon['id']]['ensembl_phase'] > 0:
										transcript_ensembl_exon_data[transcript_protein_exon['id']]["translated_sequence"]  = '-' + utilities_codons.translate_oligo(transcript_ensembl_exon_data[transcript_protein_exon['id']]['dna_sequence'][3-transcript_ensembl_exon_data[transcript_protein_exon['id']]['ensembl_phase']:])
										transcript_ensembl_exon_data[transcript_protein_exon['id']]["cds_sequence"] = transcript_ensembl_exon_data[transcript_protein_exon['id']]['dna_sequence']
									elif transcript_ensembl_exon_data[transcript_protein_exon['id']]['ensembl_phase'] < 0:
										transcript_ensembl_exon_data[transcript_protein_exon['id']]["full_protein_coding_exon"] = False

										for i in range(0,3):
											translated_sequence = utilities_codons.translate_oligo(transcript_ensembl_exon_data[transcript_protein_exon['id']]['dna_sequence'][i:])
											
											if len(transcript_ensembl_exon_data[transcript_protein_exon['id']]["sequence"]) > 1:
												split_translated_sequence = translated_sequence.split(transcript_ensembl_exon_data[transcript_protein_exon['id']]["sequence"][:-1])
											else:
												split_translated_sequence = []

											if len(split_translated_sequence) > 1 or translated_sequence.count(transcript_ensembl_exon_data[transcript_protein_exon['id']]["sequence"][:-1]) == 1:
												transcript_ensembl_exon_data[transcript_protein_exon['id']]["protein_coding_start"] = len(split_translated_sequence[0])*3 + i
												transcript_ensembl_exon_data[transcript_protein_exon['id']]["cds_sequence"] = transcript_ensembl_exon_data[transcript_protein_exon['id']]['dna_sequence'][transcript_ensembl_exon_data[transcript_protein_exon['id']]["protein_coding_start"]:]
												translated_sequence = utilities_codons.translate_oligo(transcript_ensembl_exon_data[transcript_protein_exon['id']]["cds_sequence"])
												transcript_ensembl_exon_data[transcript_protein_exon['id']]["translated_sequence"] = translated_sequence
											
											if  translated_sequence[-len(transcript_ensembl_exon_data[transcript_protein_exon['id']]["sequence"]):] == transcript_ensembl_exon_data[transcript_protein_exon['id']]["sequence"]:
											#	print(i,translated_sequence)
												transcript_ensembl_exon_data[transcript_protein_exon['id']]["protein_coding_start"] = (len(transcript_ensembl_exon_data[transcript_protein_exon['id']]['dna_sequence']) - len(transcript_ensembl_exon_data[transcript_protein_exon['id']]["sequence"])*3) + transcript_ensembl_exon_data[transcript_protein_exon['id']]['ensembl_end_phase']
												transcript_ensembl_exon_data[transcript_protein_exon['id']]["cds_sequence"] = transcript_ensembl_exon_data[transcript_protein_exon['id']]['dna_sequence'][transcript_ensembl_exon_data[transcript_protein_exon['id']]["protein_coding_start"]:]
												translated_sequence = utilities_codons.translate_oligo(transcript_ensembl_exon_data[transcript_protein_exon['id']]["cds_sequence"])
												transcript_ensembl_exon_data[transcript_protein_exon['id']]["translated_sequence"] = translated_sequence
									
									if 'cds_sequence' not in transcript_ensembl_exon_data[transcript_protein_exon['id']]:
										transcript_ensembl_exon_data[transcript_protein_exon['id']]["cds_sequence"] = transcript_ensembl_exon_data[transcript_protein_exon['id']]['dna_sequence']
										logging.error("cds_sequence not found making assumption of dna_sequence")
								except:
									logger.error(str(transcript_protein_exon_iter) + " " + str(transcript_protein_exons[transcript_protein_exon_iter]))
									utilities_error.printError()

							if transcript_protein_exon['id'] in transcript_ensembl_exon_data:
								if "translated_sequence" in transcript_ensembl_exon_data[transcript_protein_exon['id']]:
									if transcript_ensembl_exon_data[transcript_protein_exon['id']]["translated_sequence"].count("*") > 0:
										transcript_ensembl_exon_data[transcript_protein_exon['id']]["protein_coding_end"] = (len(transcript_protein_exon['sequence'])*3 - transcript_ensembl_exon_data[transcript_protein_exon['id']]["ensembl_phase"]) + 3

						#	print(transcript_protein_exon['id'])
						#	split_translated_sequence = transcript_ensembl_exon_data[transcript_protein_exon['id']]["translated_sequence"].split(transcript_ensembl_exon_data[transcript_protein_exon['id']]["sequence"][1:-1])

							logger.debug("\t".join([
								uniprot_accession_primary_isoform,
								self.options['accession'][0],
								self.options['ensembl_gene_id'],
								self.options['ensembl_transcript_id'],
								self.options['ensembl_protein_id'],
								str(len(transcript_protein_exons)),
								str(len(transcript_protein['seq'])),
								uniprot_accession
							]))
							
							if uniprot_accession != "":
								if uniprot_accession_primary_isoform in isoforms_data:
									if uniprot_accession in isoforms_data[uniprot_accession_primary_isoform]['isoform_details']: 
										transcript_ensembl_protein_id_info["details"] = isoforms_data[uniprot_accession_primary_isoform]['isoform_details'][uniprot_accession]
							
							if self.options['transcript_centric']:
								transcript_ensembl_protein_data[self.options['ensembl_transcript_id']] = transcript_ensembl_protein_id_info
							else:
								if uniprot_accession.split("-")[-1] == "1":
									transcript_ensembl_protein_data[uniprot_accession.split("-")[0]] = transcript_ensembl_protein_id_info
								else:
									transcript_ensembl_protein_data[uniprot_accession] = transcript_ensembl_protein_id_info
			else:
				logger.debug("#Missing Ensembl gene ids for " + self.options['accession'][0] + "")

		except:
			logger.error("Error with " + self.options['accession'][0])
			utilities_error.printError()

		protein_list = list(set(protein_list))
		
		logging.debug("UniProt isoform exon match found: " + str(uniprot_isoform_exon_match_found))

		for protein in transcript_ensembl_protein_data:
			protein_exon_info = {}
			for exon in transcript_ensembl_protein_data[protein]['transcript_protein_exons']:	
				if exon['id'] not in transcript_ensembl_exon_data: continue
				
				transcript_ensembl_exon_data[exon['id']].update(exon)

				for tag in ['feature_type']:
					if tag in transcript_ensembl_exon_data[exon['id']]:
						del transcript_ensembl_exon_data[exon['id']][tag]

				if 'constitutive' not in transcript_ensembl_exon_data[exon['id']]:
					if 'proteins' in transcript_ensembl_exon_data[exon['id']]:
						transcript_ensembl_exon_data[exon['id']]['constitutive'] = len(protein_list) == len(transcript_ensembl_exon_data[exon['id']]['proteins'])
					else:
						transcript_ensembl_exon_data[exon['id']]['constitutive'] = False
				
				for protein_accession in protein_list:
					if 'proteins' in transcript_ensembl_exon_data[exon['id']]:
						if protein_accession not in transcript_ensembl_exon_data[exon['id']]['proteins']:
							transcript_ensembl_exon_data[exon['id']]['proteins'][protein_accession] = False
					else:
						transcript_ensembl_exon_data[exon['id']]['proteins'] = {}
						transcript_ensembl_exon_data[exon['id']]['proteins'][protein_accession] = False

				protein_exon_info[exon['rank']] = transcript_ensembl_exon_data[exon['id']]

			transcript_ensembl_protein_data[protein]['transcript_protein_exons'] = protein_exon_info 
	
		if self.options['get_exon_data']:
			return transcript_ensembl_exon_data
		elif self.options['get_protein_data']:
			return transcript_ensembl_protein_data

	def get_exon_protein_data(self):
		self.options['get_protein_data'] = True
		return self.process_isoform_exon_data()

	def get_exon_data(self):
		self.options['get_exon_data'] = True
		return self.process_isoform_exon_data()

	def grab_human_isoform_exon_data(self):
		taxon_id = "9606"
		response_taxon_accessions = queryRunner.queryRunner("uniprot","parse_uniprot_accession_taxa",{"taxon_id":taxon_id,"structure_required":False,"reviewed":True}).run()
		counter = 0
		for taxon_accesion in response_taxon_accessions['data']:
			counter += 1
			logger.debug("#"*20 + str(counter) + " " + taxon_accesion)
			try:
				self.options['accession'] = [taxon_accesion]
				self.grab_isoform_exon_data()
			except:
				raise
	
	def grab_ensembl_gene_id(self):
		response_db_xref = queryRunner.queryRunner("uniprot","parse_db_xref",{"accession":self.options['accession']}).run()
		#https://rest.ensembl.org/sequence/id/ENST00000676264?content-type=application/json;type=protein
		#https://rest.ensembl.org/overlap/translation/ENSP00000288602?feature=translation_exon;content-type=application/json
		#pprint.pprint(response_db_xref['data']['Ensembl'])

	######### END EXON FUNCTIONS ###############
	######### END EXON FUNCTIONS ###############
	######### END EXON FUNCTIONS ###############
	######### END EXON FUNCTIONS ###############

	def grab_dna_regions(self):
		logger.debug("Grabbing: " + self.options['dna_region_coordinates'])

		if self.options['species'] == "homo_sapiens":
			url = "https://rest.ensembl.org/sequence/region/human"
			out_path = os.path.join(self.options['ensembl_path'], 'dna_region_coordinates',self.options['dna_region_coordinates'].replace(":","_") + ".dna_region.json")
		
		else:
			url = "https://rest.ensembl.org/sequence/region/" + self.options['species'] 
			out_path = os.path.join(self.options['ensembl_path'], 'dna_region_coordinates',self.options['dna_region_coordinates'].replace(":","_") + "." + self.options['species'] + ".dna_region.json")
		

		headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
		data='{ "regions" : ["' + self.options['dna_region_coordinates']  + '"]}'
		
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="POST",data=data,headers=headers,remake_age=self.options['remake_age'])

		if os.path.exists(out_path):
			try:
				with open(out_path) as outfile:
					json_data = json.load(outfile)
				
				logger.debug("Reading" + out_path)
				return json_data
			except:
				print("Error reading: " + out_path)
				raise
		else:
			logger.error(out_path)
			logger.error(status)
			return {}

	def grab_transcript_protein(self):
		logger.debug("Grabbing: " + self.options['ensembl_transcript_id'])
		url = "https://rest.ensembl.org/sequence/id/"  + self.options['ensembl_transcript_id'] + "?content-type=application/json;type=protein"
		out_path = os.path.join(self.options['ensembl_path'],'ensembl_protein',self.options['ensembl_transcript_id'] + ".ensembl_protein.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])
	
		if os.path.exists(out_path):
			try:
				with open(out_path) as outfile:
					json_data = json.load(outfile)
				
				logger.debug("Reading" + out_path)
				return json_data
			except:
				print("Error reading: " + out_path)
				raise
		else:
			logger.error(out_path)
			logger.error(status)
			return {}
		

	def grab_transcript_protein_exons(self):
		logger.debug("Grabbing: " + self.options['ensembl_protein_id'])
		url = "https://rest.ensembl.org/overlap/translation/"  + self.options['ensembl_protein_id'] + "?feature=translation_exon;content-type=application/json"  
		out_path = os.path.join(self.options['ensembl_path'],'ensembl_protein_exons',self.options['ensembl_protein_id'] + ".ensembl_protein_exons.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])
		
		if os.path.exists(out_path):
			try:
				with open(out_path) as outfile:
					json_data = json.load(outfile)
				
				logger.debug("Reading" + out_path)
				return json_data
			except:
				print("Error reading: " + out_path)
				raise
		else:
			logger.error(out_path)
			logger.error(status)
			return {}

	
	def grab_transcript_gene_exons(self):
		logger.debug("Grabbing: " + self.options['ensembl_gene_id'])
		url = "https://rest.ensembl.org/overlap/id/"  + self.options['ensembl_gene_id'].split(".")[0]  + "?content-type=application/json;feature=exon"  
		out_path = os.path.join(self.options['ensembl_path'],'ensembl_gene_exons',self.options['ensembl_gene_id'] + ".ensembl_gene_exons.json")
		
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])
		
		if os.path.exists(out_path):
			try:
				with open(out_path) as outfile:
					json_data = json.load(outfile)
				
				logger.debug("Reading" + out_path)
				return json_data
			except:
				print("Error reading: " + out_path)
				raise
		else:
			logger.error(out_path)
			logger.error(status)
			return {}
		
	def grab_gene_transcripts(self):
		logger.debug("Grabbing: " + self.options['ensembl_gene_id'])
		url = "https://rest.ensembl.org/overlap/id/"  + self.options['ensembl_gene_id'].split(".")[0]  + "?content-type=application/json;feature=transcript"  
		out_path = os.path.join(self.options['ensembl_path'],'ensembl_gene_exons',self.options['ensembl_gene_id'] + ".ensembl_gene_transcript.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])
		
		if os.path.exists(out_path):
			try:
				with open(out_path) as outfile:
					json_data = json.load(outfile)

				return json_data
			except:
				print("Error reading: " + out_path)
				raise
		else:
			logger.error(out_path)
			logger.error(status)
			return {}
		

	def grab_gene_genename(self):
		logger.debug("Grabbing: " + self.options['gene_name'])
		
		url = "https://rest.ensembl.org/xrefs/symbol/" + self.options['species']  + "/" + self.options['gene_name'] + "?db_type=core&content-type=application/json"

		if self.options['species'] == "homo_sapiens":
			out_path = os.path.join(self.options['ensembl_path'],'gene_name',self.options['gene_name'].replace("/","_") + ".gene_genename.json")
		else:
			out_path = os.path.join(self.options['ensembl_path'],'gene_name',self.options['gene_name'].replace("/","_") + "." + self.options['species'] + ".gene_genename.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])
		
		if os.path.exists(out_path):
			try:
				with open(out_path) as outfile:
					json_data = json.load(outfile)
				
				logger.debug("Reading" + out_path)
				
			except:
				print("Error reading: " + out_path)
				raise
		else:
			logger.error(out_path)
			logger.error(status)
			return []

		ensembl_gene_ids = []
		for value in json_data:
			if self.options['species'] == "homo_sapiens":
				if value['id'][0:4] == 'ENSG':
					ensembl_gene_ids.append(value['id'])
			else:
				if value['id'][0:3] == 'ENS':
					ensembl_gene_ids.append(value['id'])

		skip_ensembl_gene_ids = []
		for ensembl_gene_id in ensembl_gene_ids:
			self.options['ensembl_gene_id'] = ensembl_gene_id
			gene_information = self.grab_gene_transcript_information()
			
			if self.options['remove_scaffolds']:
				if 'seq_region_name' in gene_information:
					if len(gene_information['seq_region_name']) > 2:
						logger.error(">>>>>>>>>>>>> Removing " + ensembl_gene_id + " from gene id list - non-canonical chromosome " + gene_information['seq_region_name'])
						skip_ensembl_gene_ids.append(ensembl_gene_id)
					else:
						logger.debug(">>>>>>>>>>>>> Keeping " + ensembl_gene_id + " from gene id list- canonical chromosome " + gene_information['seq_region_name'])
				else:
					logger.error(">>>>>>>>>>>>> Removing " + ensembl_gene_id + " from gene id list -  no chromosome information")
					
			if 'display_name' in gene_information:
				if gene_information['display_name'] != self.options['gene_name']:
					logger.error(">>>>>>>>>>>>> Removing " + ensembl_gene_id + " from gene id list - display_name " + gene_information['display_name'] + " != " + self.options['gene_name'])
					skip_ensembl_gene_ids.append(ensembl_gene_id)
			else:
				logger.error(">>>>>>>>>>>>> Removing " + ensembl_gene_id + " from gene id list - no display_name")
				skip_ensembl_gene_ids.append(ensembl_gene_id)
		
		ensembl_gene_ids = list(set(ensembl_gene_ids).difference(skip_ensembl_gene_ids))
		return ensembl_gene_ids
		
	def grab_gene_uniprot_xref(self):
		logger.debug("Grabbing: " + self.options['ensembl_gene_id'])
		url = "https://rest.ensembl.org/xrefs/id/" + self.options['ensembl_gene_id'] + "?external_db=Uniprot_gn&content-type=application/json"
		out_path = os.path.join(self.options['ensembl_path'],'uniprot_xref',self.options['ensembl_gene_id'] + ".ensembl_gene.uniprot_xref.json")
		
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])
		
		if os.path.exists(out_path):
			try:
				with open(out_path) as outfile:
					json_data = json.load(outfile)
				
				logger.debug("Reading" + out_path)
				return json_data
			except:
				print("Error reading: " + out_path)
				raise
		else:
			logger.error(out_path)
			logger.error(status)
			return {}


	def grab_transcript_uniprot_xref(self):
		logger.debug("Grabbing: " + self.options['ensembl_transcript_id'])
		url = "https://rest.ensembl.org/xrefs/id/" + self.options['ensembl_transcript_id'] + "?all_levels=1;external_db=Uniprot/SPTREMBL;content-type=application/json"
		out_path = os.path.join(self.options['ensembl_path'],'uniprot_xref',self.options['ensembl_transcript_id'] + ".ensembl_transcript.uniprot_xref.json")
		
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])
		
		if os.path.exists(out_path):
			try:
				with open(out_path) as outfile:
					json_data = json.load(outfile)
				
				logger.debug("Reading" + out_path)
				return json_data
			except:
				print("Error reading: " + out_path)
				raise
		else:
			logger.error(out_path)
			logger.error(status)
			return {}

	def grab_transcript_swissprot_xref(self):
		logger.debug("Grabbing: " + self.options['ensembl_transcript_id'])
		url = "https://rest.ensembl.org/xrefs/id/" + self.options['ensembl_transcript_id'] + "?all_levels=1;external_db=Uniprot/SWISSPROT;content-type=application/json"
		out_path = os.path.join(self.options['ensembl_path'],'uniprot_xref',self.options['ensembl_transcript_id'] + ".ensembl_transcript.swissprotxref.json")
		
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])
		
		if os.path.exists(out_path):
			try:
				with open(out_path) as outfile:
					json_data = json.load(outfile)
				
				logger.debug("Reading" + out_path)
				return json_data
			except:
				print("Error reading: " + out_path)
				raise
		else:
			logger.error(out_path)
			logger.error(status)
			return {}
	
	def grab_gene_transcript_information(self):
		logger.debug("Grabbing: " + self.options['ensembl_gene_id'])
		url = "http://rest.ensembl.org/lookup/id/" + self.options['ensembl_gene_id'].split(".")[0] + "?content-type=application/json;expand=1"
		
		out_path = os.path.join(self.options['ensembl_path'],'ensembl_gene',self.options['ensembl_gene_id'] + ".ensembl_gene.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])
		
		if os.path.exists(out_path):
			try:
				with open(out_path) as outfile:
					json_data = json.load(outfile)
				
				logger.debug("Reading" + out_path)
				return json_data
			except:
				print("Error reading: " + out_path)
				utilities_error.printError()
				raise
		else:
			logger.error(out_path)
			logger.error(status)
			return {}

	def grab_gene_transcript_ids(self):
		logger.debug("Grabbing: " + self.options['ensembl_gene_id'])
		url = "http://rest.ensembl.org/lookup/id/" + self.options['ensembl_gene_id'].split(".")[0] + "?content-type=application/json;expand=1"
		
		out_path = os.path.join(self.options['ensembl_path'],'ensembl_gene',self.options['ensembl_gene_id'] + ".ensembl_gene.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])
		
		if "error" in status:
			if status["error"] == "error":
				logging.error(status)
				return []

		transcript_ensembl_protein_ids = []
		try:
			with open(out_path) as outfile:
				json_data = json.load(outfile)
			
			for transcript in json_data['Transcript']:
				if transcript['biotype'] == "protein_coding":
					transcript_ensembl_protein_ids.append(transcript['id'])
		except:
			logger.error("Error reading: " + out_path)
			#sys.exit()

		return transcript_ensembl_protein_ids

	def grab_genomic_region_genes(self):
		logger.debug("Grabbing: " + str(self.options['chromosome']) + ":" + str(min(self.options['region_start'],self.options['region_end'])) + "-" + str(max(self.options['region_start'],self.options['region_end'])))

		region_id = "_".join(
			[
				str(self.options['chromosome']),
				str(min(self.options['region_start'],self.options['region_end'])),
				str(max(self.options['region_start'],self.options['region_end']))
			])

		logger.debug("Grabbing " + region_id)
		
		url = "https://rest.ensembl.org/overlap/region/human/" + str(self.options['chromosome']) + ":" + str(min(self.options['region_start'],self.options['region_end'])) + "-" + str(max(self.options['region_start'],self.options['region_end'])) + "?feature=gene;content-type=application/json"
		out_path = os.path.join(self.options['genomes_region_path'],region_id + ".region_features.genes.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])
	
		with open(out_path) as outfile:
			json_data = json.load(outfile)

		return json_data

	def grab_genomic_region_features(self):
		region_id = "_".join(
			[
				str(self.options['chromosome']),
				str(min(self.options['region_start'],self.options['region_end'])),
				str(max(self.options['region_start'],self.options['region_end']))
			])

		logger.debug("Grabbing " + region_id)

		feature_options = []
		for feature_type in self.options['feature_type']:
			feature_options.append("feature=" + feature_type)
		feature_options_str = ";".join(feature_options)

		search_region_id = str(self.options['chromosome']) + ":" + str(min(self.options['region_start'],self.options['region_end'])) + "-" + str(max(self.options['region_start'],self.options['region_end']))
		url = "https://rest.ensembl.org/overlap/region/human/" + search_region_id + "?" + feature_options_str + ";content-type=application/json"
		out_path = os.path.join(self.options['genomes_region_path'],region_id + ".region_features." + "_".join(self.options['feature_type']) + ".json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])

		logger.debug("Reading " + out_path)

		if not os.path.exists(out_path):
			logger.error("File does not exist " + out_path)
		else:
			with open(out_path) as outfile:
				json_data = json.load(outfile)

			logger.debug("Size " + str(len(json_data)))

		return json_data
		
	def grab_transcript_sequence_information(self,sequence_information_type="protein"):
		logger.debug("Grabbing: " + self.options['ensembl_transcript_id'])
		
		if self.options['sequence_information_type'] != "protein":
			sequence_information_type = self.options['sequence_information_type']
		
		url = "https://rest.ensembl.org/sequence/id/"  + self.options['ensembl_transcript_id'] + "?content-type=application/json;type=" + sequence_information_type
		out_path = os.path.join(self.options['ensembl_path'],"sequence_information",self.options['ensembl_transcript_id'] + ".ensembl_" + sequence_information_type + ".json")


		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,JSON=True,method="GET",remake_age=self.options['remake_age'])

		if not os.path.exists(out_path):
			logger.error("File does not exist " + out_path)
			json_data = {}
		else:
			with open(out_path) as outfile:
				json_data = json.load(outfile)

		return json_data
	
	"""
	def grab_ensembl_gene_id_by_uniprot_xref(self):
		logger.debug("Grabbing: " + self.options['accession'][0])
		response_db_xref = queryRunner.queryRunner("uniprot","parse_db_xref",{"accession":self.options['accession'][0].split("-")[0]}).run()
		
		gene_ids = []
		if 'Ensembl' in response_db_xref['data']:
			for xref in response_db_xref['data']['Ensembl']:
				gene_ids.append(response_db_xref['data']['Ensembl'][xref]['gene ID'])
		return list(set(gene_ids))
	"""
	
	def grab_ensembl_transcript_id_by_uniprot_xref(self):
		logger.debug("Grabbing: " + self.options['accession'][0])
		response_db_xref = queryRunner.queryRunner("uniprot","parse_db_xref",{"accession":self.options['accession'][0].split("-")[0]}).run()
		
		
		if 'Ensembl' in response_db_xref['data']:
			for transcript_id in response_db_xref['data']['Ensembl']:
				if 'molecule' in (response_db_xref['data']['Ensembl'][transcript_id]):
					if self.options['accession'][0] == response_db_xref['data']['Ensembl'][transcript_id]['molecule']:
						return transcript_id

					if self.options['accession'][0].count("-") == 0:
						if self.options['accession'][0] + '-1' == response_db_xref['data']['Ensembl'][transcript_id]['molecule']:
							return transcript_id
		
		return ""
	
	def grab_ensembl_gene_id_by_uniprot_xref(self,as_list=True):
		logger.debug("Grabbing: " + self.options['accession'][0])
		response_db_xref = queryRunner.queryRunner("uniprot","parse_db_xref",{"accession":self.options['accession'][0].split("-")[0]}).run()
		
		gene_ids = []
		if 'Ensembl' in response_db_xref['data']:
			for transcript_id in response_db_xref['data']['Ensembl']:
				if 'molecule' in (response_db_xref['data']['Ensembl'][transcript_id]):
					if self.options['accession'][0] == response_db_xref['data']['Ensembl'][transcript_id]['molecule']:
						if as_list:
							return [response_db_xref['data']['Ensembl'][transcript_id]['gene ID']]
						else:
							return response_db_xref['data']['Ensembl'][transcript_id]['gene ID']

					if self.options['accession'][0].count("-") == 0:
						if self.options['accession'][0] + '-1' == response_db_xref['data']['Ensembl'][transcript_id]['molecule']:
							if as_list:
								return [response_db_xref['data']['Ensembl'][transcript_id]['gene ID']]
							else:
								return response_db_xref['data']['Ensembl'][transcript_id]['gene ID']

				else:
					if 'Ensembl' in response_db_xref['data']:
						for xref in response_db_xref['data']['Ensembl']:
							gene_ids.append(response_db_xref['data']['Ensembl'][xref]['gene ID'])

		return gene_ids

	def make_ensembl_transcript_uniprot_mapping(self,transcript_protein_alignment,uniprot_protein_alignment):
		
		if len(transcript_protein_alignment) == 0: return {}
		if len(uniprot_protein_alignment) == 0: return {}
		
		mapping = {
			"n_terminal_extension":False,
			"c_terminal_extension":False,
			"sequence_mismatches":[],
			"sequence_identical":transcript_protein_alignment==uniprot_protein_alignment,
			"ensembl_uniprot":{},
			"uniprot_ensembl":{},
		}

		transcript_protein_alignment_offet = 0
		uniprot_protein_alignment_offet = 0
		for i in range(0,len(transcript_protein_alignment)):
			if transcript_protein_alignment[i] != "-":
				transcript_protein_alignment_offet += 1
			if uniprot_protein_alignment[i] != "-":
				uniprot_protein_alignment_offet += 1
			
			if uniprot_protein_alignment[i] != transcript_protein_alignment[i] and transcript_protein_alignment[i] != "-" and uniprot_protein_alignment[i] != "-":
				mapping["sequence_mismatches"].append({
					"uniprot_offset":uniprot_protein_alignment_offet,
					"ensembl_offset":transcript_protein_alignment_offet,
					"ensembl_residue":transcript_protein_alignment[i],
					"uniprot_residue":uniprot_protein_alignment[i] 
					}
				)

			if uniprot_protein_alignment[i] != transcript_protein_alignment[i] and uniprot_protein_alignment[i] == "-":
				mapping["deletion"] = True
			
			if uniprot_protein_alignment[i] != transcript_protein_alignment[i] and uniprot_protein_alignment[i] == "-" and uniprot_protein_alignment[:i].count("-") == len(uniprot_protein_alignment[:i]):
				mapping["n_terminal_extension"] = True

			elif uniprot_protein_alignment[i] != transcript_protein_alignment[i] and uniprot_protein_alignment[i] == "-" and uniprot_protein_alignment[i:].count("-") == len(uniprot_protein_alignment[i:]):
				mapping["c_terminal_extension"] = True
			
			elif uniprot_protein_alignment[i] != transcript_protein_alignment[i] and transcript_protein_alignment[i] == "-":
				mapping["insertion"] = True
			
			
			mapping["ensembl_uniprot"][transcript_protein_alignment_offet] = {
				"uniprot_offset":uniprot_protein_alignment_offet,
				"aligned":transcript_protein_alignment[i]==uniprot_protein_alignment[i],
				"uniprot_gap":uniprot_protein_alignment[i] == "-",
				"ensembl_residue":transcript_protein_alignment[i],
				"uniprot_residue":uniprot_protein_alignment[i] 
			}
			mapping["uniprot_ensembl"][uniprot_protein_alignment_offet] = {
				"ensembl_offset":transcript_protein_alignment_offet,
				"aligned":transcript_protein_alignment[i]==uniprot_protein_alignment[i],
				"ensembl_gap":transcript_protein_alignment[i] == "-",
				"ensembl_residue":transcript_protein_alignment[i],
				"uniprot_residue":uniprot_protein_alignment[i] 
			}
		
		#mapping["uniprot_start"]
		#mapping["uniprot_end"]
		#mapping["ensembl_start"]
		#mapping["ensembl_end"]
		#print(min(mapping["ensembl_uniprot"]),mapping["ensembl_uniprot"][min(mapping["ensembl_uniprot"])]['uniprot_offset'])
		#print(max(mapping["ensembl_uniprot"]),mapping["ensembl_uniprot"][max(mapping["ensembl_uniprot"])]['uniprot_offset'])
		#print(min(mapping["uniprot_ensembl"]),mapping["uniprot_ensembl"][min(mapping["uniprot_ensembl"])]['ensembl_offset'])
		#print(max(mapping["uniprot_ensembl"]),mapping["uniprot_ensembl"][max(mapping["uniprot_ensembl"])]['ensembl_offset'])
	#	pprint.pprint(mapping)

	def grab_ensembl_transcript_id_by_genename(self):
		gene_ids = self.grab_ensembl_gene_id_by_genename()
		
		logger.debug(gene_ids)
		
		transcript_id_length_match = {
			"transcript_id":"",
			"gene_id":"",
			"gene_ids":[],
			"uniprot_mapped_proportion":0,
			"identity_count":0,
			"transcript_protein_unmapped_residues_count":0,
			"uniprot_protein_unmapped_residues_count":0,
			"transcript_protein_length":0,
			"uniprot_protein_length":0,
			"transcript_protein":False,
			"transcript_protein_alignment":"",
			"uniprot_protein_alignment":"",
			"identity_alignment":"",
			"match_type":"unknown"
		}
		transcript_details = {}
		gene_ids.sort()
		
		##
		##   Can match the wrong transcript when the protein is encoded by more than one gene / transcript
		##

		self.options['gene_ids'] = gene_ids
		for gene_id in gene_ids:
			self.options['ensembl_gene_id'] = gene_id
			transcript_ids = self.grab_gene_transcript_ids()
			
			if len(transcript_ids) == 0:
				logger.error("No transcript_ids for " + gene_id + " " + ",".join(gene_ids))
				continue

			self.options['sequence'] = self.grab_uniprot_sequence()
			transcript_ids.sort()

			for transcript_id in transcript_ids:
				try:
					self.options['ensembl_transcript_id'] = transcript_id
					transcript_protein = self.grab_transcript_sequence_information("protein")
					transcript_genomic = self.grab_transcript_sequence_information("genomic")

					match_oligos_check = None
					if 'match_oligos' in self.options:
						match_oligos_check = True
						if len(self.options['match_oligos']) == 0:
							match_oligos_check = False
						else:
							for match_oligo in self.options['match_oligos']:
								match_oligo_complementry = utilities_codons.get_complementry_oligo(match_oligo)[::-1]
								if transcript_genomic['seq'].count(match_oligo) == 0 and transcript_genomic['seq'].count(match_oligo_complementry) == 0:
									match_oligos_check = False

					if transcript_protein["seq"] == self.options['sequence'] and match_oligos_check in [True,None]:
						if "detailed" in self.options:
							if self.options["detailed"]:
								return {
									"transcript_id":transcript_id,
									"gene_id":gene_id,
									"gene_ids":gene_ids,
									"uniprot_mapped_proportion":1,
									"identity_count":len(transcript_protein["seq"]),
									"transcript_protein_unmapped_residues_count":0,
									"uniprot_protein_unmapped_residues_count":0,
									"transcript_protein_length":len(transcript_protein['seq']),
									"uniprot_protein_length":len(self.options['sequence']),
									"transcript_protein_alignment":transcript_protein['seq'],
									"uniprot_protein_alignment":self.options['sequence'],
									"identity_alignment":'+'*len(self.options['sequence']),
									"transcript_protein":transcript_protein['seq'],
									"match_type":"exact",
									"match_oligos_check":match_oligos_check
								}
						
						return transcript_id
					
					alignment = utilities_alignment.alignPeptides(transcript_protein['seq'],self.options['sequence'])
					
					if isinstance(alignment[0][0],list):
						transcript_cds_translated_align = "".join(alignment[0][0])
						uniprot_sequence_align = "".join(alignment[0][1])
					else:
						transcript_cds_translated_align = alignment[0][0]
						uniprot_sequence_align = alignment[0][1]
						
					identity = utilities_alignment.sequenceIdentity(transcript_cds_translated_align,uniprot_sequence_align)
					uniprot_mapped_proportion = identity.count("+")/len(self.options['sequence'])
					
					logger.debug(transcript_id + " length:" + str(transcript_id_length_match['uniprot_mapped_proportion']) + " identity:" + "%1.2f"%uniprot_mapped_proportion)
					
					transcript_details[transcript_id] = {
						"uniprot_mapped_proportion":uniprot_mapped_proportion,
						"identity_count":identity.count("+"),
						"transcript_protein_unmapped_residues_count":uniprot_sequence_align.count("-"),
						"uniprot_protein_unmapped_residues_count":transcript_cds_translated_align.count("-"),
					}

					match_details = {
						"transcript_id":transcript_id,
						"gene_id":gene_id,
						"gene_ids":gene_ids,
						"uniprot_mapped_proportion":uniprot_mapped_proportion,
						"identity_count":identity.count("+"),
						"transcript_protein_unmapped_residues_count":uniprot_sequence_align.count("-"),
						"uniprot_protein_unmapped_residues_count":transcript_cds_translated_align.count("-"),
						"transcript_protein_length":len(transcript_protein['seq']),
						"uniprot_protein_length":len(self.options['sequence']),
						"transcript_protein":transcript_protein['seq'],
						"transcript_protein_alignment":transcript_cds_translated_align,
						"uniprot_protein_alignment":uniprot_sequence_align,
						"identity_alignment":identity,
						"match_oligos_check":match_oligos_check
					}


					if match_details['transcript_protein_unmapped_residues_count'] == 0:
						match_details["match_type"] = "alignment_complete_ensembl"
					elif match_details['uniprot_protein_unmapped_residues_count'] == 0:
						match_details["match_type"] = "alignment_complete_uniprot"
					else:
						match_details["match_type"] = "alignment_incomplete_uniprot_incomplete_ensembl"
					
					"""
					print("#" + transcript_id)
					print("E","".join(alignment[0][0]))
					print("U","".join(alignment[0][1]))
					print(uniprot_mapped_proportion)
					"""
					if identity.count("+") >= transcript_id_length_match['identity_count']:
						if transcript_protein['seq'].count(self.options['sequence'][1:]) == 1:
							match_details['match_type'] = "subsequence_complete_uniprot_incomplete_ensembl"
							transcript_id_length_match = match_details
						
						if self.options['sequence'].count(transcript_protein['seq'][1:].replace("*","")) == 1:
							match_details['match_type'] = "subsequence_incomplete_uniprot_complete_ensembl"
						
					if uniprot_mapped_proportion == transcript_id_length_match['uniprot_mapped_proportion'] and uniprot_mapped_proportion > 0.5:
						if identity.count("+") >= transcript_id_length_match['identity_count'] and uniprot_sequence_align.count("-") < transcript_id_length_match['transcript_protein_unmapped_residues_count']:
							transcript_id_length_match = match_details
					elif uniprot_mapped_proportion > transcript_id_length_match['uniprot_mapped_proportion'] and uniprot_mapped_proportion > 0.5:
						if identity.count("+") > transcript_id_length_match['identity_count'] :
							transcript_id_length_match = match_details
					else:
						pass

					#print("\t".join([self.options['accession'][0],gene_id,transcript_id,str(len(self.options['sequence'])),"%1.2f"%uniprot_mapped_proportion,str(uniprot_sequence_align.count("-")),str(transcript_cds_translated_align.count("-")),match_details['match_type'],str(transcript_id == transcript_id_length_match['transcript_id']),transcript_id_length_match['transcript_id']]))		
				except:
					logging.debug("Error:" + transcript_id)
			
		if transcript_id_length_match['transcript_protein'] != False:
			logging.debug("\t".join([
				"Updating protein sequence",
				self.options['accession'][0],
				transcript_id_length_match['gene_id'],
				transcript_id_length_match['transcript_id'],
				"%1.2f"%transcript_id_length_match['uniprot_mapped_proportion'],
				str(transcript_id_length_match['uniprot_protein_unmapped_residues_count']),
				str(transcript_id_length_match['transcript_protein_unmapped_residues_count']),
				transcript_id_length_match['match_type']
			]))

			self.options['updated_uniprot_sequence'] = transcript_id_length_match['transcript_protein']
		else:
			logging.error("\t".join([
				"No match",
				self.options['accession'][0],
				transcript_id_length_match['gene_id'],
				transcript_id_length_match['transcript_id'],
				"%1.2f"%transcript_id_length_match['uniprot_mapped_proportion'],
				str(transcript_id_length_match['uniprot_protein_unmapped_residues_count']),
				str(transcript_id_length_match['transcript_protein_unmapped_residues_count']),
				transcript_id_length_match['match_type']
			]))
		
		if "detailed" in self.options:
			if self.options["detailed"]:
				self.make_ensembl_transcript_uniprot_mapping(transcript_id_length_match['transcript_protein_alignment'],transcript_id_length_match['uniprot_protein_alignment'])
				return transcript_id_length_match
			
		return transcript_id_length_match['transcript_id']
	
	def grab_ensembl_transcript_id_by_genename_peptide(self):
		gene_ids = self.grab_ensembl_gene_id_by_genename()

		mappings = {}
		for gene_id in gene_ids:
			self.options['ensembl_gene_id'] = gene_id

			transcript_ids = self.grab_gene_transcript_ids()
			self.options['sequence'] = self.grab_uniprot_sequence()

			for transcript_id in transcript_ids:
				self.options['ensembl_transcript_id'] = transcript_id
				transcript_protein = self.grab_transcript_sequence_information("protein")

				if transcript_protein['seq'].count(self.options['peptide_sequence']) > 0: 	
					if self.options['detailed']:
						if len(transcript_protein['seq'].split(self.options['peptide_sequence'])) == 2 or transcript_protein['seq'] == self.options['peptide_sequence']:
							mappings[transcript_id] = {"transcript_id":transcript_id,"gene_id":gene_id,"peptide_start":len(transcript_protein['seq'].split(self.options['peptide_sequence'])[0])}
					else:
						if len(transcript_protein['seq'].split(self.options['peptide_sequence'])) == 2:
							return {"transcript_id":transcript_id,"gene_id":gene_id,"peptide_start":len(transcript_protein['seq'].split(self.options['peptide_sequence'])[0])}
						else:
							return {"transcript_id":transcript_id,"gene_id":gene_id}

		return mappings
	"""
	def grab_ensembl_gene_id_by_genename(self):
		gene_ids = self.grab_ensembl_gene_id_by_genename()
		gene_id_length_match = "" 
		
		for gene_id in gene_ids:
			self.options['ensembl_gene_id'] = gene_id

			transcript_ids = self.grab_gene_transcript_ids()
			self.options['sequence'] = self.grab_uniprot_sequence()

			for transcript_id in transcript_ids:
				self.options['ensembl_transcript_id'] = transcript_id
				transcript_protein = self.grab_transcript_sequence_information("protein")

				if transcript_protein['seq'] == self.options['sequence']: 	
					#print("transcript_protein match",transcript_id,len(transcript_protein['seq']),len(self.options['sequence']))
					return gene_id
				elif len(transcript_protein['seq']) == len(self.options['sequence']):
					#print("transcript_protein length match",transcript_id,len(transcript_protein['seq']),len(self.options['sequence']))
					gene_id_length_match = gene_id	
				else:
					pass#print("transcript_protein mismatch",transcript_id,len(transcript_protein['seq']),len(self.options['sequence']))
			
		return gene_id_length_match
	"""

	def grab_ensembl_gene_id_by_genename(self):
		logger.debug("Grabbing: " + self.options['accession'][0])
		response_uniprot = queryRunner.queryRunner("uniprot","parse_basic",{"accession":self.options['accession'][0]}).run()
		if 'gene_name' in response_uniprot['data']:
			self.options['gene_name'] = response_uniprot['data']['gene_name']
			if self.options['gene_name'] == "":
				return {}
			else:
				return self.grab_gene_genename()
		else:
			return {}
		
	def grab_uniprot_sequence(self):
		logger.debug("Grabbing: " + self.options['accession'][0])
		response_uniprot = queryRunner.queryRunner("uniprot","parse_basic",{"accession":self.options['accession'][0]}).run()
		self.options['sequence'] = response_uniprot['data']['sequence']
		return self.options['sequence']

	def align_exons_to_genomic_transcript(self,transcript_id,transcript_genomic,transcript_cds):
		cds_genomic_alignment_path = os.path.join(self.options["cds_genomic_alignment_path"],transcript_id + ".exon_mapping.v2.json")
		self.options['remake'] = True

		if not os.path.exists(cds_genomic_alignment_path) or self.options['remake']:
			self.options['transcript_centric'] = True
			exon_data = self.get_exon_protein_data()
			self.options['transcript_centric'] = False
				
			if transcript_id in exon_data:
				transcript_protein_exons = exon_data[transcript_id]
			else:
				logger.error(self.options['accession'][0].split("-")[0] + " not in " + transcript_id + " transcript_protein_exons " + str(exon_data.keys()) + " " + cds_genomic_alignment_path)
				return ["",""]
			
			alignment_holder = ["-"]*len(transcript_genomic['seq'])
			exons_length_summer = 0

			for transcript_protein_exon in transcript_protein_exons['transcript_protein_exons']: 
				try:
					dna_sequence = transcript_protein_exons['transcript_protein_exons'][transcript_protein_exon]['cds_sequence']
					exon_start = transcript_protein_exons['transcript_protein_exons'][transcript_protein_exon]['protein_coding_start']
					exon_end = transcript_protein_exons['transcript_protein_exons'][transcript_protein_exon]['protein_coding_end']

					exons_length_summer += len(dna_sequence)

					logger.debug("Adding exon " + transcript_id + "-" + str(transcript_protein_exon) + " " +  str(exon_start)  + " " +  str(exon_end) + " " + str(len(dna_sequence)) +  " " + str(len(transcript_protein_exons['transcript_protein_exons'])) + " " + str(exons_length_summer))
					
				
					if transcript_cds['seq'].count(dna_sequence) == 0:
						if len(transcript_protein_exons['transcript_protein_exons']) == 1:
							logger.debug("Single exon transcript")
							dna_sequence = transcript_cds['seq']
						else:
							dna_sequence = dna_sequence[exon_start:exon_end]

					try:
						transcript_genomic_exon_start = len(transcript_genomic['seq'].split(dna_sequence)[0])
						transcript_genomic_exon_end = transcript_genomic_exon_start + len(dna_sequence)

						for i in range(0,len(dna_sequence)):
							alignment_holder[transcript_genomic_exon_start + i] = dna_sequence[i]

					except:
						logger.error("Error in exon mapping: " + transcript_id + "-" + str(transcript_protein_exon))
						logger.error(transcript_protein_exons['transcript_protein_exons'][transcript_protein_exon])
						utilities_error.printError()
				except:
					logger.error("Error interating through exons: " + transcript_id + "-" + str(transcript_protein_exon))

			alignment = [transcript_genomic['seq'],"".join(alignment_holder)]
			
			logger.debug("Making long gene protein - " + cds_genomic_alignment_path)
			with open(cds_genomic_alignment_path,'w') as outfile:
				json.dump(alignment, outfile)
		else:
			logger.debug("Reading " + cds_genomic_alignment_path)
			with open(cds_genomic_alignment_path) as outfile:
				alignment = json.load(outfile)
		
		return alignment

	def get_ensembl_data_for_gene(self):		
		transcript_data = {}

		if 'transcript_id' in self.options and 'gene_id' in self.options:
			logger.info("transcript_id/gene_id provided: " + self.options['transcript_id'] + " " + self.options['gene_id'] )
			self.options['ensembl_transcript_id'] = self.options['transcript_id'] 
			self.options['ensembl_gene_id'] = self.options['gene_id'] 
			gene_ids = [self.options['gene_id']]
		else:
			try:
				
				self.options["detailed"] = True
				try:
					if "peptide_sequence" in self.options:
						if len(self.options["peptide_sequence"]) > 0:
							self.options["detailed"] = False
							transcript_data = self.grab_ensembl_transcript_id_by_genename_peptide()
							self.options["detailed"] = True
							self.options['ensembl_transcript_id'] = transcript_data['transcript_id']
							self.options['ensembl_gene_id'] = transcript_data['gene_id']

							transcript_data = self.grab_ensembl_transcript_id_by_genename()
						else:
							transcript_data = self.grab_ensembl_transcript_id_by_genename()
							
					else:
						transcript_data = self.grab_ensembl_transcript_id_by_genename()
				except:
					transcript_data = self.grab_ensembl_transcript_id_by_genename()
					
				self.options['ensembl_transcript_id'] = transcript_data['transcript_id']
				self.options['ensembl_gene_id'] = transcript_data['gene_id']

				gene_ids = transcript_data['gene_ids']
			except:
				utilities_error.printError()
				return
		
		logging.debug("#"*40)
		logging.debug("Finding guides for " + ",".join(self.options['accession']))

		logger.debug("Using " + self.options['ensembl_gene_id'] + " " +  self.options['ensembl_transcript_id'])

		#----#
		
		transcript_protein = self.grab_transcript_sequence_information("protein")
		transcript_cds = self.grab_transcript_sequence_information("cds")
		transcript_genomic = self.grab_transcript_sequence_information("genomic")

		alignment = self.align_exons_to_genomic_transcript(self.options['ensembl_transcript_id'],transcript_genomic,transcript_cds)
		

		if len(transcript_protein) == 0:
			return {
				"ensembl_transcript_id":self.options['ensembl_transcript_id'],
				"ensembl_gene_id":self.options['ensembl_gene_id'],
				"ensembl_gene_ids":gene_ids
			}
		else:
			return {
				"ensembl_transcript_id":self.options['ensembl_transcript_id'],
				"ensembl_gene_id":self.options['ensembl_gene_id'],
				"ensembl_gene_ids":gene_ids,
				"transcript_protein":transcript_protein['seq'],
				"transcript_cds":transcript_cds['seq'],
				"transcript_genomic":transcript_genomic['seq'],
				"transcript_cds_genomic_alignment":alignment[1],
				"transcript_mapping":transcript_data
			}
	
	def get_ensembl_data_for_region(self):

		transcript_id = ""
		check_transcript_start = False
		
		if 'transcript_id' in self.options and 'gene_id' in self.options:
			logger.debug("transcript_id/gene_id provided")
			transcript_id = self.options['transcript_id'] 
			gene_id = self.options['gene_id'] 
		else:
			try:
				transcript_id = self.grab_ensembl_transcript_id_by_genename()
				
				if len(transcript_id) != 0:
					gene_id = self.grab_ensembl_gene_id_by_genename()
				
				logger.debug("transcript_id/gene_id provided by grab_ensembl_transcript_id_by_genename")
			except:
				logger.debug("transcript_id by_genename mapping failed")
			
			try:
				if len(transcript_id) == 0:
					transcript_id = self.grab_ensembl_transcript_id_by_uniprot_xref()
					gene_id = self.grab_ensembl_gene_id_by_uniprot_xref()
					logger.debug("transcript_id/gene_id provided by grab_ensembl_gene_id_by_uniprot_xref")
			except:
				logger.debug("transcript_id by_uniprot_xref mapping failed")
				utilities_error.printError()
			
			try:
				if len(transcript_id) == 0:
					transcript_id_response = self.grab_ensembl_transcript_id_by_genename_peptide()
					
					if len(transcript_id_response) > 0:
						transcript_id = transcript_id_response['transcript_id']
						
						logger.debug("transcript_id/gene_id provided by grab_ensembl_transcript_id_by_genename_peptide")
						if 'peptide_start' in transcript_id_response:
							self.options['peptide_start'] = transcript_id_response['peptide_start']+1 
							self.options['peptide_end'] = transcript_id_response['peptide_start'] + len(self.options['peptide_sequence'])
						
						check_transcript_start = True
						gene_id = self.grab_ensembl_gene_id_by_genename()
			except:
				logger.debug("transcript_id by_genename_peptide mapping failed")
				utilities_error.printError()
				
			if len(transcript_id) == 0:
				logger.error("No transcript matching found")
				return {"status":"error","error_type":"No transcript matching found"}

			transcript_id = transcript_id.split(".")[0]

		try:
			self.options['ensembl_transcript_id'] = transcript_id
			self.options['ensembl_gene_id'] = gene_id
			
			#----#
			
			transcript_protein = self.grab_transcript_sequence_information("protein")
			transcript_cds = self.grab_transcript_sequence_information("cds")
			transcript_genomic = self.grab_transcript_sequence_information("genomic")

			transcript_cds_translated = utilities_codons.translate_oligo(transcript_cds['seq']).replace("*","")

			if self.options['peptide_sequence'] != transcript_cds_translated[self.options['peptide_start']:self.options['peptide_end']]:
				check_transcript_start = True

			uniprot_sequence = queryRunner.queryRunner("uniprot","parse_basic",{"accession":self.options['accession'][0]}).run()['data']['sequence']
			
			if "updated_uniprot_sequence" in self.options:
				logging.error("UniProt Transcript mismatch - updating to transcript sequence")
				uniprot_sequence = self.options['updated_uniprot_sequence']
				self.options['peptide_sequence'] = uniprot_sequence
				self.options['peptide_start'] = 0
				self.options['peptide_end'] = len(self.options['updated_uniprot_sequence'])
				check_transcript_start = False
				
			if check_transcript_start:
				if transcript_cds_translated.count(uniprot_sequence) == 0:
					alignment = utilities_alignment.alignPeptides(transcript_cds_translated,uniprot_sequence)

					if isinstance(alignment[0][0],list):
						transcript_cds_translated_align = "".join(alignment[0][0])
						uniprot_sequence_align = "".join(alignment[0][1])
					else:
						transcript_cds_translated_align = alignment[0][0]
						uniprot_sequence_align = alignment[0][1]

					identity = utilities_alignment.sequenceIdentity(transcript_cds_translated_align,uniprot_sequence_align)
			

					if identity.count("+")/len(uniprot_sequence) > 0.90:
						if uniprot_sequence_align[0] == "-":
							self.options['peptide_start'] = self.options['peptide_start'] + len(uniprot_sequence_align) - len(uniprot_sequence_align.lstrip('-'))
							self.options['peptide_start'] = self.options['peptide_end'] + len(uniprot_sequence_align) - len(uniprot_sequence_align.lstrip('-'))
						if transcript_cds_translated_align[0] == "-":
							self.options['peptide_start'] = self.options['peptide_start'] + -(len(transcript_cds_translated_align) - len(transcript_cds_translated_align.lstrip('-')))
							self.options['peptide_end'] = self.options['peptide_end'] + -(len(transcript_cds_translated_align) - len(transcript_cds_translated_align.lstrip('-')))
					else:
						logging.error("Sequence mismatch " + transcript_id + " " + self.options['accession'][0])
						logging.error(uniprot_sequence)
						logging.error(transcript_cds_translated)

						logging.error("Unretreivable sequence mismatch " + transcript_id + " " + self.options['accession'][0])

						return {"status":"error","error_type":"Unretreivable sequence mismatch " + transcript_id + " " + self.options['accession'][0]}
				else:
					self.options['peptide_start'] = self.options['peptide_start']+len(transcript_cds_translated.split(uniprot_sequence)[0])
					self.options['peptide_end'] = self.options['peptide_end']+len(transcript_cds_translated.split(uniprot_sequence[0]))
				
			#----#

			region_start = max(0,self.options['peptide_start'] - 1)
			region_end = self.options['peptide_end']
		
			#----#
			if self.options['peptide_sequence'] != transcript_cds_translated[region_start:region_end] and self.options['peptide_sequence'][1:] != transcript_cds_translated[region_start:region_end][1:]:
				logger.debug(self.options['peptide_sequence'] + " not correct position in translated transcript "  + transcript_cds_translated[region_start:region_end] + ':' + str(region_start) + "-"  +str(region_start))
				
				if transcript_cds_translated.count(self.options['peptide_sequence']) > 0:
					region_start = len(transcript_cds_translated.split(self.options['peptide_sequence'])[0])
					region_end = region_start + len(self.options['peptide_sequence'])
					logger.debug("Positions updated to " + str(region_start) + " " + str(region_start))
				elif len(transcript_cds_translated[:-1]) == len(self.options['peptide_sequence']):
					region_start = 0
					region_end = len(self.options['peptide_sequence'])
					logger.debug("Positions updated to " + str(region_start) + " " + str(region_start))
				else:		
					logger.error("No matching peptide found " + str(self.options['ensembl_gene_id']) + " " + str(self.options['accession']))
					logger.error(transcript_cds_translated)
					logger.error(self.options['peptide_sequence'])
					return {"status":"error","error_type":"No matching peptide found"}
			
			peptide = transcript_protein['seq'][region_start:region_end]
			peptide_cds = transcript_cds['seq'][region_start*3:region_end*3]

			if len(transcript_cds['seq'].split(peptide_cds)) == 2:
				peptide_cds_start = len(transcript_cds['seq'].split(peptide_cds)[0])
				peptide_cds_end = len(transcript_cds['seq'].split(peptide_cds)[0]) + len(peptide_cds)
			else:
				logger.debug("Found oligo (" + peptide_cds + ") multiple times in transcript_cds " + str(len(transcript_cds['seq'].split(peptide_cds))))
				peptide_cds_start = region_start*3#len(transcript_cds['seq'].split(peptide_cds)[0])
				peptide_cds_end = region_end*3#len(transcript_cds['seq'].split(peptide_cds)[0]) + len(peptide_cds)

			peptide_cds_complementary = utilities_codons.get_complementry_oligo(peptide_cds)
			peptide_cds_translated = utilities_codons.translate_oligo(peptide_cds)

			peptide_translation_check = peptide_cds_translated == peptide
				
			#######################################################################################################

			exon_boundary_overlap = False
			peptide_genomic_check = False
			peptide_cds_gapped = ""
			peptide_genomic_translated = "x"
			peptide_cds_aligned_region = ""

			if len(transcript_genomic['seq'].split(peptide_cds)) >= 2:
				genomic_coding_region_start = False
				genomic_coding_region_end = False
				
				if len(transcript_genomic['seq'].split(peptide_cds)) == 2:
					genomic_coding_region_start = len(transcript_genomic['seq'].split(peptide_cds)[0])
					genomic_coding_region_end = len(transcript_genomic['seq'].split(peptide_cds)[0]) + len(peptide_cds)
				else:
					logger.debug("Found oligo (" + peptide_cds + ") multiple times (" + str(len(transcript_genomic['seq'].split(peptide_cds)))  + ") in transcript_genomic of " + self.options['accession'][0] + " - " + str(len(transcript_genomic['seq'].split(peptide_cds))) + " - " + str(len(transcript_genomic['seq'])))
					alignment = self.align_exons_to_genomic_transcript(transcript_id,transcript_genomic,transcript_cds)
					
					transcript_cds_alignment_offset = 0

					for i in range(0,len(alignment[0])):
						if alignment[1][i] != "-":
							if transcript_cds_alignment_offset == peptide_cds_start and genomic_coding_region_start == False:
								genomic_coding_region_start = i
							if transcript_cds_alignment_offset == peptide_cds_end and genomic_coding_region_end == False:
								genomic_coding_region_end = i

							transcript_cds_alignment_offset += 1
					
					if genomic_coding_region_start == False and genomic_coding_region_end == False:
						pass
					else:
						peptide_cds_aligned_region = alignment[1][genomic_coding_region_start:genomic_coding_region_end]

				peptide_genomic_translated = utilities_codons.translate_oligo(transcript_genomic['seq'][genomic_coding_region_start:genomic_coding_region_end])
				peptide_cds_gapped = peptide_cds
				
				if peptide_cds_aligned_region.count("-") > 0: 
					exon_boundary_overlap = True
				elif genomic_coding_region_start == False or genomic_coding_region_end == False:
					logger.error("Possible exon boundary overlap - No genomic_coding_region_start/genomic_coding_region_end  - " + transcript_id + " - " + self.options['accession'][0] + " - " + peptide_cds + " " + transcript_genomic['seq'][genomic_coding_region_start:genomic_coding_region_end] + " " + str(self.options['peptide_start']) + ":" + str(self.options['peptide_end']))
					peptide_genomic_check = False
				elif peptide_cds != transcript_genomic['seq'][genomic_coding_region_start:genomic_coding_region_end]:
					logger.error("Possible exon boundary overlap - cds transcript and genomic transcript do not match - " + transcript_id + " - " + self.options['accession'][0] + " - " + peptide_cds + " " + peptide_cds_aligned_region)
			
					peptide_genomic_check = False
				else:
					peptide_genomic_check = True
			else:
				logger.debug("Found oligo 0 times in transcript_genomic: " + peptide_cds + " " + str(len(transcript_genomic['seq'].split(peptide_cds))))
				exon_boundary_overlap = True
			
			#######################################################################################################
		
			if exon_boundary_overlap:
				logger.debug("Can't find oligo in genome - aligning to find intron-exon boundary")
				
				exon_boundary_overlap = True
				peptide_genomic_translated = ''
				peptide_genomic_coding_oligo = "" 
				peptide_genomic_oligo = "" 
				cds_counter = 0
				genomic_counter = 0
				genomic_coding_region_start = False
				genomic_coding_region_end = False

				###---------------------------------###

				logger.debug("Aligned - " + str(len(transcript_genomic['seq'])) + " " + str(len(transcript_cds['seq'])))

				if len(transcript_genomic['seq']) > 5000:
					cds_genomic_alignment_path = os.path.join(self.options["cds_genomic_alignment_path"],transcript_id + ".exon_mapping.json")

					if not os.path.exists(cds_genomic_alignment_path):
						alignment = self.align_exons_to_genomic_transcript(transcript_id,transcript_genomic,transcript_cds)
					else:
						logger.debug("Reading align_exons_to_genomic_transcript - " + cds_genomic_alignment_path)
						with open(cds_genomic_alignment_path) as outfile:
							alignment = json.load(outfile)
				else:
					cds_genomic_alignment_path = os.path.join(self.options["cds_genomic_alignment_path"],transcript_id + ".json")
					if not os.path.exists(cds_genomic_alignment_path):
						logger.debug("Making " + cds_genomic_alignment_path)
						
						try:
							alignment = pairwise2.align.globalms(transcript_genomic['seq'], transcript_cds['seq'], 2, -1, -.5, -.1)
							alignment = alignment[0]
						except:
							logging.error("Alignment failed using align_exons_to_genomic_transcript")
							alignment = self.align_exons_to_genomic_transcript(transcript_id,transcript_genomic,transcript_cds)
						
						with open(cds_genomic_alignment_path,'w') as outfile:
							json.dump(alignment, outfile)
					else:
						logger.debug("Reading pairwise2 - " + cds_genomic_alignment_path)
						with open(cds_genomic_alignment_path) as outfile:
							alignment = json.load(outfile)
				
				###---------------------------------###

				peptide_cds_gapped = ""
				
				for i in range(0,len(alignment[0])):
					if cds_counter >= peptide_cds_start and alignment[1][i] != "-" and genomic_coding_region_start == False:
						genomic_coding_region_start = genomic_counter
					
					if cds_counter == peptide_cds_end:
						genomic_coding_region_end = genomic_counter

					if genomic_coding_region_start != False and genomic_coding_region_end == False:
						peptide_genomic_oligo += alignment[0][i]
						peptide_cds_gapped += alignment[1][i]
						if alignment[1][i] != "-":
							peptide_genomic_coding_oligo += alignment[0][i]
							if len(peptide_genomic_coding_oligo)%3 == 0:
								peptide_genomic_translated += utilities_codons.translate_oligo(peptide_genomic_coding_oligo[-3:])
						else:
							if i %3 == 0:
								peptide_genomic_translated += "x"

					if alignment[1][i] != "-":
						cds_counter += 1
					
					if alignment[0][i] != "-":
						genomic_counter += 1

				

				logger.debug("Aligned - " + str(len(transcript_genomic['seq'])) + " " + str(len(transcript_cds['seq'])))
				logger.debug("Aligned - gaps in genomic=" + str(alignment[0].count("-")))

				del alignment


			#######################################################################################################
			
			if len(peptide_cds_gapped) > 0:
				length_adjust_3 = 0
				length_adjust_5 = 0
				
				peptide_cds_gapped_flank_gap_stripped = peptide_cds_gapped.rstrip('-').lstrip('-')
				
				if peptide_cds_gapped[0] == "-":
					length_adjust_3 = len(peptide_cds_gapped) - len(peptide_cds_gapped.lstrip('-'))
					peptide_cds_gapped = peptide_cds_gapped.lstrip('-')
				
				if peptide_cds_gapped[-1] == "-":
					length_adjust_5 = len(peptide_cds_gapped) - len(peptide_cds_gapped.rstrip('-'))
					peptide_cds_gapped = peptide_cds_gapped.rstrip('-')
				
				peptide_genomic = transcript_genomic['seq'][genomic_coding_region_start:genomic_coding_region_end]
				
				if genomic_coding_region_end == False:
					peptide_genomic = transcript_genomic['seq'][:len(peptide_cds_gapped)]
				
				process_start = max(0,genomic_coding_region_start-self.options['flanks'])
						
				if len(peptide_genomic) != len(peptide_cds_gapped):
					peptide_genomic = peptide_genomic[:len(peptide_cds_gapped)]
					transcript_genomic_seq = transcript_genomic['seq']
					transcript_genomic_seq = transcript_genomic_seq + "-"*self.options['flanks']  ## TO FIX ISSUES WITH TEMINAL MOTIFS
					peptide_genomic_flanks = transcript_genomic_seq[process_start:genomic_coding_region_start + len(peptide_genomic) + self.options['flanks']]
				else:
					transcript_genomic_seq = transcript_genomic['seq']
					transcript_genomic_seq = transcript_genomic_seq + "-"*self.options['flanks'] 
					peptide_genomic_flanks = transcript_genomic_seq[process_start:genomic_coding_region_end + self.options['flanks']]
		
				if genomic_coding_region_start-self.options['flanks'] < 0:
					peptide_genomic_flanks = "N"*abs(genomic_coding_region_start-self.options['flanks']) + peptide_genomic_flanks

				peptide_genomic_flanks_complementary = utilities_codons.get_complementry_oligo(peptide_genomic_flanks)

				peptide_genomic_translated_check = peptide_cds_translated == peptide_genomic_translated.replace("x","")
				peptide_check = peptide == self.options['peptide_sequence']
			else:
				logger.error("peptide_cds_gapped length is 0 " + str(self.options['peptide_start']) + ":" + str(self.options['peptide_end']))

			#######################################################################################################
			
			return {
				"status":"success",
				"transcript_id":transcript_id,
				"gene_id":self.options['ensembl_gene_id'],
				"peptide":peptide,
				"peptide_cds_translated":peptide_cds_translated,
				"peptide_genomic_translated":peptide_genomic_translated,
				"peptide_cds":peptide_cds,
				"peptide_cds_start":peptide_cds_start,
				"peptide_cds_end":peptide_cds_end,
				"peptide_cds_gapped":peptide_cds_gapped,
				"peptide_genomic":peptide_genomic,
				"peptide_genomic_flanks":peptide_genomic_flanks,
				"peptide_genomic_flanks_complementary":peptide_genomic_flanks_complementary,
				"genomic_coding_region_start":genomic_coding_region_start,
				"genomic_coding_region_end":genomic_coding_region_end,
				"peptide_cds_complementary":peptide_cds_complementary,
				"peptide_length":len(transcript_protein['seq']),
				"transcript_cds_length":len(transcript_cds['seq']),			
				"transcript_genomic_length":len(transcript_genomic['seq']),
				"peptide_check":peptide_check,
				"peptide_translation_check":peptide_translation_check,
				"peptide_genomic_check":peptide_genomic_check,
				#"peptide_genomic_complementary_check":peptide_genomic_complementary_check,
				"peptide_genomic_translated_check":peptide_genomic_translated_check,
				"exon_boundary_overlap":exon_boundary_overlap
			}
		except:
			utilities_error.printError()
			return {"status":"error","error_type":"get_ensembl_data_for_region data parsing error"}


	########### OLIGO SEARCH FUNCTIONS ########################## OLIGO SEARCH FUNCTIONS ########################## OLIGO SEARCH FUNCTIONS ###############
	########### OLIGO SEARCH FUNCTIONS ########################## OLIGO SEARCH FUNCTIONS ########################## OLIGO SEARCH FUNCTIONS ###############
	########### OLIGO SEARCH FUNCTIONS ########################## OLIGO SEARCH FUNCTIONS ########################## OLIGO SEARCH FUNCTIONS ###############
	########### OLIGO SEARCH FUNCTIONS ########################## OLIGO SEARCH FUNCTIONS ########################## OLIGO SEARCH FUNCTIONS ###############

	def map_oligo_to_gene(self):
		oligo_data = self.find_sequence_human_genome()

		oligo_mapping = {}
		for oligo in oligo_data:
			if 'matched_genes' in oligo_data[oligo]:
				#pprint.pprint(oligo_data[oligo]['matched_genes'])
				oligo_mapping[oligo] = list(oligo_data[oligo]['matched_genes'].keys())
			else:
				oligo_mapping[oligo] = []

		return oligo_mapping
	
	def find_genome_region_sequence(self):

		region_id = "_".join(
			[
				str(self.options['chromosome']),
				self.options['region_strand'],
				str(min(self.options['region_start'],self.options['region_end'])),
				str(max(self.options['region_start'],self.options['region_end']))
			])

		logger.debug("Grabbing " + region_id)
		out_path = os.path.join(self.options['genomes_region_path'],region_id + ".region_sequence.tdt")

		if not os.path.exists(out_path):
			try:
				blast_results = subprocess.run([
					"blastdbcmd",
					"-db",self.options['blast_db_path'],
					"-dbtype","nucl",
					"-strand",self.options['region_strand'],
					"-entry",self.options['chromosome'],
					"-range",str(self.options['region_start']) + "-" + str(self.options['region_end']), 
					], 
					stdout=subprocess.PIPE, 
					stderr=subprocess.PIPE
				)

				region_sequence = blast_results.stdout.decode().strip().split("\n")[-1]
				logger.debug(region_id + " -> " + region_sequence)

				logger.debug("Writing " + out_path)
				open(out_path,'w').write(region_sequence)
				return region_sequence
			except:
				utilities_error.printError()
				logger.error("blastdbcmd failed")
				return ""
		else:
			logger.debug("Reading " + out_path)
			return open(out_path).read()

	def find_sequence_human_genome(self):
		self.options['blast_db_path'] = os.path.join(self.options['data_path'],"genomes","GRCh38_dna_primary_assembly")
		
		# makeblastdb -dbtype nucl  -in "/home/data/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa" -out GRCh38_dna_primary_assembly -parse_seqids
		# blastn -task megablast -evalue 1e-4 -num_threads 1 -dust yes -db GRCh38_dna_primary_assembly -query /home/data/genomes/gRNA/gRNA_counts/gRNAs.fasta -out /home/data/genomes/gRNA/gRNA_counts/gRNAs.blast.out -outfmt 7
		
		tag = ""
		if self.options['use_filtering'] == False:
			tag = ".no_complexity_filtering"

		logger.debug("Searching " + str(len(self.options['oligo_sequence'])) + " guides")
	
		guide_list = []
		for guide in self.options['oligo_sequence']:
			guides_counts_path = os.path.join(self.options['genomes_guides_blast_path'],"counts", guide + tag + ".counts.json")
			guides_counts_coordinates_path = os.path.join(self.options['genomes_guides_blast_path'], "counts_coordinates",guide + tag + ".counts_coordinates.json")
			genome_guides_mapping_path = os.path.join(self.options['genomes_guides_mapping_path'], guide + tag + ".gRNAs.blast.json")
		
			if self.options['remake']:
				guide_list.append(guide)
			elif self.options["counts_only"] and self.options["counts_output_type"] == "coordinates" and not os.path.exists(guides_counts_coordinates_path):
				guide_list.append(guide)
			elif self.options["counts_only"] and self.options["counts_output_type"] != "coordinates" and not os.path.exists(guides_counts_path):
				guide_list.append(guide)
			elif not self.options["counts_only"] and not os.path.exists(genome_guides_mapping_path) or self.options['remake']:
				guide_list.append(guide)
			#else:
			#	logger.debug("Skipping " + guide + " data exists")
		
		logger.debug("Processing " + str(len(guide_list)) + " guides unprocessed")
		logger.debug("Skipping " + str(len(self.options['oligo_sequence']) - len(guide_list)) + " guides")
		
		add_reverse_complementry_oligo = False

		if len(guide_list) > 0:
			counter = 0
			chunks = self.options['chunks_size']
			
			for i in range(0,len(guide_list),chunks):
				try:
					logger.debug("chunk " + str(i) + " of " + str(len(guide_list)))

					if i+chunks >= len(guide_list):
						guide_list_chunk = guide_list[i:]
					else:
						guide_list_chunk = guide_list[i:i+chunks]
					
					hash = utilities_basic.params_to_hash({"sequences":guide_list_chunk})

					tag = ""
					if self.options['use_filtering'] == False:
						tag = ".no_complexity_filtering"
					
					guides_fasta_out_path = os.path.join(self.options['genomes_guides_blast_path'],'raw', hash + tag + ".gRNAs.blast.out")
					guides_fasta_out_zipped_path = os.path.join(self.options['genomes_guides_blast_path'], 'raw',hash + tag + ".gRNAs.blast.out.zip")
					guides_fasta_path = os.path.join(self.options['genomes_guides_fasta_path'], hash + tag + ".gRNAs.fasta")

					if not os.path.exists(guides_fasta_out_path) and not os.path.exists(guides_fasta_out_zipped_path):
						fasta_out = ""
						counter = 0
						for oligo in guide_list_chunk:
							counter += 1
							
							fasta_out += ">gRNA_" + oligo + "_for\n"
							fasta_out += oligo + "\n"

							if add_reverse_complementry_oligo:
								reverse_oligo = utilities_codons.get_complementry_oligo(oligo)[::-1]
								fasta_out += ">gRNA_" + reverse_oligo + "_rev\n"
								fasta_out += reverse_oligo + "\n"
						

						logger.debug("Writing " + guides_fasta_path)
						open(guides_fasta_path,"w").write(fasta_out)

						# qacc sacc pident length mismatch gapopen qstart, qend, sstart, send, evalue, bitscore
						# qacc sacc evalue qstart qend sstart send sseq

						if self.options['use_filtering']:
							filtering_options = []
						else:
							filtering_options = [
								"-dust","no",
								"-soft_masking","false"
							]

						blast_results = subprocess.run([
							"blastn",
							"-task","megablast", #"blastn-short",#
							"-num_threads","1",
							"-db",self.options['blast_db_path'],
							"-query", guides_fasta_path,
							"-out",guides_fasta_out_path,
							#23"-outfmt",'"7 delim=@ qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq sseq"',
							"-outfmt",'7 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq sseq',
							#"-outfmt",'"7 delim=@ qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq sseq"',
							"-word_size","7",
							"-evalue",self.options['evalue'],
							"-ungapped"] + filtering_options, 
							stdout=subprocess.PIPE, 
							stderr=subprocess.PIPE
						)
						
						logger.debug(blast_results)
						logger.debug("\n".join(blast_results.stdout.decode().strip().split("\n")))
						logger.debug(blast_results.returncode)
						guides_counts = [len(blast_results.stdout.decode().strip().split("\n"))]
						logger.debug(guides_counts)

						if len(blast_results.stdout.decode().strip().split("\n")) < 2:
							logger.error('Check blastn version - outfmt options changes and may need to be updated - "7 delim=@ qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq sseq"')

					else:
						logger.debug("Blast data completed exists")
					
					####--------------------------------–####
					if not os.path.exists(guides_fasta_out_zipped_path) and os.path.exists(guides_fasta_out_path):
						logger.debug("Zipping " + guides_fasta_out_path + "->" + guides_fasta_out_zipped_path)
						with zipfile.ZipFile(guides_fasta_out_zipped_path, 'w', zipfile.ZIP_DEFLATED) as myzip:
							myzip.write(guides_fasta_out_path)
	
					if os.path.exists(guides_fasta_out_zipped_path) and os.path.exists(guides_fasta_out_path):
						logger.debug("Deleting zipped file source " + guides_fasta_out_path)
						os.remove(guides_fasta_out_path)
					
					if self.options["counts_only"]:
						self.read_blast_output_counts(guides_fasta_out_zipped_path)
					else:
						self.read_blast_output(guides_fasta_out_zipped_path)
				except:
					utilities_error.printError()
					logger.error("Failed: " + "chunk " + str(i) + " of " + str(len(guide_list)) + " " + guides_fasta_out_path)
					raise

		###-----------------------------------###

		if self.options["counts_only"]:
			mismatch_counts = {}

			

			for guide in self.options['oligo_sequence']:
				if self.options["counts_output_type"] == "coordinates":
					if not os.path.exists(os.path.join(self.options['genomes_guides_blast_path'], "counts_coordinates")):
						os.mkdir(os.path.join(self.options['genomes_guides_blast_path'], "counts_coordinates"))

					guides_counts_path = os.path.join(self.options['genomes_guides_blast_path'], "counts_coordinates",guide + tag + ".counts_coordinates.json")
				else:
					if not os.path.exists(os.path.join(self.options['genomes_guides_blast_path'], "counts")):
						os.mkdir(os.path.join(self.options['genomes_guides_blast_path'], "counts"))

					guides_counts_path = os.path.join(self.options['genomes_guides_blast_path'],"counts", guide + tag + ".counts.json")

				if os.path.exists(guides_counts_path):
					with open(guides_counts_path) as outfile:
						try:
							mismatch_counts[guide] = json.load(outfile)
						except:
							logger.error(guides_counts_path + " does not parse correctly")
				else:
					logger.debug(guides_counts_path + " not found")

			return mismatch_counts
		
		###-----------------------------------###

		guide_mapping = {}
		for guide in self.options['oligo_sequence']:
			genomes_guides_path = os.path.join(self.options['genomes_guides_mapping_path'], guide + ".gRNAs.blast.json")
			
			if os.path.exists(genomes_guides_path):
				with open(genomes_guides_path) as outfile:
					guide_mapping[guide] = json.load(outfile)
			else:
				guide_mapping[guide] = {}
			
		if self.options['annotate_all_exact_matches']:
			try:
				guide_mapping = self.annotate_oligos(guide_mapping)
			except:
				pass
			
			try:
				guide_mapping = self.read_annotate_blast_output(guide_mapping)
			except:
				pass
			
		return guide_mapping

	###-------------------------------------------------###

	def read_blast_output_counts(self,guides_fasta_out_zipped_path):
		blast_headers = ['query','chromosome','identity','length','mismatches','gaps','query_start','query_end','match_start','match_end','e_value','bit_score','sstrand','qseq','sseq']
		
		logger.debug("Reading " + guides_fasta_out_zipped_path)

		tag = ""
		if self.options['use_filtering'] == False:
			tag = ".no_complexity_filtering"

		with zipfile.ZipFile(guides_fasta_out_zipped_path, 'r') as zip_ref:
			# Get the list of files in the zip archive
			file_list = zip_ref.namelist()
			
			lines = zip_ref.read(file_list[0]).decode().strip().split("\n")
			
			if self.options["counts_only"]:
				mismatch_counts = {}
				for line in lines:
					try:
						if len(line) == 0: continue
						if line[0] != "#":
							line_bits = line.split()
							try:
								guide = line_bits[0].split("_")[1]
							except:
								continue
							
							if guide not in mismatch_counts:
								mismatch_counts[guide] = {}

							match_detail = {}
							for blast_header_iter in range(0,len(blast_headers)):#,len(blast_headers))
								if blast_header_iter < len(line_bits):
									match_detail[blast_headers[blast_header_iter]] = line_bits[blast_header_iter]
								else:
									logger.error("Error: " + line + " " + str(blast_header_iter) + " " + guides_fasta_out_zipped_path)
									logger.error("Malformed BLAST output")
									os.remove(guides_fasta_out_zipped_path)
									logger.error("Deleting " + guides_fasta_out_zipped_path)
									return {}
									
							match_detail['guide_match'] = "-"*(int(match_detail['query_start'])-1) + match_detail['sseq'] + "-"*(len(guide)-int(match_detail['query_end']))
							mismatches = [match_detail['guide_match'][i] == guide[i] for i in range(0,len(guide))].count(False)
							
							if mismatches < 4:
								if self.options["counts_output_type"] == "coordinates":
									guide_coordinates = match_detail['chromosome'] + '.' + match_detail['match_start'] + ':' + match_detail['match_end']
									if mismatches not in mismatch_counts[guide]:
										mismatch_counts[guide][mismatches] = []

									mismatch_counts[guide][mismatches].append(guide_coordinates)
								else:
									if mismatches not in mismatch_counts[guide]:
										mismatch_counts[guide][mismatches] = 0

									mismatch_counts[guide][mismatches] += 1
								
					except:
						logger.debug("Unprocessed match: " + line)
						
				for guide in mismatch_counts:
					if self.options["counts_output_type"] == "coordinates":
						guides_counts_path = os.path.join(self.options['genomes_guides_blast_path'],"counts_coordinates",guide + tag + ".counts_coordinates.json")
					else:
						guides_counts_path = os.path.join(self.options['genomes_guides_blast_path'],"counts", guide + tag + ".counts.json")
						logger.debug(guide + "\t" + "\t".join([str(mismatch_counts[guide][count]) if count in mismatch_counts[guide] else "0" for count in [0,1,2,3]]) + "\t" + guides_counts_path)
						
					with open(guides_counts_path,'w') as outfile:
						json.dump(mismatch_counts[guide], outfile)

				return mismatch_counts
			
	def read_blast_output(self,guides_fasta_out_zipped_path):
		match_details = {}
		blast_headers = ['query','chromosome','identity','length','mismatches','gaps','query_start','query_end','match_start','match_end','e_value','bit_score','sstrand','qseq','sseq']
		
		logger.debug("Reading " + guides_fasta_out_zipped_path)

		with zipfile.ZipFile(guides_fasta_out_zipped_path, 'r') as zip_ref:
			# Get the list of files in the zip archive
			file_list = zip_ref.namelist()
			
			lines = zip_ref.read(file_list[0]).decode().strip().split("\n")
			
			for line in lines:
				if len(line) == 0: continue
				if line[0] != "#":
					try:
						line_bits = line.split()
						guide = line_bits[0].split("_")[1]
						
						if guide not in match_details:
							match_details[guide] = {
								"matches":{},
								"exact_matches":0,
								"exact_matches_guide_coordinates":[],
								"mismatches":{},
								"off_target_mapping":{},
								"off_target_type":{},
								"matched_genes":{}
							}

						match_detail = {}
						for blast_header_iter in range(0,len(blast_headers)):#,len(blast_headers))
							match_detail[blast_headers[blast_header_iter]] = line_bits[blast_header_iter]
						
						guide_coordinates = match_detail['chromosome'] + '.' + match_detail['match_start'] + ':' + match_detail['match_end']
						
						if guide_coordinates in match_details[guide]["matches"]:
							pass
						else:
							match_detail['guide_match'] = "-"*(int(match_detail['query_start'])-1) + match_detail['sseq'] + "-"*(len(guide)-int(match_detail['query_end']))

							pad_5 = int(match_detail['query_start'])-1
							pad_3 = len(guide)-int(match_detail['query_end'])
							
							if match_detail['sstrand'] == "plus":
								updated_match_start =  int(match_detail['match_start']) - pad_5
								updated_match_end = int(match_detail['match_end']) + pad_3
							else:
								updated_match_start =  int(match_detail['match_start']) + pad_5
								updated_match_end = int(match_detail['match_end']) - pad_3 

							self.options['region_strand'] = match_detail['sstrand']
							self.options['chromosome'] = match_detail['chromosome']
							self.options['region_start'] = min(updated_match_start, updated_match_end)
							self.options['region_end'] = max(updated_match_start, updated_match_end)

							mismatches = [match_detail['guide_match'][i] == guide[i] for i in range(0,len(guide))].count(False)
							
							if self.options['exact_matches_only']:
								if mismatches > 0:
									continue

							if mismatches == 0:
								logger.debug("0 mismatches: " + line)

								if match_detail['guide_match'] == guide:
									aligned_sequence = guide
								elif len(match_detail['sseq']) == len(guide):
									aligned_sequence = match_detail['sseq']
								else:
									logger.debug("Aligning to genome: " + match_detail['sseq'] + " - " + str(len(match_details)) + " of " + str(len(lines)))	
									aligned_sequence = self.find_genome_region_sequence()

								match_details[guide]["aligned_sequence"] = aligned_sequence

								if aligned_sequence[0:20] == guide[0:20] and aligned_sequence[-2] == guide[-2]:# and match_detail['length'] ==  str(len(guide)):
									match_details[guide]["exact_matches"] += 1
									match_details[guide]["exact_matches_guide_coordinates"].append(guide_coordinates)
									match_detail["guide_pam_match"] = True

								match_detail["mismatches"] = [aligned_sequence[i] == guide[i] for i in range(0,len(guide))].count(False)
								match_detail["guide_mismatches"] = [aligned_sequence[i] == guide[i] for i in range(0,len(guide[0:20]))].count(False)
								match_details[guide]["mismatches"][guide_coordinates] = [aligned_sequence[i] == guide[i] for i in range(0,len(guide))].count(False)

							else:
								match_detail["mismatches"] = [match_detail['guide_match'][i] == guide[i] for i in range(0,len(guide))].count(False)
								match_detail["guide_mismatches"] = [match_detail['guide_match'][i] == guide[i] for i in range(0,len(guide[0:20]))].count(False)
								match_details[guide]["mismatches"][guide_coordinates] = [match_detail['guide_match'][i] == guide[i] for i in range(0,len(guide))].count(False)
								
							match_details[guide]["matches"][guide_coordinates] = match_detail
					except:
						utilities_error.printError()
						
		###-------------------------------------------------###
		
		match_details = self.read_annotate_blast_output(match_details)
		match_details = self.annotate_oligos(match_details)

		###-------------------------------------------------###

		for guide in match_details:
			genome_guides_mapping_path = os.path.join(self.options['genomes_guides_mapping_path'], guide + ".gRNAs.blast.json")
			logger.debug("Writing " + genome_guides_mapping_path)
			with open(genome_guides_mapping_path,'w') as outfile:
				json.dump(match_details[guide], outfile)

	###-------------------------------------------------###

	def annotate_oligos(self,match_details):
		if "sequence_flanks_5'" not in self.options:
			self.options["sequence_flanks_5'"] = 0
		if "sequence_flanks_3'" not in self.options:
			self.options["sequence_flanks_3'"] = 0
			
		add_off_target_effects = True

		if len(match_details) > 1000:
			logger.error("Too many matches to annotate for :" + str(len(match_details)))
		else:
			counter = 0
			for guide in match_details:
				counter += 1
				logger.debug("Annotating oligo: " + guide + " " + str(counter) + "/" + str(len(match_details)))
					
				if 'off_target_mapping' in match_details[guide]:
					if len(match_details[guide]['off_target_mapping']) > 0:
						continue

				if add_off_target_effects:
					match_details[guide]['off_target_mapping'] = {}
					match_details[guide]['off_target_type'] = {}

					if match_details[guide]['exact_matches'] > 1 or self.options['annotate_all_exact_matches']:
						exact_matches_guide_coordinates_counter = 0
						if len(match_details[guide]['exact_matches_guide_coordinates']) > 1000:
							logger.error("Too many matches to annotate for :" + str(len(match_details[guide]['exact_matches_guide_coordinates'])))
						else:
							for match in match_details[guide]['exact_matches_guide_coordinates']:
								exact_matches_guide_coordinates_counter += 1
								match_details[guide]['off_target_mapping'][match] = []
								coordinate_match_details = match_details[guide]['matches'][match]
								
								options = {
									"chromosome":coordinate_match_details["chromosome"],
									"region_start":str(int(coordinate_match_details["match_start"]) - self.options["sequence_flanks_5'"]),
									"region_end":str(int(coordinate_match_details["match_end"]) + self.options["sequence_flanks_3'"]),
									"feature_type":["gene","exon","transcript"]
								}
								
								logger.debug("Annotating grab genomic region features: " + guide + " " + str(exact_matches_guide_coordinates_counter) + "/" + str(len( match_details[guide]['exact_matches_guide_coordinates'])))

								ensembl_response = queryRunner.queryRunner("ensembl","grab_genomic_region_features",options).run()
								
								if 'data' in ensembl_response:
									for feature in ensembl_response['data']:
										if feature['feature_type'] == "gene":
											external_name = ""
											biotype = ""
											
											if 'external_name' in feature:
												external_name = feature['external_name']
											if 'biotype' in feature:
												biotype = feature['biotype']

											match_details[guide]['off_target_mapping'][match].append({"name":external_name,"type":biotype})

											if feature['biotype'] not in match_details[guide]['off_target_type']:
												match_details[guide]['off_target_type'][feature['biotype']] = 0

											match_details[guide]['off_target_type'][feature['biotype']] += 1
								else:
									pprint.pprint(ensembl_response)
		
		return match_details

	###-------------------------------------------------###

	def read_annotate_blast_output(self,match_details):
		for guide in match_details:
			try:
				if "sequence_flanks_5'" not in self.options:
					self.options["sequence_flanks_5'"] = 0
				if "sequence_flanks_3'" not in self.options:
					self.options["sequence_flanks_3'"] = 0

				if match_details[guide]["exact_matches"] > 1 or self.options['annotate_all_exact_matches']:
					exons_parents = {}
					transcript_mapping = {}

					guide_coordinates_list = match_details[guide]["exact_matches_guide_coordinates"]  + list(match_details[guide]['mismatches'].keys())

					if len(guide_coordinates_list) > 2000:
						logger.error("Too many matches to annotate for " + guide + ":" + str(len(guide_coordinates_list)))
					else:	
						guide_coordinates_counter = 0
						for guide_coordinates in guide_coordinates_list:
							try:
								guide_coordinates_counter += 1
								logger.debug("Annotating " + guide + " " + str(guide_coordinates_counter) + "/" + str(len(guide_coordinates_list)) + " " + guide_coordinates)
								match_details[guide]['off_target_mapping'][guide_coordinates] = []
								match_detail = match_details[guide]["matches"][guide_coordinates]

								options = {
									"chromosome":match_detail["chromosome"],
									"region_start":str(int(match_detail["match_start"]) - self.options["sequence_flanks_5'"]),
									"region_end":str(int(match_detail["match_end"]) + self.options["sequence_flanks_3'"]),
									"feature_type":["gene","exon","transcript"]
								}
								ensembl_response = queryRunner.queryRunner("ensembl","grab_genomic_region_features",options).run()
								
								if 'data' not in ensembl_response:
									if 'error' not in match_details[guide]['off_target_type']:
										match_details[guide]['off_target_type']['error'] = 0

									match_details[guide]['off_target_type']['error'] += 1

								elif len(ensembl_response['data']) == 0:
									if 'intergenic' not in match_details[guide]['off_target_type']:
										match_details[guide]['off_target_type']['intergenic'] = 0

									match_details[guide]['off_target_type']['intergenic'] += 1
								else:
									for feature in ensembl_response['data']:
										if feature['feature_type'] == "transcript":
											transcript_mapping[feature['id']] = {"gene_id":feature['Parent'],"biotype":feature['biotype']}

									for feature in ensembl_response['data']:
										if feature['feature_type'] == "exon":
											if transcript_mapping[feature['Parent']]["gene_id"] not in exons_parents:
												exons_parents[transcript_mapping[feature['Parent']]["gene_id"]] = {
													"coding":False,
													"exons":{},
													"exons_parent_transcript":{}
													}
											
											if transcript_mapping[feature['Parent']]["biotype"] == "protein_coding":
												exons_parents[transcript_mapping[feature['Parent']]["gene_id"]]['coding'] = True
											
											exons_parents[transcript_mapping[feature['Parent']]["gene_id"]]['exons'][feature['id']] = feature
											exons_parents[transcript_mapping[feature['Parent']]["gene_id"]]['exons_parent_transcript'][feature['id']] = transcript_mapping[feature['Parent']]

									for feature in ensembl_response['data']:	
										if feature['feature_type'] == "gene":
											external_name = ""
											biotype = ""
										
											if 'external_name' in feature:
												external_name = feature['external_name']
											if 'biotype' in feature:
												biotype = feature['biotype']

											exon_overlap = None
											if feature['id'] in exons_parents:
												exon_overlap = exons_parents[feature['id']]['coding']

											if exon_overlap:
												if external_name not in match_details[guide]["matched_genes"]:
													match_details[guide]["matched_genes"][external_name] = 0
												
												match_details[guide]["matched_genes"][external_name] +=1

											match_details[guide]['off_target_mapping'][guide_coordinates].append({
												"name":external_name,
												"id":feature['id'],
												"type":biotype,
												"exon_overlap":exon_overlap
												}
											)

											#if guide_coordinates in match_details[guide]['mismatches']:
											#	logger.debug(match_details[guide]['mismatches'][guide_coordinates],"\t",match_details[guide]['matches'][guide_coordinates]['guide_match'],"\t",guide,"\t",feature['id'],"\t",external_name,"\t",feature['biotype'])
											
											if feature['biotype'] not in match_details[guide]['off_target_type']:
												match_details[guide]['off_target_type'][feature['biotype']] = 0

											match_details[guide]['off_target_type'][feature['biotype']] += 1
							except:
								utilities_error.printError()

								if 'error' not in match_details[guide]['off_target_type']:
									match_details[guide]['off_target_type']['error'] = 0

								match_details[guide]['off_target_type']['error'] += 1

						match_details[guide]['gene_exon_mapping'] = exons_parents
						match_details[guide]['gene_exon_matches'] = 0
						for gene_id in match_details[guide]['gene_exon_mapping']:
							if  match_details[guide]['gene_exon_mapping'][gene_id]['coding']:
								match_details[guide]['gene_exon_matches'] += 1
			except:
				pass

		return match_details