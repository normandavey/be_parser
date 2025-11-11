import sys, re, os, json, pprint, inspect, copy
		

file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../guide_design/"))
import guideAnnotation 

sys.path.append(os.path.join(file_path,"../utilities/"))

import utilities_codons
import utilities_basic
import utilities_pssm
import option_reader
import utilities_error
import utilities_gzip


sys.path.append(os.path.join(file_path,"../data_management/"))
import queryRunner

#-----
import logging
logging.basicConfig(level=logging.INFO)
#-----

##################################################################################################
##################################################################################################
##################################################################################################

__docs___ = """
Basic approach is:
- grab the paired reads for each insert
- find the designed flanks, discard anything without a designed flank
- find the gRNA, discard anything without a designed gRNA
- find the reporter barcode, discard anything without a barcode
- extract the target segment of reporter
- count the distinct target segment edits for each gRNA
next,
- calculate the frequency of each distinct target segment edits for each gRNA
- compare the target segment to the gRNA and calculate differences
- map difference to the protein level to define the protein edits
- merge the protein edits across the genomic target segment edits and calculate their frequency
"""

class baseEditingReporterParser():
	def __init__(self):
		self.guides_data = {}
		self.options = {
			'calculate_reporter_metrics':True,
			'map_unedited_reptorter_to_similar_reporter':False,
			"robust_guide_count":False,
			"use_full_reporter":False,
			"guide_annotation_file":"",		
			"experiment_design":{	
				"reporter":{
					"pair_sequencing_file":"",
					"oligo_format":{
						"insert_oligo_5'_flank":"GCGAATTCAAAC",
						"insert_oligo_3'_flank":"",
						"insert_oligo_5'_flank_requirement":"",
						"insert_oligo_3'_flank_requirement":"AA"
					},
					"insert_length":41,
					"insert_trim":{"start":11,"end":31},
					"add_barcode":True,
					"barcode_length":4,
					"use_complement":True,
					"sample_information":{
						"group":"1",
						"timepoint":"1",
						"replicate":"1"
					}
				},
				"guide":{
					"pair_sequencing_file":"",
					"oligo_format":{
						"insert_oligo_5'_flank":"GAAAGGACGAAACACCG",
						"insert_oligo_3'_flank":"",
						"insert_oligo_5'_flank_requirement":"",
						"insert_oligo_3'_flank_requirement":""
					},
					"insert_length":20,
					"insert_trim":{},
					"use_complement":False,
					"add_barcode":False,
					"sample_information":{
						"group":"1",
						"timepoint":"1",
						"replicate":"1"
					}
				}
			},
			"task":"count",
			"analysis_name":"undefined",
			"data_directory":"./",
			"rhAmpSeq_analysis_design_file":"",
			"remake":False,
			"subsample":False,
			"subsample_cutoff":1000000,
			"reporter_count_cutoff":0.01,
			"barcode_length":4,
			"subset_tag":"",
			"reporter_skip_no_edit_annotation":False,
			"guides_data_json_file":None
		}

		self.options.update(option_reader.load_commandline_options(self.options,{}))
		self.read_metrics = {}

	##################################################################################################

	def setup_analysis(self):	
		self.options["data_directory"] = os.path.join(self.options["data_directory"],"reporter")

		if not os.path.exists(self.options["data_directory"]):
			os.mkdir(self.options["data_directory"] )

		self.read_metrics = {
			"global_guide_edit_percentage":None,
			"not_designed_guide":0,
			"no_barcode":0,
			"incorrect_barcode":0,
			"no_guide":0,
			"no_reporter":0,
			"unfiltered":0,
			"processed":0
		}


	##################################################################################################

	def setup_rhAmpSeq_analysis(self):	
		self.read_metrics = {
			"global_guide_edit_percentage":None,
			"not_designed_guide":0,
			"no_barcode":0,
			"incorrect_barcode":0,
			"no_guide":0,
			"no_reporter":0,
			"unfiltered":0,
			"processed":0
		}

		if not os.path.exists(self.options["rhAmpSeq_analysis_design_file"]):
			logging.error(self.options["rhAmpSeq_analysis_design_file"]  + " does not exist")
			sys.exit()

		self.options ["experiment_design"]['amplicon'] = {}
		import pandas as pd

		delimiter = "\t"
		key_column = "guide"
		value_columns = ["guide","gene",'accession',"chrom","guide start","guide end","amplicon start","amplicon end","amplicon length (nt)"]

		df = pd.read_csv(self.options ["rhAmpSeq_analysis_design_file"], delimiter=delimiter)

		data_dict = {
			row[key_column]: {col: row[col] for col in value_columns}
			for _, row in df.iterrows()
		}

		self.options["amplicon_flanks_size"] = 0
		self.options['reporter_count_cutoff'] = 0.001
		
		for guide in data_dict:
			self.options["experiment_design"]['amplicon'][guide] = {
				"accession":data_dict[guide] ['accession'],
				"gene":data_dict[guide]['gene'],
				"oligo_format":{
					"insert_oligo_5'_flank":"",
					"insert_oligo_3'_flank":"",
					"insert_oligo_5'_flank_requirement":"",
					"insert_oligo_3'_flank_requirement":""
				},
				"insert_length":20,
				"insert_trim":{},
				"use_complement":False,
				"add_barcode":False,
				"sample_information":{
					"group":"1",
					"timepoint":"1",
					"replicate":"1"
				},
				"strand":""
			}

			dna_region_coordinate = "chr" + str(data_dict[guide]['chrom']) + ":" + str(data_dict[guide]['guide start']) + "-" + str(data_dict[guide]['guide end'])
			guide_region_sequence = queryRunner.queryRunner("ensembl","grab_dna_regions",{"dna_region_coordinates":dna_region_coordinate}).run()
			dna_region_coordinate = "chr" + str(data_dict[guide]['chrom']) + ":" + str(data_dict[guide]['amplicon start']) + "-" + str(data_dict[guide]['amplicon end'])
			region_sequence_data = queryRunner.queryRunner("ensembl","grab_dna_regions",{"dna_region_coordinates":dna_region_coordinate}).run()
			region_sequence = region_sequence_data['data'][0]['seq']
			
			flanks = {
				"insert_oligo_5'_flank":"",
				"insert_oligo_3'_flank":"",
			}

			if guide_region_sequence['data'][0]['seq'].count(guide) == 1:
				self.options["experiment_design"]['amplicon'][guide]['strand'] = "positive"
				if self.options["amplicon_flanks_size"] > 0:
					self.options["experiment_design"]['amplicon'][guide]['oligo_format']["insert_oligo_5'_flank"] = region_sequence.split(guide)[0][-10 - self.options["amplicon_flanks_size"]:-self.options["amplicon_flanks_size"]]
				else:
					self.options["experiment_design"]['amplicon'][guide]['oligo_format']["insert_oligo_5'_flank"] = region_sequence.split(guide)[0][-10 - self.options["amplicon_flanks_size"]:]

				self.options["experiment_design"]['amplicon'][guide]['oligo_format']["insert_oligo_3'_flank"] = region_sequence.split(guide)[1][self.options["amplicon_flanks_size"]:self.options["amplicon_flanks_size"] + 10]
				self.options["experiment_design"]['amplicon'][guide]['reporter_sequence'] =region_sequence.split(self.options["experiment_design"]['amplicon'][guide]['oligo_format']["insert_oligo_3'_flank"])[0].split(self.options["experiment_design"]['amplicon'][guide]['oligo_format']["insert_oligo_5'_flank"])[1]
			elif guide_region_sequence['data'][0]['seq'].count(utilities_codons.get_reverse_complementry_oligo(guide)) == 1:
				self.options["experiment_design"]['amplicon'][guide]['strand'] = "negative"
				self.options["experiment_design"]['amplicon'][guide]['use_complement'] = True
				if self.options["amplicon_flanks_size"] > 0:
					self.options["experiment_design"]['amplicon'][guide]['oligo_format']["insert_oligo_5'_flank"] = region_sequence.split(utilities_codons.get_reverse_complementry_oligo(guide))[0][-10 - self.options["amplicon_flanks_size"]:-self.options["amplicon_flanks_size"]]	
				else:
					self.options["experiment_design"]['amplicon'][guide]['oligo_format']["insert_oligo_5'_flank"] = region_sequence.split(utilities_codons.get_reverse_complementry_oligo(guide))[0][-10 - self.options["amplicon_flanks_size"]:]	
				
				self.options["experiment_design"]['amplicon'][guide]['oligo_format']["insert_oligo_3'_flank"] = region_sequence.split(utilities_codons.get_reverse_complementry_oligo(guide))[1][self.options["amplicon_flanks_size"]:self.options["amplicon_flanks_size"]+10]
				self.options["experiment_design"]['amplicon'][guide]['reporter_sequence'] = region_sequence.split(self.options["experiment_design"]['amplicon'][guide]['oligo_format']["insert_oligo_3'_flank"])[0].split(self.options["experiment_design"]['amplicon'][guide]['oligo_format']["insert_oligo_5'_flank"])[1]
			else:
				logging.error(guide)
			
			self.options["experiment_design"]['amplicon'][guide]['insert_length'] = len(self.options["experiment_design"]['amplicon'][guide]['reporter_sequence'] )


	##################################################################################################
	
	def is_high_quality(self,quality_str, threshold):
		"""Check if all bases meet the minimum Phred score."""
		return all(ord(q) - 33 >= threshold for q in quality_str)

	def read_fastq_gzipped_files(self):
		logging.debug("Reading zipped file")
		processed_oligos = 1
		
		import gzip
		
		self.oligo_list[self.options['pair_type']] = [] 

		self.read_metrics[self.options['pair_type']] = {
			"unique_counts":0,
			"total":0,
			"filtered_phred":0,
			"skipped_reporter":0,
			"no_flanks":0,
			"incorrect_oligo_insert_length":0,
			"incorrect_oligo_insert_length_distribution":{},
			"processed":0
		}

		if isinstance(self.options["fastq_zipped_file"],str):
			fastq_gzipped_file = self.options["fastq_zipped_file"]
		elif isinstance(self.options["fastq_zipped_file"],list):
			fastq_gzipped_file = self.options["fastq_zipped_file"][0]
				
		if self.options["subsample"]:
			data_json_path = os.path.join(fastq_gzipped_file.replace(".fastq.gz","") + '.' + str(self.options["subsample_cutoff"]) + self.options["subset_tag"] + "." + self.options['pair_type'] + '_list.reporter.jsonc')
		else:
			data_json_path = os.path.join(fastq_gzipped_file.replace(".fastq.gz","")  + self.options["subset_tag"] + '.' + self.options['pair_type'] + '_list.reporter.jsonc')
		
		self.options['re_split_reporter'] = False
		if sum([self.options['sequencing_data_design']["oligo_format"]["insert_oligo_5'_flank"].count(nuc) for nuc in ['A','T','C','G']]) != len(self.options['sequencing_data_design']["oligo_format"]["insert_oligo_5'_flank"]):
			self.options['re_split_reporter'] = True
		if sum([self.options['sequencing_data_design']["oligo_format"]["insert_oligo_3'_flank"].count(nuc) for nuc in ['A','T','C','G']]) != len(self.options['sequencing_data_design']["oligo_format"]["insert_oligo_3'_flank"]):
			self.options['re_split_reporter'] = True

		if self.options['re_split_reporter']:
			self.options['insert_oligo_5_flank_re'] = re.compile(self.options['sequencing_data_design']["oligo_format"]["insert_oligo_5'_flank"])
			self.options['insert_oligo_3_flank_re'] = re.compile(self.options['sequencing_data_design']["oligo_format"]["insert_oligo_3'_flank"])
			
		if not os.path.exists(data_json_path) or self.options['remake'] == True:	
			
			logging.info("Making: " + data_json_path)
			logging.info("Reading: " + fastq_gzipped_file)

			guides_list = list(self.guides_data.keys())
			with utilities_gzip.open(fastq_gzipped_file,'rt') as f:
				head = []
				line = ""
				oligo_counter = 0

				try:
					for line in f:
						try:
							if self.options["subsample"]:
								if oligo_counter >= self.options["subsample_cutoff"]:
									break
								
							head.append(line)
							
							###-----------------###

							if len(head) == 4:
								try:
									is_high_quality = True
									if self.options['pair_type'] == "reporter":
										if oligo_counter >= len(self.oligo_list['guide']):
											logging.error("oligo_counter >= len(self.oligo_list['guide'])")
											break
											
									if is_high_quality == False:
										self.read_metrics[self.options['pair_type']]["filtered_phred"] += 1
									else:
										oligo = head[1].strip()

										if self.options['pair_type'] == 'guide' and self.options['robust_guide_count']:
											insert_oligo = None
											oligo_counter += 1

											for i in range(0,len(oligo) - self.options["experiment_design"]["guide"]["insert_length"]):
												if oligo[i:i+self.options["experiment_design"]["guide"]["insert_length"]] in guides_list:
													insert_oligo = oligo[i:i+self.options["experiment_design"]["guide"]["insert_length"]]
													break
													
											self.oligo_list[self.options['pair_type']].append(insert_oligo)
											
											if oligo_counter < 10:
												print(self.options['pair_type'],"robust_guide_count",oligo,insert_oligo)

											if oligo_counter%100000 == 0:
												logging.info("Chunk " + str(self.read_metrics[self.options['pair_type']]["processed"] ) + "/" + str(len(self.oligo_list[self.options['pair_type']])) + " of " + self.options['analysis_name'] + ' ' + self.options['pair_type'] + " " + oligo)
												
										else:
											self.read_metrics[self.options['pair_type']]["total"] += 1
											skip = False
											
											if self.options['pair_type'] == "reporter":
												if len(self.oligo_list[self.options['pair_type']]) < 10:
													print(head)

												if self.oligo_list['guide'][oligo_counter] == None:
													self.read_metrics[self.options['pair_type']]["skipped_reporter"] += 1
													skip = True
											
											oligo_counter += 1
											
											if oligo_counter%100000 == 0:
												logging.info("Chunk " + str(len(self.oligo_list[self.options['pair_type']])) + " of " + self.options['analysis_name'] + ' ' + self.options['pair_type'] + " " + oligo)

											if skip == False:
												processed_oligos += 1
												
												insert_oligo = self.process_oligos(oligo)
												
												if oligo_counter< 10:
													print(self.options['pair_type'],oligo,insert_oligo)

												self.oligo_list[self.options['pair_type']].append(insert_oligo)
											else:
												self.oligo_list[self.options['pair_type']].append(None)
												if len(self.oligo_list[self.options['pair_type']]) < 10:
													print("Skipping")
								except:
									logging.error(utilities_error.getError())
									
								head = []
						except:
							logging.error(line)	
							logging.error(utilities_error.getError())
				except:
					logging.error(utilities_error.getError())

			#####------------------------#####
			
			logging.info("Writing " + data_json_path)
			pprint.pprint(self.read_metrics)
			
			if self.options['pair_type'] == "guide" or self.options['pair_type'] == "amplicon":
				logging.info("Unique oligos: " + str(len(set(self.oligo_list[self.options['pair_type']]))))
				logging.info("Non-unique oligos: " + str(len(self.oligo_list[self.options['pair_type']])))
				self.read_metrics[self.options['pair_type']]["unique_counts"] = len(set(self.oligo_list[self.options['pair_type']]))
				
			utilities_basic.write_to_json(data_json_path,self.oligo_list[self.options['pair_type']],zipped=True)
			logging.info("Finished Making " + data_json_path + " " + str(len(self.oligo_list[self.options['pair_type']])) + " reads")
		else:
			logging.info("Reading " + data_json_path)
			self.oligo_list[self.options['pair_type']] = utilities_basic.read_from_json(data_json_path,zipped=True)
			logging.info("Finished Reading " + data_json_path + " " + str(len(self.oligo_list[self.options['pair_type']])) + " reads")
		
	##################################################################################################

	def process_oligos(self,oligo):		
		insert_found = False

		if self.options['re_split_reporter']:
			insert_oligo_split = []
			if self.options['sequencing_data_design']["oligo_format"]["insert_oligo_5'_flank"] != "" and self.options['sequencing_data_design']["oligo_format"]["insert_oligo_3'_flank"] != "":
				insert_oligo_split_5 = re.split(self.options['insert_oligo_5_flank_re'],oligo)[-1]
				insert_oligo_split = re.split(self.options['insert_oligo_3_flank_re'],insert_oligo_split)[0]
				if insert_oligo_split_5 != oligo and insert_oligo_split != oligo:
					insert_found = True

			elif self.options['sequencing_data_design']["oligo_format"]["insert_oligo_5'_flank"] != "":
				insert_oligo_split = re.split(self.options['insert_oligo_5_flank_re'],oligo)[-1]
				if insert_oligo_split != oligo:
					insert_found = True

			elif self.options['sequencing_data_design']["oligo_format"]["insert_oligo_3'_flank"] != "":
				insert_oligo_split = re.split(self.options['insert_oligo_3_flank_re'],oligo)[0]
				if insert_oligo_split != oligo:
					insert_found = True
		else:
			insert_oligo_split = []

			if self.options['sequencing_data_design']["oligo_format"]["insert_oligo_5'_flank"] != "" and self.options['sequencing_data_design']["oligo_format"]["insert_oligo_3'_flank"] != "":
				insert_oligo_split_5 = oligo.split(self.options['sequencing_data_design']["oligo_format"]["insert_oligo_5'_flank"])[-1]
				insert_oligo_split = insert_oligo_split_5.split(self.options['sequencing_data_design']["oligo_format"]["insert_oligo_3'_flank"])[0]
			
				if insert_oligo_split_5 != oligo and insert_oligo_split != oligo:
					insert_found = True
	
			elif self.options['sequencing_data_design']["oligo_format"]["insert_oligo_5'_flank"] != "":
				insert_oligo_split = oligo.split(self.options['sequencing_data_design']["oligo_format"]["insert_oligo_5'_flank"])[-1]
				if insert_oligo_split != oligo:
					insert_found = True
					
			elif self.options['sequencing_data_design']["oligo_format"]["insert_oligo_3'_flank"] != "":
				insert_oligo_split = oligo.split(self.options['sequencing_data_design']["oligo_format"]["insert_oligo_5'_flank"])[0]	
				if insert_oligo_split != oligo:
					insert_found = True
		
		#-------#

		if not insert_found:
			self.read_metrics[self.options['pair_type']]["no_flanks"] += 1
			return
		else:
			insert_oligo = insert_oligo_split[:self.options['insert_length']]

		insert_oligo_5_flank_requirement = self.options['sequencing_data_design']["oligo_format"]["insert_oligo_5'_flank_requirement"]
		insert_oligo_3_flank_requirement = self.options['sequencing_data_design']["oligo_format"]["insert_oligo_3'_flank_requirement"]


		if len(insert_oligo) not in self.read_metrics[self.options['pair_type']]["incorrect_oligo_insert_length_distribution"]:
			self.read_metrics[self.options['pair_type']]["incorrect_oligo_insert_length_distribution"][len(insert_oligo)] = 0
		
		self.read_metrics[self.options['pair_type']]["incorrect_oligo_insert_length_distribution"][len(insert_oligo)] += 1
		
		if len(insert_oligo) == 0:
			return
		
		#-------#
		
		if insert_oligo_3_flank_requirement != "":
			insert_flanks = oligo.split(insert_oligo)
			if insert_flanks[-1][0:len(insert_oligo_3_flank_requirement)] != insert_oligo_3_flank_requirement:
				self.read_metrics[self.options['pair_type']]["no_flanks"] += 1
				return

		if len(insert_oligo) < self.options['insert_length']:
			self.read_metrics[self.options['pair_type']]["incorrect_oligo_insert_length"] += 1
			return
		
		#-------#
		
		if insert_oligo_5_flank_requirement != "":
			insert_flanks = oligo.split(insert_oligo)
			
		if 'use_complement' in self.options:
			if self.options['use_complement'] == True:	
				insert_oligo = utilities_codons.get_reverse_complementry_oligo(insert_oligo)

		barcode = ""
		if 'add_barcode' in self.options:
			if self.options['add_barcode'] == True:	
				barcode = insert_oligo[-self.options["barcode_length"]:]
				insert_oligo = insert_oligo[:-self.options["barcode_length"]]
		
		untrimmed_insert = insert_oligo
		if 'insert_trim' in self.options:
			if 'start' in self.options['insert_trim'] and 'end' in self.options['insert_trim']:
				insert_oligo = insert_oligo[self.options['insert_trim']['start']:self.options['insert_trim']['end']]
			elif 'start' in self.options['insert_trim']:
				insert_oligo = insert_oligo[self.options['insert_trim']['start']:]
			elif 'end' in self.options['insert_trim']:
				insert_oligo = insert_oligo[:self.options['insert_trim']['end']]
	
		#-------#

		self.read_metrics[self.options['pair_type']]["processed"] += 1
		
		if barcode != "":
			return [insert_oligo,barcode,untrimmed_insert]
		else:
			return insert_oligo
	
	##################################################################################################

	def calculate_reporter_edit_frequency(self):
		for guide in self.guides_data:
			self.guides_data[guide]["reporter_edit_frequency"] = {}
			unedited_reporter = guide
			
			if unedited_reporter in self.guides_data[guide]["ngs_data"]:
				self.guides_data[guide]["genomic_mutant_frequency"] = 1 - self.guides_data[guide]["ngs_data"][unedited_reporter]/self.guides_data[guide]['ngs_total_count_filtered']
			else:
				self.guides_data[guide]["genomic_mutant_frequency"] = 0

			###------------------###

			self.guides_data[guide]["reporter_edit_frequency"] = {}
			for reporter in self.guides_data[guide]["ngs_data"]:
				self.guides_data[guide]["reporter_edit_frequency"][reporter] = self.guides_data[guide]["ngs_data"][reporter]/self.guides_data[guide]["ngs_total_count_filtered"]
	
	##################################################################################################
	
	def map_designed_guide_edits(self):
		protein_gRNAs = {}
		for guide in self.guides_data:
			if self.guides_data[guide]['accession'] not in protein_gRNAs:
				protein_gRNAs[self.guides_data[guide]['accession']] = []
			
			protein_gRNAs[self.guides_data[guide]['accession']].append(guide)

		guideAnnotationObj = guideAnnotation.guideAnnotation()

		guideAnnotationObj.options['species'] = "9606"
		guideAnnotationObj.options['editors'] = ["ABE","CBE"]
		guideAnnotationObj.options['behive_base_editor'] = "ABE"
		guideAnnotationObj.options['protein_gRNAs'] = protein_gRNAs

		designed_guide_edits = guideAnnotationObj.annotate_protein_guides()

		for guide in designed_guide_edits:
			self.guides_data[guide]['designed_gRNAs_expected_edits'] = designed_guide_edits[guide]
			
		del guideAnnotationObj

	##################################################################################################

	def map_reporter_edits(self):
		protein_gRNAs_reporter = {}
		protein_gRNAs_reporter_low_frequency = {}

		for guide in self.guides_data:
			
			unedited_reporter = guide
			total_reporter_count = sum(self.guides_data[guide]["ngs_data"].values())
			for reporter in self.guides_data[guide]["ngs_data"]:
				trimmed_reporter = reporter
				
				accession = self.guides_data[guide]['accession']
				if accession not in protein_gRNAs_reporter:
					protein_gRNAs_reporter[accession] = {}
					protein_gRNAs_reporter_low_frequency[accession] = {}

				if unedited_reporter not in protein_gRNAs_reporter[accession]:
					protein_gRNAs_reporter[accession][unedited_reporter] = {}
					protein_gRNAs_reporter_low_frequency[accession][unedited_reporter] = {}

				if unedited_reporter == trimmed_reporter:
					protein_gRNAs_reporter[accession][unedited_reporter][trimmed_reporter] = self.guides_data[guide]["ngs_data"][trimmed_reporter]/total_reporter_count
				elif self.guides_data[guide]["ngs_data"][reporter]/total_reporter_count <= self.options['reporter_count_cutoff']: 
					protein_gRNAs_reporter_low_frequency[accession][unedited_reporter][trimmed_reporter] = self.guides_data[guide]["ngs_data"][trimmed_reporter]/total_reporter_count
				else:
					protein_gRNAs_reporter[accession][unedited_reporter][trimmed_reporter] = self.guides_data[guide]["ngs_data"][trimmed_reporter]/total_reporter_count
	
		###---------------###

		accession_counter = 0
		for accession in protein_gRNAs_reporter:
			try:
				accession_counter += 1
				
				if accession == None: continue
				logging.info(">"*5 + str(accession) + " " + str(accession_counter)  + " of " + str(len(protein_gRNAs_reporter)))
				guideAnnotationObj = guideAnnotation.guideAnnotation()
				guideAnnotationObj.options['species'] = "9606"
				guideAnnotationObj.options['accession'] = accession
				guideAnnotationObj.options['reporter_edited_gRNAs'] = protein_gRNAs_reporter[accession]
				reporter_edits_data = guideAnnotationObj.annotate_protein_reporter_edited_gRNAs()

				for guide in self.guides_data:
					if self.guides_data[guide]['accession'] != accession: continue

					unedited_reporter = guide

					if unedited_reporter in reporter_edits_data:
						self.guides_data[guide]['reporter_edits'] = reporter_edits_data[unedited_reporter]

					if 'reporter_edits' not in self.guides_data[guide]:
						continue
				
					compare_report_list = copy.deepcopy(list(protein_gRNAs_reporter[accession][guide].keys()))
					
					reporter_edits_update = {}
					
					if self.options['map_unedited_reptorter_to_similar_reporter']:
						for reporter in protein_gRNAs_reporter_low_frequency[accession][unedited_reporter]:
							if reporter in self.guides_data[guide]['reporter_edits']: continue
							
							most_similar_reporter = self.find_most_similar(unedited_reporter,reporter, compare_report_list)
						
							min_dissimilarity = min(list(most_similar_reporter.values()))
							sorted_reporters_scores = dict(sorted(protein_gRNAs_reporter[accession][guide].items(), key=lambda item: item[1],reverse=True))
							sorted_reporters = list(sorted_reporters_scores.keys())

							mapped_reporter = None
							for sorted_reporter in sorted_reporters:
								if sorted_reporter in most_similar_reporter:
									if mapped_reporter == None and most_similar_reporter[sorted_reporter] == min_dissimilarity:
										mapped_reporter = sorted_reporter
										break
							
							if mapped_reporter != None:
								reporter_edits_update[reporter] = {} 
								reporter_edits_update[reporter]['mutated_genomic_region'] = reporter
								reporter_edits_update[reporter]['collapsed_from'] = mapped_reporter

								for tag in ['changes', 'difference','edits','genome_edited', 'proteome_edited']:
									if mapped_reporter in self.guides_data[guide]['reporter_edits']:
										if tag in self.guides_data[guide]['reporter_edits'][mapped_reporter]:
											reporter_edits_update[reporter][tag] = copy.deepcopy(self.guides_data[guide]['reporter_edits'][mapped_reporter][tag])
										else:
											print("Error: " + tag + " not in reporter_edits: " + mapped_reporter)
									else:
										pass
										print("Error: " + mapped_reporter + " not in reporter_edits")
										print(mapped_reporter == unedited_reporter)
							else:
								print(
									"ERROR",
									accession,"\t",
									min_dissimilarity,"\t",
									guide,
									reporter,"\t",
									mapped_reporter,"\t",
									unedited_reporter == mapped_reporter
								)

					self.guides_data[guide]['reporter_edits'].update(reporter_edits_update)
					
				del guideAnnotationObj
			except:
				logging.error(str(accession) + " not processed by annotate_protein_reporter_edited_gRNAs")
				logging.error(utilities_error.getError())
				raise
			
		del protein_gRNAs_reporter
						

	##################################################################################################

	def hamming_distance(self, unedited_reporter: str, seq1: str, seq2: str) -> int:
		"""Calculate the Hamming distance between two equal-length strings."""

		if unedited_reporter == seq2:
			distance = sum(seq1[i] != seq2[i] for i in range(0,len(unedited_reporter)))
		else:
			distance = sum(seq1[i] != seq2[i] and unedited_reporter[i] != seq2[i] for i in range(0,len(unedited_reporter)))
			
		return distance

	def find_most_similar(self,unedited_reporter: str, target: str, sequences: list) -> str:
		"""Find the most similar 10-mer string from a list based on Hamming distance."""
		if not sequences:
			return None
		
		hamming_distances = {}
		for seq in sequences:
			hamming_distances[seq] = self.hamming_distance(unedited_reporter,target, seq)
		return hamming_distances 

	##################################################################################################

	def calculate_reporter_splice_site_edit_frequency(self):
		for guide in self.guides_data:
			try:
				self.guides_data[guide]["splice_site_edit"] = {}
				
				self.guides_data[guide]["splice_site_mutant_frequency"] = 0
				if 'reporter_edits' not in self.guides_data[guide]:
					continue

				for reporter in self.guides_data[guide]['reporter_edits']:
					if reporter not in self.guides_data[guide]["reporter_edit_frequency"]:
						logging.error(reporter + " not in reporter_edit_frequency data")
						continue

					splice_site_edits = ""

					if 'splice_site' in self.guides_data[guide]['reporter_edits'][reporter]:
						if self.guides_data[guide]['reporter_edits'][reporter]['splice_site']:
							splice_site_edits = self.guides_data[guide]['reporter_edits'][reporter]['splice_site_info']['splice_site_type']
							self.guides_data[guide]["splice_site_mutant_frequency"] += self.guides_data[guide]["reporter_edit_frequency"][reporter]
						
						
							if splice_site_edits not in self.guides_data[guide]["splice_site_edit"]:
								self.guides_data[guide]["splice_site_edit"][splice_site_edits] = 0

							self.guides_data[guide]["splice_site_edit"][splice_site_edits] += self.guides_data[guide]["reporter_edit_frequency"][reporter]
			except:
				logging.error(utilities_error.getError())

	##################################################################################################
	
	def calculate_reporter_protein_edit_frequency(self):
		for guide in self.guides_data:
			reporter_edit_frequency_sum = 0
			try:
				self.guides_data[guide]["protein_edit"] = {
					"single_edit":{},
					"multiple_edit":{}
				}
				
				if self.guides_data[guide]['target_type'] == "targets": 
					self.guides_data[guide]["protein_mutant_frequency"] = 0
					if 'reporter_edits' not in self.guides_data[guide]:
						logging.error("'reporter_edits' not in guide " + guide)
						continue
					
					for reporter in self.guides_data[guide]['reporter_edits']:
					
						if reporter not in self.guides_data[guide]["reporter_edit_frequency"]:
							logging.error(reporter + " not in reporter_edit_frequency data")
							continue
						
						if self.options['reporter_skip_no_edit_annotation']:
							if len(self.guides_data[guide]['reporter_edits'][reporter]['changes']) == 0:
								continue

						reporter_edit_frequency_sum += self.guides_data[guide]["reporter_edit_frequency"][reporter]
						
						if 'changes' not in self.guides_data[guide]['reporter_edits'][reporter]:
							changes = [""]
						else:
							if len(self.guides_data[guide]['reporter_edits'][reporter]['changes']) > 0:
								self.guides_data[guide]["protein_mutant_frequency"] += self.guides_data[guide]["reporter_edit_frequency"][reporter]
						
							changes = self.guides_data[guide]['reporter_edits'][reporter]['changes']
						
							if len(changes) == 0: changes = [""]

						for mutation in changes:
							if mutation not in self.guides_data[guide]["protein_edit"]["single_edit"]:
								self.guides_data[guide]["protein_edit"]["single_edit"][mutation] = 0

							self.guides_data[guide]["protein_edit"]["single_edit"][mutation] += self.guides_data[guide]["reporter_edit_frequency"][reporter]
						
						mutation = ",".join(changes)
						if mutation not in self.guides_data[guide]["protein_edit"]["multiple_edit"]:
							self.guides_data[guide]["protein_edit"]["multiple_edit"][mutation] = 0

						self.guides_data[guide]["protein_edit"]["multiple_edit"][mutation] += self.guides_data[guide]["reporter_edit_frequency"][reporter]
				else:
					self.guides_data[guide]["protein_mutant_frequency"] = 0

				###------------------###
					
				self.guides_data[guide]["protein_edit"]["single_edit"] = dict(sorted(self.guides_data[guide]["protein_edit"]["single_edit"].items(), key=lambda item: item[1],reverse=True))
				self.guides_data[guide]["protein_edit"]["multiple_edit"] = dict(sorted(self.guides_data[guide]["protein_edit"]["multiple_edit"].items(), key=lambda item: item[1],reverse=True))
			
				single_edit_best = list(self.guides_data[guide]["protein_edit"]["single_edit"].keys())[0] if len(self.guides_data[guide]["protein_edit"]["multiple_edit"]) > 0 else "-"
				self.guides_data[guide]["most_common_protein_single_edit"] = single_edit_best + ":" + "%1.2f"%self.guides_data[guide]["protein_edit"]["single_edit"][single_edit_best] if len(self.guides_data[guide]["protein_edit"]["single_edit"]) > 0 else "-"
				multiple_edit_best = list(self.guides_data[guide]["protein_edit"]["multiple_edit"].keys())[0] if len(self.guides_data[guide]["protein_edit"]["multiple_edit"]) > 0 else "-"
				self.guides_data[guide]["most_common_protein_multiple_edit"] = multiple_edit_best + ":" + "%1.2f"%self.guides_data[guide]["protein_edit"]["multiple_edit"][multiple_edit_best] if len(self.guides_data[guide]["protein_edit"]["multiple_edit"]) > 0 else "-"
			except:
				logging.error(utilities_error.getError())

	##################################################################################################

	def process_reporter_screen_output(self):
		reporter_rows = []
		for guide in self.guides_data:
			for reporter in self.guides_data[guide]['ngs_data']:
				reporter_row = []

				try:
					reporter_row += [
						str(self.guides_data[guide]['accession']),
						str(self.guides_data[guide]['annotation']),
						guide,
						reporter,
						"".join([guide[i] if guide[i] != reporter[i] else "-" for i in  range(0,len(guide))]),
						"".join([reporter[i] if guide[i] != reporter[i] else "-" for i in  range(0,len(guide))]),
						str(guide!=reporter),
						str(self.guides_data[guide]['ngs_total_count']),
						str(self.guides_data[guide]['ngs_total_count_filtered']),
						"%1.3f"%self.guides_data[guide]["genomic_mutant_frequency"],
						"%1.3f"%self.guides_data[guide]["protein_mutant_frequency"],
						str(self.guides_data[guide]['ngs_data'][reporter]),
						"%1.3f"%(self.guides_data[guide]['ngs_data'][reporter]/self.guides_data[guide]['ngs_total_count_filtered']),
					]

					if 'reporter_edits' in self.guides_data[guide] and reporter in self.guides_data[guide]['reporter_edits']:
						try:
							reporter_row += [str(len(self.guides_data[guide]['reporter_edits'][reporter]['changes'])) if 'changes' in self.guides_data[guide]['reporter_edits'][reporter] else "-"]
						except:
							reporter_row += ["-"]
						
						try:
							reporter_row += [str(self.guides_data[guide]['reporter_edits'][reporter]['edit_counts']) if 'edit_counts' in self.guides_data[guide]['reporter_edits'][reporter] else "-"]
						except:
							reporter_row += ["-"]
						
						try:
							reporter_row += [str(self.guides_data[guide]['reporter_edits'][reporter]["wildtype_region"]) if 'wildtype_region' in self.guides_data[guide]['reporter_edits'][reporter] else "-"]
						except:
							reporter_row += ["-"]

						try:
							reporter_row += [str(self.guides_data[guide]['reporter_edits'][reporter]["mutated_region"]) if 'mutated_region' in self.guides_data[guide]['reporter_edits'][reporter] else "-"]
						except:
							reporter_row += ["-"]
						
						try:
							reporter_row += [",".join(self.guides_data[guide]['reporter_edits'][reporter]['changes']) if 'changes' in self.guides_data[guide]['reporter_edits'][reporter] else "-"]
						except:
							reporter_row += ["-"]
						
						try:
							reporter_row += [self.guides_data[guide]['target_type']]
						except:
							reporter_row += ["-"]
				except:
					reporter_row += [
						"failed"
					]

				reporter_rows.append("\t".join(reporter_row)) 

		header = [
			"accession",
			"annotation",
			"guide",
			"reporter",
			"guide_edited_nucleotides",
			"reporter_edited_nucleotides",
			"reporter_edited",
			"total_reporter_ngs_count",
			"total_reporter_ngs_count_filtered",
			"total_reporter_genomic_mutant_frequency",
			"total_report_protein_mutant_frequency",
			"reporter_ngs_count",
			"reporter_edit_frequency",
			"reporter_protein_changes_count",
			"reporter_genomic_changes_count",
			"reporter_wildtype_region",
			"reporter_mutated_region",
			"reporter_protein_changes",
			"target_type"
		]

		if self.options["subsample"]:
			data_tdt_path = os.path.join(self.options["data_directory"],self.options['analysis_name'] + '.' + str(self.options["subsample_cutoff"]) + '.detailed.reporter.tdt')
		else:
			data_tdt_path = os.path.join(self.options["data_directory"],self.options['analysis_name'] + '.detailed.reporter.tdt')

		logging.info("Writing " + data_tdt_path)
		open(data_tdt_path,"w").write("\t".join(header) + "\n" + "\n".join(reporter_rows))

	##################################################################################################

	def process_reporter_counts(self):	
		for guide in self.guides_data:
			try:
				unedited_reporter = guide
				
				if self.read_metrics["global_guide_edit_percentage"] == None or len(self.read_metrics["global_guide_edit_percentage"]) == 0:
					self.read_metrics["global_guide_edit_percentage"] = [0]*len(unedited_reporter)

				guide_count = [0]*len(unedited_reporter)
				
				for reporter in self.guides_data[guide]['ngs_data']:
					for i in range(0,len(unedited_reporter)):
						if reporter[i] == unedited_reporter[i]:
							guide_count[i] += self.guides_data[guide]['ngs_data'][reporter]/self.guides_data[guide]['ngs_total_count_filtered']
							self.read_metrics["global_guide_edit_percentage"][i] +=(self.guides_data[guide]['ngs_data'][reporter]/self.guides_data[guide]['ngs_total_count_filtered'])/len(self.guides_data)
				
				self.guides_data[guide]['guide_edit_percentage'] = ["%3.3f"%(abs(1-x)*100) for x in guide_count]
			except:
				logging.error(utilities_error.getError())
			
	##################################################################################################

	def process_guide_screen_output(self):		
		
		rows = []
		for guide in self.guides_data:
			try:
				unedited_reporter = guide

				if unedited_reporter in self.guides_data[guide]['ngs_data']:
					###------------------###

					single_edit_proportion = [sorted_reporter+ ":" + "%1.3f"%self.guides_data[guide]["protein_edit"]["single_edit"] [sorted_reporter] for sorted_reporter in self.guides_data[guide]["protein_edit"]["single_edit"]  if self.guides_data[guide]["protein_edit"]["single_edit"] [sorted_reporter] >= self.options['reporter_count_cutoff']]
					multiple_edit_proportion = [sorted_reporter + ":" + "%1.3f"%self.guides_data[guide]["protein_edit"]["multiple_edit"] [sorted_reporter] for sorted_reporter in self.guides_data[guide]["protein_edit"]["multiple_edit"]  if self.guides_data[guide]["protein_edit"]["multiple_edit"] [sorted_reporter] >= self.options['reporter_count_cutoff']]
					
					###------------------###
					
					splice_site_edit_proportion = [splice_site_type + ":" + "%1.3f"%self.guides_data[guide]["splice_site_edit"][splice_site_type] for splice_site_type in self.guides_data[guide]["splice_site_edit"]]
					
					###------------------###
					
					row = [
						str(self.options['analysis_name']),
						str(self.guides_data[guide]['accession']),
						str(self.guides_data[guide]['annotation']),
						guide,
						guide,
						str(self.guides_data[guide]['ngs_data'][unedited_reporter]),
						str(self.guides_data[guide]['ngs_total_count'] if 'ngs_total_count' in self.guides_data[guide] else "-"),
						"%1.3f"%self.guides_data[guide]["genomic_mutant_frequency"] if 'genomic_mutant_frequency' in self.guides_data[guide] else "-",
						"%1.3f"%self.guides_data[guide]["protein_mutant_frequency"] if 'protein_mutant_frequency' in self.guides_data[guide] else "-",
						single_edit_proportion[0] if len(single_edit_proportion) > 0 else "-",
						multiple_edit_proportion[0] if len(multiple_edit_proportion) > 0 else "-",
						"|".join(single_edit_proportion),
						"|".join(multiple_edit_proportion),
						"|".join(splice_site_edit_proportion)
					]

					row += self.guides_data[guide]['guide_edit_percentage']
					row += [self.guides_data[guide]['target_type']]
					
					rows.append("\t".join(row))
			except:
				row = [
					str(self.guides_data[guide]['accession']),
					str(self.guides_data[guide]['annotation']),
					guide,
					str(utilities_error.getError())
				]
				
				rows.append("\t".join(row))

		###------------------###

		header = [
			"sample",
			"accession",
			"annotation",
			"guide",
			"reporter",
			"wt NGS count",
			"total NGS count",
			"genomic_mutant_frequency",
			"protein_mutant_frequency",
			"most_common_protein_single_edit",
			"most_common_protein_multiple_edit",
			"protein_single_edits",
			"protein_multiple_edits",
			"splice_site_edits"
		] + [
			"gRNA p" + str(x+1) + " mut %" for x in range(0,len(list(self.guides_data[guide]['guide_edit_percentage']  if 'guide_edit_percentage' in self.guides_data[guide] else [])))
		] + [
			"target_type",
		]

		###------------------###

		if self.options["subsample"]:
			data_tdt_path = os.path.join(self.options["data_directory"],self.options['analysis_name'] + '.'  + str(self.options["subsample_cutoff"]) + '.edits.reporter.tdt')
		else:
			data_tdt_path = os.path.join(self.options["data_directory"],self.options['analysis_name'] + '.edits.reporter.tdt')
		
		logging.debug("Writing " + data_tdt_path)
		open(data_tdt_path,"w").write("\t".join(header) + "\n" + "\n".join(rows))

		###------------------###

		if self.options["subsample"]:
			data_tdt_path = os.path.join(self.options["data_directory"],self.options['analysis_name'] + '.'  + str(self.options["subsample_cutoff"]) + '.edits.reporter.global.tdt')
		else:
			data_tdt_path = os.path.join(self.options["data_directory"],self.options['analysis_name'] + '.edits.reporter.global.tdt')
		
		try:
			global_tdt_content = "\t".join(["sample","unfiltered reads"] + ["gRNA p" + str(x+1) + " mut %" for x in range(0,len(self.read_metrics["global_guide_edit_percentage"]))]) + "\n" + "\t".join([self.options['analysis_name'],str(self.read_metrics["unfiltered"])] + ["%3.3f"%(abs(1-x)*100) for x in self.read_metrics["global_guide_edit_percentage"]])
		
			logging.debug("Writing " + data_tdt_path)
			open(data_tdt_path,"w").write(global_tdt_content)
		except:
			logging.error("Failed writing " + data_tdt_path)
			
		###------------------###
		if self.options["subsample"]:
			data_stats_json_path = os.path.join(self.options["data_directory"],self.options['analysis_name'] + '.'  + str(self.options["subsample_cutoff"]) + '.stats.reporter.json')
		
		else:
			data_stats_json_path = os.path.join(self.options["data_directory"],self.options['analysis_name'] +  '.stats.reporter.json')
		
		data_stats_json_path = os.path.join(self.options["data_directory"],self.options['analysis_name'] +  '.stats.reporter.json')
		logging.debug("Writing " + data_stats_json_path)
		utilities_basic.write_to_json(data_stats_json_path,self.read_metrics,zipped=False)

	##################################################################################################

	def plot_reporter_distribution(self):
		import pandas as pd
		import seaborn as sns
		import matplotlib.pyplot as plt

		plt.set_loglevel (level = 'warning')

		plot_data = []
		for guide in self.guides_data:
			try:
				unedited_reporter = guide

				guide_start = len(unedited_reporter.split(guide)[0])
				for x in range(0,len(unedited_reporter)):
					plot_data.append({
						"guide":guide,
						"position":"p" + str(x + 1 - guide_start),
						"guide_edit_percentage":float(self.guides_data[guide]['guide_edit_percentage'][x]),
						"nucleotide":unedited_reporter[x]
						}
					)
			except:
				logging.error(guide + "  failed")
		
		df = pd.DataFrame(plot_data)
		
		df.set_index('guide', inplace=True)

		plt.figure(figsize=(40,6))
		ax = sns.boxplot(
			data=df,
			y='guide_edit_percentage',
			x='position',
			hue='nucleotide',
			hue_order=['A','C','G','T'],
			fliersize=0.5
		)
		ax.set_ylim([0, 100]) 
		ax.set_xlabel("position")
		ax.set_ylabel("guide_edit_percentage")
		
		###------------------###

		if self.options["subsample"]:
			data_boxplot_path = os.path.join(self.options["data_directory"],self.options['analysis_name']  + '.'  + str(self.options["subsample_cutoff"]) + '.reporter.boxplot.png')
		else:
			data_boxplot_path = os.path.join(self.options["data_directory"],self.options['analysis_name'] + '.reporter.boxplot.png')
		
		plt.savefig(data_boxplot_path)

	##################################################################################################

	def process_guide_annotation_file(self):
		logging.info("Processing annotation file: " + self.options["guide_annotation_file"])

		self.protein_reporter_edits = {}
		self.accession_mapping = {}

		file_content = open(self.options["guide_annotation_file"]).read().strip().split("\n")
		
		file_content_header = file_content[0]
		file_content_data = file_content[1:]

		file_content_header = file_content_header.split("\t")
		
		guide_index = file_content_header.index("guide")
		barcode_index = file_content_header.index("bar_code")
		reporter_index =  file_content_header.index("reporter")
		annotation_index = file_content_header.index("annotation")
		target_type_index = file_content_header.index("type")

		line_counter = 0

		accession_check_dict = {}
		for line in file_content_data:
			line_counter += 1

			if len(line.split("\t")) < 2: continue

			identifier = line.split("\t")[annotation_index].split("|")[0]
			
			if line.split("\t")[target_type_index] == "target":
				line.split("\t")[target_type_index] = "targets"

			if line.split("\t")[target_type_index] in ["targets",'target']:
				if identifier not in self.protein_reporter_edits: 
					self.protein_reporter_edits[identifier] = {}
				self.protein_reporter_edits[identifier][line.split("\t")[guide_index]] = {}

			accession = None
			if identifier not in accession_check_dict:
				logging.info("Adding" + identifier)
				accession_check_dict[identifier] = queryRunner.queryRunner("uniprot","check_accession",{"accession":identifier}).run()
			
			if accession_check_dict[identifier] ['status'] == "Error":
				if identifier not in self.accession_mapping:
					gene_name_check = queryRunner.queryRunner("uniprot","check_gene_name",{"identifier":identifier}).run()
					
					if gene_name_check['data']:
						gene_name_accession = queryRunner.queryRunner("uniprot","get_uniprot_mapping",{"identifier":[identifier],"mapping_from":"GENENAME","mapping_to":"UniProtKB"}).run()
						
						try:
							self.accession_mapping[identifier] = gene_name_accession['data'][identifier]['Entry']
							accession = self.accession_mapping[identifier]
						except:
							self.accession_mapping[identifier] = None
				else:
					accession = self.accession_mapping[identifier]
			else:
				if accession_check_dict[identifier]['data']:
					accession = identifier

			self.guides_data[line.split("\t")[guide_index]] = {
				"barcode":line.split("\t")[barcode_index],
				"reporter":line.split("\t")[reporter_index],
				"target_type":line.split("\t")[target_type_index],
				"accession":accession,
				"annotation":line.split("\t")[annotation_index],
				"ngs_total_count":0,
				"ngs_total_count_filtered":0,
				"ngs_total_count_edited":0,
				"ngs_data":{}
			}

		logging.info("guide_annotation_file: " + str(len(self.guides_data)))

		###------------------###
	
	def process_rhAmpSeq_ngs_files(self):		
		self.setup_rhAmpSeq_analysis()

		for guide in self.options["experiment_design"]["amplicon"]:
			print("#"*100)
			print("#",guide)
			print( self.options["experiment_design"]["amplicon"][guide]['strand'])

			if guide == "": continue

			self.guides_data[guide] = {
				"barcode":"",
				"reporter":self.options["experiment_design"]['amplicon'][guide]['reporter_sequence'] ,
				"target_type":"targets",
				"accession":self.options["experiment_design"]['amplicon'][guide]['accession'],
				"annotation":self.options["experiment_design"]['amplicon'][guide]['gene'],
				"ngs_total_count":0,
				"ngs_total_count_filtered":0,
				"ngs_total_count_edited":0,
				"ngs_data":{}
			}

			self.options['sequencing_data_design'] = self.options["experiment_design"]["amplicon"][guide]
			self.options['pair_type'] = "amplicon"
			
			self.options["subset_tag"] = "." + guide  + "." + str(self.options["amplicon_flanks_size"])
			self.oligo_list = {}
			
			for tag in ['insert_length','sample_information','use_complement','insert_trim','add_barcode','barcode_length']:
				if tag in self.options["experiment_design"]["amplicon"][guide]:
					self.options[tag] = self.options["experiment_design"]["amplicon"][guide][tag]

			self.read_fastq_gzipped_files()

			edited = {True:0,False:0}
			for i in range(0,len(self.oligo_list['amplicon'])):
				reporter = self.oligo_list['amplicon'][i]

				self.guides_data[guide]["ngs_total_count"] += 1

				if self.oligo_list['amplicon'][i] == None: 
					continue

				if reporter not in self.guides_data[guide]["ngs_data"]:
					self.guides_data[guide]["ngs_data"][reporter] = 0

				self.guides_data[guide]["ngs_data"][reporter] += 1
				self.guides_data[guide]["ngs_total_count_filtered"] += 1

				if reporter == self.guides_data[guide]['reporter']:
					self.guides_data[guide]["ngs_total_count_edited"] += 1

		logging.info("Running map_reporter_edits")
		self.map_reporter_edits()
		logging.info("Completed map_reporter_edits")


		logging.info("Running calculate_reporter_edit_frequency")
		self.calculate_reporter_edit_frequency()
		logging.info("Running process_reporter_counts")
		self.process_reporter_counts()

		###------------------###

		logging.info("Running calculate_reporter_protein_edit_frequency")
		self.calculate_reporter_protein_edit_frequency()
		logging.info("Running calculate_reporter_splice_site_edit_frequency")
		self.calculate_reporter_splice_site_edit_frequency()

		###------------------###

		logging.info("Running process_guide_screen_output")
		self.process_guide_screen_output()
		logging.info("Running process_reporter_screen_output")
		self.process_reporter_screen_output()

		###------------------###

		logging.info("Running plot_reporter_distribution")
		self.plot_reporter_distribution()
		logging.info("Running make_logo")
		self.make_logo()
		

	def process_ngs_file(self):		
		if self.options["subsample"]:
			ngs_data_json_path = os.path.join(self.options['experiment_design']['guide']['pair_sequencing_file'].replace("_R1_001.fastq.gz","") + '.' + str(self.options["subsample_cutoff"]) + self.options["subset_tag"] + ".ngs_data.reporter.jsonc")
		else:
			ngs_data_json_path = os.path.join(self.options['experiment_design']['guide']['pair_sequencing_file'].replace("_R1_001.fastq.gz","") + self.options["subset_tag"] + '.ngs_data.reporter.jsonc')
			
		if self.options["remake"] or not os.path.exists(ngs_data_json_path):
			logging.info("Making: " + ngs_data_json_path)

			self.process_guide_annotation_file()

			self.oligo_list = {}
			for pair_type in ["guide","reporter"]:#self.options['experiment_design']:
				logging.info("Processing " + pair_type + " from " + self.options['experiment_design'][pair_type]['pair_sequencing_file'])
				self.options['fastq_zipped_file'] = self.options['experiment_design'][pair_type]['pair_sequencing_file']
				self.options['sequencing_data_design'] = self.options["experiment_design"][pair_type]
				self.options['pair_type'] = pair_type

				for tag in ['insert_length','sample_information','use_complement','insert_trim','add_barcode','barcode_length']:
					if tag in self.options["experiment_design"][pair_type]:
						self.options[tag] = self.options["experiment_design"][pair_type][tag]

				self.read_fastq_gzipped_files()

			logging.info("read- reporter:" + str(len(self.oligo_list['reporter'])) + " guide:" + str(len(self.oligo_list['guide'])))
			logging.info("processed - reporter:" + str(self.read_metrics["reporter"]["processed"]) + " guide:" + str(self.read_metrics["guide"]["processed"]))
			
			if len(self.oligo_list['reporter']) != len(self.oligo_list['guide']):
				logging.error("NGS length mismatch")
				logging.error(self.options['experiment_design']["guide"]['pair_sequencing_file'])
				logging.error(self.options['experiment_design']["reporter"]['pair_sequencing_file'])
				sys.exit()
				
			###------------------###
				
			for i in range(0,len(self.oligo_list['guide'])):
				try:
					guide = self.oligo_list['guide'][i]
					[reporter,barcode,full_reporter] = self.oligo_list['reporter'][i] if self.oligo_list['reporter'][i] != None else [None,None,None]

					if guide == None:
						self.read_metrics["no_guide"] += 1 
					if reporter == None: 
						self.read_metrics["no_reporter"] += 1 
					if guide not in self.guides_data: 
						self.read_metrics["not_designed_guide"] += 1 
					
					if guide == None:
						continue

					if guide not in self.guides_data: 
						continue
					
					self.guides_data[guide]["ngs_total_count"] += 1

					if reporter == None: 
						continue

					if barcode != self.guides_data[guide]['barcode']: 
						self.read_metrics["no_barcode"] += 1 
						continue

					if barcode != self.guides_data[guide]['barcode']: 
						self.read_metrics["incorrect_barcode"] += 1 
						continue

					self.read_metrics["unfiltered"] += 1

					if reporter not in self.guides_data[guide]["ngs_data"]:
						self.guides_data[guide]["ngs_data"][reporter] = 0

					self.guides_data[guide]["ngs_data"][reporter] += 1
					self.guides_data[guide]["ngs_total_count_filtered"] += 1

					if reporter == self.guides_data[guide]['reporter']:
						self.guides_data[guide]["ngs_total_count_edited"] += 1
				except:
					logging.error(i)
					logging.error(self.oligo_list['guide'][i])
					logging.error(self.oligo_list['reporter'][i])
			
			counter = 0
			for guide in self.guides_data:
				counter += 1
				
			###---------------------###
			
			utilities_basic.write_to_json(ngs_data_json_path,self.guides_data,zipped=True)

		else:
			logging.debug("Reading: " + ngs_data_json_path)
			self.guides_data  = utilities_basic.read_from_json(ngs_data_json_path,zipped=True)
		
		self.oligo_list = {}
		return self.guides_data

	def process_screen(self,return_data=True):		
		self.setup_analysis()
		
		if self.options["subsample"]:
			data_json_path = os.path.join(self.options['experiment_design']['guide']['pair_sequencing_file'].replace(".fastq.gz","") + "." + str(self.options["subsample_cutoff"]) + "." + str(self.options["reporter_count_cutoff"]) + '.data.reporter.jsonc')
			metrics_json_path = os.path.join(self.options['experiment_design']['guide']['pair_sequencing_file'].replace(".fastq.gz","") + "." + str(self.options["subsample_cutoff"]) + "." + str(self.options["reporter_count_cutoff"])+ '.metrics.reporter.json')
		else:
			data_json_path = os.path.join(self.options['experiment_design']['guide']['pair_sequencing_file'].replace(".fastq.gz","") + "." + str(self.options["reporter_count_cutoff"]) + '.data.reporter.jsonc')
			metrics_json_path = os.path.join(self.options['experiment_design']['guide']['pair_sequencing_file'].replace(".fastq.gz","") + "." + str(self.options["reporter_count_cutoff"]) + '.metrics.reporter.json')
		
		logging.info("Processing " + data_json_path)
		if not os.path.exists(data_json_path) or self.options['remake']:
			logging.info("Making " + data_json_path)
		
			###------------------###

			logging.info("Running process_ngs_file")
			self.process_ngs_file()
			logging.info("Completed process_ngs_file")
			
			###------------------###
			
			logging.info("Running map_reporter_edits")
			self.map_reporter_edits()
			logging.info("Completed map_reporter_edits")

			###------------------###

			logging.debug("Writing " + data_json_path)
			utilities_basic.write_to_json(data_json_path,self.guides_data,zipped=True)

			logging.debug("Writing " + metrics_json_path)
			utilities_basic.write_to_json(metrics_json_path,self.read_metrics,zipped=False)
		else:
			if return_data:
				logging.debug("Reading " + data_json_path)
				self.guides_data = utilities_basic.read_from_json(data_json_path,zipped=True)
			
			self.read_metrics = utilities_basic.read_from_json(metrics_json_path,zipped=False)

		###------------------###
		
		if self.options['calculate_reporter_metrics']:
			logging.info("Running calculate_reporter_edit_frequency")
			self.calculate_reporter_edit_frequency()
			logging.info("Running process_reporter_counts")
			self.process_reporter_counts()

			###------------------###

			logging.info("Running calculate_reporter_protein_edit_frequency")
			self.calculate_reporter_protein_edit_frequency()
			logging.info("Running calculate_reporter_splice_site_edit_frequency")
			self.calculate_reporter_splice_site_edit_frequency()

			###------------------###

			logging.info("Running process_guide_screen_output")
			self.process_guide_screen_output()
			logging.info("Running process_reporter_screen_output")
			self.process_reporter_screen_output()

			###------------------###

			logging.info("Running plot_reporter_distribution")
			self.plot_reporter_distribution()
			logging.info("Running make_logo")
			self.make_logo()
		
		if self.options["guides_data_json_file"] != None:
			utilities_basic.write_to_json(self.options["guides_data_json_file"],self.guides_dat,zipped=False)
		return self.guides_data 
	
	##################################################################################################

	def make_position_logo(self):	
		try:	
			logo_data = {
				'reporter':{},
				'edited_reporter':{},
				'guide':{}
			}

			print("#######",self.options['analysis_name'])
			for guide in self.guides_data.keys():
				for i in range(3,9):
					if self.guides_data[guide]["reporter"][i] not in ['A','C']: continue
					flanks = self.guides_data[guide]["reporter"][i] + "-" + self.guides_data[guide]["reporter"][i-1] + "-" + self.guides_data[guide]["reporter"][i] + "-" + self.guides_data[guide]["reporter"][i+1] + "-" 
					flanks = flanks + self.guides_data[guide]["reporter"][i-1] + self.guides_data[guide]["reporter"][i] + self.guides_data[guide]["reporter"][i+1] + "-" 
					flanks = flanks + str(i)

					reporter_edit_counts = {'A':[], 'C':[], 'T':[], 'G':[], 'N':[]}
					reporter_edit_counts_total = 0
					for reporter in self.guides_data[guide]["ngs_data"]:
						reporter_edit_counts_total += self.guides_data[guide]["ngs_data"][reporter]
						reporter_edit_counts[reporter[i]].append(self.guides_data[guide]["ngs_data"][reporter])
					
					if reporter_edit_counts_total == 0: continue
					
					if flanks not in logo_data['reporter']:
						logo_data['reporter'][flanks] = {
							'A':sum(reporter_edit_counts['A'])/reporter_edit_counts_total, 
							'C':sum(reporter_edit_counts['C'])/reporter_edit_counts_total, 
							'T':sum(reporter_edit_counts['T'])/reporter_edit_counts_total, 
							'G':sum(reporter_edit_counts['G'])/reporter_edit_counts_total, 
							'N':sum(reporter_edit_counts['N'])/reporter_edit_counts_total,
							'reporter_edit_counts_total':reporter_edit_counts_total

						}
			flanks_list = list(logo_data['reporter'].keys())
			flanks_list.sort()
		
		except:
			print(utilities_error.getError())
				

	##################################################################################################

	def make_logo(self):	
		try:	
			logo_data = {
				'reporter':{},
				'edited_reporter':{},
				'guide':{}
			}
			nucs = ['A', 'C', 'T', 'G', 'N']

			print("#######",self.options['analysis_name'])
			for guide in self.guides_data.keys():
				if len(guide) == 0: continue
				if len(logo_data['guide']) == 0:
					for nuc in nucs:
						logo_data['guide'][nuc] = []

						for i in range(0,len(guide)):
							logo_data['guide'][nuc].append(0)
				
				for i in range(0,len(guide)):
					logo_data['guide'][guide[i]][i] += 1

				##---------------------------------##

				if len(logo_data['reporter']) == 0:
					for nuc in nucs:
						logo_data['reporter'][nuc] = []

						for i in range(0,len(guide)):
							logo_data['reporter'][nuc].append(0)
				
				for i in range(0,len(guide)):
					logo_data['reporter'][guide[i]][i] += 1

				##---------------------------------##
				
				for reporter in self.guides_data[guide]["ngs_data"]:
					if len(logo_data['edited_reporter']) == 0:
						for nuc in nucs:
							logo_data['edited_reporter'][nuc] = []

							for i in range(0,len(reporter)):
								logo_data['edited_reporter'][nuc].append(0)

					for i in range(0,len(reporter)):
						logo_data['edited_reporter'][reporter[i]][i] += self.guides_data[guide]["ngs_data"][reporter]/self.guides_data[guide]["ngs_total_count_filtered"]
			
			for pair_type in ['guide','reporter','edited_reporter']:
				print("###",pair_type)
				ngs_count = sum([logo_data[pair_type][nuc][0] for nuc in nucs])
				logo_data_normalised = utilities_pssm.normalise_division(logo_data[pair_type],ngs_count,aas=list(logo_data[pair_type].keys()))
				utilities_pssm.print_pssm(logo_data_normalised,aas=list(logo_data[pair_type].keys()))

			logo_data_normalised_reporter = utilities_pssm.pssm_normalise_column_sum_to_one_fold_change(logo_data['reporter'],logged=False,aas=['A','C','G','T'])
			logo_data_normalised_edited_reporter = utilities_pssm.pssm_normalise_column_sum_to_one_fold_change(logo_data['edited_reporter'],logged=False,aas=['A','C','G','T'])
			logo_data_divided = utilities_pssm.divide_pssms(logo_data_normalised_edited_reporter,logo_data_normalised_reporter,aas=['A','C','G','T'])

			print("#logo_data__reporter")
			utilities_pssm.print_pssm(logo_data['reporter'],aas=['A','C','G','T'])
			print("#logo_data_normalised_reporter")
			utilities_pssm.print_pssm(logo_data_normalised_reporter,aas=['A','C','G','T'])
			print("#logo_data_normalised_edited_reporter")
			utilities_pssm.print_pssm(logo_data_normalised_edited_reporter,aas=['A','C','G','T'])
			print("#logo_data_divided")
			utilities_pssm.print_pssm(logo_data_divided,aas=['A','C','G','T'])

			if self.options["subsample"]:
				plot_file = os.path.join(self.options["data_directory"],self.options['analysis_name'] + "." + str(self.options["subsample_cutoff"]) + '.change.logo.png')
			else:
				plot_file = os.path.join(self.options["data_directory"],self.options['analysis_name'] + '.change.logo.png')

			print("Writing: " + plot_file)
			utilities_pssm.plot_pssm_as_logo(logo_data_divided,outfile=plot_file,output="png")
		except:
			print(utilities_error.getError())
		
	##################################################################################################
