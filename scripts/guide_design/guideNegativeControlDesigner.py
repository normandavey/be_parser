import os, sys, inspect, json, pprint, random, subprocess, time, copy,re
				
file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../"))
import config_reader
import option_reader

sys.path.append(os.path.join(file_path,"../utilities"))
import utilities_error
import utilities_codons
import utilities_alignment
import utilities_basic
import utilities_oligos

sys.path.append(os.path.join(file_path,"../data_management/"))
import queryRunner

#-----
import logging
logger = logging.getLogger(__name__)
#-----

import guideDesigner

##------------------------------------------------------------------##
##------------------------------------------------------------------##
##------------------------------------------------------------------##

class guideNegativeControlDesigner():

	def __init__(self):
		self.options = {
			"intergenic_guide_sequence_length":37,
			"intergenic_guide_sequence_samples":20,
			"intergenic_guide_sequence_requirements":{
				"positions":{32:["G"]},
				"mismatch":{"start":11,"end":31}
		}}

	def make_negative_control_guides(self):

		self.options["intergenic_guide_sequence_length"] = (self.options['guide_length'] - self.options['pam_length']) + self.options["guide_reporter_flank_5'"] + self.options["guide_reporter_flank_3'"]
	
		selected_sequences = {}
		selected_sequences['intergenic'] = {}
		selected_sequences['non_targetting'] = {}
		iteration = 0

		random.seed(1)	

		guideDesignerObj = guideDesigner.guideDesigner()
		guideDesignerObj.options = self.options

		intergenic_guide_sequence_samples = self.options["intergenic_guide_sequence_samples"]+(self.options["intergenic_guide_sequence_samples"]*0.05)
		while len(selected_sequences['intergenic']) <  intergenic_guide_sequence_samples and len(selected_sequences['non_targetting']) < intergenic_guide_sequence_samples:
			tested_sequences = {}
		
			self.options['blast_db_path'] = os.path.join(self.options['data_path'],"genomes","Homo_sapiens.GRCh38.dna.primary_assembly.fa")
			length = os.stat(self.options['blast_db_path']).st_size

			iteration += 1
			logger.info(str(iteration) + ' iteration - ' + str(len(selected_sequences['intergenic'])) + " of " + str(self.options["intergenic_guide_sequence_samples"]))
			
			new_core_sequence = []
			new_shuffled_core_sequence = []

			chunks = 100

			with open(self.options['blast_db_path']) as input:
				while True:
					try:
						random_offset = random.randint(0, length - self.options["intergenic_guide_sequence_length"])
						input.seek(random_offset)
						sequence = input.read(self.options["intergenic_guide_sequence_length"]+1)
						sequence = sequence.replace("\n","")[:self.options["intergenic_guide_sequence_length"]]
						
						use_sequence = True
						if len(self.options["intergenic_guide_sequence_requirements"]['positions']) > 0:
							use_sequence = False
							for position in self.options["intergenic_guide_sequence_requirements"]['positions']:
								if sequence[int(position)] in self.options["intergenic_guide_sequence_requirements"]['positions'][position]:
									use_sequence = True

						if use_sequence == False:
							#logger.debug(sequence + " discarded - intergenic_guide_sequence_requirements")
							continue
						

						core_sequence = sequence[self.options["intergenic_guide_sequence_requirements"]['mismatch']['start']:self.options["intergenic_guide_sequence_requirements"]['mismatch']['end']]
						construct = guideDesignerObj.format_guide(core_sequence,second_insert=sequence)
						sequence_check_data = utilities_oligos.check_oligo_sequence_features(construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])
						if sequence_check_data['check_status']['all'] == False:
							logger.debug(construct + " discarded - oligo_sequence_features")
							continue
						
						if "intergenic_guide_sequence_requirements" in self.options:
							edit_window_requirements = []
							for edit_window_requirement in self.options["intergenic_guide_sequence_requirements"]['edit_window']:
								edit_window_requirements.append([core_sequence[x-1] for x in self.options['edit_window']].count(edit_window_requirement)>0)
							
							if edit_window_requirements.count(False) > 0:
								logger.debug(core_sequence +  " " + shuffled_core_sequence +  " discarded - intergenic_guide_sequence_requirements edit_window")
								continue

						##------------------------------------------------------------------------##
						use_shuffle = False
						use_shuffle_count = 0
						
						while use_shuffle == False and use_shuffle_count < 10:
							
							use_shuffle_count += 1
							use_shuffle_fail_reason = ""
							shuffled_core_sequence = ''.join(random.sample(core_sequence, len(core_sequence)))
							shuffled_core_sequence_extended = sequence.replace(core_sequence,shuffled_core_sequence)
							use_shuffle = True

							construct = guideDesignerObj.format_guide(shuffled_core_sequence,second_insert=shuffled_core_sequence_extended)
							shuffled_sequence_check_data = utilities_oligos.check_oligo_sequence_features(construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])
							if shuffled_sequence_check_data['check_status']['all'] == False:
								use_shuffle_fail_reason = construct + " discarded - shuffled oligo_seqeunce_features"
								use_shuffle = False
							
							if "intergenic_guide_sequence_requirements" in self.options:
								edit_window_requirements = []
								for edit_window_requirement in self.options["intergenic_guide_sequence_requirements"]['edit_window']:
									edit_window_requirements.append([shuffled_core_sequence[x-1] for x in self.options['edit_window']].count(edit_window_requirement)>0)
								
								if edit_window_requirements.count(False) > 0:
									use_shuffle_fail_reason = core_sequence +  " " + shuffled_core_sequence +  " discarded - intergenic_guide_sequence_requirements edit_window" + str([shuffled_core_sequence[x-1] for x in self.options['edit_window']])
									use_shuffle = False

							#logger.debug("Making shuffle " + core_sequence + " " + shuffled_core_sequence + " " + str(use_shuffle_count) + " " + str(use_shuffle) + " " + use_shuffle_fail_reason)

						if use_shuffle == False:
							logger.debug("Skipping " + core_sequence + " " + shuffled_core_sequence + " no shuffle with correct sequence properties found")
							continue

						if sequence not in tested_sequences:
							
							tested_sequences[core_sequence] = {
								"core_sequence":core_sequence,
								"core_sequence_extended":sequence,
								"shuffled_core_sequence":shuffled_core_sequence,
								"shuffled_core_sequence_extended":shuffled_core_sequence_extended
							}
							
							new_core_sequence.append(core_sequence)
							new_shuffled_core_sequence.append(shuffled_core_sequence)

		
						if len(new_core_sequence) > chunks:
							break
					except:
						pass
			
			options = {}
			options['remake'] = self.options['remake']
			options['oligo_sequence'] =  new_shuffled_core_sequence + new_core_sequence
			options["counts_only"] = True
			options["evalue"] = "1"
			options["annotate_all_exact_matches"] = False

			try:

				hash = utilities_basic.params_to_hash(options)
				guides_json_out_path = os.path.join(self.options['genomes_guides_fasta_path'], hash + ".gRNAs.blast.jsonc")
				if not os.path.exists(guides_json_out_path):
					logger.debug("Writing " + guides_json_out_path)
					random_sequences_counts = queryRunner.queryRunner("ensembl","find_sequence_human_genome",options).run()['data']
					utilities_basic.write_to_json(guides_json_out_path,random_sequences_counts,zipped=True,normalise_json=True)
				else:
					logger.debug("Reading " + guides_json_out_path)
					random_sequences_counts = utilities_basic.read_from_json(guides_json_out_path,zipped=True)

				for tested_sequence in new_core_sequence:
					#logger.debug("Checking similarity for " + tested_sequence)
					if tested_sequences[tested_sequence]['shuffled_core_sequence'] not in random_sequences_counts:
						random_sequences_counts[tested_sequences[tested_sequence]['shuffled_core_sequence']] = {}
				
					similar_guides = []
					for i in ["0","1","2","3","4","5"]:
						if tested_sequence in random_sequences_counts:
							if i in random_sequences_counts[tested_sequence]:
								similar_guides.append(random_sequences_counts[tested_sequence][i])
							else:
								similar_guides.append(0)
						else:
							similar_guides.append(0)

					similar_guides_shuffled = []

					for i in ["0","1","2","3","4","5"]:
						if tested_sequences[tested_sequence]['shuffled_core_sequence'] in random_sequences_counts:
							if i in random_sequences_counts[tested_sequences[tested_sequence]['shuffled_core_sequence']]:
								similar_guides_shuffled.append(random_sequences_counts[tested_sequences[tested_sequence]['shuffled_core_sequence']][i])
							else:
								similar_guides_shuffled.append(0)
						else:
							similar_guides_shuffled.append(0)
							
					if (similar_guides[0] == 1 and sum(similar_guides[1:4]) == 0) and sum(similar_guides_shuffled[0:4]) == 0:
						selected_sequences['intergenic'][tested_sequence] = {
							'sequence':tested_sequences[tested_sequence]['core_sequence'],
							'partner':tested_sequences[tested_sequence]['shuffled_core_sequence'],
							'sequence_extended':tested_sequences[tested_sequence]['core_sequence_extended'],
							'similarity':random_sequences_counts[tested_sequence]
						}

						selected_sequences['non_targetting'][tested_sequences[tested_sequence]['shuffled_core_sequence']] = {
							'sequence':tested_sequences[tested_sequence]['shuffled_core_sequence'],
							'partner':tested_sequence,
							'sequence_extended':tested_sequences[tested_sequence]['shuffled_core_sequence_extended'],
							'similarity':random_sequences_counts[tested_sequences[tested_sequence]['shuffled_core_sequence']]
						}
			except:
				logger.error(options)

			oligo_sequence_for_mapping = []
			
			for source in ['intergenic','non_targetting']:
				for sequence in list(selected_sequences[source].keys()):
					if 'gene_mapping' not in selected_sequences[source][sequence]:
						oligo_sequence_for_mapping.append(sequence)

			options = {}
			options['oligo_sequence'] = oligo_sequence_for_mapping
			options["counts_only"] = False
			options["evalue"] = "1"
			options["annotate_all_exact_matches"] = True
		
			random_sequences_annotated = queryRunner.queryRunner("ensembl","find_sequence_human_genome",options).run()
			
			filter_sequences_annotated = []
			for source in ['intergenic','non_targetting']: 
				for sequence in selected_sequences[source]:
					if sequence in oligo_sequence_for_mapping:
						selected_sequences[source][sequence]['gene_mapping'] = {}
						try:
							if len(random_sequences_annotated['data'][sequence]) != 0:
								selected_sequences[source][sequence]['gene_mapping'] = random_sequences_annotated['data'][sequence]

								for locus in selected_sequences[source][sequence]['gene_mapping']['off_target_mapping']:
									if len(selected_sequences[source][sequence]['gene_mapping']['off_target_mapping'][locus]) > 0:
										#if selected_sequences[source][sequence]['gene_mapping']['off_target_mapping'][locus][0]['exon_overlap']:
										print("Filtering",locus,sequence,selected_sequences[source][sequence]['gene_mapping']['off_target_mapping'][locus][0])
										filter_sequences_annotated.append(sequence)
							else:
								logger.debug(sequence + " not in random_sequences_annotated - " + source)
						except:
							logger.error([sequence,source])
			
			print("Filtering " + str(filter_sequences_annotated))
			for sequence in filter_sequences_annotated:
				try:
					if sequence in selected_sequences['intergenic']:
						del selected_sequences['non_targetting'][selected_sequences['intergenic'][sequence]['partner']]
						del selected_sequences['intergenic'][sequence]	
					if sequence in selected_sequences['non_targetting']:
						del selected_sequences['non_targetting'][sequence]
						del selected_sequences['intergenic'][selected_sequences['non_targetting'][sequence]['partner']]	
				except:
					logging.error(utilities_error.getError())

		return selected_sequences
