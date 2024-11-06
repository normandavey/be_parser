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

import guideNegativeControlDesigner
import guideSpliceSiteDesigner
import guideProteinDesigner
import guideAnnotation

#-----
import logging
logger = logging.getLogger(__name__)
#-----


class guideDesigner():
	##------------------------------------------------------------------##

	def __init__(self):

		###################################################

		self.options = {
			"job_file":"",
			"username":"",
			"password":"",
			"is_superuser":True,
			"library_name":"library",
			"testing":False,
			"species":"homo_sapiens",
			'skip_stops':True,
			'only_stops':False,
			'only_splice_sites':False,
			'include_designed_controls':False,
			'include_designed_guides':True,
			'include_splice_sites':False,
			'include_intergenic_sites':False,
			'include_nontargetting_sites':False,
			'include_default_controls':False,
			"default_controls_file":"./static/be_guide_controls_reporter_default.500_250_250.tdt"
		}

		self.options.update(config_reader.load_configeration_options(sections=["general"]))
		self.options.update(option_reader.load_commandline_options(self.options,{}))
		
		###################################################
		
		self.options['flanks'] = 20
		self.options["bar_code_length"] = 4
		self.options['background_source'] = "genomic"
		self.options['annotate_guide'] = True
		self.options['strand_types'] = ['positive','negative']
		self.options['upload_guides'] = False
		self.options['reset_remote_dataset']  = True

		###################################################

		self.options['guide_length'] = 23
		self.options["guide_reporter_flank_5'"] = 11
		self.options["guide_reporter_flank_3'"] = 6
		self.options['pam_length'] = 3
		self.options['editors'] = ["CBE","ABE","Dual"]
		self.options['output_editors'] = ["CBE","ABE","Dual"]
		self.options['edit_window'] = [4,5,6,7,8]

		self.options['target_nucleotides'] = {
			"ABE":{"G":"g","T":"t","C":"c","A":"A"},
			"CBE":{"G":"g","T":"t","C":"C","A":"a"},
			"Dual":{"G":"g","T":"t","C":"C","A":"A"}
		}

		self.options['editable_nucleotides'] = {
			"ABE":{"C":"c","A":"G","T":"t","G":"g"},
			"CBE":{"C":"T","A":"a","T":"t","G":"g"},
			"Dual":{"C":"T","A":"G","T":"t","G":"g"}
		}
		
		###################################################

		"""
		In summary we need to exclude guides with following sites:
		A) CGTCTC (Esp3I/BsmBI site)
		B) GAGACG (Esp3I/BsmBI site)
		C) TTTT (termination signal for RNA pol III)
		"""
		self.options['disallowed_sequences'] = ['CGTCTC','GAGACG','TTTTT']
		
		###################################################
		"""
		The oligonucleotide design is as follows (Esp3I/BsmBI sites underlined):
		5’-[Forward Primer]CGTCTCACACCG[sgRNA, 20 nt]GTTTCGAGACG[Reverse Primer]-3’

		5’- AGGCACTTGCTCGTACGACGCGTCTCACACCGNNNNNNNNNNNNNNNNNNNNGTTTCGAGACGTTAAGGTGCCGGGCCCACAT -3’


		Reporter
		[Fw_primer,20nt]CGTCTCACACCG[gRNA,20nt]GTTTGAGAGCTAGAAATAGCAAGTTCAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT[flank,11nt][gRNA,20nt][flank,10nt]GTTTCGAGACG[Rv_primer,20nt]

		The forward and reverse primer sequences can be the same as those we have used before, e.g.:
		Fw_primer:
		5'- AGGCACTTGCTCGTACGACG -3'

		Rv_primer:
		5'- TTAAGGTGCCGGGCCCACAT -3'

		"""

		report = {
			"all":{
				"5'":"AGGCACTTGCTCGTACGACGCGTCTCACACCG",
				"internal":"GTTTGAGAGCTAGAAATAGCAAGTTCAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT",
				"3'":"GTTTCGAGACGTTAAGGTGCCGGGCCCACAT"
			}
		}

		self.options['primers'] = {
			"NGG":{
				"5'":"AGGCACTTGCTCGTACGACGCGTCTCACACCG",
				"3'":"GTTTCGAGACGTTAAGGTGCCGGGCCCACAT"
			},
			"NGN":{
				"5'":"AGGCACTTGCTCGTACGACGCGTCTCACACCG",
				"3'":"GTTTCGAGACGTTAAGGTGCCGGGCCCACAT"
			},
			"all":{
				"5'":"AGGCACTTGCTCGTACGACGCGTCTCACACCG",
				"3'":"GTTTCGAGACGTTAAGGTGCCGGGCCCACAT"
			}
		}

		###################################################

		self.options['genomes_guides_path'] = os.path.join(self.options['data_path'],"genomes","gRNA")
		self.options['genomes_region_path'] = os.path.join(self.options['data_path'],"genomes","regions")
		self.options['genomes_guides_libraries_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","gRNA_libraries")
		self.options['genomes_guides_mapping_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","gRNA_mapping")
		self.options['genomes_guides_fasta_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","gRNA_fasta")
		self.options['genomes_guides_library_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","gRNA_libraries",self.options['library_name'])
		self.options['genomes_guides_library_regions_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","gRNA_libraries",self.options['library_name'],"gRNA_regions")
		self.options['genomes_guides_counts_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","gRNA_counts")
		self.options['genomes_guides_library_guides_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","gRNA_libraries",self.options['library_name'],"guides")
		
		if self.options['background_source'] == "genomic":
			self.options['genome_path'] = os.path.join(self.options['data_path'],"genomes","Homo_sapiens.GRCh38.dna.primary_assembly.collapsed.fa")

		if self.options['background_source'] == "cds":
			self.options['genome_path'] = os.path.join(self.options['data_path'],"genomes","Homo_sapiens.GRCh38.cds.all.collapsed.fa")
	
		if self.options['background_source'] == "cdna":
			self.options['genome_path'] = os.path.join(self.options['data_path'],"genomes","Homo_sapiens.GRCh38.cdna.all.collapsed_lines.fa")
	
		if not os.path.exists(self.options['genomes_guides_path']):
			os.mkdir(self.options['genomes_guides_path'])
		
		if not os.path.exists(self.options['genomes_region_path']):
			os.mkdir(self.options['genomes_region_path'])
		
		if not os.path.exists(self.options['genomes_guides_libraries_path']):
			os.mkdir(self.options['genomes_guides_libraries_path'])

		if not os.path.exists(self.options['genomes_guides_mapping_path']):
			os.mkdir(self.options['genomes_guides_mapping_path'])

		if not os.path.exists(self.options['genomes_guides_fasta_path']):
			os.mkdir(self.options['genomes_guides_fasta_path'])
		
		if not os.path.exists(self.options['genomes_guides_library_path']):
			os.mkdir(self.options['genomes_guides_library_path'])
	
		if not os.path.exists(self.options['genomes_guides_library_regions_path']):
			os.mkdir(self.options['genomes_guides_library_regions_path'])

		if not os.path.exists(self.options['genomes_guides_counts_path']):
			os.mkdir(self.options['genomes_guides_counts_path'])
			
		if not os.path.exists(self.options['genomes_guides_library_guides_path']):
			os.mkdir(self.options['genomes_guides_library_guides_path'])

		self.blosum_matrix = utilities_alignment.read_blosum_matrix()
		
		self.control_gRNAs = {}
		self.test_gRNAs = {}
		self.genes = {}
		self.region_guides_mapping = {}
			
	###-----------------_-------------_-_--_-_---_---_-------_------###
	
	def make_negative_control_guides(self):
		logger.info("Finding negative control guides")
		guideNegativeControlDesignerObj = guideNegativeControlDesigner.guideNegativeControlDesigner()
		guideNegativeControlDesignerObj.options = self.options
		return guideNegativeControlDesignerObj.make_negative_control_guides()

	###-----------------_-------------_-_--_-_---_---_-------_------###
	
	def make_splice_site_guides(self):
		
		logger.info("Finding splice sites for " + str(len(self.options['accession'])) + " proteins")
		guideSpliceSiteDesignerObj = guideSpliceSiteDesigner.guideSpliceSiteDesigner()
		guideSpliceSiteDesignerObj.options = self.options
		guideSpliceSiteDesignerObj.options['accessions_splice_site'] = self.options['accession']
		guideSpliceSiteDesignerObj.options['splice_site_set_type'] = "test"
		guideSpliceSiteDesignerObj.options['editors'] = self.options['splice_site_editors'] 
		return guideSpliceSiteDesignerObj.make_splice_site_guides()
	
	###-----------------_-------------_-_--_-_---_---_-------_------###

	def make_control_splice_site_guides(self):
		guideSpliceSiteDesignerObj = guideSpliceSiteDesigner.guideSpliceSiteDesigner()
		guideSpliceSiteDesignerObj.options = self.options

		if 'accessions_control_splice_site' in self.options:
			logger.info("Finding control splice sites from " + str(len(self.options['accessions_control_splice_site'])) + " proteins")
			guideSpliceSiteDesignerObj.options['accessions_splice_site'] = self.options['accessions_control_splice_site']
			guideSpliceSiteDesignerObj.options['splice_site_set_type'] = "control"
			guideSpliceSiteDesignerObj.options['editors'] = self.options['splice_site_editors'] 
			return guideSpliceSiteDesignerObj.make_splice_site_guides()
		else:
			logger.error("accessions_control_splice_site is not set")
			return {}

	###-----------------_-------------_-_--_-_---_---_-------_------###

	def make_guides(self):
		if self.options['include_designed_guides']:
			guide_designs_peptide_response = self.make_designed_guides()
		elif self.options['include_designed_controls']:
			guide_designs_peptide_response = self.make_designed_guides()

		if self.options['include_intergenic_sites'] or self.options['include_nontargetting_sites']:
			guide_designs_peptide_response = self.make_negative_control_guides()

			if self.options['include_intergenic_sites']:
				self.control_gRNAs['intergenic_control'] = guide_designs_peptide_response['intergenic']
			if self.options['include_nontargetting_sites']:
				self.control_gRNAs['non_targetting_control'] = guide_designs_peptide_response['non_targetting']

		
		if self.options['include_control_splice_sites']:
			guide_designs_peptide_response = self.make_control_splice_site_guides()
			self.control_gRNAs['control_splice_sites'] = guide_designs_peptide_response

		if self.options['include_splice_sites']:
			guide_designs_peptide_response = self.make_splice_site_guides()
			self.control_gRNAs['splice_sites'] = guide_designs_peptide_response
	
		if self.options['upload_guides']:
			self.uploadGuides()

		self.make_genescript_list()
	
	###-----------------_-------------_-_--_-_---_---_-------_------###
	
	def make_designed_guides(self):
		errors = {}

		guide_designs = {}
		mutation_designs = {}
		protein_designs = {}
		peptide_guide_candidates = {}
		no_guides = {}
		no_protein_data = {}

		out_path = os.path.join(self.options['genomes_guides_libraries_path'], self.options['library_name'] + ".regions.json")
		logger.debug("Writing: " + out_path )
		with open(out_path,'w') as outfile:
			json.dump(self.options['regions'], outfile)

		##---------------------------------##

		if len(self.options['residues']) > 0:
			protein_offsets = {}
			for accession in self.options['residues']:
				if accession not in protein_offsets:
					protein_offsets[accession] = []

				protein_offsets[accession] = []
				for x in list(self.options['residues'][accession].keys()):
					try:
						protein_offsets[accession].append(int(x))
					except:
						logger.error("make_designed_guides: offet - " + accession + " " + x)
		else:
			protein_offsets = {}
			for region in self.options['regions']:
				if region['accession'] not in protein_offsets:
					protein_offsets[region['accession']] = []
				
				protein_offsets[region['accession']] += list(range(int(region['peptide_start']),int(region['peptide_end']+1)))

		counter = 0

		logger.info("Processing " + str(len(protein_offsets.keys())) + " accessions")
		logger.info(",".join(list(protein_offsets.keys())))
		
		for accession in protein_offsets:
			try:
				#if accession != "P53990-2": continue
				counter += 1

				no_guides[accession] = copy.deepcopy(protein_offsets[accession])

				uniprot_basic_data = queryRunner.queryRunner("uniprot","parse_basic",{"accession":accession}).run()
				
				try:
					gene_name = uniprot_basic_data['data']['gene_name']
					protein_name = uniprot_basic_data['data']['protein_name']
					sequence = uniprot_basic_data['data']['sequence']
				except:
					gene_name = ""
					protein_name = ""
					sequence = ""
				
				logger.info("#"*30)
				logger.info("#" + str(counter) + "/" + str(len(protein_offsets)) + " - total guides:" + str(len(guide_designs)) + " - " + accession + ' ' + gene_name + ' ' + protein_name + ' ' + str(len(sequence)))
				logger.info("#"*30)
				sequence_peptides = [
					sequence[i].lower()
					if i+1 in protein_offsets[accession] else "-" 
					for i in range(0,len(sequence))
				]
				
				guideProteinDesignerObj = guideProteinDesigner.guideProteinDesigner()
				guideProteinDesignerObj.options.update(self.options)
				guideProteinDesignerObj.options['accession'] = [accession.replace("-1","")]
				guideProteinDesignerObj.options['mutant_residue_offsets'] = protein_offsets[accession]
				protein_guides = guideProteinDesignerObj.find_guides()
				
				protein_designs[accession] = {
					"gene_name":gene_name,
					"protein_name":protein_name,
					"sequence":sequence,
					"guides":[],
					"sequence_peptides":sequence_peptides,
					"sequence":sequence,
					"annotation":"-"
				}
					
				for gRNA in protein_guides['guide_designs']:
					try:
						for editor in protein_guides['guide_designs'][gRNA]:

							protein_guides['guide_designs'][gRNA][editor]['gene_name'] = gene_name
							if 'pos' not in  protein_guides['guide_designs'][gRNA][editor]: 
								logging.error("pos not in protein_guides['guide_designs'][gRNA][editor]")
								continue
							
							for offset in protein_guides['guide_designs'][gRNA][editor]['pos']:
								if offset in no_guides[accession]:
									no_guides[accession].remove(offset)

								if accession in self.options['residues']:
									if gRNA not in protein_designs[accession]['guides']: protein_designs[accession]['guides'].append(gRNA)
									if offset in self.options['residues'][accession]:
										if self.options['residues'][accession][offset] not in self.region_guides_mapping:
											self.region_guides_mapping[self.options['residues'][accession][offset]] = {}

										if gRNA not in self.region_guides_mapping[self.options['residues'][accession][offset]]:
											self.region_guides_mapping[self.options['residues'][accession][offset]][gRNA] = protein_guides['guide_designs'][gRNA][editor]['mutations']
									else:
										pass
								else:
									if gRNA not in protein_designs[accession]['guides']: protein_designs[accession]['guides'].append(gRNA)
									#logger.debug("No region_guides_mapping info")
									
					except:
						#logger.error(utilities_error.getError())
						logger.error([accession,gRNA,len(protein_guides['guide_designs']),utilities_error.getError()])
				
				##---------------------##
				
				guide_designs.update(protein_guides['guide_designs'])
				peptide_guide_candidates.update(protein_guides['peptide_guide_candidates'])
				mutation_designs.update(protein_guides['mutation_designs'])

				if accession in protein_guides['mutation_designs']:
					no_guides[accession] = list(set(protein_offsets[accession]).difference(set(list(protein_guides['mutation_designs'][accession].keys()))))
				else:
					no_guides[accession] = list(set(protein_offsets[accession]))

				logger.info(accession + ' - designed: ' + str(len(protein_offsets[accession])))
				logger.info(accession + ' - guides: ' + str(len(protein_guides['guide_designs'])))
				logger.info(accession + ' - not edited: ' + str(len(no_guides[accession])))

				if len(protein_guides['mutation_designs']):
					no_protein_data[accession] = {
						"offsets_length":protein_offsets[accession],
						"protein_guides":len(protein_guides['guide_designs'])
					}
			except:
				logger.error("# find_guides failed for " + accession)

		###-----------------_-------------_-_--_-_---_---_-------_------#

		single_guide_regions = {}
		for region in self.region_guides_mapping:
			if len(self.region_guides_mapping[region]) not in single_guide_regions:
				single_guide_regions[len(self.region_guides_mapping[region])] = 0

			single_guide_regions[len(self.region_guides_mapping[region])] += 1

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".guide_region_counts.json")
		logger.debug("Writing: " + out_path )
		with open(out_path,'w') as outfile:
			json.dump(single_guide_regions, outfile)

		###-----------------_-------------_-_--_-_---_---_-------_------#

		out_path = os.path.join(self.options['genomes_guides_libraries_path'], self.options['library_name'] + ".protein_designs.json")
		logger.debug("Writing: " + out_path )
		with open(out_path,'w') as outfile:
			json.dump(protein_designs, outfile)

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".mutation_designs.json")
		logger.debug("Writing: " + out_path )
		with open(out_path,'w') as outfile:
			json.dump(mutation_designs, outfile)

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".guides.json")
		logger.debug("Writing: " + out_path )
		with open(out_path,'w') as outfile:
			json.dump(guide_designs, outfile)
	
		###-----------------_-------------_-_--_-_---_---_-------_------#

		self.make_guides_guide_centric_list(protein_designs,mutation_designs,guide_designs)
		self.make_guides_list(guide_designs)
		self.make_guides_protein_centric_list(protein_designs,mutation_designs)

		###-----------------_-------------_-_--_-_---_---_-------_------#

		logger.info("#"*100)
		logger.info("# " + str(len(errors)) + " ERRORs")
		logger.info("#"*100)

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".errors.json")
		logger.debug("Writing: " + out_path )
		logger.debug(list(errors.keys()))
		with open(out_path,'w') as outfile:
			json.dump(errors, outfile)

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".no_guides.json")
		logger.debug("Writing: " + out_path )
		with open(out_path,'w') as outfile:
			json.dump(no_guides, outfile)

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".skipped_guides.json")
		logger.debug("Writing: " + out_path )
		with open(out_path,'w') as outfile:
			json.dump(peptide_guide_candidates, outfile)

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".no_protein_data.json")
		logger.debug("Writing: " + out_path )
		with open(out_path,'w') as outfile:
			json.dump(no_protein_data, outfile)
			
		###-----------------_-------------_-_--_-_---_---_-------_------#

		self.make_guides_stats(guide_designs)

		self.test_gRNAs["protein"] = protein_designs
		self.test_gRNAs["mutation"] = mutation_designs
		self.test_gRNAs["guide"] = guide_designs

		del protein_designs
		del mutation_designs
		del guide_designs

	###-----------------_-------------_-_--_-_---_---_-------_------###
	###------------- OUTPUTS OUTPUTS OUTPUTS  OUTPUTS ---_-----_----###
	###-----------------_-------------_-_--_-_---_---_-------_------###

	def annotate_guide_off_target(self,oligo_list,retry=False):
		guide_off_target_json_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name']  + ".guide_off_target.cache.jsonc")
		guide_off_target_no_complexity_filtering_json_path = os.path.join(self.options['genomes_guides_libraries_path'], self.options['library_name'] + ".no_complexity_filtering.guide_off_target.cache.jsonc")
		guide_off_target = {}
		
		if not os.path.exists(guide_off_target_json_path):
			logging.debug("Writing: " + guide_off_target_json_path)
			guide_off_target= queryRunner.queryRunner("ensembl","find_sequence_human_genome",{"oligo_sequence":oligo_list,"counts_only":True,"remake":False}).run()
			if 'data' not in guide_off_target:
				pprint.pprint(guide_off_target)

			guide_off_target = guide_off_target['data']
			utilities_basic.write_to_json(guide_off_target_json_path,guide_off_target,zipped=True)
			
		else:
			logging.debug("Reading: " + guide_off_target_json_path)
			guide_off_target = utilities_basic.read_from_json(guide_off_target_json_path,zipped=True)
			
		if len(guide_off_target) != len(oligo_list):
			logging.error([guide_off_target_json_path,len(guide_off_target),len(oligo_list)])
			logging.error("Deleting: Guides missing - " + guide_off_target_json_path)

			if not os.path.exists(guide_off_target_no_complexity_filtering_json_path):
				logging.debug("Writing off target with no complextity filtering: " + guide_off_target_no_complexity_filtering_json_path)
				missing_oligo_list = list(set(oligo_list).difference(guide_off_target.keys()))
				guide_off_target_no_complexity_filtering = queryRunner.queryRunner("ensembl","find_sequence_human_genome",{"oligo_sequence":missing_oligo_list,"counts_only":True,"use_filtering":False,"chunks_size":10,"remake":True}).run()
				
				guide_off_target.update(guide_off_target_no_complexity_filtering['data'])
				utilities_basic.write_to_json(guide_off_target_no_complexity_filtering_json_path, guide_off_target_no_complexity_filtering,zipped=True)
			else:
				logging.debug("Reading: " + guide_off_target_no_complexity_filtering_json_path)
				guide_off_target_no_complexity_filtering = utilities_basic.read_from_json(guide_off_target_no_complexity_filtering_json_path,zipped=True)
				guide_off_target.update(guide_off_target_no_complexity_filtering['data'])
			
			if retry == False: 
				pass
			else:
				logging.error("Guides missing - " + guide_off_target_json_path)
				logging.error("Retry failed")
		
		guide_off_target_stats = {"input":len(oligo_list),"output":len(guide_off_target)}
		for guide in guide_off_target:
			try:
				if guide_off_target[guide]["0"] not in guide_off_target_stats:
					guide_off_target_stats[guide_off_target[guide]["0"] ] = 0
			
				guide_off_target_stats[guide_off_target[guide]["0"]] += 1
			except:
				logging.error("Guides '0' missing - " + guide_off_target_json_path)

		guide_off_target_stats_json_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name']  + ".guide_off_target.stat.json")
		utilities_basic.write_to_json(guide_off_target_stats_json_path, guide_off_target_stats)

		logging.debug("#### OFF TARGET Statistics ####")
		pprint.pprint(guide_off_target_stats)

		return guide_off_target

	def make_guides_list(self,guide_designs):
		guide_rows = []

		common_tags = [
			"peptide_sequence",
			"accession",
			"gene_name",
			"peptide_start",
			"peptide_end",
			"strand_type",
			"exact_matches_guide_coordinates",
			"peptide_genomic_flanks",
			"exact_matches",
			"exact_matches_mappings"
		]

		tags = [
			"peptide_sequence_translated",
			"mutated_oligo",
			"mutations",
			"codon_wt",
			"codon_mut",
			"pos",
			"wt",
			"mut",
			"wildtype_region",
			"mutated_region",
			"changes_region",
			"stop_region",
			"exon_boundary_overlap"
		]

	
		headers = ["guide","pam"] + common_tags
		for output_editor in self.options['output_editors']:    
			headers += [output_editor + "_"  + tag for tag in tags]
		
		guide_rows.append("\t".join(headers))

		counter = 0
		guide_match_counter = {}
		guide_exon_match_counter = {}
		guide_off_target_data = {}
		if self.options['annotate_guide']:
			self.guideAnnotationObj = guideAnnotation.guideAnnotation()
			self.guideAnnotationObj.options['species'] = "9606"
			guide_off_target_data = self.annotate_guide_off_target(list(guide_designs.keys()))
		
		for guide in guide_designs:
			counter += 1
			#logger.debug("Adding " + guide + " " + str(counter) + " of " + str(len(guide_designs)))
			#self.options['sequence'] = [guide]		
			#logger.debug("Adding " + guide)		
		
			try:
				if self.options['annotate_guide']:
					exact_matches = "-"
					exact_matches_mappings = "-"
					exact_matches_guide_coordinates = "-"

					if guide in guide_off_target_data:
						
						exact_matches = str(guide_off_target_data[guide]["0"] if "0" in guide_off_target_data[guide] else 0)
						exact_matches_mappings = "|".join([str(i) + ":" + str(guide_off_target_data[guide][str(i)]) if str(i) in guide_off_target_data[guide] else str(i) + ":" + "0" for i in range(0,4)])

				else:
					exact_matches = "-"
					exact_matches_mappings = "-"
					exact_matches_guide_coordinates = "-"
			except:
				exact_matches = "?"
				exact_matches_mappings = "?"
				exact_matches_guide_coordinates = "?"
			
			row_data = {
				'common':[],
			}

			for editor in self.options['editors']:
				row_data[editor] = editor
			
			use_row = False 
			for editor in self.options['output_editors']:
				if editor in guide_designs[guide]:
					use_row = True
					guide_designs[guide][editor]['guide'] = guide
					guide_designs[guide][editor]['exact_matches'] = exact_matches
					guide_designs[guide][editor]['exact_matches_mappings'] = exact_matches_mappings
					guide_designs[guide][editor]['exact_matches_guide_coordinates'] = exact_matches_guide_coordinates
					row_data['common'] = [guide[:20],guide[20:]]
					
					for tag in common_tags:
						if tag in guide_designs[guide][editor]:
							row_data['common'].append(str(guide_designs[guide][editor][tag]))
						else:
							print(tag)
							row_data['common'].append("")
							logger.error(guide + " " + editor + " " + tag)

					row_data[editor] = []
					for tag in tags:
						if tag in guide_designs[guide][editor]:
							try:
								if tag in ["peptide_sequence_translated","mutated_oligo","wildtype_region","mutated_region"]:
									row_data[editor].append(guide_designs[guide][editor][tag])
								elif tag in ["pos"]:
									row_data[editor].append(",".join([str(x) for x in guide_designs[guide][editor][tag]]))
								elif tag in ["stop_region","exon_boundary_overlap"]:
									row_data[editor].append(str(guide_designs[guide][editor][tag]))
								else:
									row_data[editor].append(",".join(guide_designs[guide][editor][tag]))
							except:
								#logging.error(guide + " " + tag) 
								row_data[editor].append("?")

				else:
					row_data[editor] = [""]*len(tags)

			if use_row:
				row_list = row_data['common']
				for output_editor in self.options['output_editors']:    
					row_list += row_data[output_editor]

				guide_rows.append("\t".join(row_list))
			
		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".guides.tdt")
		logger.debug("Writing: " + out_path )
		open(out_path,'w').write("\n".join(guide_rows))

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".guides.match.stats.json")
		logger.debug("Writing: " + out_path )
		with open(out_path,'w') as outfile:
			json.dump(guide_match_counter, outfile)

		print("#guide_match_counter")
		logger.debug(guide_match_counter)

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".guides.exon_match.stats.json")
		logger.debug("Writing: " + out_path )
		with open(out_path,'w') as outfile:
			json.dump(guide_exon_match_counter, outfile)

		print("#guide_exon_match_counter")
		logger.debug(list(guide_exon_match_counter.keys()))

	def make_guides_protein_centric_list(self,protein_designs,mutation_designs):
		region_counts = 0
		guide_counts = 0
		accession_counts = 0
		targeted_residues = 0
		mutated_residues = 0
		
		protein_rows = []

		for accession in protein_designs:
			try:
				accession_counts += 1
			
				#regions = [protein_designs[accession]['sequence_peptides'][i] if str(i) in list(mutation_designs[accession].keys()) else "x" for i in range(0,len(protein_designs[accession]['sequence_peptides']))]
				#region_counts += len(regions)

				if accession in mutation_designs:
					regions = [
						protein_designs[accession]['sequence'][i].upper() 
						if i+1 in list(mutation_designs[accession].keys()) else "-"#protein_designs[accession]['sequence_peptides'][i].lower()
						for i in range(0,len(protein_designs[accession]['sequence_peptides']))
					]
				else:
					regions = [
						protein_designs[accession]['sequence_peptides'][i].upper() 
						for i in range(0,len(protein_designs[accession]['sequence_peptides']))
					]

				
				region_counts += len(regions)
				

				targeted_residues += len(protein_designs[accession]['sequence_peptides']) - protein_designs[accession]['sequence_peptides'].count("-")
				mutated_residues += len(mutation_designs[accession].keys())
				guide_counts += len(protein_designs[accession]['guides'])

				edit_string = "".join(regions).replace("x","-")
				target_string = "".join(protein_designs[accession]['sequence_peptides']) 

				#print(accession,"\t",protein_designs[accession]['gene_name'],"\t",len(protein_designs[accession]['guides']))
				
				proviz_link = 'https://slim.icr.ac.uk/proviz/proviz.php?uniprot_acc=' + accession + '&tools=alphafold&tracks=peptides,edited residues,FB4,' + edit_string + ',0,' + str(len(edit_string)) + ';peptides,targeted residues,f0d9b1,' + target_string + ',0,' + str(len(target_string))
				edit_string_proviz_link = '=hyperlink("' + proviz_link + '","targeted residues")'
				
				depmap_score = ""
				try:
					depmap_response = queryRunner.queryRunner("depmap","get_gene_dependency",{"accession":accession}).run()
					depmap_score = "%1.2f"%depmap_response['data']['CRISPR']['mean_gene_effect']
				except:
					pass
				
				coverage = "0"
				try:
					coverage ="%1.2f"%(len(mutation_designs[accession].keys())/len(regions))
				except:
					pass

				protein_row = [	
					str(len(protein_designs[accession]['guides'])),
					#str(guide_counts),
					#str(accession_counts),
					#str(region_counts),
					str(len(regions)),
					str(len(protein_designs[accession]['sequence_peptides']) - protein_designs[accession]['sequence_peptides'].count("-")),
					str(len(mutation_designs[accession].keys())),
					coverage,
					protein_designs[accession]['annotation'],
					accession,
					protein_designs[accession]['gene_name'],
					protein_designs[accession]['protein_name'],
					depmap_score,
					edit_string,
					edit_string_proviz_link,
					#",".join(protein_designs[accession]['guides']),
					",".join([str(offset) for offset in mutation_designs[accession].keys()])
					]
				
				protein_rows.append("\t".join(protein_row))
			except:
				logger.error(accession)
				logger.error(utilities_error.getError())

			#print(len(protein_designs))
			
		header = [
			"guides count",
			"protein regions",
			"targeted residues total",
			"mutated residues total",
			"targetted/mutated",
			"annotation",
			"accession",
			"gene",
			"protein",
			"depmap_score",
			"region sequence",
			"mutated_residues"
			]
		
		logger.info("#"*100)
		logger.info("# PROTEIN INFO")
		logger.info("#"*100)
		logger.info("\t".join(header))
		
		logger.info("\n".join(protein_rows))

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".protein.tdt")

		logger.debug("Writing: " + out_path )
		open(out_path,'w').write("\n".join(["\t".join(header)] + protein_rows))

		logger.info("#"*100)
		logger.info("# PROTEIN STATS")
		logger.info("#"*100)

		protein_stats = {
			"guide_counts":str(guide_counts),
			"accession_counts":str(accession_counts),
			"region_counts":str(region_counts),
			"targeted_residues":str(targeted_residues),
			"mutated_residues":str(mutated_residues),
		}

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".protein_stats.json")
		logger.debug("Writing: " + out_path)

		with open(out_path,'w') as outfile:
			json.dump(protein_stats, outfile)

	def format_guide(self,guide_sequence,pam_type="all",second_insert=""):
		if "internal" not in self.options['primers'][pam_type]:
			return self.options['primers'][pam_type]["5'"].lower() + guide_sequence + self.options['primers'][pam_type]["3'"].lower() 
		else:
			return self.options['primers'][pam_type]["5'"].lower()  + guide_sequence + self.options['primers'][pam_type]["internal"].lower()  + second_insert + self.options['primers'][pam_type]["3'"].lower() 
	
	def make_bar_code(self,guide='',second_insert=''):
		bar_code = ""
		bar_code_created = False
		bar_code_iteration_counter = 0
		if 'barcode_guide_construct' in self.options:
			if self.options['barcode_guide_construct']:
				while bar_code_created == False: 	
					bar_code_iteration_counter += 1
					
					bar_code = "".join([["A","T","C","G"][int(random.random()*self.options["bar_code_length"])] for x in range(0,self.options["bar_code_length"])])

					guide_construct = self.format_guide(guide,second_insert=second_insert + bar_code)
					sequence_check_data = utilities_oligos.check_oligo_sequence_features(guide_construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])
					
					if sequence_check_data['check_status']['all'] == False:
						bar_code_created = False
					elif [bar_code.count(nuc) == self.options["bar_code_length"] for nuc in ["A","T","C","G"]].count(True) == 1:
						bar_code_created = False
					else:
						bar_code_created = True

					if bar_code_iteration_counter > 10:
						bar_code_created = True

		return bar_code
		
	def make_genescript_list(self):
		genescript_table = ["\t".join(["*Sequence(5'to3')","No.","guide","reporter","bar_code","annotation","passes_sequence_feature_checks","construct_length","type"]) + "\n"]
		genescript_repeat_table = ["\t".join(["*Sequence(5'to3')","No.","guide","reporter","bar_code","annotation","passes_sequence_feature_checks","construct_length","type"]) + "\n"]
		
		counter = 0
		guide_constructs = []
		splice_sites = []
		control_splice_sites = []
		target_guides = []
		intergenic_control = []
		non_targetting_control = []
		default_controls = []
		
		if self.options['include_designed_guides'] or self.options['include_designed_controls']:
			protein_designs = self.test_gRNAs["protein"] 
			guide_designs = self.test_gRNAs["guide"] 

			for guide in guide_designs:
				try:
					pam_type = ""
					guide_candidate_reporter = ""
					mutations = {}

					for editor in self.options['editors']:
						mutations[editor] = ""

					mutations_editors = []
					mutations_editors_mutation = []
					
					for editor in self.options['editors']:
						if editor in guide_designs[guide] and editor in self.options['output_editors']:
							pam_type = guide_designs[guide][editor]['pam_type'] if 'pam_type' in guide_designs[guide][editor] else "-"
							guide_candidate_reporter = guide_designs[guide][editor]['guide_candidate_reporter']
							guide_sequence = guide_designs[guide][editor]['guide']
							accession = guide_designs[guide][editor]['accession']
							peptide_start =  guide_designs[guide][editor]['peptide_start']
							peptide_end =  guide_designs[guide][editor]['peptide_end']
							mutations[editor] = ":".join([editor, ",".join(guide_designs[guide][editor]['changes_region'])])
							mutations_editors.append(editor)

							mutations_editors_mutation.append(mutations[editor])
						else:
							mutations_editors_mutation.append("-")
					
					if len((set(self.options['output_editors']).intersection(mutations_editors))) > 0:
						annotation = ''
						annotation = "|".join([accession,protein_designs[accession]['gene_name'],",".join(mutations_editors_mutation)])

						bar_code = self.make_bar_code(guide=guide_sequence[:20],second_insert=guide_candidate_reporter)
						guide_construct = self.format_guide(guide_sequence[:20],second_insert=guide_candidate_reporter + bar_code)
						guide_candidate_reporter = guide_candidate_reporter + bar_code
						
						sequence_check_data = utilities_oligos.check_oligo_sequence_features(guide_construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])
						row = [guide_construct,str(counter),guide_sequence[:20],guide_candidate_reporter,bar_code,annotation,",".join(sequence_check_data['check_status']['failed_check_type']),str(len(guide_construct)),"targets"]
						
						if sequence_check_data['check_status']['all']:
							counter += 1
							#print("\t".join(row) + "\n" )
							
							if guide_sequence[:20] not in guide_constructs:
								genescript_table.append("\t".join(row))
							else:
								genescript_repeat_table.append("\t".join(row))
								
							guide_constructs.append(guide_sequence[:20])
							target_guides.append(guide_construct)
				except:
					logger.error(guide + " not added to genscript oligo list")
					logger.error(utilities_error.getError())

		if self.options['include_intergenic_sites']:
			for guide in self.control_gRNAs['intergenic_control'].keys():
				
				guide_sequence = self.control_gRNAs['intergenic_control'][guide]['sequence']
				guide_sequence_extended = self.control_gRNAs['intergenic_control'][guide]['sequence_extended']

				bar_code = self.make_bar_code(guide=guide_sequence[:20],second_insert=guide_sequence_extended)

				guide_construct = self.format_guide(guide_sequence[:20],second_insert=guide_sequence_extended + bar_code)
				guide_candidate_reporter = guide_sequence_extended + bar_code
				pam_type = "?"
				
				try:
					annotation = str(self.control_gRNAs['intergenic_control'][guide]['gene_mapping']['off_target_mapping'])
				except:
					annotation = ""

				sequence_check_data = utilities_oligos.check_oligo_sequence_features(guide_construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])
				row = [guide_construct,str(counter),guide_sequence[:20],guide_candidate_reporter,bar_code,annotation,",".join(sequence_check_data['check_status']['failed_check_type']),str(len(guide_construct)),"intergenic_control"]
				
				if sequence_check_data['check_status']['all'] and len(intergenic_control) < self.options["intergenic_guide_sequence_samples"]:
					counter += 1
					#print("\t".join(row) + "\n" )
					intergenic_control.append(guide_construct)
					
					if guide_sequence[:20] not in guide_constructs:
						genescript_table.append("\t".join(row))
					else:
						genescript_repeat_table.append("\t".join(row))

					guide_constructs.append(guide_sequence[:20])
			
		if self.options['include_nontargetting_sites']:
			for guide in list(self.control_gRNAs['non_targetting_control'].keys()):
				
				guide_sequence = self.control_gRNAs['non_targetting_control'][guide]['sequence']
				guide_sequence_extended = self.control_gRNAs['non_targetting_control'][guide]['sequence_extended']
				bar_code = self.make_bar_code(guide=guide_sequence[:20],second_insert=guide_sequence_extended)

				guide_construct = self.format_guide(guide_sequence[:20],second_insert=guide_sequence_extended + bar_code)
				guide_candidate_reporter = guide_sequence_extended + bar_code
				pam_type = "?"

				try:
					annotation = str(self.control_gRNAs['non_targetting_control'][guide]['gene_mapping']['off_target_mapping'])
				except:
					annotation = ""

				
				sequence_check_data = utilities_oligos.check_oligo_sequence_features(guide_construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])
				row = [guide_construct,str(counter),guide_sequence[:20],guide_candidate_reporter,bar_code,annotation,",".join(sequence_check_data['check_status']['failed_check_type']),str(len(guide_construct)),"non_targetting_control"]
				
				if sequence_check_data['check_status']['all'] and len(non_targetting_control) < self.options["intergenic_guide_sequence_samples"]:
					counter += 1
					#print("\t".join(row) + "\n" )
					
					if guide_sequence[:20] not in guide_constructs:
						genescript_table.append("\t".join(row))
					else:
						genescript_repeat_table.append("\t".join(row))

					guide_constructs.append(guide_sequence[:20])
					non_targetting_control.append(guide_construct)

		if self.options['include_splice_sites']:
			for guide in self.control_gRNAs['splice_sites'].keys():
				guide_sequence = guide
				for editor in self.control_gRNAs['splice_sites'][guide]:

					guide_sequence_extended = self.control_gRNAs['splice_sites'][guide][editor]['guide_candidate_extended']
					bar_code = self.make_bar_code(guide=guide_sequence[:20],second_insert=guide_sequence_extended)
					guide_construct = self.format_guide(guide_sequence[:20],second_insert=guide_sequence_extended + bar_code)

					guide_candidate_reporter = guide_sequence_extended + bar_code
					pam_type = self.control_gRNAs['splice_sites'][guide][editor]['pam_sequence']
					annotation = "|".join([
						self.control_gRNAs['splice_sites'][guide][editor]['accession'] if "accession" in self.control_gRNAs['splice_sites'][guide][editor] else "-",
						self.control_gRNAs['splice_sites'][guide][editor]['gene_name'],
						self.control_gRNAs['splice_sites'][guide][editor]['strand'],
						self.control_gRNAs['splice_sites'][guide][editor]['splice_site_type'],
						"intron_" + str(self.control_gRNAs['splice_sites'][guide][editor]['intron_count'])
					])

				sequence_check_data = utilities_oligos.check_oligo_sequence_features(guide_construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])

				row = [guide_construct,str(counter),guide_sequence[:20],guide_candidate_reporter,bar_code,annotation,",".join(sequence_check_data['check_status']['failed_check_type']),str(len(guide_construct)),"splice_sites"]
				
				if sequence_check_data['check_status']['all']:
					counter += 1
					#print("\t".join(row) + "\n" )
					
					if guide_sequence[:20] not in guide_constructs:
						genescript_table.append("\t".join(row))
					else:
						genescript_repeat_table.append("\t".join(row))

					guide_constructs.append(guide_sequence[:20])
					splice_sites.append(guide_construct)

		if self.options['include_control_splice_sites']:
			for guide in self.control_gRNAs['control_splice_sites'].keys():
				guide_sequence = guide
				for editor in self.control_gRNAs['control_splice_sites'][guide]:

					guide_sequence_extended = self.control_gRNAs['control_splice_sites'][guide][editor]['guide_candidate_extended']
					bar_code = self.make_bar_code(guide=guide_sequence[:20],second_insert=guide_sequence_extended)

					guide_construct = self.format_guide(guide_sequence[:20],second_insert=guide_sequence_extended + bar_code)
					guide_candidate_reporter = guide_sequence_extended + bar_code
					pam_type = self.control_gRNAs['control_splice_sites'][guide][editor]['pam_sequence']
					annotation = "|".join([
						self.control_gRNAs['control_splice_sites'][guide][editor]['gene_name'],
						self.control_gRNAs['control_splice_sites'][guide][editor]['strand'],
						self.control_gRNAs['control_splice_sites'][guide][editor]['splice_site_type'],
						"intron_" + str(self.control_gRNAs['control_splice_sites'][guide][editor]['intron_count'])
					])

				sequence_check_data = utilities_oligos.check_oligo_sequence_features(guide_construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])
				row = [guide_construct,str(counter),guide_sequence[:20],guide_candidate_reporter,bar_code,annotation,",".join(sequence_check_data['check_status']['failed_check_type']),str(len(guide_construct)),"control_splice_sites"]
				
				if sequence_check_data['check_status']['all']:
					counter += 1
					#print("\t".join(row) + "\n" )
					control_splice_sites.append(guide_construct)
					
					if guide_sequence[:20] not in guide_constructs:
						genescript_table.append("\t".join(row))
					else:
						genescript_repeat_table.append("\t".join(row))

					guide_constructs.append(guide_sequence[:20])

		if self.options['include_default_controls']:
			for line in open(self.options["default_controls_file"]).read().split("\n")[1:]:
				counter += 1
	
				line_bits = line.split("\t")
					
				guide_sequence = line_bits[0]
				guide_sequence_extended = line_bits[1]
				bar_code = line_bits[2]
				control_type = line_bits[4]
				annotation = line_bits[3]

				guide_construct = self.format_guide(guide_sequence,second_insert=guide_sequence_extended + bar_code)
				guide_candidate_reporter = guide_sequence_extended + bar_code
				
				sequence_check_data = utilities_oligos.check_oligo_sequence_features(guide_construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])
				row = [guide_construct,str(counter),guide_sequence,guide_candidate_reporter,bar_code,annotation,",".join(sequence_check_data['check_status']['failed_check_type']),str(len(guide_construct)),control_type]
				
				if sequence_check_data['check_status']['all']:
					counter += 1
					default_controls.append(guide_construct)
					
					if guide_sequence not in guide_constructs:
						genescript_table.append("\t".join(row))
					else:
						genescript_repeat_table.append("\t".join(row))

					guide_constructs.append(guide_sequence)

					if control_type == "control_splice_sites":	
						control_splice_sites.append(guide_construct)
					
					if control_type == "intergenic_control":	
						intergenic_control.append(guide_construct)

					if control_type == "non_targetting_control":	
						non_targetting_control.append(guide_construct)

		#----##----##----##----##----##----##----#

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".genescript.gRNA.tdt")
		open(out_path,'w').write("\n".join(genescript_table))

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".genescript.repeated.gRNA.tdt")
		open(out_path,'w').write("\n".join(genescript_repeat_table))

		print(str(len(genescript_table) - 1) + " gRNAs")
		print(str(len(genescript_repeat_table) - 1) + " repeated gRNAs")
		print(str(len(set(guide_constructs))) + " unique gRNAs")
		print(str(len(target_guides)) + " target guides gRNAs")
		print(str(len(splice_sites)) + " target splice sites gRNAs")
		print(str(len(control_splice_sites)) + " control splice sites gRNAs")
		print(str(len(intergenic_control)) + " intergenic control gRNAs")
		print(str(len(non_targetting_control)) + " non targetting control gRNAs")
		print(str(len(default_controls)) + " default control set gRNAs")

		logger.debug("Writing: " + out_path )

	def make_guides_guide_centric_list(self,protein_designs,mutation_designs,guide_designs):
		guide_rows = {
			'all':{}
		}

		for editor in self.options['editors']:
			guide_rows[editor] = []

		for accession in protein_designs:
			if accession not in mutation_designs:
				logger.error(accession + " not in mutation_designs")
				continue

			offsets = list(mutation_designs[accession].keys())
			offsets.sort()

			for offset in offsets:
				for editor in self.options['editors']:
					if editor in mutation_designs[accession][offset]['guides']:
						for strand in ['negative','positive']:
							for guide in mutation_designs[accession][offset]['guides'][editor][strand]:
							
								try:
									guide_mutation_blosum_effect = str(min(mutation_designs[accession][offset]['guides'][editor][strand][guide]['mutation_blosum_effect']))
								except:
									guide_mutation_blosum_effect = "?"

								try:
									offset_mutation_blosum_effect = str(min(mutation_designs[accession][offset]['mutation_blosum_effect']))
								except:
									offset_mutation_blosum_effect = "?"

								try:
									peptide_genomic_flanks = guide_designs[guide][editor]["peptide_genomic_flanks"]
								except:
									peptide_genomic_flanks = "?"

								guide_sequence = guide[0:-3]
								pam_sequence = guide[-3:]
								pam_type = "NGN"
								
								if pam_sequence[-2:] == "GG":
									pam_type = "NGG"

								guide_row = [ 
									editor,
									strand,
									guide[:20],
									guide[20:],
									peptide_genomic_flanks,
									pam_sequence,
									pam_type,
									#self.options['primers']["pam_type"]["5'"] + guide_sequence + self.options['primers'][pam_type]["3'"],
									self.format_guide(guide_sequence),
									"%1.2f"%(float(guide.count('G') + guide.count('C'))/len(guide)),
									str(float(guide.count('G') + guide.count('C'))/len(guide) > 0.4 and float(guide.count('G') + guide.count('C'))/len(guide) < 0.75),
									#str(guide.count('CGTCTC') > 0),
									#str(guide.count('TTTT') > 0),
									mutation_designs[accession][offset]['wildtype_residue'],
									",".join(mutation_designs[accession][offset]['guides'][editor][strand][guide]['mutations']),
									guide_mutation_blosum_effect,
									offset_mutation_blosum_effect,
									offset,
									accession,
									protein_designs[accession]['gene_name'],
									protein_designs[accession]['protein_name']
									]
								
								guide_rows[editor].append("\t".join(guide_row))

		
		header_guide_row = [ 
			"editor",
			"strand",
			"guide_targetting",
			"guide_pam",
			"guide_genomic_flanks",
			"guide_pam_sequence",
			"guide_pam_type",
			"guide_formatted"
			"guide_gc_content",
			"guide_gc_content_check",
			"guide_wildtype_residue",
			"guide_mutated_residue",
			"guide_mutation_blosum_effect",
			"offset_mutation_blosum_effect",
			"offset",
			"accession",
			"gene_name",
			"protein_name"
		]

		for editor in self.options['editors']:
			out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + "." + editor + ".gRNA.tdt")
			logger.debug("Writing: " + out_path )
			open(out_path,'w').write("\n".join(["\t".join(header_guide_row)] + guide_rows[editor]))
		
		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".all.gRNA.tdt")
		logger.debug("Writing: " + out_path )

		all_gRNA = ["\t".join(header_guide_row)]
		for editor in self.options['editors']:
			all_gRNA += guide_rows[editor]

		open(out_path,'w').write("\n".join(all_gRNA))

	def make_guides_stats(self,guide_designs):

		guide_design_specificity = {"Both":0,"total":0}

		for editor in self.options['editors']:
			guide_design_specificity[editor] = 0
			
		for guide_design in guide_designs:
			if "ABE" in guide_designs[guide_design] and "CBE" in guide_designs[guide_design]:
				try:
					if len(guide_designs[guide_design]['ABE']['mutations']) > 0 and len(guide_designs[guide_design]['CBE']['mutations']) > 0:
						guide_design_specificity["Both"] += 1
				except:
					pass

			for editor in self.options['editors']:
				if editor not in guide_designs[guide_design]: continue

				if 'mutations' in guide_designs[guide_design][editor]:
					if len(guide_designs[guide_design][editor]['mutations']) > 0:
						guide_design_specificity[editor] += 1

			guide_design_specificity["total"] += 1
			
		logger.info("#"*100)
		logger.info("# GUIDE STATS")
		logger.info("#"*100)

		logger.info(guide_design_specificity)

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + ".editor_stats.json")
		logger.debug("Writing: " + out_path)

		with open(out_path,'w') as outfile:
			json.dump(guide_design_specificity, outfile)


	###-----------------_-------------_-_--_-_---_---_-------_------###
	###------------- UPLOAD  UPLOAD  UPLOAD   UPLOAD  ---_-----_----###
	###-----------------_-------------_-_--_-_---_---_-------_------###
	def makeGuidesRemoteDataDataset(self):
		return self.makeGuidesRemoteDataset("Motif Mutation Guides","motif_mutation_guides")

	def makeGuidesRemoteIndexDataset(self):
		return self.makeGuidesRemoteDataset("Motif Mutation Guides Index","motif_mutation_guides_index")

	def makeGuidesRemoteDataset(self,dataset_name,dataset_id):

		# Get current dataset
		logger.info("Get current dataset:" + dataset_name)
		options = {
			"database_name":"datasets",
			"dataset_name":dataset_name,
			"username":self.options['username'],
			"password":self.options['password'],
			"is_superuser":self.options['is_superuser']
		}

		response = queryRunner.queryRunner("database_access","get_dataset_details_by_name",options).run()
		dataset_ids = response

		if self.options['reset_remote_dataset'] or len(dataset_ids['data']) == 0:
			for dataset in dataset_ids['data']:
				delete_dataset_id = dataset['dataset_id']

				options = {
					"database_name":"datasets",
					"dataset_id":delete_dataset_id,
					"dataset_name":dataset_name,
					"username":self.options['username'],
					"password":self.options['password'],
					"is_superuser":self.options['is_superuser']
				}
				response = queryRunner.queryRunner("database_access","delete_dataset",options).run()
				logger.debug(response)
				
			guides_dataset_definition =  {
				"name": dataset_name,
				"dataset_id":dataset_id,
				"description": "Guides dataset",
				"is_public": True,
				"groups": [],
				"owner": "ndavey",
				"details": {},
				"items": [],
				"is_superuser":self.options['is_superuser']		
			}

			options = {
				"database_name":"datasets",
				"dataset_name":dataset_name,
				"dataset_id":dataset_id,
				"id":"slimprints",
				"dataset_instance":guides_dataset_definition,
				"username":self.options['username'],
				"password":self.options['password'],
				"is_superuser":self.options['is_superuser']
			}
			
			response = queryRunner.queryRunner("database_access","add_dataset",options).run()
			logger.debug(response)
			dataset_id = response['data']['dataset_id']
		else:
			dataset_id = dataset_ids['data'][0]['dataset_id']

		return dataset_id

	def uploadGuides(self,protein_designs,mutation_designs):

		protein_designs = self.test_gRNAs["protein"] 
		mutation_designs = self.test_gRNAs["mutation"] 
	
		dataset_id_data = self.makeGuidesRemoteDataDataset()
		dataset_id_index = self.makeGuidesRemoteIndexDataset()

		logger.debug("Making Guides dataset -  dataset id: " + str(dataset_id_data))
		logger.debug("Making Guides dataset index -  dataset id: "  + str(dataset_id_index))
		
		items = []
		index_items = []
		
		try:
			for accession in protein_designs:
				protein_offsets = list(mutation_designs[accession].keys())
				protein_offsets.sort()

				from itertools import groupby
				from operator import itemgetter

				ranges = []
				regions = [i for i in range(0,len(protein_designs[accession]['sequence_peptides'])) if protein_designs[accession]['sequence_peptides'][i] != "-"]
				
				for k,g in groupby(enumerate(regions),lambda x:x[0]-x[1]):
					group = (map(itemgetter(1),g))
					group = list(map(int,group))
					ranges.append([group[0],group[-1]])

				index_item_details = {
					"mutations":0,
					"residues":[],
					"guides_counts":0,
					"regions":ranges
				}

				tags = ["gRNA","SLiMs"]

				for peptide_range in ranges:
					protein_offset_info = {}
					peptide_range_guides = 0
					for protein_offset in range(peptide_range[0],peptide_range[1]+2):
						
						protein_offset = str(protein_offset)

						if protein_offset in protein_offsets:
							guides = {}
							for editor in self.options['editors']:
								if editor in mutation_designs[accession][protein_offset]['guides']:
									for strand in ['negative','positive']:
										for guide in mutation_designs[accession][protein_offset]['guides'][editor][strand]:
											guides[guide] = mutation_designs[accession][protein_offset]['guides'][editor][strand][guide]
											guides[guide]['guide'] = guide 
											guides[guide]['strand'] = strand 
											guides[guide]['editor'] = editor 

							try:
								mutation_blosum_effect = mutation_designs[accession][protein_offset]['mutation_blosum_effect']
							except:
								mutation_blosum_effect = []

							protein_offset_info[protein_offset] = {
								"mutation_blosum_effect":mutation_blosum_effect,
								"mutations":mutation_designs[accession][protein_offset]['mutations'],
								"wildtype_residue":mutation_designs[accession][protein_offset]['wildtype_residue'],
								"guides_counts":len(guides),
								"guides":guides
								}
							
							peptide_range_guides += len(guides)
							index_item_details["mutations"] = index_item_details['mutations']+len(mutation_designs[accession][protein_offset]['mutations'])
							index_item_details["residues"] = index_item_details['residues'] + [protein_offset]
							index_item_details["guides_counts"] = index_item_details['guides_counts']+len(guides)
							

					details = {
						"peptide":"".join(protein_designs[accession]['sequence_peptides'][peptide_range[0]:peptide_range[1]+1]),
						"protein_offset_info":protein_offset_info
					}

					items.append(
						{
						"protein_acc": accession,
						"protein_start": peptide_range[0]+1,
						"protein_end": peptide_range[1]+1,
						"type": "gRNA",
						"details": details,
						"source": "gRNA",
						"owner": self.options['username'],
						"tags":tags
						}
					)

				index_items.append(
					{
					"protein_acc": accession,
					"type": "gRNA",
					"details": index_item_details,
					"source": "gRNA",
					"owner": self.options['username'],
					"tags":tags
					}
				)

				###############################################################################################
				###############################################################################################
				logger.debug("Uploading " + str(len(items)) + " items")
				if len(items) > 100:
					options = {
						"database_name":"datasets",
						"dataset_id":dataset_id_data,
						"dataset_instance_data":items,
						"username":self.options['username'],
						"password":self.options['password'],
						"debug":True
					}

					response = queryRunner.queryRunner("database_access","add_dataset_instances",options).run()
					items = []

				if len(index_items) > 100:
					options = {
						"database_name":"datasets",
						"dataset_id":dataset_id_index,
						"dataset_instance_data":index_items,
						"username":self.options['username'],
						"password":self.options['password'],
						"debug":True
					}
				
					response = queryRunner.queryRunner("database_access","add_dataset_instances",options).run()
					index_items = []
					
				###############################################################################################
				###############################################################################################

			###############################################################################################
			###############################################################################################

			if len(items) > 0:
				options = {
					"database_name":"datasets",
					"dataset_id":dataset_id_data,
					"dataset_instance_data":items,
					"username":self.options['username'],
					"password":self.options['password'],
					"debug":True
				}

				response = queryRunner.queryRunner("database_access","add_dataset_instances",options).run()

			if len(index_items) > 0:
				options = {
					"database_name":"datasets",
					"dataset_id":dataset_id_index,
					"dataset_instance_data":index_items,
					"username":self.options['username'],
					"password":self.options['password'],
					"debug":True
				}
			
				response = queryRunner.queryRunner("database_access","add_dataset_instances",options).run()
		
			###############################################################################################
			###############################################################################################
			# 	
		except:
			utilities_error.printError()

if __name__ == "__main__":
	pass
