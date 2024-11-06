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

import guideAnnotation
import guideDesigner

#-----
import logging
logger = logging.getLogger(__name__)
#-----

##------------------------------------------------------------------##
##------------------------------------------------------------------##
##------------------------------------------------------------------##

class guideProteinDesigner():

	def __init__(self):

		self.options = config_reader.load_configeration_options(sections=["general"])
		self.options.update({
			"splice_site_set_type":"unknown",
			"intergenic_guide_sequence_length":41,
			"intergenic_guide_sequence_samples":2000,
			"skip_stops":True,
			"only_stops":False,
			"intergenic_guide_sequence_requirements":{
				"positions":{32:["G"]},
				"mismatch":{"start":11,"end":31}
			},
			'strand_types':['positive','negative'],
			"mutant_residue_offsets":[]
		}
		)

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

		###-----------------_-------------_-_--_-_---_---_-------_------###

	def find_guides(self):

		hash = utilities_basic.params_to_hash({
			"editors":self.options['editors'],
			"edit_window":self.options['edit_window'],
			"target_nucleotides":self.options['target_nucleotides'],
			"editable_nucleotides":self.options['editable_nucleotides'],
			'skip_stops':self.options["skip_stops"],
			'only_stops':self.options["only_stops"]
		})

		genomes_guides_path = os.path.join(self.options['data_path'],"genomes","gRNA","gRNA_libraries",self.options['library_name'] + hash + ".guides.gRNA.json")
		  
		if not os.path.exists(genomes_guides_path) or self.options['remake'] or True:
			guide_designs = {}
			mutation_designs = {}
			peptide_guide_candidates = {}
			errors = {}

			guideDesignerObj = guideDesigner.guideDesigner()
			guideDesignerObj.options = self.options

			for editor in self.options['editors']:
				peptide_guide_candidates[editor] = {}

			grna_list = {}
			
			for accession in self.options['accession']:
				uniprot_basic_data = queryRunner.queryRunner("uniprot","parse_basic",{"accession":accession}).run()

				if 'gene_name' in  uniprot_basic_data['data']:
					gene_name = uniprot_basic_data['data']['gene_name']
				else:
					gene_name = accession

				sequence = uniprot_basic_data['data']['sequence']
				
				
				transcript_list = []
				if 'transcript_selection' in self.options:
					if accession in self.options['transcript_selection']:
						options = {}
						options['transcript_id'] = self.options['transcript_selection'][accession]['transcript_id']
						options['gene_id'] = self.options['transcript_selection'][accession]['gene_id']

						transcript_list.append(options) 
				else:
					options = {
						"accession":accession,
						"peptide_start":1,
						"peptide_end":len(sequence),
						"flanks":self.options['flanks'],
						"peptide_sequence":sequence,
						"species":self.options['species'],
						"detailed":True
					}
					response = queryRunner.queryRunner("ensembl","grab_ensembl_transcript_id_by_genename_peptide",options).run()
					transcript_list = response['data'].values()

				logger.info(str(len(transcript_list)) + " transcripts matching protein")
				for transcript in transcript_list:
					transcript_id = transcript['transcript_id']
					gene_id = transcript['gene_id']

					options = {
						"accession":accession,
						"peptide_start":1,
						"peptide_end":len(sequence),
						"flanks":self.options['flanks'],
						"peptide_sequence":sequence,
						"species":self.options['species'],
						'transcript_id':transcript_id,
						'gene_id':gene_id
					}

					logging.info("#"*40)
					logging.info("Finding guides for " + accession + " " + gene_name + " " + gene_id + " " + transcript_id)

					response = queryRunner.queryRunner("ensembl","get_ensembl_data_for_gene",options).run()
					
					if response['status'] == "error":
						continue
					
					if 'data' not in response:
						pprint.pprint(options)
						pprint.pprint(response)
						continue

					if response['data']['transcript_protein'] != sequence:
						logging.error("transcript_protein != sequence")

					for strand_type in self.options['strand_types']:
						try:
							peptide_genomic = response['data']['transcript_genomic']
							peptide_cds_gapped = response['data']['transcript_cds_genomic_alignment']

							
							if strand_type == "negative":
								peptide_genomic = utilities_codons.get_complementry_oligo(peptide_genomic)[::-1]
								peptide_cds_gapped = utilities_codons.get_complementry_oligo(peptide_cds_gapped)[::-1]

							for offset in range(0,len(peptide_cds_gapped)):
								
								protein_coding = False
								try:
									guide_candidate_start = offset
									guide_candidate_end = guide_candidate_start + self.options['guide_length']
									
									guide_candidate = ""
									guide_candidate_reporter = ""
									guide_candidate_gapped = ""
									guide_candidate_gapped_extended = ""
									guide_candidate_formatted = ""
								
									transcript_end = False
									
									for i in range(guide_candidate_start-self.options["guide_reporter_flank_5'"],guide_candidate_end+self.options["guide_reporter_flank_3'"] - self.options['pam_length']):
										if i >= len(peptide_genomic):
											transcript_end = True

									if transcript_end:
										continue
									
									for i in range(guide_candidate_start-self.options["guide_reporter_flank_5'"],guide_candidate_end+self.options["guide_reporter_flank_3'"] - self.options['pam_length']):
										guide_candidate_reporter += peptide_genomic[i]
										guide_candidate_gapped_extended += peptide_cds_gapped[i]

									iter_position = 0
									for i in range(guide_candidate_start,guide_candidate_end):
										iter_position += 1

										guide_candidate += peptide_genomic[i]
										guide_candidate_gapped += peptide_cds_gapped[i]

										if peptide_cds_gapped[i] != "-":
											if iter_position in self.options['edit_window']:
												protein_coding = True

										if peptide_cds_gapped[i] == "-":
											if peptide_cds_gapped[i-1] != "-":
												guide_candidate_formatted += "|"

											guide_candidate_formatted += peptide_genomic[i].lower()
										else:
											if peptide_cds_gapped[i-1] == "-":
												guide_candidate_formatted += "|"
											guide_candidate_formatted += peptide_genomic[i]
									
									if strand_type == "positive":
										strand_type_formatted = "+"
									else:
										strand_type_formatted = "-"
								
									if guide_candidate.count("-"):
										continue

									#-------------------------#
									if protein_coding == False:
										continue

									if len(re.findall(self.options['pam_motif'],guide_candidate[-3:])) == 0:
										continue

									construct = guideDesignerObj.format_guide(guide_candidate[:20],second_insert=guide_candidate_reporter + "XXXX")
									
									sequence_check_data = utilities_oligos.check_oligo_sequence_features(construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])
									
									if sequence_check_data['check_status']['all'] == False:
										logger.debug("Skipping avoided_oligo_sequences: " + guide_candidate_reporter + " " + guide_candidate + " " + ",".join(sequence_check_data['check_status']['failed_check_type']))
										
										for editor in self.options['editors']:
											peptide_guide_candidates[editor][guide_candidate] = "Failed sequence checks"

										continue
									
									sequence_check_data = utilities_oligos.check_oligo_sequence_features(construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])

									#-------------------------#

									grna_list[guide_candidate] = {
										"strand":strand_type,
										"strand_type_formatted":strand_type_formatted,
										"pam_sequence":guide_candidate[-3:],
										"guide_candidate":guide_candidate[:-3],
										"guide_candidate_start":guide_candidate_start,
										"guide_candidate_end":guide_candidate_end,
										"guide_candidate_reporter":guide_candidate_reporter,
										"guide_candidate_gapped_extended":guide_candidate_gapped_extended,
										"guide_candidate_formatted":guide_candidate_formatted
									}
								except:
									if offset + 20 < len(peptide_genomic):
										logging.error(["region error for guide definition",offset,len(peptide_genomic)])
									if offset - 20 < 0:
										logging.error(["region error for guide definition",offset,len(peptide_genomic)])

									utilities_error.printError()
						except:
							utilities_error.printError()

					if len(grna_list) == 0:
						logging.error("No guides for protein: " + accession + " " + gene_id + " " + transcript_id)
					else:
						logging.info(str(len(grna_list)) + " guides for protein: " + accession + " " + gene_id + " " + transcript_id)
						self.options['gene_id'] = gene_id
						self.options['transcript_id'] = transcript_id
					
			##-------------------------------------------------##
		
			logger.info(">>>>> guides discovery step complete <<<<<<<")
			logger.info("Unfiltered guide count for " + accession + ": " + str(len(grna_list)))
			
			##-------------------------------------------------##

			guideAnnotationObj = guideAnnotation.guideAnnotation()
			guideAnnotationObj.options.update(self.options)
			guideAnnotationObj.options['accession'] = accession
			guideAnnotationObj.options['editors'] = self.options['editors']
			guideAnnotationObj.options['gene_id'] = self.options['gene_id']
			guideAnnotationObj.options['transcript_id'] = self.options['transcript_id']
			guideAnnotationObj.options['gRNA'] = list(grna_list.keys())

			annotated_guides = guideAnnotationObj.annotate_guide()
			
			edit_positions_possible = []
			edit_positions_used = []

			for editor in annotated_guides:
				for guide_candidate in annotated_guides[editor]:
					try:
						changes = self.format_edit_changes(annotated_guides[editor][guide_candidate]['difference'])
						edit_positions_possible += changes["pos"]
						if len(self.options['mutant_residue_offsets']) > 0:
							if len(set(self.options['mutant_residue_offsets']).intersection(changes["pos"])) == 0:
								continue

						edit_positions_used += changes["pos"]

						if annotated_guides[editor][guide_candidate]['splice_site'] == False and len(annotated_guides[editor][guide_candidate]['changes']) == 0:	
							peptide_guide_candidates[editor][guide_candidate] = "No functional edits"
							continue

						if '*' in changes['mut']:
							if self.options['skip_stops']:
								logging.error("Skipping stop:" + str(changes['mut']))
								continue

						if self.options['only_stops']:	
							if '*' not in changes['mut']:
								continue

						##------------------------------##

						if guide_candidate not in guide_designs:
							guide_designs[guide_candidate] = {} 
						
						if editor not in guide_designs[guide_candidate]:
							guide_designs[guide_candidate][editor] = {}

						##------------------------------##

						if accession not in mutation_designs:
							mutation_designs[accession] = {}

						for pos in changes['pos']:
							if pos not in mutation_designs[accession]:
								try:
									if pos > len(sequence):
										mutation_designs[accession][pos] = {
											"wildtype_residue":">",
											"guides":{},
											"mutations":[],
											"mutation_blosum_effect":[]
										}
									else:
										mutation_designs[accession][pos] = {
											"wildtype_residue":sequence[pos-1],
											"guides":{},
											"mutations":[],
											"mutation_blosum_effect":[]
										}
								except:
									logging.error(str(pos) + " not in sequence (" + str(len(sequence)) + ")")
								
							mutation_designs[accession][pos]['guides'][guide_candidate] = {
								"pos":changes["pos"],
								"mutations":changes["change"]
							}

							mutation_designs[accession][pos]['mutations'] += changes["change"]
							
						##------------------------------##
						
						guide_designs[guide_candidate][editor] = {
							"guide_candidate_reporter":grna_list[guide_candidate]['guide_candidate_reporter'],
							"peptide_genomic_flanks":"-",
							"peptide_genomic_flanks_mutated":annotated_guides[editor][guide_candidate]['mutated_genomic_region'],
							"peptide_sequence":annotated_guides[editor][guide_candidate]['wildtype_region'],
							"mutated_oligo":annotated_guides[editor][guide_candidate]['mutated_genomic_region'],
							"mutated_oligo_sub_edits":annotated_guides[editor][guide_candidate]['mutated_genomic_region_sub_edits'],
							"peptide_sequence_translated":annotated_guides[editor][guide_candidate]['wildtype_region'],
							"mutations":changes["change"],
							"codon_wt":changes["codon_wt"],
							"codon_mut":changes["codon_mut"],
							"pos":changes["pos"],
							"wt":changes["wt"],
							"mut":changes["mut"],
							"edit_window": annotated_guides[editor][guide_candidate]["edit_window"],
							"edit_window_mutated": annotated_guides[editor][guide_candidate]["edit_window_mutated"],
							"splice_site_info":annotated_guides[editor][guide_candidate]['splice_site_info'],
							"accession":accession,
							"peptide_start":changes["peptide_start"],
							"peptide_end":changes["peptide_end"],
							"strand_type":grna_list[guide_candidate]['strand'],
							"pam_sequence":grna_list[guide_candidate]['pam_sequence'],
							"pam_type":grna_list[guide_candidate]['pam_sequence'],
							"wildtype_region":annotated_guides[editor][guide_candidate]['wildtype_region'],
							"mutated_region":annotated_guides[editor][guide_candidate]['mutated_region'],
							"difference_region":annotated_guides[editor][guide_candidate]['difference'],
							"difference_region_edits":annotated_guides[editor][guide_candidate]['edits'],
							"changes_region":annotated_guides[editor][guide_candidate]['changes'],
							"stop_region":annotated_guides[editor][guide_candidate]['stop'],
							"exon_boundary_overlap":annotated_guides[editor][guide_candidate]['splice_site']
						}
					except:
						logger.debug(editor + " " + guide_candidate + " " + str(utilities_error.getError()))

			##-------------------------------------##
					
			alphafold_response = queryRunner.queryRunner("alphafold","get_alphafold_structural_classification",{"accession":accession}).run()
	
			##-------------------------------------##

			try:
				logging.info("".join(alphafold_response['data']['structural_module_class_classification']))
			except:
				pass

			logging.info("".join([
				sequence[o]
				if o+1 in list(set(edit_positions_possible)) else sequence[o].lower()
				for o in range(0,len(sequence)) 
			]) + " edits possible - " + str(len(list(set(edit_positions_possible)))))

			logging.info("".join([
				sequence[o]
				if o+1 in list(set(self.options['mutant_residue_offsets'])) else "x"
				for o in range(0,len(sequence)) 
			]) + " edits designed - " + str(len(list(set(self.options['mutant_residue_offsets'])))))

			logging.info("".join([
				sequence[o] 
				if o+1 in list(set(edit_positions_used)) else "x"
				for o in range(0,len(sequence)) 
			]) + " edits made - " + str(len(list(set(edit_positions_used)))))
		
			##-------------------------------------##
			
			response = {
				"guide_designs":guide_designs,
				"mutation_designs":mutation_designs,
				"peptide_guide_candidates":peptide_guide_candidates,
				"errors":errors
			}
	
			logger.debug("Writing " + genomes_guides_path)
			utilities_basic.write_to_json(genomes_guides_path,response,zipped=False)
		else:
			logger.debug("Reading " + genomes_guides_path)

			utilities_basic.read_from_json(genomes_guides_path,zipped=True)
			response = json.loads(open(genomes_guides_path, 'r').read())

		return response
	
	###-------------------------------###
		
	def format_edit_changes(self,edit_changes):
		changes = {
			"peptide_start":"-",
			"peptide_end":"-",
			"change":[],
			"codon_wt":[],
			"codon_mut":[],
			"pos":[],
			"wt":[],
			"mut":[]
		}

		for edit_pos in edit_changes:
			changes["pos"].append(edit_pos)
			changes["change"].append(edit_changes[edit_pos]["change"])
			changes["codon_wt"].append(edit_changes[edit_pos]["wildtype_codon"])
			changes["codon_mut"].append(edit_changes[edit_pos]["mutation_codon"])
			changes["mut"].append(edit_changes[edit_pos]["mutation"])
			changes["wt"].append(edit_changes[edit_pos]["wildtype"])
		
		if len(edit_changes) > 0:
			changes["peptide_start"] = min(edit_changes.keys())
			changes["peptide_end"] = max(edit_changes.keys())

		return changes
		
if __name__ == "__main__":
	guideProteinDesignerObj = guideProteinDesigner()
	guideProteinDesignerObj.options['accession'] = ["O60216"]#,"Q14683","Q8N3U4","Q9UQE7"]
	#guideProteinDesignerObj.options['mutant_residue_offsets'] = [2,3,4,5,6,7,8,9,52,53,54,55,56,57,58,59]#,"Q14683","Q8N3U4","Q9UQE7"]
	guideProteinDesignerObj.options['editors'] = ["Dual"]
	guideProteinDesignerObj.options['edit_window'] = [3,4,5,6,7,8,9]
	guideProteinDesignerObj.options['flanks'] = 20
	guideProteinDesignerObj.options['guide_length'] = 23
	guideProteinDesignerObj.options['species'] = "9606"
	guideProteinDesignerObj.options['library_name'] = "test_protein"
	guideProteinDesignerObj.options["guide_reporter_flank_5'"] = 11
	guideProteinDesignerObj.options["guide_reporter_flank_3'"] = 6
	guideProteinDesignerObj.options['pam_length'] = 3
	guideProteinDesignerObj.options['pam_motif'] = ".G."

	guideProteinDesignerObj.find_guides()
	