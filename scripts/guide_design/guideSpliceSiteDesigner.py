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

class guideSpliceSiteDesigner():

	def __init__(self):
		self.options = {
			"splice_site_set_type":"unknown",
			"intergenic_guide_sequence_length":41,
			"intergenic_guide_sequence_samples":2000,
			"intergenic_guide_sequence_requirements":{
				"positions":{32:["G"]},
				"mismatch":{"start":11,"end":31}
		}}

		###-----------------_-------------_-_--_-_---_---_-------_------###


	def find_splice_sites(self,strand_type,peptide_genomic,peptide_cds_gapped):
		splice_site_offsets_set = []				
	
		boundary_counter = 0
		exon_count = len(re.split(r'\-+', peptide_cds_gapped)) - 2
		
		for offset in range(0,len(peptide_genomic)):
			splice_site = False
			splice_site_type = "-"
			splice_site_offsets = []

			if peptide_cds_gapped[offset] != "-" and peptide_cds_gapped[offset-1] == "-": 
				if strand_type == "negative": 
					boundary_counter_protein = exon_count-boundary_counter
				else:
					boundary_counter_protein = boundary_counter

				if boundary_counter_protein == 0:
					splice_site_type = "Kozak"
					splice_site_offsets = [offset-5,offset-4,offset-3,offset-2,offset-1]
					splice_site = peptide_genomic[offset-5:offset]
					splice_site_reverse = utilities_codons.get_complementry_oligo(peptide_genomic[offset-5:offset])[::-1]
					splice_site_match = [splice_site,splice_site_reverse]
					splice_site_match.sort()
				else:
					if strand_type == "positive": 
						splice_site_type = "3'"
					else:
						splice_site_type = "5'"

					splice_site = peptide_genomic[offset-2:offset]
					splice_site_offsets = [offset-2,offset-1]
					splice_site_reverse = utilities_codons.get_complementry_oligo(peptide_genomic[offset-2:offset])[::-1]
					splice_site_match = [splice_site,splice_site_reverse]
					splice_site_match.sort()

				splice_site_offsets = {
					"offsets":splice_site_offsets,
					"splice_site_type":splice_site_type,
					"splice_site_match":splice_site_match,
					"boundary_counter":boundary_counter_protein
				}

				splice_site_offsets_set.append(splice_site_offsets)
			
			if peptide_cds_gapped[offset-1] != "-" and peptide_cds_gapped[offset] == "-":
				boundary_counter += 1
				if strand_type == "negative": 
					boundary_counter_protein = exon_count-boundary_counter
				else:
					boundary_counter_protein = boundary_counter

				if boundary_counter_protein == 0:
					splice_site_type = "Kozak"
					splice_site = peptide_genomic[offset:offset+5]
					splice_site_offsets = [offset,offset+1,offset+2,offset+3,offset+4]
					splice_site_reverse = utilities_codons.get_complementry_oligo(peptide_genomic[offset:offset+5])[::-1]
					splice_site_match = [splice_site,splice_site_reverse]
					splice_site_match.sort()
				else:
					if strand_type == "positive": 
						splice_site_type = "5'"
					else:
						splice_site_type = "3'"

					splice_site = peptide_genomic[offset:offset+2]
					splice_site_offsets = [offset,offset+1]
					splice_site_reverse = utilities_codons.get_complementry_oligo(peptide_genomic[offset:offset+2])[::-1]
					splice_site_match = [splice_site,splice_site_reverse]
					splice_site_match.sort()

				splice_site_offsets = {
					"offsets":splice_site_offsets,
					"splice_site_type":splice_site_type,
					"splice_site_match":splice_site_match,
					"boundary_counter":boundary_counter_protein
				}

				splice_site_offsets_set.append(splice_site_offsets)

		return splice_site_offsets_set

	def make_splice_site_guides(self):
	
		self.splice_site_guides = {}
		self.splice_site_guides_protein_mapping = {}

		peptide_guide_candidates = {}
		peptide_guide_candidates_type = {
			"+|3'":[],
			"+|5'":[],
			"+|Kozak":[],
			"-|3'":[],
			"-|5'":[],
			"-|Kozak":[]
		} 

		validated_splice_site_guides = ["CTCACCATGGTGGCGGATGG","TTACTCACCATGGTGGCGGA","CGTACCATGATGTTGCCCAT","CCGTACCATGATGTTGCCCA","CTCACCTCTGGATTGGCTTC","GCCCACCTGCACCATGGGGT","CCTCACCTTGATGGGGATGC","ACCTCACCTGGATCTCCACA","CTCACCAAAGGACTCGTTGA","TCTTACCTCCCGCCTGCTGC","TACTTACACCCTGCTGCCAC","ATACTTACTGCAAAAAGAGC","GCTCACCTCTGAAGCTTCAA","AGTTTCACCTGTATAACTAG","GTTACCTGTCTGTAGTCTCT","CCTTTACCTGTTGGAGCTGA","CCCTTTACCTGTTGGAGCTG","CACACCTGAGCATGAGCCCT","GTCACCTGTAAGGCCTCTCC","CTCACCTTAAGGGCATCAGC","CTCACCATGGCTGTTGCGCG","CCGCTTACCTTCGTGAGTCT","ATCACCTGTTCTGTTGACTG","AGAAGTACCTCTCTCTGTTT","GTACTACCTTATATCCTTCT","TGTACTACCTTATATCCTTC","GGACTTACCTTTCTCATTGG","GAACAAACCTTATGCACATA","GATTTACCTTTCCTCTGTGT","CATACCTATGACATTTACAA","CCGTACCTGTGATTGTTCCA","GAATACCTCTAGCCATTCAG","ACTCACGAGGTGCTGGGAGG","CATACCTCTGAATGATCAAT","CACTCACTGTTGCACTCTGC","AACTCACCCCAAAGTCCTGA","GCTCACCTGGCAGAAGCTGC","ACGTACCTGTTGAAAAGTGC","TAAAATACCTCCATTGATTC","CTTACCTGGCTCTTCTTGGT","CCTCACCTTTGAAGAGAAAG","TCATTACCTGGTTGTCAAGC","CCTACCTGAAACACTGACCT","CACTCACCTCTTGGAGCTGG","ACTACCTTATGATAGAGTTC","TATACCCAGGCTATAGTCCC","CACTGACCTTGTCCAGCTGC","CACTCACATTCTTAAAGAAG","TACTCACCTTGAGTCAAACC","ACTTACCTGTCTGGGACTGA","ATTACCTTTGGTGAATATCT","CTTACAAGTTCAACCGGAGG","CTCACCAGCTCCTCTATTGT","TTTAATACCTGGAGAAGATA","CATACCTCTGTGACATTGAA","CCTACCTGTGGCTAGCAGAG","GTCTACCTCCAATTCCCGCA","ACTTACTGAATCCTCAGCCT","CCAACCTCGTCCACCCGCCC","AGCTCACCTTTATAAATAAG","TAGTTACCTTTATTTGACCA","CCAACCTGGCAAAACAGAAT"]
		
		rows = []
		counter = 0
		for accession in self.options['accessions_splice_site']:
			counter += 1
			logging.info("Finding splice sites for " + accession + " - " + str(counter) + "/" + str(len(self.options['accessions_splice_site'])))
			
			self.splice_site_guides_protein_mapping[accession] = {
				"+|3'":[],
				"+|5'":[],
				"+|Kozak":[],
				"-|3'":[],
				"-|5'":[],
				"-|Kozak":[]
			}

			uniprot_basic_data = queryRunner.queryRunner("uniprot","parse_basic",{"accession":accession}).run()
			gene_name = uniprot_basic_data['data']['gene_name']
			sequence = uniprot_basic_data['data']['sequence']

			hash = utilities_basic.params_to_hash({
				"editors":self.options['editors'],
				"edit_window":self.options['edit_window'],
				"target_nucleotides":self.options['target_nucleotides'],
				"editable_nucleotides":self.options['editable_nucleotides'],
				"include_splice_sites":self.options["include_splice_sites"],
				"only_splice_sites":self.options["only_splice_sites"],
				'skip_stops':self.options["skip_stops"],
				'only_stops':self.options["only_stops"]
			})

			options = {
				"accession":accession,
				"peptide_start":0,
				"peptide_end":len(sequence),
				"flanks":self.options['flanks'],
				"peptide_sequence":sequence,
				"species":self.options['species']
			}

			if 'transcript_selection' in self.options:
				if accession in self.options['transcript_selection']:
					options['transcript_id'] = self.options['transcript_selection'][accession]['transcript_id']
					options['gene_id'] = self.options['transcript_selection'][accession]['gene_id']

			genomes_guides_path = os.path.join(self.options['genomes_guides_library_regions_path'],accession + "." + hash + ".splice_site.gRNA.json")
			
			#("#",accession,genomes_guides_path,os.path.exists(genomes_guides_path),self.options['remake'])
			
			if not os.path.exists(genomes_guides_path) or self.options['remake']:
				response = queryRunner.queryRunner("ensembl","get_ensembl_data_for_region",options).run()
				
				if response['data']['status'] == "error":
					pprint.pprint(response)
					continue

				grna_list = []
				for editor in self.options['editors']:
					for strand_type in self.options['strand_types']:
						
						#peptide_genomic = response['data']['peptide_genomic']
						#peptide_cds_gapped = response['data']['peptide_cds_gapped']
						peptide_genomic = response['data']['peptide_genomic_flanks']
						peptide_cds_gapped = "-"*self.options['flanks'] + response['data']['peptide_cds_gapped'] + "-"*self.options['flanks']
					
						protein_codon_mapping = {}
						codon = []
						codon_offsets = []
						protein_offset = []
						protein_offset = 0
						protein_sequence = ""
						
						for offset in range(0,len(peptide_cds_gapped)):
							if peptide_cds_gapped[offset] != "-":
								codon.append(peptide_cds_gapped[offset])
								codon_offsets.append(offset)
						
								if len(codon) == 3:
									protein_offset += 1
									codon_aa = utilities_codons.translate_oligo("".join(codon))
									protein_sequence += codon_aa

									for codon_offset in codon_offsets:
										protein_codon_mapping[codon_offset] = {
											"protein_offset":protein_offset,
											"residue":codon_aa,
											"first_nuclotide":codon_offsets[0] == codon_offset,
											"codon_offsets":copy.deepcopy(codon_offsets),
											"codon_contiguous":int(sum(codon_offsets)/3)==codon_offsets[1],
											"exon":True
										}

									codon = []
									codon_offsets = []
							else:
								protein_codon_mapping[offset] = {
									"residue":"x",
									"split_codon":len(codon) > 0,
									"codon_offsets":copy.deepcopy(codon_offsets),
									"exon":False
								}
					
						if strand_type == "negative":
							peptide_genomic = utilities_codons.get_complementry_oligo(peptide_genomic)[::-1]
							peptide_cds_gapped = utilities_codons.get_complementry_oligo(peptide_cds_gapped)[::-1]

						splice_sites = self.find_splice_sites(strand_type,peptide_genomic,peptide_cds_gapped)
						
						for splice_site in splice_sites:
							splice_site_offsets = splice_site["offsets"]
							splice_site_type = splice_site["splice_site_type"]
							splice_site_match = splice_site["splice_site_match"]
							boundary_counter = splice_site["boundary_counter"]
							
							for splice_site_offset in splice_site_offsets:
								if splice_site_type == "Kozak":
									print(splice_site_offset)
								for edit_window_position in self.options['edit_window']:
									try:
										window_position_edited  = []
										guide_candidate_start = splice_site_offset-edit_window_position
										guide_candidate_end = guide_candidate_start + self.options['guide_length']
										
										guide_candidate = ""
										guide_candidate_extended = ""
										guide_candidate_gapped = ""
										guide_candidate_gapped_extended = ""
										guide_candidate_formatted = ""
										splice_sites = ""
										edits = ""
										edit_window = ""

										iter_position = 0
										splice_site_edited = False
										transcript_end = False

										for i in range(guide_candidate_start-self.options["guide_reporter_flank_5'"],guide_candidate_end - self.options['pam_length'] + self.options["guide_reporter_flank_3'"]):
											if i >= len(peptide_genomic):
												transcript_end = True

										if transcript_end:
											continue
										
										for i in range(guide_candidate_start-self.options["guide_reporter_flank_5'"],guide_candidate_end - self.options['pam_length'] + self.options["guide_reporter_flank_3'"]):
											guide_candidate_extended += peptide_genomic[i]
											guide_candidate_gapped_extended += peptide_cds_gapped[i]
									
										for i in range(guide_candidate_start,guide_candidate_end):
											iter_position += 1

											guide_candidate += peptide_genomic[i]
											guide_candidate_gapped += peptide_cds_gapped[i]

											if peptide_cds_gapped[i] == "-":
												if peptide_cds_gapped[i-1] != "-":
													guide_candidate_formatted += "|"

												guide_candidate_formatted += peptide_genomic[i].lower()
											else:
												if peptide_cds_gapped[i-1] == "-":
													guide_candidate_formatted += "|"
												guide_candidate_formatted += peptide_genomic[i]
											
											edited = False
										
											if peptide_genomic[i] not in self.options['target_nucleotides'][editor]: continue #N etc

											if self.options['target_nucleotides'][editor][peptide_genomic[i]] == peptide_genomic[i] and iter_position in self.options['edit_window']:
												edited = True
												edits += self.options['target_nucleotides'][editor][peptide_genomic[i]]
											else:
												edits += "-"

											if i in splice_site_offsets:
												if edited == False:
													splice_sites += "x"
												else:
													splice_site_edited = True
													splice_sites += self.options['target_nucleotides'][editor][peptide_genomic[i]]

													if iter_position not in window_position_edited:
														window_position_edited.append(iter_position)
											else:
												splice_sites += "-"

											if iter_position in self.options['edit_window']:
												edit_window += "e"
											else:
												edit_window += "-"

										#---------------------------------––------#

										guideDesignerObj = guideDesigner.guideDesigner()
										guideDesignerObj.options = self.options
										construct = guideDesignerObj.format_guide(guide_candidate[:-3],second_insert=guide_candidate_extended)
										
										sequence_check_data = utilities_oligos.check_oligo_sequence_features(construct.upper(),avoided_sequences=self.options['avoided_oligo_construct'])
										
										if sequence_check_data['check_status']['all'] == False:
											#logger.debug(construct + " discarded - oligo_sequence_features " + ",".join(sequence_check_data['check_status']['failed_check_type']))
											continue

										#---------------------------------––------#
										
										pam_sequence = guide_candidate[-3:]
										pam_type = "None"
										if pam_sequence[-2] == "G":
											pam_type = "NGN"
										
										if pam_sequence[-2:] == "GG":
											pam_type = "NGG"
										
										if strand_type == "positive":
											strand_type_formatted = "+"
										else:
											strand_type_formatted = "-"
									
										if pam_type in ["NGN","NGG"] and splice_site_edited:
											grna_list.append(guide_candidate)
											
											###----------------------------------------------------------------------###
											###                  PEPTIDE OVERLAPPING SPLICE SITE                     ### 
											###----------------------------------------------------------------------###
											last_type = False
											peptide_formatted = []
											protein_offsets = []
											
											for offset_iter in range(guide_candidate_start-self.options["guide_reporter_flank_5'"],guide_candidate_end - self.options['pam_length'] + self.options["guide_reporter_flank_3'"]):
												try:
													if strand_type == "negative":
														offset_iter = len(peptide_cds_gapped) - offset_iter

													if offset_iter not in protein_codon_mapping:
														logging.debug(["offset not in protein_codon_mapping",offset_iter,len(peptide_genomic)])
														continue
													#if guide_candidate[:-3] == "GGCTCACCTGTTTGTTGTGG":print(peptide_formatted)
													#if guide_candidate[:-3] == "GGCTCACCTGTTTGTTGTGG":print(protein_codon_mapping[offset_iter])

													if protein_codon_mapping[offset_iter]['exon'] == False:
														if last_type == "exon":
															if protein_codon_mapping[offset_iter]["split_codon"]:
																last_codon_offset = protein_codon_mapping[offset_iter]["codon_offsets"][0]
																peptide_formatted.append(protein_codon_mapping[last_codon_offset]['residue'].lower())

															peptide_formatted.append("|")
															peptide_formatted_splice_site_position = len(peptide_formatted)
															last_type = "intron"

														if last_type == False: last_type = "intron"
														
														if offset_iter%3 == 0:
															peptide_formatted.append("x")
													else:
														if 'residue' in protein_codon_mapping[offset_iter]:
															if last_type == "intron":
																protein_offsets.append(protein_codon_mapping[offset_iter]['protein_offset'])
																peptide_formatted.append("|")
																peptide_formatted_splice_site_position = len(peptide_formatted)

																if protein_codon_mapping[offset_iter]['codon_contiguous']:
																	peptide_formatted.append(protein_codon_mapping[offset_iter]['residue'])
																else:
																	peptide_formatted.append(protein_codon_mapping[offset_iter]['residue'].lower())
																
																last_type = "exon"

															elif protein_codon_mapping[offset_iter]['first_nuclotide'] and protein_codon_mapping[offset_iter]['protein_offset'] not in protein_offsets:
																peptide_formatted.append(protein_codon_mapping[offset_iter]['residue'])
															
															if last_type == False: last_type = "exon"
												except:
													logging.error(["Error: offset not in protein_codon_mapping",offset_iter,len(peptide_genomic)])


											peptide_formatted = peptide_formatted[peptide_formatted_splice_site_position-6:peptide_formatted_splice_site_position+5]
											if strand_type == "negative":
												peptide_formatted.reverse()
												splice_site_transcript_position = len(peptide_cds_gapped) - splice_site_offsets[1]
												splice_site_protein_position = int(len(peptide_cds_gapped[splice_site_offsets[1]:].replace("-",""))/3)
												boundary_counter_formatted =boundary_counter  
											else:
												splice_site_transcript_position = splice_site_offsets[1]
												boundary_counter_formatted = boundary_counter
												splice_site_protein_position = int(len(peptide_cds_gapped[:splice_site_transcript_position].replace("-",""))/3)
											
											peptide_formatted = "".join(peptide_formatted)
											
											###----------------------------------------------------------------------###
											###----------------------------------------------------------------------###
											###----------------------------------------------------------------------###
										
											row = "\t".join([
												str(len(peptide_guide_candidates)),
												guide_candidate[:-3],
												#str(splice_site_transcript_position),
												accession,
												gene_name,
												editor,
												#pam_type,
												pam_sequence,
												strand_type_formatted,
												splice_site_type,
												str(splice_site_transcript_position-self.options['flanks']),
												str(splice_site_protein_position),
												str(boundary_counter_formatted),
												",".join(splice_site_match),
												peptide_formatted,
												guide_candidate_formatted,
												splice_sites,
												#edit_window,
												edits,
												",".join([str(x) for x in window_position_edited]),
												str(guide_candidate[:-3] in validated_splice_site_guides),
												",".join(sequence_check_data['check_status']['failed_check_type'])
											])

											if row not in rows:
												logger.debug(row)
												rows.append(row)

											splice_site_strand_type = strand_type_formatted + '|' + splice_site_type

											use_guide = True
											if self.options['splice_site_set_type'] == "control":
												if "splice_sites_guide_sequence_requirement" in self.options:
													if "required_edits" in self.options['splice_sites_guide_sequence_requirement']:
														if splice_site_strand_type in self.options['splice_sites_guide_sequence_requirement']['required_edits']:
															if splice_sites.count(self.options['splice_sites_guide_sequence_requirement']['required_edits'][splice_site_strand_type]) == 0:
																use_guide = False
															
											if use_guide:
												if guide_candidate not in peptide_guide_candidates:
													peptide_guide_candidates[guide_candidate] = {}

												peptide_guide_candidates[guide_candidate][editor] = {
													"guide":guide_candidate[:-3],
													'accession':accession,
													'guide_candidate_extended':	guide_candidate_extended,
													'guide_candidate_gapped_extended':guide_candidate_gapped_extended,
													'guide_candidate_extended_flanks':guide_candidate_extended.split(guide_candidate[:-3]),
													'gene_name':gene_name,
													'splice_site_set_type':self.options["splice_site_set_type"],
													'editor':editor,
													'pam_type':pam_type,
													'pam_sequence':pam_sequence,
													'strand':strand_type_formatted,
													'splice_site_type':splice_site_type,
													'intron_count':boundary_counter,
													'splice_site_protein_position':splice_site_protein_position,
													'splice_site_transcript_position':splice_site_transcript_position - self.options['flanks'],
													'splice_site_match':splice_site_match,
													'guide_exon_boundary':guide_candidate_formatted,
													'guide_exon_gapped':guide_candidate_gapped,
													'splice_site_edited':splice_sites,
													'edit_window':edit_window,
													'all_site_edits':edits,
													'edit_window_position_splice_site_edits':window_position_edited,
													'validated':guide_candidate[:-3] in validated_splice_site_guides
												}

												if guide_candidate not in self.splice_site_guides_protein_mapping[accession][splice_site_strand_type]:
													self.splice_site_guides_protein_mapping[accession][splice_site_strand_type].append(guide_candidate)

												if guide_candidate not in peptide_guide_candidates_type[splice_site_strand_type]:
													peptide_guide_candidates_type[splice_site_strand_type].append(guide_candidate)
											else:
												logger.debug("Skipping: " + guide_candidate)
												
									except:
										if splice_site_offset + 20 < len(peptide_genomic):
											logging.error(["region error for splice sites guide definition",splice_site_offset,edit_window_position,len(peptide_genomic)])
										if splice_site_offset - 20 < 0:
											logging.error(["region error for splice sites guide definition",splice_site_offset,edit_window_position,len(peptide_genomic)])

										raise


				#response = queryRunner.queryRunner("ensembl","map_oligo_to_gene",{"oligo_sequence":grna_list,"annotate_all_exact_matches":True,"exact_matches_only":True,"remake":True}).run()

		##----------------------------------------##

		if self.options['splice_site_set_type'] == "control":
			random.seed(1)
			filtered_peptide_guide_candidates = {}
			keep_set = []

			if "control_splice_sites_guide_sequence_requirement" in self.options:
				if "splice_site_types" in self.options['control_splice_sites_guide_sequence_requirement']:
					for splice_site_type in self.options['control_splice_sites_guide_sequence_requirement']['splice_site_types']:
						try:
							random_keep_set = random.sample(peptide_guide_candidates_type[splice_site_type], self.options['control_splice_sites_guide_sequence_samples'])
							keep_set += random_keep_set
						except:
							pass
		
			for guide_candidate in keep_set:
				filtered_peptide_guide_candidates[guide_candidate] = peptide_guide_candidates[guide_candidate]

			peptide_guide_candidates = filtered_peptide_guide_candidates
		else:
			if 'target_splice_sites_guide_sequence_samples' in self.options:
				filtered_peptide_guide_candidates = {}
				keep_set = []

				for accession in self.splice_site_guides_protein_mapping:
					search_set = self.splice_site_guides_protein_mapping[accession]["-|5'"]
					random_keep_set = random.sample(search_set, min(len(search_set), self.options['target_splice_sites_guide_sequence_samples']))
					if len(random_keep_set) < self.options['target_splice_sites_guide_sequence_samples']:
						search_set = self.splice_site_guides_protein_mapping[accession]["-|3'"] + self.splice_site_guides_protein_mapping[accession]["+|3'"]
						random_keep_set += random.sample(search_set, min(len(search_set), self.options['target_splice_sites_guide_sequence_samples'] - len(random_keep_set)))

					keep_set += random_keep_set

				for guide_candidate in keep_set:
					filtered_peptide_guide_candidates[guide_candidate] = peptide_guide_candidates[guide_candidate]

				peptide_guide_candidates = filtered_peptide_guide_candidates

		del self.splice_site_guides_protein_mapping

		##----------------------------------------##
			
		headers = [
			'guide',
			'construct',
			'accession',
			'gene_name',
			'editor',
			#'pam_type',
			'pam_sequence',
			'strand',
			'splice_site_type',
			'splice_site_transcript_position',
			'splice_site_protein_position',
			'intron_count',
			'splice_site_match',
			'guide_exon_boundary',
			'guide_exon_boundary_peptide',
			'splice_site_edited',
			#'splice_site_edits',
			'all_site_edits',
			'edit_window_position_splice_site_edits',
			'validated'
		]

		#rows.sort()
		#print('\t'.join(headers) + '\n' + '\n'.join(rows))
		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + "." + self.options["splice_site_set_type"] + ".splice_sites.tdt")
		logger.debug("Writing: " + out_path )
		open(out_path,'w').write('\t'.join(headers) + '\n' + '\n'.join(rows))

		out_path = os.path.join(self.options['genomes_guides_libraries_path'],self.options['library_name'] + "." + self.options["splice_site_set_type"] + ".splice_sites.json")
		logger.debug("Writing: " + out_path )
		with open(out_path,'w') as outfile:
			json.dump(peptide_guide_candidates, outfile)
		
		return peptide_guide_candidates
	
		
		
		
		