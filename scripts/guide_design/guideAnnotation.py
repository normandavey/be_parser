import os, sys, inspect, json, pprint, random, subprocess, time, copy,re
				
file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../"))
import config_reader

sys.path.append(os.path.join(file_path,"../utilities"))
import utilities_error
import utilities_codons
import utilities_basic

sys.path.append(os.path.join(file_path,"../data_management/"))
import queryRunner

#-----
import guideSpliceSiteDesigner

#-----

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#-----

##------------------------------------------------------------------##
##------------------------------------------------------------------##
##------------------------------------------------------------------##

class guideAnnotation():

	def __init__(self):

		self.options = config_reader.load_configeration_options(sections=['general'])
		self.options.update({
			'guide_length':23,
			'flanks':20,
			'guides_efficiency_flanks':50,
			'behive_flank_5':20,
			'behive_flank_3':10,
			'behive_base_editor':'',
			'behive_celltype':'',
			'pam_length':3,
			'editors':["ABE"],
			'edit_window':[4,5,6,7,8],
			'target_nucleotides':{
				"ABE":{"G":"g","T":"t","C":"c","A":"A"},
				"CBE":{"G":"g","T":"t","C":"C","A":"a"},
				"Dual":{"G":"g","T":"t","C":"C","A":"A"}
			},
			'editable_nucleotides':{
				"ABE":{"C":"c","A":"G","T":"t","G":"g"},
				"CBE":{"C":"T","A":"a","T":"t","G":"g"},
				"Dual":{"C":"T","A":"G","T":"t","G":"g"}
			}
		})

		self.options['genomes_path'] = os.path.join(self.options['data_path'],"genomes")
		self.options['genomes_guides_path'] = os.path.join(self.options['data_path'],"genomes","gRNA")
		self.options['genomes_guides_mutation_annotation_path'] = os.path.join(self.options['data_path'],"genomes","gRNA","gRNA_mutation_annotation")

		if not os.path.exists(self.options['genomes_guides_path']):
			os.mkdir(self.options['genomes_guides_path'])

		if not os.path.exists(self.options['genomes_guides_path']):
			os.mkdir(self.options['genomes_guides_path'])

		if not os.path.exists(self.options['genomes_guides_mutation_annotation_path']):
			os.mkdir(self.options['genomes_guides_mutation_annotation_path'])

	###########################################
	
	def calculate_gRNA_edits(self,gRNA,editor):
		"""
		Calculate the gRNA edits based on the given gRNA sequence and editor.
		@param gRNA - the gRNA sequence
		@param editor - the editor used for editing
		@return editor_edit_dict - a dictionary containing information about the edits made
		"""
		tmp_seq_multi_mutation = list(gRNA.lower())

		editor_edit_dict = {}
		edit_dict = {}
		mutated_genomic_region_sub_edits = []
		edit_check = False

		for edit_window_position in self.options['edit_window']:
			edit_window_position = edit_window_position - 1
			tmp_seq = list(gRNA.lower())
			wt_nucleotide = gRNA[edit_window_position]
			if wt_nucleotide not in self.options['editable_nucleotides'][editor]: continue

			tmp_seq[edit_window_position] =  self.options['editable_nucleotides'][editor][tmp_seq[edit_window_position].upper()]
			tmp_seq_multi_mutation[edit_window_position] =  self.options['editable_nucleotides'][editor][tmp_seq_multi_mutation[edit_window_position].upper()]

			if gRNA[edit_window_position].upper() != self.options['editable_nucleotides'][editor][tmp_seq_multi_mutation[edit_window_position].upper()].upper() :
				edit_dict[edit_window_position+1] = {
					"wt":gRNA[edit_window_position],
					"mut":self.options['editable_nucleotides'][editor][tmp_seq_multi_mutation[edit_window_position].upper()]
				}

				mutated_genomic_region_sub_edits.append("".join(tmp_seq).upper())
				edit_check = True

		editor_edit_dict['edits'] = edit_dict
		editor_edit_dict['edit_counts'] = len(edit_dict)
		editor_edit_dict['genome_edited'] = edit_check
		editor_edit_dict['edit_window'] = gRNA[min(self.options['edit_window'])-1:max(self.options['edit_window'])]
		editor_edit_dict['edit_window_mutated'] = "".join(tmp_seq_multi_mutation[min(self.options['edit_window'])-1:max(self.options['edit_window'])])
		editor_edit_dict['mutated_genomic_region'] = "".join(tmp_seq_multi_mutation).upper()
		editor_edit_dict['mutated_genomic_region_sub_edits'] = mutated_genomic_region_sub_edits
		
		return editor_edit_dict
		
	###########################################
	
	def calculate_gRNA_edits_behive(self,gRNA,gRNA_edited):
    
		editor_edit_dict = {}
		edit_dict = {}
		edit_check = False

		for i in range(0,len(gRNA)):
			if gRNA[i] != gRNA_edited[i]:
				edit_dict[i - self.options['behive_flank_5']] = {
					"wt":gRNA[i],
					"mut":gRNA_edited[i]
				}

				edit_check = True

		editor_edit_dict['edits'] = edit_dict
		editor_edit_dict['edit_counts'] = len(edit_dict)
		editor_edit_dict['genome_edited'] = edit_check
		editor_edit_dict['edit_window'] = gRNA
		editor_edit_dict['edit_window_mutated'] = gRNA_edited
		editor_edit_dict['mutated_genomic_region'] = gRNA_edited
		editor_edit_dict['mutated_genomic_region_sub_edits'] = [gRNA_edited]
		
		return editor_edit_dict
		
	###########################################
	
	def find_guide_position(self,gRNA,gRNA_edits,peptide_genomic_flanks):
		"""
		This method finds the position of a guide RNA in a given genomic context.
		@param gRNA - the guide RNA sequence
		@param gRNA_edits - dictionary containing information about the guide RNA edits
		@param peptide_genomic_flanks - the genomic context where the guide RNA is located
		@return None
		"""
		gRNA_complementry = utilities_codons.get_complementry_oligo(gRNA)[::-1]

		gRNA_start = -1
		if peptide_genomic_flanks.count(gRNA) > 0:
			gRNA_edits['strand'] = 'positive'
			gRNA_start = len(peptide_genomic_flanks.split(gRNA)[0])
		
		if peptide_genomic_flanks.count(gRNA) > 1:
			gRNA_start = len(peptide_genomic_flanks.split(gRNA)[0])
			logger.error(gRNA + " multiple matches " + str(peptide_genomic_flanks.count(gRNA)) + "@" + str(gRNA_start))

		if peptide_genomic_flanks.count(gRNA_complementry) > 0:
			gRNA_edits['strand'] = 'negative'
			gRNA_start = len(peptide_genomic_flanks.split(gRNA_complementry)[0])
	
		if peptide_genomic_flanks.count(gRNA_complementry) > 1:
			gRNA_start = len(peptide_genomic_flanks.split(gRNA_complementry)[0])
			logger.error(gRNA + " multiple matches " + str(peptide_genomic_flanks.count(gRNA_complementry)) + "@" + str(gRNA_start))
			
		if peptide_genomic_flanks.count(gRNA_complementry) == 0 and peptide_genomic_flanks.count(gRNA) == 0:
			logger.error(gRNA + " not in peptide_genomic_flanks")

		gRNA_edits['start'] = gRNA_start
		return gRNA_edits
	
	###########################################
	
	def find_splice_sites(self,gRNA_edits,splice_sites):
		"""
		Find splice sites in the given gRNA edits based on the provided splice sites information.
		@param gRNA_edits - the gRNA edits to be checked for splice sites
		@param splice_sites - the list of splice sites to compare against
		@return gRNA_edits with updated splice site information
		"""
		splice_site_check = False
		splice_site_info = {}
		for splice_site in splice_sites:
			if splice_site['boundary_counter'] == 0: continue
			for edit_window_position in gRNA_edits['edits']:
				if (edit_window_position + gRNA_edits['start']) in splice_site['offsets']:
					splice_site_check = True
					splice_site_info = splice_site
		
		gRNA_edits['splice_site'] = splice_site_check
		gRNA_edits['splice_site_info'] = splice_site_info

		return gRNA_edits
	
	###########################################
	
	def get_mutated_dna(self,gRNA,gRNA_edits,peptide_genomic_flanks,editor):
		"""
		This method takes a gRNA, gRNA edits, peptide genomic flanks, and an editor to mutate the DNA sequence based on the provided information.
		@param gRNA - the guide RNA sequence
		@param gRNA_edits - the edited guide RNA sequence
		@param peptide_genomic_flanks - the genomic flanks of the peptide
		@param editor - the editor for making the edits
		@return The mutated peptide genomic flanks
		"""
		gRNA_complementry = utilities_codons.get_complementry_oligo(gRNA)[::-1]

		mutated_peptide_genomic_flanks = False
		if peptide_genomic_flanks.count(gRNA) == 0 and peptide_genomic_flanks.count(gRNA_complementry) == 0:
			logger.error("Guide not found: " + gRNA + " " + str(peptide_genomic_flanks.count(gRNA_complementry)) + " " + str(self.options['accession']))

		if peptide_genomic_flanks.count(gRNA) > 0:
			mutated_peptide_genomic_flanks = peptide_genomic_flanks.replace(gRNA,gRNA_edits['mutated_genomic_region'])
			
		if peptide_genomic_flanks.count(gRNA_complementry) > 0:
			mutated_gRNA_complementry = utilities_codons.get_complementry_oligo(gRNA_edits['mutated_genomic_region'])[::-1]
			mutated_peptide_genomic_flanks = peptide_genomic_flanks.replace(gRNA_complementry,mutated_gRNA_complementry)
			
		return mutated_peptide_genomic_flanks
	
	###########################################
	
	def get_mutated_cds(self,peptide_cds_gapped_flanks,mutated_peptide_genomic_flanks):
		"""
		Generate mutated coding sequences based on the input peptide coding sequences with gapped flanks and mutated peptide genomic flanks.
		@param self - the object instance
		@param peptide_cds_gapped_flanks - the peptide coding sequences with gapped flanks
		@param mutated_peptide_genomic_flanks - the mutated peptide genomic flanks
		@return The mutated peptide coding sequences with gapped flanks.
		"""
		mutated_peptide_cds_gapped_flanks = []
		
		for i in range(0,len(peptide_cds_gapped_flanks)):
			if peptide_cds_gapped_flanks[i] == "-":
				mutated_peptide_cds_gapped_flanks.append("-")
			else:
				mutated_peptide_cds_gapped_flanks.append(mutated_peptide_genomic_flanks[i])
			
		return "".join(mutated_peptide_cds_gapped_flanks)
	
	###########################################
	
	def get_mutated_protein(self,gRNA,gRNA_edits,peptide_genomic_flanks,peptide_cds_gapped_flanks,editor):
		"""
		This method takes a gRNA, gRNA edits, peptide genomic flanks, peptide CDS gapped flanks, and an editor to mutate a protein.
		@param gRNA - the gRNA sequence
		@param gRNA_edits - the edits made to the gRNA
		@param peptide_genomic_flanks - the genomic flanks of the peptide
		@param peptide_cds_gapped_flanks - the CDS gapped flanks of the peptide
		@param editor - the editor used for mutation
		@return None
		"""
		mutated_peptide_genomic_flanks = self.get_mutated_dna(gRNA,gRNA_edits,peptide_genomic_flanks,editor)

		mutated_peptide_cds_gapped_flanks = self.get_mutated_cds(peptide_cds_gapped_flanks,mutated_peptide_genomic_flanks)
				
		wildtype_oligo = peptide_cds_gapped_flanks.replace("-","")
		mutation_oligo = mutated_peptide_cds_gapped_flanks.replace("-","")
		
		wildtype_oligo_translated = utilities_codons.translate_oligo(wildtype_oligo)
		mutation_oligo_translated = utilities_codons.translate_oligo(mutation_oligo)

		difference = {}
		changes = []
		stop = False
		proteome_edited = False

		for i in range(0,len(wildtype_oligo_translated)):
			if wildtype_oligo_translated[i] != mutation_oligo_translated[i]:
				change = wildtype_oligo_translated[i] + ":" + str(i + 1) + "->" + mutation_oligo_translated[i]
				difference[i + 1] = {
					"wildtype":wildtype_oligo_translated[i], 
					"mutation":mutation_oligo_translated[i],
					"wildtype_codon":wildtype_oligo[i*3:(i+1)*3],
					"mutation_codon":mutation_oligo[i*3:(i+1)*3],
					"change":change
				}
				changes.append(change)
				proteome_edited = True

			
			if mutation_oligo_translated[i] == "*" and wildtype_oligo_translated[i] != "*":
				stop = True
				#logger.debug("Stop found:" + wildtype_oligo_translated + " " + str(i) + "->" + mutation_oligo_translated)

		wildtype_region = "" 
		mutation_region = ""
		
		if len(difference) > 0:
			region_centre = int(sum(list(difference.keys()))/len(difference))
			for i in range(max(0,region_centre-6),min(region_centre+6,len(wildtype_oligo_translated))):
				try:
					if i + 1 in difference:
						wildtype_region += wildtype_oligo_translated[i]
						mutation_region += mutation_oligo_translated[i]
					else:
						wildtype_region += wildtype_oligo_translated[i].lower()
						mutation_region += mutation_oligo_translated[i].lower()
				except:
					logger.error("get_mutated_protein differences: " + str(difference))
					
		gRNA_edits.update({
			"proteome_edited":proteome_edited,
			"wildtype_region":wildtype_region,
			"mutated_region":mutation_region,
			"difference":difference,
			"changes":changes,
			"stop":stop
		})

		return gRNA_edits
	
	###########################################

	def annotate_guide(self):
		"""
		Annotate the guide using UniProt data and create a guide annotation table.
		@param self - the object itself
		@return annotated guide with gRNA edits and annotations
		"""
  
		uniprot_basic_data = queryRunner.queryRunner("uniprot","parse_basic",{"accession":self.options['accession']}).run()
		if len(uniprot_basic_data['data']) == 0:
			uniprot_basic_data = queryRunner.queryRunner("uniprot","check_uniprot_entry",{"accession":self.options['accession']}).run()
			logging.error("skipping no data " + self.options['accession'])
			logging.error(uniprot_basic_data['data'])
			return {}

		if len(uniprot_basic_data['data']['sequence']) > 10000: 
			logging.info("skipping long proteins " + self.options['accession'] + " " + str(len(uniprot_basic_data['data']['sequence'])))
			return {}  # Skip long proteins

		options = {
			"accession":self.options['accession'],
			"peptide_start":0,
			"peptide_end":len(uniprot_basic_data['data']['sequence']),
			"flanks":self.options['flanks'],
			"peptide_sequence":uniprot_basic_data['data']['sequence'],
			"species":self.options['species'],
			"match_oligos":self.options['gRNA'],
			"debug":True
		}

		if 'transcript_selection' in self.options:
			if self.options['accession'] in self.options['transcript_selection']:
				options['transcript_id'] = self.options['transcript_selection'][self.options['accession']]['transcript_id']
				options['gene_id'] = self.options['transcript_selection'][self.options['accession']]['gene_id']
		
		if 'gene_id' in self.options and 'transcript_id' in self.options:
			options['gene_id'] = self.options['gene_id']
			options['transcript_id'] = self.options['transcript_id']

		annotated_gRNAs = {}
		
		try:
			response = queryRunner.queryRunner("ensembl","get_ensembl_data_for_gene",options).run()

			if "error_type" in response:
				pprint.pprint(response)
			else:				
				transcript_cds_genomic_alignment =  response['data']['transcript_cds_genomic_alignment'] 
				transcript_genomic = response['data']['transcript_genomic']
				logging.info("Annotating " + self.options['accession'] + " " + response['data']['ensembl_gene_id'] + " gRNAs:" + str(len(self.options['gRNA'])))
		

				guideSpliceSiteDesignerObj = guideSpliceSiteDesigner.guideSpliceSiteDesigner()
				guideSpliceSiteDesignerObj.options = self.options
				splice_sites = guideSpliceSiteDesignerObj.find_splice_sites("positive",transcript_genomic,transcript_cds_genomic_alignment)
				
				for editor in self.options["editors"]:
					annotated_gRNAs[editor] = {}
					
					for gRNA in self.options['gRNA']:
						try:
							gRNA_edits = self.calculate_gRNA_edits(gRNA,editor)
							gRNA_edits = self.find_guide_position(gRNA,gRNA_edits,transcript_genomic)
							gRNA_edits = self.find_splice_sites(gRNA_edits,splice_sites)
							gRNA_edits = self.get_mutated_protein(gRNA,gRNA_edits,transcript_genomic,transcript_cds_genomic_alignment,editor)

							gRNA_edits.update({"accession":self.options['accession']})

							annotated_gRNAs[editor][gRNA] = gRNA_edits
						except:
							logger.error(self.options['accession'] + " " + gRNA + " " + str(len(transcript_genomic)) + " " + str(len(transcript_cds_genomic_alignment)))
					
				self.make_guide_annotation_table(annotated_gRNAs)
			
				self.options['genomes_guides_mutation_annotation_path']

			return annotated_gRNAs
		except:
			logger.error(utilities_error.getError())

	###########################################
	
	def annotate_protein_guides(self):
		"""
		Annotate protein guides based on the provided options and protein gRNAs.
		@param self - the object instance
		@return annotated_protein_gRNAs - a dictionary of annotated protein gRNAs
		"""
		annotated_protein_gRNAs = {}
		counter = 0
		protein_centric = True
		if "protein_centric" in self.options:
			protein_centric = self.options['protein_centric']
		
		guide_counter = 0

		accessions = list(self.options['protein_gRNAs'].keys())
		
		for accession in accessions:

			logger.debug("#" + accession)

			counter += 1
			try:
				self.options['accession'] = accession
				self.options['gRNA'] = self.options['protein_gRNAs'][accession]

				logger.debug("#" + accession + "\t" + str(len(self.options['protein_gRNAs'])) + "\t" + str(counter))
				guide_counter += len(self.options['protein_gRNAs'][accession])

				if protein_centric:
					annotated_protein_gRNAs[accession] = self.annotate_guide()
				else:
					annotation_response = self.annotate_guide()
					
					for editor in annotation_response:
						if editor not in annotated_protein_gRNAs:
							annotated_protein_gRNAs[editor] = {}

						annotated_protein_gRNAs[editor].update(annotation_response[editor])
						logger.debug("annotate_protein_guides: " + str([len(self.options['gRNA']),len(annotation_response[editor]),len(annotated_protein_gRNAs[editor]),guide_counter]))
			except:
				logger.error(accession)
				logger.error(utilities_error.getError())
		
		return annotated_protein_gRNAs

	###########################################
	
	def annotate_guide_conservation(self,orth_db="metazoa", annotated_gRNAs=None):
		"""
		Annotate the conservation information for the guide RNAs based on the provided orthology database.
		@param self - the instance of the class
		@param orth_db - the orthology database to use (default is "metazoa")
		@param annotated_gRNAs - a dictionary of annotated guide RNAs (default is None)
		@return A dictionary of annotated guide conservation information
		"""
		annotated_guide_protein_centric = {}
		annotated_guide_conservation = {}

		if annotated_gRNAs == None:
			annotated_gRNAs =  self.annotate_protein_guides()
		
		for editor in annotated_gRNAs:
			for gRNA in annotated_gRNAs[editor]:
				if 'accession' in  annotated_gRNAs[editor][gRNA]:
					accession = annotated_gRNAs[editor][gRNA]['accession']
					if accession not in annotated_guide_protein_centric: annotated_guide_protein_centric[accession] = {}
					if editor not in annotated_guide_protein_centric[accession]: annotated_guide_protein_centric[accession][editor] = {}
					annotated_guide_protein_centric[accession][editor][gRNA] = [int(o) for o in annotated_gRNAs[editor][gRNA]['difference'].keys()]
				else:
					print(annotated_gRNAs[editor][gRNA])
		
		for accession in annotated_guide_protein_centric:
			data_access_response = queryRunner.queryRunner("database_access","get_attributes",{"accession":accession,"database_name":"database"}).run()
			for editor in annotated_guide_protein_centric[accession]:
				if editor not in annotated_guide_conservation: annotated_guide_conservation[editor] = {}
				for gRNA in  annotated_guide_protein_centric[accession][editor]:
					gRNA_annotation = {}
					for uniprot_offset in annotated_guide_protein_centric[accession][editor][gRNA]:
						try:
							gRNA_annotation[uniprot_offset] = {
								"conservation_rlc":data_access_response['data']['Conservation'][orth_db][str(uniprot_offset)],
								"conservation_raw":data_access_response['data']['Conservation_raw'][orth_db][str(uniprot_offset)]
							}
						except:
							gRNA_annotation[uniprot_offset] = {
								"conservation_rlc":-1,
								"conservation_raw":-1
							}
					
					annotated_guide_conservation[editor][gRNA] = gRNA_annotation

		return annotated_guide_conservation
	
	###########################################
 
	def annotate_guide_features(self,annotated_gRNAs=None):
		"""
		Annotate guide features based on annotated gRNAs and protein-centric information.
		@param self - the object instance
		@param annotated_gRNAs - a dictionary of annotated gRNAs (default is None)
		@return annotated_guide_features - a dictionary of annotated guide features
		"""
		
		annotated_guide_protein_centric = {}
		annotated_guide_features = {}

		if annotated_gRNAs == None:
			annotated_gRNAs =  self.annotate_protein_guides()
		
		for editor in annotated_gRNAs:
			for gRNA in annotated_gRNAs[editor]:
				if 'accession' in  annotated_gRNAs[editor][gRNA]:
					accession = annotated_gRNAs[editor][gRNA]['accession']
					if accession not in annotated_guide_protein_centric: annotated_guide_protein_centric[accession] = {}
					if editor not in annotated_guide_protein_centric[accession]: annotated_guide_protein_centric[accession][editor] = {}
					annotated_guide_protein_centric[accession][editor][gRNA] = [int(o) for o in annotated_gRNAs[editor][gRNA]['difference'].keys()]
				else:
					print(annotated_gRNAs[editor][gRNA])
		
		for accession in annotated_guide_protein_centric:
			data_access_response = queryRunner.queryRunner("database_access","get_features",{"accession":accession,"residue_centric":True,"database_name":"database"}).run()
			

			for editor in annotated_guide_protein_centric[accession]:
				if editor not in annotated_guide_features: annotated_guide_features[editor] = {}

				if data_access_response['status'] == 'Error':
					logger.debug("Skipping: database_access - features " + accession)
					continue

				for gRNA in  annotated_guide_protein_centric[accession][editor]:
					gRNA_annotation = {}
					for uniprot_offset in annotated_guide_protein_centric[accession][editor][gRNA]:
						if uniprot_offset in data_access_response['data']:
								gRNA_annotation.update(data_access_response['data'][uniprot_offset])
						
					annotated_guide_features[editor][gRNA] = gRNA_annotation
		
		return annotated_guide_features
	
	###########################################
 
	def annotate_guide_motifs(self,annotated_gRNAs=None):
		"""
		Annotate guide motifs based on the provided annotated gRNAs.
		@param annotated_gRNAs - the annotated guide RNAs
		@return annotated_guide_motif - a dictionary containing annotated guide motifs
		"""

		annotated_guide_protein_centric = {}
		annotated_guide_motif = {}

		self.options['momap_base_path'] = os.path.join(os.path.join(self.options['data_path'],"momap"))
		if not os.path.exists(self.options['momap_base_path']):
			logger.error("Making: " + self.options['momap_base_path'])
			os.mkdir(self.options['momap_base_path'])

		self.options['momap_data_path'] = os.path.join(os.path.join(self.options['data_path'],"momap","motif_containing"))
		if not os.path.exists(self.options['momap_data_path']):
			logger.error("Making: " + self.options['momap_data_path'])
			os.mkdir(self.options['momap_data_path'])

		if annotated_gRNAs == None:
			annotated_gRNAs =  self.annotate_protein_guides()
		
		for editor in annotated_gRNAs:
			for gRNA in annotated_gRNAs[editor]:
				if 'accession' in  annotated_gRNAs[editor][gRNA]:
					accession = annotated_gRNAs[editor][gRNA]['accession']
					if accession not in annotated_guide_protein_centric: annotated_guide_protein_centric[accession] = {}
					if editor not in annotated_guide_protein_centric[accession]: annotated_guide_protein_centric[accession][editor] = {}
					annotated_guide_protein_centric[accession][editor][gRNA] = [int(o) for o in annotated_gRNAs[editor][gRNA]['difference'].keys()]
				else:
					print(annotated_gRNAs[editor][gRNA])
		
		for accession in annotated_guide_protein_centric:
			out_path = os.path.join(self.options['momap_data_path'], accession + ".motif_containing.jsonc")
			
			if not os.path.exists(out_path):
				data_access_response = queryRunner.queryRunner("database_access","get_instances_by_motif_accession",{"accession":accession,"database_name":"motifs"}).run()
				utilities_basic.write_to_json(out_path,data_access_response,zipped=True)
			else:
				data_access_response = utilities_basic.read_from_json(out_path,zipped=True)
			
			if "data" not in data_access_response:
				logger.error("Deleting " + out_path)
				os.remove(out_path)
				logger.error(data_access_response)
				
				for editor in annotated_guide_protein_centric[accession]:
					if editor not in annotated_guide_motif: annotated_guide_motif[editor] = {}

				continue

			for editor in annotated_guide_protein_centric[accession]:
				if editor not in annotated_guide_motif: annotated_guide_motif[editor] = {}
				for gRNA in  annotated_guide_protein_centric[accession][editor]:
					gRNA_annotation = {}
					for uniprot_offset in annotated_guide_protein_centric[accession][editor][gRNA]:
						for motif in data_access_response['data']:
							if uniprot_offset in range(motif['motif_start_residue_number'],motif['motif_end_residue_number']+1):

								specificity_classes = []

								if motif['specificity_classes'] != None:
									specificity_classes = [specificity_class['specificity_acc'] + ":" + specificity_class['specificity_id'] + ":" + specificity_class['name'] if specificity_class != None else "" for specificity_class in motif['specificity_classes']]
								
								domain_proteins = []
								if motif['domain_protein_relation'] != None:
									for domain_protein in motif['domain_protein_relation']:
										domain_protein_data = []
										for tag in [					
											'domain_pfam_accession',
											'domain_pfam_long_name',
											'domain_pfam_name',
											'domain_protein_uniprot_accession',
											'domain_protein_gene_name'
										]:
											if domain_protein[tag] == None:
												domain_protein_data.append("-")
											else:
												domain_protein_data.append(str(domain_protein[tag]))	
										
										if ":".join(domain_protein_data) not in domain_proteins:
											domain_proteins.append(":".join(domain_protein_data))

								pmids = []
								if 'references' in motif:
									if motif['references'] == None: continue
									for reference in motif['references']:
										if 'pmid' in reference:
											pmids.append(reference['pmid'])
								
								motif_filtered = {
									"motif_sequence":motif['motif_sequence'],
									"motif_start_residue_number":motif['motif_start_residue_number'],
									"motif_end_residue_number":motif['motif_end_residue_number'],
									"interaction_types":motif['interaction_types'],
									"pmid":pmids,
									"binding_proteins":domain_proteins,
									"specificity_classes":specificity_classes
								}

								gRNA_annotation[motif['group_by_id']] = motif_filtered
					
					annotated_guide_motif[editor][gRNA] = list(gRNA_annotation.values())
		
		return annotated_guide_motif

	###########################################

	def annotate_guide_accessiblity(self,annotated_gRNAs=None):
		"""
		Annotate the accessibility of guide RNAs based on protein-centric information.
		@param self - the object instance
		@param annotated_gRNAs - a dictionary of annotated guide RNAs
		@return a dictionary of annotated guide accessibility
		"""
		annotated_guide_protein_centric = {}
		annotated_guide_accessibility = {}

		if annotated_gRNAs == None:
			annotated_gRNAs =  self.annotate_protein_guides()

		for editor in annotated_gRNAs:
			for gRNA in annotated_gRNAs[editor]:
				if 'accession' in  annotated_gRNAs[editor][gRNA]:
					accession = annotated_gRNAs[editor][gRNA]['accession']
					if accession not in annotated_guide_protein_centric: annotated_guide_protein_centric[accession] = {}
					if editor not in annotated_guide_protein_centric[accession]: annotated_guide_protein_centric[accession][editor] = {}
					annotated_guide_protein_centric[accession][editor][gRNA] = [int(o) for o in annotated_gRNAs[editor][gRNA]['difference'].keys()]
				else:
					print(annotated_gRNAs[editor][gRNA])

		for accession in annotated_guide_protein_centric:
			tessellation_accessibility_response = queryRunner.queryRunner("alphafold","get_tessellation_accessibility",{"accession":accession}).run()
			uniprot_topology_response = queryRunner.queryRunner("uniprot","get_protein_topology_string",{"accession":accession}).run()
			alphafold_structural_classification_response = queryRunner.queryRunner("alphafold","get_alphafold_structural_classification",{"accession":accession}).run()
			alphafold_accessibility_response = queryRunner.queryRunner("alphafold","get_accessibility",{"accession":accession}).run()
			alphafold_pLDDT_response = queryRunner.queryRunner("alphafold","get_pLDDT",{"accession":accession}).run()
			alphafold_disorder_response = queryRunner.queryRunner("alphafold","get_accessibility_windowed_scores",{"accession":accession}).run()
			
			for editor in annotated_guide_protein_centric[accession]:
				if editor not in annotated_guide_accessibility: annotated_guide_accessibility[editor] = {}
				for gRNA in  annotated_guide_protein_centric[accession][editor]:
					gRNA_annotation = {}
					for uniprot_offset in annotated_guide_protein_centric[accession][editor][gRNA]:
						try:
							tessellation_accessibility_sidechain_atom = tessellation_accessibility_response['data']['chain_details']['A']['accessible_residues'][str(uniprot_offset)]['side_chain_score']
							tessellation_accessibility_any_atom = tessellation_accessibility_response['data']['chain_details']['A']['accessible_residues'][str(uniprot_offset)]['any_atm_score']
							tessellation_accessibility = tessellation_accessibility_sidechain_atom > 0 or tessellation_accessibility_any_atom > 0
							tessellation_accessibility_type = ""

							if tessellation_accessibility_sidechain_atom + tessellation_accessibility_any_atom == 0:
								tessellation_accessibility_type = "buried"
							elif tessellation_accessibility_sidechain_atom == 0:
								tessellation_accessibility_type = "sidechain_buried"
							elif tessellation_accessibility_sidechain_atom > 0:
								tessellation_accessibility_type = "sidechain_accessible"

							gRNA_annotation[uniprot_offset] = {
								"topology":uniprot_topology_response['data'][uniprot_offset-1],
								"structural_module_class":alphafold_structural_classification_response['data']['structural_module_class_classification'][uniprot_offset-1],
								"alphafold_accessibility":alphafold_accessibility_response['data'][uniprot_offset-1],
								"alphafold_pLDDT":alphafold_pLDDT_response['data'][uniprot_offset-1],
								"alphafold_disorder":alphafold_disorder_response['data'][uniprot_offset-1],
								"tessellation_accessibility":tessellation_accessibility_response['data']['chain_details']['A']['accessible_residues'][str(uniprot_offset)],
								"tessellation_accessibility_sidechain_atom":tessellation_accessibility_response['data']['chain_details']['A']['accessible_residues'][str(uniprot_offset)]['side_chain_score'],
								"tessellation_accessibility_any_atom":tessellation_accessibility_response['data']['chain_details']['A']['accessible_residues'][str(uniprot_offset)]['any_atm_score'],
								"tessellation_accessibility":tessellation_accessibility,
								"tessellation_accessibility_type":tessellation_accessibility_type
							}
						except:
							pass
					
					annotated_guide_accessibility[editor][gRNA] = gRNA_annotation
		
		return annotated_guide_accessibility
	
	###########################################

	def annotate_guide_snps(self,annotated_gRNAs=None):
		"""
		Annotate guide SNPs based on annotated gRNAs and protein-centric annotations.
		@param self - the object instance
		@param annotated_gRNAs - a dictionary of annotated gRNAs
		@return annotated_guide_snps - a dictionary of annotated guide SNPs
		"""

		annotated_guide_protein_centric = {}
		annotated_guide_snps = {}
	
		if annotated_gRNAs == None:
			annotated_gRNAs =  self.annotate_protein_guides()

		for editor in annotated_gRNAs:
			for gRNA in annotated_gRNAs[editor]:
				accession = annotated_gRNAs[editor][gRNA]['accession']
				if accession not in annotated_guide_protein_centric: annotated_guide_protein_centric[accession] = {}
				if editor not in annotated_guide_protein_centric[accession]: annotated_guide_protein_centric[accession][editor] = {}
				annotated_guide_protein_centric[accession][editor][gRNA] = [int(o) for o in annotated_gRNAs[editor][gRNA]['difference'].keys()]

		for accession in annotated_guide_protein_centric:
			snp_data = queryRunner.queryRunner("mutations","parse_snps",{"accession":accession}).run()
			snp_positions = {}

			try:
				for snp in snp_data['data']:
					for i in range(snp['start'],snp['end']+1):
						if i not in snp_positions:
							snp_positions[i] = []
						
						snp_positions[i].append(snp)
			except:
				pass

			for editor in annotated_guide_protein_centric[accession]:
				if editor not in annotated_guide_snps: annotated_guide_snps[editor] = {}
				for gRNA in  annotated_guide_protein_centric[accession][editor]:
					for uniprot_offset in annotated_guide_protein_centric[accession][editor][gRNA]:
						if uniprot_offset in snp_positions:
							if gRNA not in annotated_guide_snps[editor]:
								annotated_guide_snps[editor][gRNA] = {}
							if uniprot_offset not in annotated_guide_snps[editor][gRNA]:
								annotated_guide_snps[editor][gRNA][uniprot_offset] = []
							
							annotated_guide_snps[editor][gRNA][uniprot_offset] += snp_positions[uniprot_offset]

		return annotated_guide_snps

	###########################################

	def annotate_guides_efficiency(self):
   
		sys.path.append(os.path.abspath(os.path.join(file_path,"../tools/be_predict")))
		try:
			from efficiency import predict as be_efficiency_model
		except:
			return {
				"status":"error",
				"error_type":"BE-hive be_predict_efficiency import failed. Requires python 3.7. Please install https://github.com/maxwshen/be_predict_efficiency",
				"expected_path":os.path.abspath(os.path.join(file_path,"../tools/be_predict"))
			}
		try:
			sys.path.append(os.path.abspath(os.path.join(file_path,"../tools/be_predict")))
			from bystander import predict as bystander_model
		except:
			return {
				"status":"error",
				"error_type":"BE-hive be_predict_bystander import failed. Requires python 3.7. Please install https://github.com/maxwshen/be_predict_bystander",
				"expected_path":os.path.abspath(os.path.join(file_path,"../tools/be_predict"))
			}
		
		annotated_guide_efficiency = {}

		be_efficiency_model.init_model(base_editor = self.options['behive_base_editor'], celltype = self.options['behive_celltype'])
		bystander_model.init_model(base_editor = self.options['behive_base_editor'], celltype = self.options['behive_celltype'])

		for accession in self.options['protein_gRNAs']:
			try:
				uniprot_basic_data = queryRunner.queryRunner("uniprot","parse_basic",{"accession":accession}).run()

				options = {
					"accession":accession,
					"peptide_start":0,
					"peptide_end":len(uniprot_basic_data['data']['sequence']),
					"flanks":self.options['guides_efficiency_flanks'],
					"peptide_sequence":uniprot_basic_data['data']['sequence'],
					"species":self.options['species']
				}
			
				
				response = queryRunner.queryRunner("ensembl","get_ensembl_data_for_region",options).run()
				
				peptide_genomic_flanks = response['data']['peptide_genomic_flanks']
				peptide_genomic_flanks_complementary = response['data']['peptide_genomic_flanks_complementary'][::-1]
				peptide_cds_gapped_flanks = "-"*self.options['guides_efficiency_flanks']+ response['data']['peptide_cds_gapped'] + "-"*self.options['guides_efficiency_flanks']
				
				for gRNA in self.options['protein_gRNAs'][accession]:
					try:
						if peptide_genomic_flanks.count(gRNA):
							flanks = peptide_genomic_flanks.split(gRNA)
							flank_5 = flanks[0][-self.options['behive_flank_5']:]
							flank_3 = flanks[1][:self.options['behive_flank_3']]
							gRNA_flanked = "".join([flank_5,gRNA,flank_3])
						elif peptide_genomic_flanks_complementary.count(gRNA):
							flanks =peptide_genomic_flanks_complementary.split(gRNA)
							flank_5 = flanks[0][-self.options['behive_flank_5']:]
							flank_3 = flanks[1][:self.options['behive_flank_3']]
							gRNA_flanked = "".join([flank_5,gRNA,flank_3])
						else:
							logger.error("No genomic match for " + gRNA)

						pred_be_efficiency = be_efficiency_model.predict(gRNA_flanked)
						pred_bystander_model, stats_bystander_model = bystander_model.predict(gRNA_flanked)

						headers = list(pred_bystander_model.columns.values)

						gRNA_edits_pooled_efficiency = {}
						gRNA_edits_pooled_efficiency_by_aa = {}
						maximum_frequency = {
							"all":{"edit":"","frequency":0},
							"single":{"edit":"","frequency":0}
						}

						uniprot_offsets = []
						behive_logit_score = 0

						if stats_bystander_model['Assumed protospacer sequence'] == gRNA:
							for index, row in pred_bystander_model.iterrows():
								gRNA_flanked_edit = list(gRNA_flanked)
								if row['Predicted frequency'] > 0.05:
									for header in headers:
										if header != 'Predicted frequency':
											if row[header] != gRNA_flanked_edit[int(header[1:])+20-1] :
												gRNA_flanked_edit[int(header[1:])+20-1] = row[header]#.lower()
												
									gRNA_edits = self.calculate_gRNA_edits_behive(gRNA_flanked,"".join(gRNA_flanked_edit))
									gRNA_edits = self.get_mutated_protein(gRNA_flanked,gRNA_edits,peptide_genomic_flanks,peptide_cds_gapped_flanks,"")
									gRNA_edits['behive_predicted_frequency'] = row['Predicted frequency']
									gRNA_edits['behive_logit_score'] = pred_be_efficiency['Predicted logit score']
									gRNA_edits['gRNA'] = gRNA

									uniprot_offsets = uniprot_offsets + list(gRNA_edits['difference'].keys())

									behive_logit_score = gRNA_edits['behive_logit_score']
									
									if ",".join(gRNA_edits['changes']) not in gRNA_edits_pooled_efficiency:
										gRNA_edits_pooled_efficiency[",".join(gRNA_edits['changes'])] = {"behive_predicted_frequency_aggregated":0,"behive_behive_logit_score_aggregated":0,"edits":{}}
									
									gRNA_edits_pooled_efficiency[",".join(gRNA_edits['changes'])]["behive_predicted_frequency_aggregated"] += row['Predicted frequency']
									gRNA_edits_pooled_efficiency[",".join(gRNA_edits['changes'])]["edits"]["".join(gRNA_flanked_edit)] = gRNA_edits

									if len(gRNA_edits['changes']) == 0:
										if "" not in gRNA_edits_pooled_efficiency_by_aa:
											gRNA_edits_pooled_efficiency_by_aa[""] = {"behive_predicted_frequency_aggregated":0}

										gRNA_edits_pooled_efficiency_by_aa[""]["behive_predicted_frequency_aggregated"] += row['Predicted frequency']
									else:
										for edit in gRNA_edits['changes']:
											if edit not in gRNA_edits_pooled_efficiency_by_aa:
												gRNA_edits_pooled_efficiency_by_aa[edit] = {"behive_predicted_frequency_aggregated":0}

											gRNA_edits_pooled_efficiency_by_aa[edit]["behive_predicted_frequency_aggregated"] += row['Predicted frequency']
						else:
							logger.error([stats_bystander_model['Assumed protospacer sequence'] == gRNA,stats_bystander_model['Assumed protospacer sequence'],gRNA])

						for edit in gRNA_edits_pooled_efficiency:
							if maximum_frequency["all"]['frequency'] < gRNA_edits_pooled_efficiency[edit]['behive_predicted_frequency_aggregated']:
								maximum_frequency["all"] = {"edit":edit,"frequency":gRNA_edits_pooled_efficiency[edit]['behive_predicted_frequency_aggregated']}

						for edit in gRNA_edits_pooled_efficiency_by_aa:
							if maximum_frequency["single"]['frequency'] < gRNA_edits_pooled_efficiency_by_aa[edit]['behive_predicted_frequency_aggregated']:
								maximum_frequency["single"] = {"edit":edit,"frequency":gRNA_edits_pooled_efficiency_by_aa[edit]['behive_predicted_frequency_aggregated']}

						annotated_guide_efficiency[gRNA] = {
							"gRNA_edits_pooled_efficiency":gRNA_edits_pooled_efficiency,
							"gRNA_edits_pooled_efficiency_by_aa":gRNA_edits_pooled_efficiency_by_aa,
							"maximum_frequency":maximum_frequency,
							"behive_base_editor":self.options['behive_base_editor'],
							"behive_celltype":self.options['behive_celltype'],
							"behive_logit_score":behive_logit_score,
							"uniprot_offsets":[int(uniprot_offset) for uniprot_offset in list(set(uniprot_offsets))]
						}
					except:
						logger.error(gRNA)
						logger.error(utilities_error.getError())
			except:
				logger.error(accession)

		return annotated_guide_efficiency
		
	###########################################

	def make_guide_annotation_table(self,annotated_gRNAs):
		rows = []
		rows = []
		editing_guides = {}
		for editor in annotated_gRNAs:
			editing_guides[editor] = 0
			for gRNA in annotated_gRNAs[editor]:
				row = [
					str(editing_guides[editor]), 
					editor,
					gRNA,
					annotated_gRNAs[editor][gRNA]['edit_window'],
					annotated_gRNAs[editor][gRNA]['edit_window_mutated'],
					str(annotated_gRNAs[editor][gRNA]['strand']),
					str(annotated_gRNAs[editor][gRNA]['genome_edited']),
					str(annotated_gRNAs[editor][gRNA]['proteome_edited']),
					str(annotated_gRNAs[editor][gRNA]['splice_site']),
					str(annotated_gRNAs[editor][gRNA]['stop']),
					annotated_gRNAs[editor][gRNA]['wildtype_region'],
					annotated_gRNAs[editor][gRNA]['mutated_region'],
					",".join(annotated_gRNAs[editor][gRNA]['changes']),
					",".join([str(x) for x in annotated_gRNAs[editor][gRNA]['difference'].keys()]),
				]

				if len(annotated_gRNAs[editor][gRNA]['changes']) > 0 and not annotated_gRNAs[editor][gRNA]['splice_site']:
					editing_guides[editor] += 1

				rows.append("\t".join(row))

		