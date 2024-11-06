import os, pprint, sys,json, logging

try:
	import svg
except:
	logging.error("base_editing_protein_plots")
	logging.error("python3 -m pip install --user svg.py")
	logging.error("python3 -m pip install --user svg.py")
	logging.error("python3 -m pip install --user svg.py")
	logging.error("python3 -m pip install --user svg.py")
	logging.error("python3 -m pip install --user svg.py")
	logging.error("python3 -m pip install --user svg.py")
	logging.error("python3 -m pip install --user svg.py")
	logging.error("python3 -m pip install --user svg.py")
	logging.error("python3 -m pip install --user svg.py")
	logging.error("python3 -m pip install --user svg.py")
	logging.error("python3 -m pip install --user svg.py")

	#raise

try:
	from cairosvg import svg2png
except:
	logging.error("base_editing_protein_plots")
	logging.error("python3 -m pip install --user cairosvg")
	logging.error("python3 -m pip install --user cairosvg")
	logging.error("python3 -m pip install --user cairosvg")
	logging.error("python3 -m pip install --user cairosvg")
	logging.error("python3 -m pip install --user cairosvg")
	logging.error("python3 -m pip install --user cairosvg")
	logging.error("python3 -m pip install --user cairosvg")
	logging.error("python3 -m pip install --user cairosvg")
	logging.error("python3 -m pip install --user cairosvg")
	logging.error("python3 -m pip install --user cairosvg")
	logging.error("python3 -m pip install --user cairosvg")

	#raise

from textwrap import dedent

sys.path.append(os.path.join(os.path.dirname(__file__),"../../data_management"))
import queryRunner


sys.path.append(os.path.join(os.path.dirname(__file__),"../../utilities"))
import utilities_basic


"""
The plots for each protein/condition are as follows:

Top bar shows features:
- grey: Domains defined by Pfam and AlphaFold

Middle plot shows attributes:
- blue plot : residue accessibility
- green plot : residue conservation

Bottom plot shows guide data:
- circle: guide
- y-axis: fold change
- circle size: p-value
- circle opacity: the specificity of the guide - lighter colours have some off target effect
"""

tags = [
	["design","accession"],
	["design","gene"],
	["design","peptide"],
	["design","mutations"],
	["design","peptide_start"],
	["peptools","Accessibility_raw"],
	["peptools","Accessibility"],
	["peptools","pLDDT_raw"],
]

guide_tags = [
	["conservation_metazoa"],
	["abe_edit_mutation"],
	["cbe_edit_mutation"],
	["limma_scores","log_fold_change"],
	["limma_scores","log_p_value"],
	["limma_scores","log_p_value_adj"],
	["limma_scores","wildtype"],
	["limma_scores","mutant"],	
]

header = [
	"sample",
	'gRNA_id',
	"gRNA_sequence",
	"accession",
	"gene",
	"position",
	"mutations",
	"domain",
	"accessibility",
	"conservation_metazoa",
	"log_fold_change",
	"log_p_value",
	"log_p_value_adj",
	"normalised counts tpm",
	"normalised counts fold change",
	"counts",
	"log2_tpm_counts"
]

class baseEditingProteinPlots():

	##################################################################################################

	def __init__(self):
		self.data = {}
		self.annotation = {}
		self.json = {}
		self.options = {}
		self.table_rows = []

	##################################################################################################
	
	def plot_samples(self):
		for sample in self.options['samples'] :
			self.options['sample'] = sample
			self.options['guides_annotation_file'] = self.options['samples'][sample]['guides_annotation_file'] #'/home/data/base_editing/results/tiled_pilot_' + self.options['sample'].split("_")[0] + '/tiled_pilot_' + self.options['sample'] + '.18.guide.annotation.json'
		
			self.plot_proteins()

		for accession in self.json:
			json_file = os.path.abspath("./data/" + accession + ".screens.json")
			status = utilities_basic.write_to_json(json_file,self.json[accession])

		tdt_file = "complete_screen.tdt"
		open(tdt_file,"w").write("\n".join(["\t".join(header)]+ self.table_rows))	

	def plot_proteins(self):
		if 'guides_annotation_file' in self.options:
			plot_guide_data = json.loads(open(self.options['guides_annotation_file']).read())
		else:
			plot_guide_data = self.options['plot_guide_data']

		for guide_id in plot_guide_data.keys():
			try:
				if plot_guide_data[guide_id]['control'] == False:
					if plot_guide_data[guide_id]['accession'] not in self.data:
						self.data[plot_guide_data[guide_id]['accession']] = {}
					
					if int(plot_guide_data[guide_id]['peptide_start']) not in self.data[plot_guide_data[guide_id]['accession']] :
						self.data[plot_guide_data[guide_id]['accession']][int(plot_guide_data[guide_id]['peptide_start'])] = {}

					self.data[plot_guide_data[guide_id]['accession']][int(plot_guide_data[guide_id]['peptide_start'])][guide_id] = plot_guide_data[guide_id]
			except:
				print("Error: " + guide_id)

		for accession in self.data.keys():
			try:
				self.options['accession'] = accession
				self.plot_protein()
			except:
				raise

	def get_protein_motifs(self):	
		regions = []
		return regions
		response = queryRunner.queryRunner("database_access","get_instances_by_motif_accession",{"database_name":"motifs","accession":self.options['accession']}).run()
		for motif in response['data']:	
			specificity_classes = {}
			if motif['specificity_classes'] != None:
				specificity_classes = ", ".join([specificity_class['specificity_id'] for specificity_class in motif['specificity_classes'] if specificity_class['specificity_id'] != None])

			if motif['interaction_types'] != None:
				interaction_types = ", ".join(motif['interaction_types'])

			pdb_ids = "" 
			if motif['pdb_ids'] != None:
				pdb_ids = ",".join(motif['pdb_ids'])

			pmids = "" 
			try:
				pmids = ",".join([str(reference['pmid']) for reference in motif['references']])
			except:
				pass

			details = {
				"interaction_types":interaction_types,
				"curation_ids": ",".join([str(x) for x in motif['curation_ids']]),
				"pmid":pmids ,
				"pdb":pdb_ids,
				"curation_source":",".join(motif['sources']),
				"specificity_classes":specificity_classes
			}
		
			###-----------------_-------------_-_--_-_---_---_-------_------#

			peptide_start = motif['motif_start_residue_number']
			peptide_end = motif['motif_end_residue_number']
			motif_sequence =  motif['motif_sequence']
			peptide_sequence = motif['motif_sequence']


			if peptide_start == 1:
				peptide_start = 2
				peptide_sequence = peptide_sequence[1:]
				motif_sequence = motif_sequence[1:]
		
			region = {
				'accession':motif['motif_protein_uniprot_accession'].replace("-1",""),
				'peptide_start':peptide_start,
				'peptide_end':peptide_end,
				'peptide_sequence':peptide_sequence,
				'motif':motif_sequence,
				'source':'didi',
				'specificity_classes':specificity_classes,
				"details":details
			}
			regions.append(region)

		return regions

	def plot_protein(self):	
		try:
			logging.debug("Making: " + self.options['accession'])
			if self.options['accession'] not in self.json:
				
				uniprot_data = queryRunner.queryRunner("uniprot","parse_uniprot",accession=self.options['accession']).run()
				if 'data' not in uniprot_data:
					return {}
				
				self.json[self.options['accession']] = {
					"screens":{},
					"info":{},
					"motifs":{},
					"structured_modules":{},
					"annotation":{}
				}

				self.json[self.options['accession']]['info'].update(uniprot_data['data'])
				self.json[self.options['accession']]['motifs'] = self.get_protein_motifs()

			self.json[self.options['accession']]['screens'][self.options['sample']] = {}

			png_file = os.path.join(self.options['protein_png_output_directory_path'],self.options['sample'] + ".png")
			svg_file = os.path.join(self.options['protein_png_output_directory_path'],self.options['sample'] + ".svg")

			logging.debug("Making png file: " + png_file)

			if os.path.exists(png_file) and False: return 

			size_multiplier = 5
			x_padding = 30
			y_padding = 40
			plot_height = 25
			y_base_line = 70

			svg_width = (len(self.json[self.options['accession']]['info']['sequence']) + x_padding*2)*size_multiplier
			svg_height = (y_padding*2 + 100)*size_multiplier

			elements = [
				svg.Style(
					text=dedent("""
						.offset { font: bold """ + str(9*size_multiplier/2) + """px Helvetica; }
						.domain_name { font: bold """ + str(9*size_multiplier/2) + """px Helvetica; font-family:"Arial, Helvetica, sans-serif";}
					"""),
				),
			]

			###-–-----------------------###
			
			if self.options['accession'] not in self.annotation:
				alphafold_domain_response = queryRunner.queryRunner("alphafold","get_domains_from_pae_matrix",{"accession":self.options['accession'],"remake":False,"parse_pfam":True}).run()
				alphafold_neigbours_response = queryRunner.queryRunner("evolution","get_alphafold_centrality_scores",{"accession":self.options['accession'],"domain":"all","parse_pfam":"True","get_active_sites":"True","get_neighbors":"True","side_chain_accessibility":True,"get_accessibility":True,}).run()
				accessibility_response = queryRunner.queryRunner("alphafold","get_accessibility",{"accession":self.options['accession'],"remake":False}).run()
				accessibility_windowed_response = queryRunner.queryRunner("alphafold","get_accessibility_windowed_scores",{"accession":self.options['accession'],"remake":False,"window_size":3}).run()
				conservation_response = queryRunner.queryRunner("evolution","get_conservation_scores",{"accession":self.options['accession'],"remake":False,"orthdb_taxon_id":"metazoa"}).run()
				
				
				self.annotation[self.options['accession']] = {
					"alphafold_domain_response":alphafold_domain_response,
					"alphafold_neigbours_response":alphafold_neigbours_response,
					"accessibility_response":accessibility_response,
					"accessibility_windowed_response":accessibility_windowed_response,
					"conservation_response":conservation_response
				}

				try:
					self.json[self.options['accession']]['annotation']['conservation'] = conservation_response['data']['WCS']
				except:
					self.json[self.options['accession']]['annotation']['conservation'] = {}

				try:
					self.json[self.options['accession']]['annotation']['accessibility'] = accessibility_windowed_response['data']
				except:
					self.json[self.options['accession']]['annotation']['accessibility'] = {}
			else:
				alphafold_domain_response = self.annotation[self.options['accession']]['alphafold_domain_response']
				alphafold_neigbours_response = self.annotation[self.options['accession']]['alphafold_neigbours_response']
				accessibility_response = self.annotation[self.options['accession']]['accessibility_response']
				accessibility_windowed_response = self.annotation[self.options['accession']]['accessibility_windowed_response']
				conservation_response = self.annotation[self.options['accession']]['conservation_response']

			###-–-----------------------###

			motifs = []
			for motif in self.json[self.options['accession']]['motifs']:
				
				motifs.append(svg.Rect(
					x=(x_padding + motif['peptide_start'])*size_multiplier,  
					y=10,
					width=(motif['peptide_end'] - motif['peptide_start'])*size_multiplier,
					height=(10)*size_multiplier,
					stroke="orange",
					fill="orange",
					stroke_width=0,
					fill_opacity="0.50",
				))

			domains = []
			domain_offset_mapping = {}

			self.json[self.options['accession']]['structured_modules'] = {}
			
			if 'data' in alphafold_domain_response:
				for region in alphafold_domain_response['data']:
					self.json[self.options['accession']]['structured_modules'][region] = {
						"name":"-"
					}
					domains.append(svg.Rect(
						x=(x_padding + alphafold_domain_response['data'][region]['start'])*size_multiplier,  
						y=10,
						width=(alphafold_domain_response['data'][region]['end'] - alphafold_domain_response['data'][region]['start'])*size_multiplier,
						height=(10)*size_multiplier,
						stroke="black",
						fill="black",
						stroke_width=0,
						fill_opacity="0.10",
					))

					domains.append(svg.Rect(
						x=(x_padding + alphafold_domain_response['data'][region]['start'])*size_multiplier, 
						y=10,
						width=(2)*size_multiplier,
						height=(10)*size_multiplier,
						stroke="black",
						fill="black",
						stroke_width=0,
						fill_opacity="0.25",
					))

					domains.append(svg.Rect(
						x=(x_padding + alphafold_domain_response['data'][region]['end'])*size_multiplier, 
						y=10,
						width=(2)*size_multiplier,
						height=(10)*size_multiplier,
						stroke="black",
						fill="black",
						stroke_width=0,
						fill_opacity="0.25",
					))

					pfam_names = []
					domain_starts = []
					domain_ends = []
					for pfam_accession in alphafold_domain_response['data'][region]['pfam']:
						
						pfam_names.append(alphafold_domain_response['data'][region]['pfam'][pfam_accession]['pfam_acc'])
						for match in alphafold_domain_response['data'][region]['pfam'][pfam_accession]['matches']:
							domain_starts.append(alphafold_domain_response['data'][region]['start'])
							domain_ends.append(alphafold_domain_response['data'][region]['end'])
							for i in range(alphafold_domain_response['data'][region]['start'] ,alphafold_domain_response['data'][region]['end'] ):
								domain_offset_mapping[i] = alphafold_domain_response['data'][region]['pfam'][pfam_accession]['pfam_acc']

						
					domains.append(svg.Text(
						x=(x_padding + (alphafold_domain_response['data'][region]['end'] + alphafold_domain_response['data'][region]['start'])/2)*size_multiplier, 
						y=10 + (10)*size_multiplier - 10,
						class_=["domain_name"], 
						text_anchor="middle",
						text="/".join(pfam_names)
					))

					self.json[self.options['accession']]['structured_modules'][region] = {
						"name":"/".join(pfam_names),
						"pfam_identifers":list(alphafold_domain_response['data'][region]['pfam'].keys()),
						"region_start":alphafold_domain_response['data'][region]['start'],
						"region_end":alphafold_domain_response['data'][region]['end'],
						#"direct_neighbors":alphafold_neigbours_response['data'][region]['direct_neighbors'],
						#"accessible_residues":alphafold_neigbours_response['data'][region]['accessible_residues']
					}
				
			###-–-----------------------###

			
			plots = []
			for i in range(1,len(self.json[self.options['accession']]['info']['sequence'])):
				try:
					plots.append(svg.Rect(
						x=(x_padding + i)*size_multiplier, 
						y=(y_padding - conservation_response['data']['WCS'][i]*plot_height)*size_multiplier,
						width=size_multiplier,
						height=(conservation_response['data']['WCS'][i]*plot_height)*size_multiplier,
						#stroke="black",
						fill="#3CB371",
						stroke_width=0,
						fill_opacity="0.5"
					))
				except:
					pass
				
				try:
					plots.append(svg.Rect(
						x=(x_padding + i)*size_multiplier, 
						y=(y_padding - accessibility_windowed_response['data'][i]*plot_height)*size_multiplier,
						width=size_multiplier,
						height=(2),
						fill="#4169E1",
						stroke_width=0
					))
				except:
					pass

			###---------------------------###

			background_line = []
			background_line.append(svg.Rect(
				x=(x_padding)*size_multiplier, 
				y=(y_padding) + y_base_line*size_multiplier,
				width=(len(self.json[self.options['accession']]['info']['sequence']))*size_multiplier,
				height=(1)*size_multiplier,
				fill="grey",
				fill_opacity="0.75",
			))


			background_intervals = []
			for i in range(1,4):
				background_intervals.append(svg.Rect(
					x=(x_padding)*size_multiplier, 
					y=y_padding + (y_base_line - 10*i)*size_multiplier,
					width=(len(self.json[self.options['accession']]['info']['sequence']))*size_multiplier,
					height=(1)*size_multiplier,
					fill="grey",
					fill_opacity="0.15",
					stroke="white",
					stroke_width=0
				))
				background_intervals.append(svg.Rect(
					x=(x_padding)*size_multiplier, 
					y=y_padding + (y_base_line + 10*i)*size_multiplier,
					width=(len(self.json[self.options['accession']]['info']['sequence']))*size_multiplier,
					height=(1)*size_multiplier,
					fill="grey",
					fill_opacity="0.15",
					stroke="white",
					stroke_width=0
				))
			
			markers = []
			for i in range(1,len(self.json[self.options['accession']]['info']['sequence']),100):
				markers.append(svg.Rect(
					x=(x_padding + i)*size_multiplier, 
					y=(y_padding)*size_multiplier - (4)*size_multiplier,
					width=(2),
					height=(4)*size_multiplier,
					stroke="black",
					fill="black",
					stroke_width=1
				))
				
				markers.append(svg.Text(
					x=(x_padding + i)*size_multiplier, 
					y=(y_padding - 5)*size_multiplier,
					class_=["offset"], 
					text_anchor="middle",
					text=str(i)
				))
			
			rows = []
			offsets =  list(self.data[self.options['accession']].keys())
			offsets.sort()

			for offset in offsets:
				self.json[self.options['accession']]['screens'][self.options['sample']][offset] = {}
				for guide_id in self.data[self.options['accession']][offset]:
					try:
							
						self.json[self.options['accession']]['screens'][self.options['sample']][offset][guide_id] = {}

						guide_data = self.data[self.options['accession']][offset][guide_id]

						guide_off_target_exact = 10

						try:
							guide_off_target_exact = int(guide_data['guide_off_target_exact'])
						except:
							pass
						
						#if min([guide_data["data"]['0']['1'],guide_data["data"]['0']['2'],guide_data["data"]['18']['1'],guide_data["data"]['18']['2']]) < 30:
						#	continue

						domain_overlap = ''
						if int(guide_data["peptide_start"])-1 in domain_offset_mapping:
							domain_overlap = domain_offset_mapping[int(guide_data["peptide_start"])-1]

						try:
							wcs = "%1.2f"%conservation_response['data']['WCS'][int(guide_data["peptide_start"])-1]
						except:
							wcs = ""
						
						try:
							try:
								elements.append(svg.Circle(
									cx=(x_padding + int(guide_data["peptide_start"]))*size_multiplier, 
									cy=y_padding + (y_base_line - guide_data["limma_scores"]["log_fold_change"]*10)*size_multiplier,
									#cy=y_padding + y_base_line - math.log2((guide_data["data"]['0']['1'] + guide_data["data"]['0']['2'])/(guide_data["data"]['18']['1'] + guide_data["data"]['18']['2']))*10,
									r=(guide_data["limma_scores"]["log_p_value"]/2)*size_multiplier,
									stroke="black",
									stroke_opacity="0.75",
									stroke_width=1,
									#fill_opacity=accessibility_response['data'][int(guide_data["peptide_start"])-1],
									fill="red",
									fill_opacity=1/guide_off_target_exact
								))
							except:
								pass

						except:
							raise
					except:
						logging.error(guide_id)

			#---------------#
			background = []
			background.append(svg.Rect(
					x=0, 
					y=0,
					width=max(svg_width,1000*size_multiplier),
					height=svg_height,
					fill="white"
				))
				
			canvas = svg.SVG(
				width=max(svg_width,1000*size_multiplier),
				height=svg_height,
				elements=[background + plots + background_line + background_intervals + domains + motifs + markers + elements]
			)
			
			open(svg_file,'w').write(str(canvas))

			self.table_rows += rows	
			svg2png(bytestring=str.encode(str(canvas)),write_to=png_file, dpi=600)
		
		except Exception as e:
			logging.error(e)
			raise
		#---------------#

	##################################################################################################
	
################################################################################
## SECTION II: MAIN PROGRAM                                                  ##
################################################################################

################################################################################
## END OF SECTION II                                                          ##
################################################################################

