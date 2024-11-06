import os
import sys

import inspect

file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"../"))
import option_reader

sys.path.append(os.path.join(file_path,"../evolution"))
import genetreeParser
import gopherParser

sys.path.append(os.path.join(file_path,"../conservation"))
import slimprints
import conservationScorer

sys.path.append(os.path.join(file_path,"../pocket_discovery"))

try:
	from alphafold_pocket_scorer import AlphaFoldPocketScorer
	from pdb_pocket_scorer import PDBPocketScorer
except:
	print("Error importing xProtCAS")

from dataManager import dataManager

class evolutionManager(dataManager):

	##------------------------------------------------------------------##
	## Inherits functions from dataManager
	##------------------------------------------------------------------##

	def setup_data(self):
		self.default_task_options = {
			"is_superuser":True,
			"accession":None,
			"orthdb_taxon_id":None,
			"verbose":False,
			"debug":False,
			"logfile":False,
			"username":None,
			"password":None,
			#"homologue_type":"orthologues",
			#"region_start":None,
			#"region_end":None,
			"conservation_score_type":"WCS_p",
			"remake_alignment_age":1800,
			"weighted_scores_only": False,
			"reset_remote_dataset":False,
			"domain":"all",
			"download_file":False,
			"side_chain_accessibility":True
		}

		self.allowed_options = [
			#"is_superuser",
			"accession",
			"task",
			"pdb_id",
			"chain_id",
			"conservation_score_type",
			"orthdb_taxon_id",
			"conservation_score_type",
			"homologue_type",
			"region_start",
			"region_end",
			"pmid",
			"remake_alignment_age",
			"proteome_type",
			"weighted_scores_only",
			"run_peptools",
			"swissprot_only",
			"slimprint_require_feature",
			"slimprint_sig_cutoff",
			"username",
			"password",
			"split_resolution",
			"split_accessibility",
			"parse_pfam",
			"iteration_number",
			"conservation_from_api",
			"get_accessibility",
			"get_neighbors",
			"get_active_sites",
			"get_known_interface_contacts",
			"get_known_interface_contacts_from_api",
			"get_modifications",
			"get_mutations",
			"tessellation_accessibility_cutoff",
			"weight_with_accessibility",
			"task"
		]

		self.allowed_options_admin = [
			"taxon_id",
			"reset_remote_dataset",	
			"outfile",
			"print_results",
			"remake"
		]

		self.task_options = [
			"find_homologues_genetree",
			"find_homologues_gopher",
			"get_searchdb_taxons",
			"get_searchdb_list",
			"run_slimprints",
			"get_conservation_scores",
			"get_conservation_score",
			"get_conservation_scores_by_pdb",
			"parse_conservation_scores",
			"get_alignment_query_sequence_gopher",
			"get_alignment_fasta",
			"get_alphafold_centrality_scores",
			"get_alphafold_centrality_pockets",
			"get_pdb_centrality_scores",
			"get_pockets_in_human_proteins",
			"get_searchdb_species_list",
			"help"
		]

		self.task_options_admin = [
			"cluster_slimprints",
			"upload_slimprints_data",
			"get_r4s_conservation_scores",
			"make_search_database"
		]

		self.required = {
			'all':[],
			'find_homologues_genetree':['accession'],
			'find_homologues_gopher':['accession','orthdb_taxon_id'],
			'get_alignment_query_sequence_gopher':['accession','orthdb_taxon_id'],
			'get_conservation_scores': ['accession','orthdb_taxon_id'],
			'parse_conservation_scores': ['accession','orthdb_taxon_id'],
			'get_alignment': ['accession', 'orthdb_taxon_id'],
			'get_searchdb_taxons': ['orthdb_taxon_id'],
			"get_searchdb_list":[],
			'run_slimprints':['accession','orthdb_taxon_id'],
			'get_alphafold_centrality_scores':['accession'],
			'get_alphafold_centrality_pockets':['accession'],
			'get_pdb_centrality_scores':['pdb_id'],
			'get_r4s_conservation_scores': ['accession','orthdb_taxon_id'],
			'get_searchdb_species_list': ['orthdb_taxon_id']
		}

		self.test_options = {
			'accession':"Q05323",
			'orthdb_taxon_id':"1570291",
			'pdb_id':"2AST",
			"chain_id":"A"
			}

		self.required_type = {
			'accession':"string",
			'orthdb_taxon_id':"string",
			'remake':"bool"
		}

		self.convert_type = {
			'orthdb_taxon_id':"string"
		}

		self.required_format = {
			'accession':"\A([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\-{0,1}[0-9]*\Z"
		}

		# delimiter to split options
		self.list_options = {
			'accession':","
		}

		self.required_valid_options_list = {
			'remake':[True,False],
			'conservation_score_type':['AWCS','WCS','ABS','ABV','WCS_p','ABV_p'],
			"homologue_type":["paralogues","orthologues","homologues"],
			'tasks':self.task_options
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

		if task_options["task"] == "find_homologues_genetree":
			findHomologuesObj = genetreeParser.findHomologues()
			findHomologuesObj.options.update(task_options)
			data = findHomologuesObj.run_job()

		if task_options["task"] == "find_homologues_gopher":
			findHomologuesObj = gopherParser.findHomologues()
			findHomologuesObj.options.update(task_options)
			data = findHomologuesObj.run_job()

		if task_options["task"] == "get_alignment_query_sequence_gopher":
			findHomologuesObj = gopherParser.findHomologues()
			findHomologuesObj.options.update(task_options)
			data = findHomologuesObj.get_alignment_query_sequence_gopher()

		if task_options["task"] == "get_searchdb_taxons":
			findHomologuesObj = gopherParser.findHomologues()
			findHomologuesObj.options.update(task_options)
			data = findHomologuesObj.get_searchdb_taxons()

		if task_options["task"] == "make_search_database":
			findHomologuesObj = gopherParser.findHomologues()
			findHomologuesObj.options.update(task_options)
			data = findHomologuesObj.make_search_database()

		if task_options["task"] == "get_searchdb_list":
			findHomologuesObj = gopherParser.findHomologues()
			findHomologuesObj.options.update(task_options)
			data = findHomologuesObj.get_searchdb_list()

		if task_options["task"] == "run_slimprints":
			slimprintsHelperObj = slimprints.SLiMPrints()
			slimprintsHelperObj.options.update(task_options)
			data = slimprintsHelperObj.runSLiMPrints()

		if task_options["task"] == "cluster_slimprints":
			slimprintsHelperObj = slimprints.SLiMPrints()
			slimprintsHelperObj.options.update(task_options)
			data = slimprintsHelperObj.clusterSLiMPrints()
			
		if task_options["task"] == "upload_slimprints_data":
			slimprintsHelperObj = slimprints.SLiMPrints()
			slimprintsHelperObj.options.update(task_options)
			data = slimprintsHelperObj.uploadSLiMPrints()

		if task_options['task'] == 'get_conservation_scores':
			conservationScorerObj = conservationScorer.conservationScorer()
			conservationScorerObj.options.update(task_options)
			#conservationScorerObj.options['accession'] = self.options['accession']
			#conservationScorerObj.options['orthdb_taxon_id'] = self.options['orthdb_taxon_id']
			data = conservationScorerObj.getConservationScores()

		if task_options['task'] == 'get_conservation_score':
			conservationScorerObj = conservationScorer.conservationScorer()
			conservationScorerObj.options.update(task_options)
			data = conservationScorerObj.getConservationScore()

		if task_options['task'] == 'get_conservation_scores_by_pdb':
			conservationScorerObj = conservationScorer.conservationScorer()
			conservationScorerObj.options.update(task_options)
			#conservationScorerObj.options['accession'] = self.options['accession']
			#conservationScorerObj.options['orthdb_taxon_id'] = self.options['orthdb_taxon_id']
			data = conservationScorerObj.getConservationScoresByPDB()

		if task_options['task'] == 'parse_conservation_scores':
			conservationScorerObj = conservationScorer.conservationScorer()
			conservationScorerObj.options.update(task_options)
			data = conservationScorerObj.parseConservationScores()

		if task_options['task'] == 'get_alignment_fasta':
			findHomologuesObj = gopherParser.findHomologues()
			findHomologuesObj.options.update(task_options)
			findHomologuesObj.run_job()
			data =  findHomologuesObj.get_alignment()

		if task_options['task'] == 'get_alphafold_centrality_scores':
			alphaFoldPocketScorerObj = AlphaFoldPocketScorer()
			alphaFoldPocketScorerObj.options.update(task_options)
			data = alphaFoldPocketScorerObj.get_pocket_scores()
		
		if task_options['task'] == 'get_alphafold_centrality_pockets':
			alphaFoldPocketScorerObj = AlphaFoldPocketScorer()
			task_options['parse_pfam'] = True
			task_options['split_accessibility'] = True
			task_options['get_accessibility'] = True
			task_options['get_neighbors'] = True
			task_options['get_active_sites'] = True
			task_options['get_known_interface_contacts'] = True
			task_options['get_modifications'] = True
			task_options['get_mutations'] = True
			alphaFoldPocketScorerObj.options.update(task_options)
			data = alphaFoldPocketScorerObj.get_pocket_scores()

		if task_options['task'] == 'get_pdb_centrality_scores':
			pdbPocketScorerObj = PDBPocketScorer()
			pdbPocketScorerObj.options.update(task_options)
			data = pdbPocketScorerObj.get_pocket_scores()
			
		if task_options['task'] == 'get_pockets_in_human_proteins':
			alphaFoldPocketScorerObj = AlphaFoldPocketScorer()
			alphaFoldPocketScorerObj.options.update(task_options)
			data = alphaFoldPocketScorerObj.get_pockets_in_human_proteins()

		if task_options['task'] == 'get_r4s_conservation_scores':
			conservationScorerObj = conservationScorer.conservationScorer()
			conservationScorerObj.options.update(task_options)
			data = conservationScorerObj.getRate4SiteConservationScores()


		if task_options['task'] == 'get_searchdb_species_list':
			findHomologuesObj = gopherParser.findHomologues()
			findHomologuesObj.options.update(task_options)
			data = findHomologuesObj.get_searchdb_species_list()

		return data

##------------------------------------------------------------------##
##------------------------   END    END  ---------------------------##
##------------------------------------------------------------------##

if __name__ == "__main__":

	dataManagerObj = evolutionManager()
	dataManagerObj.main()
