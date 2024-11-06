import os,json,pprint,sys,inspect,math,statistics,re,random,copy

import pandas as pd

file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../"))
sys.path.append(os.path.join(file_path,"../utilities/"))
sys.path.append(os.path.join(file_path,"../data_management/"))
sys.path.append(os.path.join(file_path,"../limma/"))

import utilities_error 
import utilities_codons 
import utilities_stats
import utilities_log
import utilities_plots
import config_reader
import option_reader
import queryRunner

import seaborn as sns
import matplotlib.pyplot as plt
import operator as op

from scipy import stats

###  python3 baseEditingScreenParser.py --counts_directory /Users/adminndavey/Documents/Work/data/base_editing/SUB22_122-124_count_normalized  --region_annotation_file /Users/adminndavey/Documents/Work/data/genomes/gRNA/gRNA_libraries/motif.regions.tdt --debug True --benchmark_control_data True


#-----
import logging
logging.basicConfig(level=logging.DEBUG)
#-----

import seaborn as sns
import matplotlib.pyplot as plt
import operator as op


################################################################################
## SECTION I: simbaParser                                                     ##
################################################################################

short_names = {
	'Essential splice site':'ESS',
	'Non-targeting control':"NT",
	'Intergenic control':"I",
	'Plasmid backbone sequence':"PB"
}


class baseEditingControlPlotter():

	##################################################################################################

	def __init__(self):

		self.data = {}

		##---------------------------------##
		self.options = {}
		self.options["debug"] = True	
		self.options["remake"] = False
		self.options.update(option_reader.load_commandline_options(self.options,{}))

		##---------------------------------##

		utilities_log.setupLogger(self.options)

		##---------------------------------##

	##################################################################################################

	def make_plot_directory(self):
		self.options["plot_dir_path"] = os.path.join(self.options["data_directory"],"plots")
		if not os.path.exists(self.options["plot_dir_path"]):
			os.mkdir(self.options["plot_dir_path"])
		
	##################################################################################################

	def plot_count_data(self):
		self.make_plot_directory()
		guide_ids = list(self.data.keys())
		guide_ids.sort()

		counts_png_path =  os.path.join(self.options["plot_dir_path"],self.options['screen_name'] + '.' + str(self.options['compare_time_point']) + '.counts.png')

		plot_data = []
		for guide_id in guide_ids:
			try:
				for replicate in self.data[guide_id]['data_t0_normalised'][self.options['compare_time_point']]:
					plot_data.append({"replicate":str(replicate),"gRNAs counts log2 cpm":self.data[guide_id]['data_log2_cpm_normalised'][self.options['compare_time_point']][replicate]})
			except:
				logging.error(guide_id)
				pass

		#######################

		fig, ax = plt.subplots()
		custom_params = {"axes.spines.right": False, "axes.spines.top": False}
		sns.set_theme(style="ticks", rc=custom_params)

		plt.rcParams["figure.figsize"] = [10,10]
		
		df = pd.DataFrame(plot_data)
		g = sns.boxenplot(data=df,x="replicate",y="gRNAs counts log2 cpm")
		g.figure.set_size_inches(5, 5)

		g.set_xlabel("replicate",fontsize=10, weight="bold")
		g.set_ylabel("gRNAs counts log2 cpm",fontsize=10, weight="bold")
		
		plt.tight_layout()
		plt.savefig(counts_png_path,dpi = 500)
		logging.debug("Writing: " +  counts_png_path)
		plt.clf()

	##################################################################################################

	def plot_control_t0_normalised_data(self):
		self.make_plot_directory()
		guide_ids = list(self.data.keys())
		guide_ids.sort()

		scatter_png_path =  os.path.join(self.options["plot_dir_path"],self.options['screen_name'] + '.' + str(self.options['compare_time_point']) + '.scatter.png')
		box_png_path =  os.path.join(self.options["plot_dir_path"],self.options['screen_name'] + '.' + str(self.options['compare_time_point']) + '.box.png')

		plot_data = []
		for guide_id in guide_ids:
			if self.data[guide_id]['control']:
				guide_type = self.data[guide_id]['control_type']
				if guide_type == "Plasmid backbone sequence":
					continue
			else:
				guide_type = "target"
				

			guide_data = {
				"guide_id":guide_id,
				"guide_type":guide_type
			}

			try:
				sizes = [] 
				guide_data["fold change"] = [] 
				for replicate in self.data[guide_id]['data_t0_normalised'][self.options['compare_time_point']]:
					guide_data["fold change r" + str(replicate)] =  self.data[guide_id]['data_t0_normalised'][self.options['compare_time_point']][replicate]
					guide_data["fold change"].append(self.data[guide_id]['data_t0_normalised'][self.options['compare_time_point']][replicate])
					try:
						sizes.append(math.log10(self.data[guide_id]['data'][0][replicate]))
					except:
						sizes.append(0)

				guide_data["fold change"] = sum(guide_data["fold change"])/3
				guide_data['size'] =  min(sizes) 
						
				plot_data.append(guide_data)
			except:
				logging.error(guide_id)
				raise
	
		#########
		#########
		
		fig, ax = plt.subplots()
		custom_params = {"axes.spines.right": False, "axes.spines.top": False}
		sns.set_theme(style="ticks", rc=custom_params)
		
		plt.rcParams["figure.figsize"] = [20,20]
		
		df = pd.DataFrame(plot_data)
		df = df.set_index("guide_id") 

		df = df.sort_values(by=['guide_type'],ascending=False)
		
		variables = ["fold change r1", "fold change r2", "fold change r3"]#,"count t0 r1", "count t0 r2", "count t0 r3"]
		g = sns.PairGrid(df, corner=True, hue="guide_type", vars=variables)
			
		#g.map_diag(sns.histplot,bins=20)
		g.map_diag(sns.kdeplot,common_norm=False)
		g.map_lower(sns.scatterplot, size=df["size"], sizes=(0.1, 20), size_norm=(2,4),alpha = 0.5)
		g.axes[1,0].set_ylim(-7,7)
		g.axes[2,0].set_ylim(-7,7)
		g.axes[1,1].set_ylim(-7,7)
		g.axes[1,0].set_xlim(-7,7)
		g.axes[2,0].set_xlim(-7,7)
		g.axes[1,1].set_xlim(-7,7)
		#g.map_upper(sns.kdeplot)
		g.add_legend(fontsize='8', title_fontsize='10',loc='upper right')
		
		logging.debug("Writing: " +  scatter_png_path)
		plt.tight_layout()
		plt.savefig(scatter_png_path,dpi = 500)
		plt.clf()

		logging.debug("Writing: " +  scatter_png_path)
		
		#########
		#########

		fig, ax = plt.subplots()
		custom_params = {"axes.spines.right": False, "axes.spines.top": False}
		sns.set_theme(style="ticks", rc=custom_params)
		
		g = sns.boxenplot(df,x="fold change", y="guide_type",showfliers = True)#, orient='h',scale='exponential')#,order=)
		g.axes.set_xlim(-5,2)
		#g.set_axis_labels("guide class", "log fold change", labelpad=10)
	#	plt.xticks(plt.xticks()[0], #)
		
		g.figure.set_size_inches(12, 3)

		g.set_xlabel("log fold change",fontsize=10, weight="bold")
		g.set_ylabel("guide class",fontsize=10, weight="bold")
		
		plt.tight_layout()
		plt.savefig(box_png_path,dpi = 500)
		plt.clf()
		#sys.exit()

	##################################################################################################
	
	def plot_control_boxplot(self,sorted_keys = ['target','Essential splice site','Non-targeting control','Intergenic control'],sorted_keys_label = ['Target', 'ESS','NTC', 'IC']):
		self.make_plot_directory()
		box_png_path =  os.path.join(self.options["plot_dir_path"],self.options['screen_name'] + '.' + str(self.options['compare_time_point']) + '.box.png')

		base_editor_variants_distribution = {
			"abe_only":[],
			"cbe_only":[],
			"both":[]
		}

		try:
			guide_ids = list(self.data.keys())
			guide_ids.sort()

			for guide_id in guide_ids:
				try:
					if self.data[guide_id]['control']:
						base_editor_variant = self.data[guide_id]['control_type']
					else:
						base_editor_variant = "target"
					
					"""
					if 'abe_edit_mutation' not in self.data[guide_id] and 'cbe_edit_mutation' not in self.data[guide_id]:
						
					elif 'abe_edit_mutation' not in self.data[guide_id] and 'cbe_edit_mutation' in self.data[guide_id]:
						base_editor_variant = "target"
					elif 'cbe_edit_mutation' not in self.data[guide_id] and 'abe_edit_mutation' in self.data[guide_id]:
						base_editor_variant = "target"
					elif self.data[guide_id]['abe_edit_mutation'] != "" and self.data[guide_id]['cbe_edit_mutation'] != "":
						base_editor_variant = "target"
					elif self.data[guide_id]['abe_edit_mutation'] != "" and self.data[guide_id]['cbe_edit_mutation'] == "":
						base_editor_variant = "target"
					elif self.data[guide_id]['cbe_edit_mutation'] != "" and self.data[guide_id]['abe_edit_mutation'] == "":
						base_editor_variant = "target"
					"""

					if base_editor_variant not in base_editor_variants_distribution:
						base_editor_variants_distribution[base_editor_variant] = []

					try:
						base_editor_variants_distribution[base_editor_variant].append(self.data[guide_id]['limma_scores']['log_fold_change'])
						
						#mean_row.append("%1.2f"%(sum(mean_proportion)/len(mean_proportion)))
					except:
						pass
				except:
					logging.error(guide_id)
					raise
			
			####################################################

			
			sorted_vals = {}
			sorted_keys_label_iter = 0
			for sorted_key in sorted_keys:
				if sorted_key in base_editor_variants_distribution:
					sorted_vals[sorted_key] = base_editor_variants_distribution[sorted_key]
					sorted_keys_label[sorted_keys_label_iter] = sorted_keys_label[sorted_keys_label_iter] + " n=" + str(len(base_editor_variants_distribution[sorted_key]))
					sorted_keys_label_iter +=1

			fig, ax = plt.subplots()
			custom_params = {"axes.spines.right": False, "axes.spines.top": False}
			sns.set_theme(style="ticks", rc=custom_params)
			
			g = sns.boxenplot(data=list(sorted_vals.values()),showfliers = True)#, orient='h',scale='exponential')#,order=)
			g.set_axis_labels("guide class", "log fold change", labelpad=10)
			plt.xticks(plt.xticks()[0], sorted_keys_label)
			
			g.figure.set_size_inches(10, 10)

			g.set_xlabel("guide class",fontsize=10, weight="bold")
			g.set_ylabel("log fold change",fontsize=10, weight="bold")
			
			plt.tight_layout()
			plt.savefig(box_png_path,dpi=300)
			
			plt.clf()
			logging.debug("Writing: " + box_png_path)

			###--------------------------###
		except:
			logging.error(box_png_path)

	###--------------------###

	def plot_control_roc_curves(self,positive_control_labels=['Essential splice site'],negative_control_labels=['Non-targeting control','Intergenic control']):
		self.make_plot_directory()

		graph_data = {}
		graph_keys = []
		graph_labels = {}
		graph_colours = {}

		colours = ['#777','#C28C69','#6F9D72','#AA6163','#8b2e57','#25e8b7'] 
		

		screen_label = self.options['screen_name'] + ' Day ' + str(self.options['compare_time_point'])
		
		
		guide_ids = list(self.data.keys())
		guide_ids.sort()

		graph_labels[screen_label] = "Essential Splice Sites"
		graph_colours[screen_label] = colours[len(graph_colours)]
		graph_data[screen_label] = {
			"status":[],
			"scores":[]
		}

		for positive_control_label in positive_control_labels:
			graph_keys.append(positive_control_label)
			graph_labels[positive_control_label] = positive_control_label
			graph_colours[positive_control_label] = colours[len(graph_colours)]
			graph_data[positive_control_label] = {
				"status":[],
				"scores":[]
			}


		for guide_id in guide_ids:
			try:
				if self.data[guide_id]['control']:
					status = False
					if self.data[guide_id]['control_type'] in positive_control_labels:
						status = True
					elif self.data[guide_id]['control_type'] in negative_control_labels:
						status = False
					else:
						continue

					fc = -self.data[guide_id]['limma_scores']['log_fold_change']

					graph_data[screen_label]['status'].append(status)
					graph_data[screen_label]["scores"].append(fc)

					if self.data[guide_id]['control_type'] in negative_control_labels:
						for positive_control_label in positive_control_labels:
							graph_data[positive_control_label]['status'].append(status)
							graph_data[positive_control_label]["scores"].append(fc)

					if self.data[guide_id]['control_type'] in positive_control_labels:
						graph_data[self.data[guide_id]['control_type']]['status'].append(status)
						graph_data[self.data[guide_id]['control_type']]["scores"].append(fc)
			except:
				pass
				
		if len(graph_data) > 2:
			graph_keys.append(screen_label)
				
		roc_path_png =  os.path.join(self.options["plot_dir_path"],self.options['screen_name'] + '.' + str(self.options['compare_time_point']) + '.roc.png')
		
		logging.debug("Writing " + roc_path_png)
		utilities_plots.plot_roc_curve(roc_path_png,graph_data,graph_keys,graph_labels,graph_colours,show_thresholds=True,show_youdens=True)
		plt.clf()

	##################################################################################################

		
			
	