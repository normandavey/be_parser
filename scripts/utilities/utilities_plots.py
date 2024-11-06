import math,sys,traceback,pprint,os,hashlib
import zlib,json, base64
import warnings
warnings.filterwarnings("ignore", message="the sets module is deprecated")

import logging
import utilities_error


###################################################################
### Errors
###################################################################

def plot_roc_curve(roc_path_png,graph_data,graph_keys,graph_labels,graph_colours,show_thresholds=False,show_youdens=False,normalise_negative_attributes_list=[]):
	"""	Args:
		roc_path_png (str): out path for the graph png
		graph_data (dict of lists): each attrbutes includes lists of scores and status (for status of a particular list attribute ["True":"False"])
		graph_keys (list): list of labels to plot from graph_data
		graph_labels (dict): label dictionary
		graph_colours (dict): colour dictionary
		show_thresholds (bool, optional): Defaults to False.
		show_youdens (bool, optional): Defaults to False.
		normalise_negative_attributes_list (list): list of negative attributes that are to be normalised by x*-1 
	"""
	try:
		from sklearn.metrics import roc_curve  
		from sklearn.metrics import roc_auc_score  
	
		import matplotlib.pyplot as plt
		import seaborn as sns
		import pandas as pd
		import numpy as np
		
		plt.rc("axes.spines", top=False, right=False)
		
		sns.set_style("ticks")

		rc={'axes.labelsize': 18}
		plt.rcParams.update(**rc)

		
		in_graph = graph_keys
		in_graph_colours = graph_colours
		use_cols = graph_keys
		
		attribute_statistics = {}
		
		for roc_output_type in ["single"]:
			
			sns.set_style("ticks")
			
			if roc_output_type == "single":
				plt.rcParams['figure.figsize'] = (5, 5)
				fig, axs = plt.subplots(1, 1)
			else:
				plt.rcParams['figure.figsize'] = (5*len(in_graph),5)
				fig, axs = plt.subplots(1, len(in_graph))
		
			plot_number = 0
			youdens_j_rows = [
					"\t".join(["AUC","j","type"])
				]

			for graph_type in use_cols:
				logging.debug("Plotting: " + graph_type)
				try:

					boxplot_df = pd.DataFrame(graph_data[graph_type])
					
					if graph_type not in attribute_statistics:
						attribute_statistics[graph_type] = {}

					data_points = boxplot_df['scores']

					if graph_type in normalise_negative_attributes_list:
						data_points = -data_points
					
					roc_labels = [] 
					for classification in boxplot_df['status']:
						roc_labels.append(classification == True)
					
					roc_pssm_scores = data_points
				
					roc_labels = np.array(roc_labels)
					roc_pssm_scores = np.array(roc_pssm_scores)

					auc = roc_auc_score(roc_labels, roc_pssm_scores)
					print(auc)
					fpr, tpr, thresholds = roc_curve(roc_labels, roc_pssm_scores,pos_label=1)

					if show_thresholds:
						rows = []
						rows.append("\t".join(["fpr","tpr","delta","threshold"]))
						for threshold_i in range(0,len(fpr)):
							rows.append("\t".join(["%1.3f"%fpr[threshold_i], "\t","%1.3f"%tpr[threshold_i],"\t","%1.3f"%(fpr[threshold_i]-tpr[threshold_i]), "\t","%1.3f"%thresholds[threshold_i],graph_type]))

						print("\n".join(rows))

				
					youdens_j = thresholds[np.argmax(tpr - fpr)]

					youdens_j_rows.append("\t".join([
						"%1.3f"%auc,
						"%1.3f"%youdens_j,
						graph_type
					]))

					attribute_statistics[graph_type]['auc'] = auc
					attribute_statistics[graph_type]['youdens_j'] = youdens_j

					if roc_output_type == "single":
						axs.plot(fpr, tpr, color=in_graph_colours[graph_type], label=graph_labels[graph_type] + ' (AUC:%.2f' % auc + ')') # - J:%.2f' % youdens_j + ')
					else:
						axs[plot_number].plot(fpr, tpr, color=in_graph_colours[graph_type], label=graph_type + ' (AUC:%.2f' % auc + ' - J:%.2f' % youdens_j + ')')
						axs[plot_number].plot([0, 1], [0, 1], color='#333333', linestyle='--',linewidth=0.5)
						axs[plot_number].set_xlabel('False Positive Rate',fontsize=14, weight="bold")
						axs[plot_number].set_ylabel('True Positive Rate',fontsize=14, weight="bold")
						axs[plot_number].set_title(graph_type,fontsize=18, weight="bold")
						axs[plot_number].legend(loc='lower right',fontsize=10)

					plot_number += 1
				except:
					utilities_error.printError()

			if roc_output_type == "single":
				axs.plot([0, 1], [0, 1], color='#333333', linestyle='--')
				axs.set_xlabel('False Positive Rate',fontsize=10, weight="bold")
				axs.set_ylabel('True Positive Rate',fontsize=10, weight="bold")
				axs.set_title("",fontsize=18, weight="bold")
				axs.legend(loc='lower right',fontsize=9)
		
			sns.despine()

			if show_youdens:
				print("\n".join(youdens_j_rows))

			plt.savefig(roc_path_png,dpi = 300)
			logging.debug("Writing Figure:" + roc_path_png)

		return attribute_statistics
	except:
		logging.error("Plotting failed")
		raise