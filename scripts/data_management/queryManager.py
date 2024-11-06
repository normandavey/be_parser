import sys

if sys.version_info[0] != 3 or sys.version_info[1] < 5:
    print("This script requires Python version 3.5 or greater")
    sys.exit(1)

import os
import inspect
file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"../"))

from queueManager import queueManager

class queryManager:

	##------------------------------------------------------------------##

	def __init__(self):
		self.options = {
			"dataset_type":"",
			"dataset_type":"",
			"is_superuser":True,
			"test_all":False,
			"rest_api":"http://disorder.icr.ac.uk/webservices/rest/get/"
			}

		self.dataset_types = [
			"cancer",
			"ptms",
			"database_access",
			"depmap",
			"exons",
			"ensembl",
			"alphafold",
			"comparimotif",
			"consensus",
			"ontology",
			"proppd",
			"druggability",
			"uniprot",
			"annotator",
			"idi_annotator",
			"pdb",
			"dssp",
			"interaction",
			"motif_alignment",
			"mutations",
			"pfam",
			"pubmed",
			"accessibility",
			"architecture",
			"cell_cycle",
			"interface",
			"elm",
			"pssm",
			"motif",
			"pocket",
			"hmm",
			"upload_access",
			"e3s",
			"evolution",
			"complex",
		]

	##------------------------------------------------------------------##

	def queryFarmer(self,task="run_main"):
		sys_argv = sys.argv
			
		if self.options["dataset_type"] == "ensembl":
			import ensemblManager
			dataManagerObj = ensemblManager.ensemblManager()

		elif self.options["dataset_type"] == "database_access":
			import databaseAccessManager
			dataManagerObj = databaseAccessManager.databaseAccessManager()

		elif self.options["dataset_type"] == "depmap":
			import depmapManager
			dataManagerObj = depmapManager.depmapManager()

		elif self.options["dataset_type"] == "accessibility":
			import accessibilityManager
			dataManagerObj = accessibilityManager.accessibiltyManager()

		elif self.options["dataset_type"] == "dssp":
			import dsspManager
			dataManagerObj = dsspManager.dsspManager()

		elif self.options["dataset_type"] == "architecture":
			import architectureManager
			dataManagerObj = architectureManager.architectureManager()

		elif self.options["dataset_type"] == "interface":
			import interfaceManager
			dataManagerObj = interfaceManager.interfaceManager()

		elif self.options["dataset_type"] == "pfam":
			import pfamManager
			dataManagerObj = pfamManager.pfamManager()

		elif self.options["dataset_type"] == "pssm":
			import pssmManager
			dataManagerObj = pssmManager.pssmManager()

		elif self.options["dataset_type"] == "uniprot":
			import uniprotManager
			dataManagerObj = uniprotManager.uniprotManager()

		elif self.options["dataset_type"] == "alphafold":
			import alphafoldManager
			dataManagerObj = alphafoldManager.alphafoldManager()
		else:
			return {"status":"Error","error_type":"Analysis type [" + self.options['dataset_type'] + "] does not exist."}

		## PURGE COMMANDLINE
		sys.argv = sys.argv[:1]

		dataManagerObj.options.update(self.options)
		dataManagerObj.options['debug'] = True

		if task == "check_options":
			response = dataManagerObj.setupTask(self.options)
		elif task == "run_main":
			response = dataManagerObj.main()
		elif task == "check_queue_requirement":
			response = self.options["task"] in dataManagerObj.task_skip_queue
		else:
			response = {"status":"Error","error_type":"Task [" + task + "] does not exist."}

		sys.argv = sys_argv

		return response

	##------------------------------------------------------------------##

	def processCommandline(self):

		### IF CALLED FROM COMMANDLINE ###

		remove_sys_argv_values = []

		if "--dataset_type" in sys.argv:
			dataset_type_index = sys.argv.index("--dataset_type")

			if dataset_type_index != -1:
				self.options['dataset_type'] = sys.argv[dataset_type_index+1]
				remove_sys_argv_values += [dataset_type_index, dataset_type_index+1]

		remove_sys_argv_values.sort()
		sys_argv_updates = []

		for i in range(0,len(sys.argv)):
			if i not in remove_sys_argv_values:
				sys_argv_updates.append(sys.argv[i])

		sys.argv = sys_argv_updates
		return sys_argv_updates

	#########################################################
	#########################################################
	#########################################################
	#########################################################

	def main(self):
		### IF CALLED FROM COMMANDLINE ###
		sys_argv_updates = self.processCommandline()

		#-------------#

		if self.options['dataset_type'] not in self.dataset_types:
			status = {
				"status":"Error",
				"error_type":"Dataset type '" + self.options['dataset_type'] + "' is not a valid option. Please set --dataset_type from " + ",".join(self.dataset_types),
				"args":sys.argv,
				"options":self.options
				}

			if 'verbose' in self.options:
				if self.options['verbose'] == True:
					print(status)
			else:
				print(status)

			return status

		#########################################################
	
		if "queue" in self.options:
			if self.options["server"]:
				queueManagerObj = queueManager()
				queueManagerObj.options.update(self.options)
				status = queueManagerObj.queue_director()
				return status
			else:
				return {
				"status":"Completed",
				"options":self.options,
				"data":self.queryFarmer(),
				}
		else:
			response = self.queryFarmer()
			return response

	#########################################################
	#########################################################
	#########################################################
	#########################################################

if __name__ == "__main__":
	queryManagerObj = queryManager()
	queryManagerObj.main()
