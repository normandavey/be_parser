import os
import sys
import re
import json
import copy

import inspect
file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"../"))
sys.path.append(os.path.join(file_path,"../utilities/"))

import utilities_error
import config_reader

#-----
import logging
logger = logging.getLogger(__name__)
#-----

class dataManager:

	##------------------------------------------------------------------##

	def __init__(self):
		self.options = {}
		self.outfile = None 
		self.default_task_options = {}
		self.allowed_options = []
		self.allowed_options_admin = []
		self.task_options = []
		self.task_options_admin = []
		self.task_skip_queue = ["help"]
		self.required = {}
		self.required_type = {}
		self.required_format= {
			'accession':"\A([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\-{0,1}[0-9]*\Z",
			'pmid':'^[0-9]+$',
			'pdb_id':'^[A-Z0-9]{4}$',
			'pfam_id':'^PF[0-9]{5}$',
			'taxon_id':'^[0-9]{1,10}$'
		}
		self.required_valid_options_list = {}
		self.descriptions = {}
		self.list_options = {}
		self.convert_type = {
			'pmid':"string",
		}

		self.options.update(config_reader.load_configeration_options(sections=["general"]))
	
		self.setup_data()

	##------------------------------------------------------------------##

	def setupLogger(self):
		if 'debug' not in self.options: self.options['debug'] = True

		if self.options['debug']:
			logging_level = logging.DEBUG
		elif self.options['verbose']:
			logging_level = logging.INFO
		else:
			logging_level = logging.ERROR
		
		if 'logfile' not in self.options: self.options['logfile'] = False
		elif self.options['logfile'] == "": self.options['logfile'] = False
		elif isinstance(self.options['logfile'],(str)): pass
		elif self.options['logfile']: self.options['logfile'] = os.path.join(file_path,'data_manager.log')
		
		# create console handler with a higher log level
		formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(funcName)s] %(message)s')
		ch = logging.StreamHandler()
		ch.setLevel(logging_level)
		ch.setFormatter(formatter)

		if self.options['logfile'] != False:
			# create file handler which logs even debug messages
			formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(funcName)s] %(message)s')
			fh = logging.FileHandler(self.options['logfile'])
			fh.setLevel(logging.DEBUG)
			fh.setFormatter(formatter)

			logging.basicConfig(
				format='%(asctime)s %(levelname)-8s [%(funcName)s] %(message)s',
				datefmt='%H:%M:%S',
				level=logging.DEBUG,
				handlers=[ch,fh]
				)
			
			logger = logging.getLogger(__name__)
			logger.info("Writing log: " + self.options['logfile'])
		
		else:
			logging.basicConfig(
				format='%(asctime)s %(levelname)-8s [%(funcName)s] %(message)s',
				datefmt='%H:%M:%S',
				level=logging_level
				)

			logger = logging.getLogger(__name__)

	##------------------------------------------------------------------##

	def main(self):
		self.setupLogger()
		
		if self.options['task'] == "test" or self.options['task'] == "test_rest_api":
			from timeit import default_timer as timer

			response = {"test_results":{}}
			tasks = self.task_options

			for task in tasks:
				start = timer()
				if task == "help": continue
				
				row = ["%-20s"%self.options["dataset_type"],"%-40s"%task]

				task_options = copy.deepcopy(self.test_options)
				task_options['task'] = task
				task_options['print'] = True
				
				try:
					if self.options['task'] == "test":
						test_query_response = self.runTask(task_options)
						
						if test_query_response['status'] == "Error":
							response['test_results'][task] = test_query_response
							row += [response['test_results'][task]['status'],"-"]
							print("ERROR - @main",response['test_results'][task]['error_type'])
						else:
							end = timer()
							time_elapsed = end - start
							if test_query_response['data'] in [True,False]:
								row += [test_query_response['status'],str(1),"%1.3fs"%time_elapsed]
								response['test_results'][task] = {
									'status':test_query_response['status'],
									'result_length':1
								}
							else:
								row += [test_query_response['status'],str(len(test_query_response['data'])),"%1.3fs"%time_elapsed]
								response['test_results'][task] = {
									'status':test_query_response['status'],
									'result_length':len(test_query_response['data'])
								}
					if self.options['task'] == "test_rest_api":
						server_url = self.options["rest_api"] 
						dataset_type = self.options["dataset_type"] 
						api_url = os.path.join(server_url,dataset_type) + "/"

						params = self.test_options
						params['task'] = task

						sys.path.append(os.path.join(file_path,"../rest/"))
						import restDownloader
						
						start = timer()
						api_response = restDownloader.run_query_with_waiter(api_url,params,request_method="GET")

						try:
							end = timer()
							time_elapsed = end - start

							size = len(json.dumps(api_response['data']).encode("utf-8"))/1024

							if api_response['data'] in [True,False]:
								row += [api_response['status'],str(1),api_url,"%1.3fs"%time_elapsed,"%1.0f"%size+ "kb"]
								response['test_results'][task] = {
									'status':api_response['status'],
									'result_length':1
								}
							else:
								row += [api_response['status'],str(len(api_response['data'])),api_url,"%1.3fs"%time_elapsed,"%1.0f"%size+ "kb"]
								response['test_results'][task] = {
									'status':api_response['status'],
									'result_length':len(api_response['data'])
								}
							
							if 'status' in api_response['data']:
								if api_response['data']['status'] == "Error":
									import pprint
									pprint.pprint(api_response['data'])

							#import pprint
							#pprint.pprint(api_response['data'])

						except:
							print("ERROR - api_response",api_response)
							raise

					print("\t".join(row))
				
				except Exception as e:
					response['test_results'][task] = {
						'status':"Error",
						'error_type':str(e)
					}

					print("\t".join(row))
					print(response['test_results'][task]['status'] + "-" + str(e))
					
			self.makeOutfile(response)

			if self.options["verbose"]:
				import pprint
				pprint.pprint(response)

				
		else:
			if 'outfile' in self.options and self.options['outfile']  != "":
				self.outfile = self.options['outfile'] 
				del self.options['outfile'] 
		
			task_options = {
				'task':self.options['task']
			}
			response = self.runTask(task_options)

		return response

	##------------------------------------------------------------------##

	def setupTask(self,task_options):

		###----######----######----###

		task_options = self.addCommandLineOption(task_options)
		task_options = self.addDefaultOption(task_options)
		task_options = self.cleanOption(task_options)
		task_options = self.convertType(task_options)
		task_options = self.checkListOptions(task_options)

		if "is_superuser" not in task_options:
			task_options["is_superuser"] = True
		
		###----######----######----###

		if "task" in task_options:
			if task_options["task"] == "help":
				help = self.formatHelp()

				if self.options["verbose"]:
					import pprint
					pprint.pprint(help)

				return {
				"status": "Sucess",
				"options":task_options,
				"data":help
				}

		###----######----######----###
		#	 check_input options	#
		###----######----######----###

		status = self.checkOptions(task_options)

		task_details = status
		task_details["options"] = task_options

		return task_details

	##------------------------------------------------------------------##

	def makeOutfile(self,content):
		if self.outfile not in [None,""]:
			with open(self.outfile ,"w") as outfile:
				json.dump(content, outfile)

	##------------------------------------------------------------------##

	def runTask(self,task_options):
		try:
			if self.options["verbose"] and self.options['task'] not in ["test","test_rest_api"]: logger.info("Task:" + task_options["task"])

			###----###

			if task_options["task"] == "help":
				help = self.formatHelp()

				if self.options["verbose"]:
					import pprint
					pprint.pprint(help)

				out_json = {
				"status": "Success",
				"options":task_options,
				"data":help
				}

				if self.options['task'] not in ["test","test_rest_api"]:
					self.makeOutfile(out_json)

				return out_json

			###----###

			if 'outfile' in task_options:
				del task_options['outfile'] 
				
			task_details = self.setupTask(task_options)

			status = task_details
			task_options = task_details['options']
		
			###----###
			if status['status'] != "Error":
				data = self.getData(status,task_options)
				
				if data == None:
					data = {}
					status ={
						'status':"Error",
						"error_type":"Nonetype returned"
						}
				elif isinstance(data,(bool)):
					pass
				else:
					if 'data' in data:
						data = data['data']

					if 'status' in data:
						if data['status'] == "Error":
							status = data

			###----###
			
			if status['status'] == "Error":
				if self.options['task'] not in ["test","test_rest_api"]:
					if self.options["verbose"]:
						import pprint
						pprint.pprint(task_details)
						pprint.pprint(status)

					self.makeOutfile(status)
				
					if status['status'] == 'Error':
						logger.error(status['error_type'])

				return status
			else:
				if self.options['task'] not in ["test","test_rest_api"]:
					if self.options["verbose"]:
						import pprint
						pprint.pprint(data,width=200, compact=True)

					self.makeOutfile(data)
				
				###----###
				
				return {
				"status": status['status'],
				"options":task_options,
				"data":data
				}
		except Exception as e:
			
			status = {
			"status":"Error",
			"error_type":utilities_error.getError()
			}

			if self.options['task'] not in ["test","test_rest_api"]:
				self.makeOutfile(status)
			
				if self.options["verbose"]:
					import pprint
					pprint.pprint(status)
			
			return status

	##------------------------------------------------------------------##

	def addDefaultOption(self,task_options):
		for default_task_option in self.default_task_options:
			if default_task_option not in task_options:
				task_options[default_task_option] = self.default_task_options[default_task_option]

		return task_options

	##------------------------------------------------------------------##

	def addCommandLineOption(self,task_options):
		for commandline_task_option in self.options:
			if commandline_task_option not in task_options and commandline_task_option in self.allowed_options:
				task_options[commandline_task_option] = self.options[commandline_task_option]

			if commandline_task_option not in task_options and commandline_task_option in self.allowed_options_admin and self.options["is_superuser"]:
				task_options[commandline_task_option] = self.options[commandline_task_option]
			
			if commandline_task_option in self.allowed_options_admin and not self.options["is_superuser"]:
				logger.debug(commandline_task_option + " option not allowed - superuser only")

		return task_options

	##------------------------------------------------------------------##

	def cleanOption(self,task_options):
		task_options_clean = {}
		for task_option in task_options:
			if task_options[task_option] != None:
				task_options_clean[task_option]= task_options[task_option]

		return task_options_clean

	##------------------------------------------------------------------##

	def checkListOptions(self,task_options):
		task_options_clean = {}

		for task_option in task_options:
			if task_option in self.list_options:
				if isinstance(task_options[task_option], (list)):
					task_options_clean[task_option] = task_options[task_option]
				else:
					task_options_clean[task_option] = str(task_options[task_option]).split(self.list_options[task_option])
			else:
				task_options_clean[task_option] = task_options[task_option]
		
		return task_options_clean

	##------------------------------------------------------------------##

	def convertType(self,task_options):
		task_options_clean = {}

		for task_option in task_options:
			if task_options[task_option] == "None" or task_options[task_option] == "none": 
				task_options_clean[task_option] = None
			elif task_options[task_option] == "False" or task_options[task_option] == "false": 
				task_options_clean[task_option] = False
			elif task_options[task_option] == "True" or task_options[task_option] == "true": 
				task_options_clean[task_option] = True
			elif task_option in self.convert_type:
				
				if isinstance(task_options[task_option], (list)):
					tmp_list = {}
					for value in task_options[task_option]:
						if self.convert_type[task_option] == "int":
							tmp_list.append(int(value))
						elif self.convert_type[task_option] == "string":
							tmp_list.append(str(value))
						else:
							tmp_list.append(value)

					task_options_clean[task_option] = tmp_list
				else:
					if self.convert_type[task_option] == "int":
						task_options_clean[task_option] = int(task_options[task_option])
					elif self.convert_type[task_option] == "string":
						task_options_clean[task_option] = str(task_options[task_option])
					else:
						task_options_clean[task_option] = task_options[task_option]
			else:
				task_options_clean[task_option] = task_options[task_option]

		return task_options_clean

	##------------------------------------------------------------------##

	def formatHelp(self):
		self.allowed_options.sort()

		help_dict = {"options":{},"tasks":{}}

		for allowed_option in self.allowed_options:

			required_tasks = []
			required_type = ""
			required_format = ""
			required_valid_options_list = []
			list_delimiter = "None"
			decription = "",
			default_options = ""
			set_options = ""

			for task_option in self.task_options:
				if task_option in self.required:
					if task_option != "all" and allowed_option in self.required[task_option]:
						required_tasks.append(task_option)

			if allowed_option in self.required_type:
				required_type = self.required_type[allowed_option]

			if allowed_option in self.required_format:
				required_format = self.required_format[allowed_option]

			if allowed_option in self.required_valid_options_list:
				required_valid_options_list = self.required_valid_options_list[allowed_option]

			if allowed_option in self.descriptions:
				decription = self.descriptions[allowed_option]

			if allowed_option in self.list_options:
				list_delimiter = self.list_options[allowed_option]

			if allowed_option in self.default_task_options:
				default_options = self.default_task_options[allowed_option]

			if allowed_option in self.options:
				set_options = self.options[allowed_option]

			help_dict["options"][allowed_option] = {
				"required_global":allowed_option in self.required['all'],
				"required_tasks":required_tasks,
				"required_type":required_type,
				"required_format":required_format,
				"required_valid_options_list":required_valid_options_list,
				"list_delimiter":list_delimiter,
				"decription":decription,
				"default_options":default_options,
				"set_options":set_options
			}

		for task in self.task_options:
			help_dict["tasks"][task] = {}

			if task in self.required:
				help_dict["tasks"][task]['required_options'] = self.required[task]


		return help_dict

	##------------------------------------------------------------------##

	def checkOptionsTypes(self,task_options):
		failed = []
		status = {"status":"Success","task_options":task_options}

		for task_option in task_options:
			if task_option in self.required_type:
				required_type = self.required_type[task_option]

				if required_type == "int":
					if not isinstance(task_options[task_option], (int)):
						try:
							status["task_options"][task_option] = int(task_options[task_option])
						except:
							failed.append(task_option)
				if required_type == "bool":
					if not isinstance(task_options[task_option], (bool)):
						try:
							status["task_options"][task_option] = bool(task_options[task_option])
						except:
							failed.append(task_option)
				if required_type == "float":
					if not isinstance(task_options[task_option], (float)):
						try:
							status["task_options"][task_option] = float(task_options[task_option])
						except:
							failed.append(task_option)
				if required_type == "str":
					if not isinstance(task_options[task_option], (str)):
						try:
							status["task_options"][task_option] = str(task_options[task_option])
						except:
							failed.append(task_option)
				if required_type == "list":
					if not isinstance(task_options[task_option], (list)):
						try:
							status["task_options"][task_option] = list(task_options[task_option])
						except:
							failed.append(task_option)
				if required_type == "dict":
					if not isinstance(task_options[task_option], (dict)):
						try:
							status["task_options"][task_option] = dict(task_options[task_option])
						except:
							failed.append(task_option)
		if len(failed) > 0:
			status["status"] = "Error"
			status["error_type"] = "Type error: " + ",".join(failed)
			status["help"] = self.required_type

		return status

	##------------------------------------------------------------------##

	def checkValidOptionLists(self,task_options):
		failed = []
		status = {"status":"Success"}

		for task_option in task_options:
			if task_option in self.required_valid_options_list:
				if task_options[task_option] not in self.required_valid_options_list[task_option]:
					failed.append(task_option)

		if len(failed) > 0:
			status["status"] = "Error"
			status["error_type"] = "Not a valid option error: " + ",".join(failed)
			status["help"] = self.required_valid_options_list

		return status

	##------------------------------------------------------------------##

	def checkOptionsFormat(self,task_options):
		failed = []
		status = {"status":"Success"}

		for task_option in task_options:
			if task_option in self.required_format:
				if task_options[task_option] not in [None,"None"]:
					if isinstance(task_options[task_option], (list)):
						
						for val in task_options[task_option]:
							if not re.match(self.required_format[task_option], val):
								failed.append(task_option)
					else:
						try:
							if not re.match(self.required_format[task_option], task_options[task_option]):
								failed.append(task_option)
						except:
							logger.error(task_option)
	
		if len(failed) > 0:
			status["status"] = "Error"
			status["error_type"] = "Format error: " + task_option + " " + ",".join(failed) + " - " + str(task_options[task_option])
			status["help"] = self.required_format

		return status

	##------------------------------------------------------------------##

	def checkTasks(self,task_options):
		status = {}
		
		if "task" not in task_options:
			status["status"] = "Error"
			status["error_type"] = "No task selected: " + ", ".join(self.task_options )
		elif task_options["task"] in self.task_options:
			status = {"status":"Success"}
		else:
			if task_options["task"] not in self.task_options + self.task_options_admin:
				status["status"] = "Error"
				status["error_type"] = "Task not in task list: " + ", ".join(self.task_options )
			else:
				status["status"] = "Error"
				status["error_type"] = "Not sure how you ended up here. Task not in default task list"
	
		if task_options["task"] in self.task_options_admin:
			if 'is_superuser' in task_options:
				if task_options['is_superuser'] :
					status = {"status":"Success","message":"Restricted task - running as superuser or local user"}
				else:
					status = {"status":"Error","error_type":"Restricted task - run as superuser or local user"}
			else:
				status = {"status":"Error","error_type":"Restricted task - run as superuser or local user"}

		return status

	##------------------------------------------------------------------##

	def checkRequiredOptions(self,task_options):
		status = {"status":"Success"}
		for task_option in self.required['all']:
			if task_option not in task_options:
				status["status"] = "Error"
				status["error_type"] = "Global option " + task_option + " not supplied"

		if task_options["task"] in self.required:
			for task_option in self.required[task_options["task"]]:
				if task_option not in task_options:
					status["status"] = "Error"
					status["error_type"] = "Required task option " + task_option + " not supplied in " + task_options["task"]

		return status

	##------------------------------------------------------------------##

	def checkOptions(self,task_options):
		status = {"status":"Success"}
		if status["status"] != "Error":
			status = self.checkTasks(task_options)

		if status["status"] != "Error":
			status = self.checkRequiredOptions(task_options)

		if status["status"] != "Error":
			status = self.checkOptionsTypes(task_options)
			task_options = status['task_options']
		
		if status["status"] != "Error":
			status = self.checkOptionsFormat(task_options)
		
		if status["status"] != "Error":
			status = self.checkValidOptionLists(task_options)

		return status
