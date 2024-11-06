import os
import sys
import re
import json
import random
import string 
#import md5
import hashlib
import json
import traceback
import time
import stat

import inspect

file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"../"))

import config_reader

class queueManager:

	##------------------------------------------------------------------##

	def __init__(self):
		self.options = {}
		self.user_jobs = {}

		self.paths = config_reader.load_configeration_options(sections=["general"])
		
		self.removed_cmd_options = [
			"remote_addr",
			"tool", 
			"queue_tasks",
			"test_all",
			"server",
			"queue",
			"cached",
			"host",
			"logged_in",
			"job_id",
			"add_header"
			]
	##------------------------------------------------------------------##

	####################################################################################
	####################################################################################

	def populate_user_jobs(self):

		if not os.path.exists(os.path.join(self.paths["data_path"],"webservices")):
			os.mkdir(os.path.join(self.paths["data_path"],"webservices"))

		for tool in os.listdir(os.path.join(self.paths["data_path"],"webservices")):

			if not os.path.isdir(os.path.join(self.paths["data_path"],"webservices",tool)): continue

			for file in os.listdir(os.path.join(self.paths["data_path"],"webservices",tool)):
				try:
					if file.split(".")[-1] == "json":
						file_age = time.time() - os.stat(os.path.join(self.paths["data_path"],"webservices",tool,file))[stat.ST_MTIME]

						#if file_age/60/60/24 > 14: continue

						with open(os.path.join(self.paths["data_path"],"webservices",tool,file)) as json_file:
							data = json.load(json_file)

							if 'remote_addr' in data:
								remote_addr = data['remote_addr']
								del data['remote_addr']
							else:
								remote_addr = "n/a"

							if remote_addr not in self.user_jobs:
								self.user_jobs[remote_addr] = {}

							job_id = data['job_id']
							del data['job_id']

							self.user_jobs[remote_addr][job_id] = data

				except Exception as e:
					print("ERROR - @ populate_user_jobs",e)

	####################################################################################

	def params_to_hash(self,params):
		param_str = ""
		param_keys = list(params.keys())
		param_keys.sort()

		for param in param_keys:
			if param not in self.removed_cmd_options:
				param_str += param+"="+str(params[param])+"\n"
		
		hash = str(hashlib.md5(param_str.encode()).hexdigest())
		return hash

	####################################################################################

	def create_jobfile(self,job_options):
		if not os.path.exists(os.path.join(self.paths["data_path"],"webservices")):
			os.mkdir(os.path.join(self.paths["data_path"],"webservices"))

		if not os.path.exists(os.path.join(self.paths["data_path"],"webservices",job_options["tool"])):
			os.mkdir(os.path.join(self.paths["data_path"],"webservices",job_options["tool"]))

		results_path = os.path.join(self.paths["data_path"],"webservices",job_options["tool"])

		pbs_file = os.path.join(results_path,job_options['job_id'] + ".sh")
		log_file = os.path.join(results_path,job_options['job_id'] + ".log")
		output_file = os.path.join(results_path,job_options['job_id'] + ".output.txt")
		error_file = os.path.join(results_path,job_options['job_id'] + ".error.txt")
		options_file = os.path.join(results_path,job_options['job_id'] + ".options.json")

		with open(options_file, 'w') as outfile:
			json.dump(job_options, outfile)

		cmd_options_str = " ".join([ "--" + option_tag + ' "' + str(job_options['cmd_options'][option_tag]) + '"' for option_tag in job_options['cmd_options']] )
		print(cmd_options_str)
		pbs = """#!/bin/bash
		#PBS -l nodes=1:ppn=1
		#PBS -o """  + output_file + """
		#PBS -e """  + error_file + """
		#PBS -d """ + results_path + """
		#PBS -N """ + job_options['job_id'] + """

		""" + job_options['cmd'] + """ """ + cmd_options_str + """ > """ + log_file + """
		"""

		pbs = pbs.replace("\t","")
		open(pbs_file,"w").write(pbs)

		job_status = os.popen("/sbin/runuser -l  submitter -c 'qsub -l walltime=100:00:00 " + pbs_file + "'").read().strip()

		job_status_output = {
		"job_options":job_options,
		"job_status":job_status,
		"cmd":job_options['cmd'] + " " + cmd_options_str
		}

		return job_status_output

	####################################################################################

	def show_jobs(self):
		qstat_flags = {"E":"Job is exiting after having run",
		"H":"Job is held",
		"Q":"Job is queued, eligible to be run or routed",
		"R":"Job is running",
		"T":"Job is being moved to new location",
		"W":"Job is waiting for its execution time",
		"S":"Job is suspend"}

		qstat_output = os.popen("runuser -l submitter -c 'qstat -f'").read()

		qstat_data = {"jobs":{},"status":{},"job_count":0}

		qstat_str = qstat_output.strip().split("Job Id: ")

		for job_data in qstat_str:
			try:
				job_data_line = job_data.strip().split("\n")

				if len(job_data_line) > 1:
					tmp_dict = {}
					for line in job_data_line[1:]:
						line_bits = line.strip().split(" = ")
						if len(line_bits) > 1:
							tmp_dict[line_bits[0].strip()] = line_bits[1].strip()

					job_id = tmp_dict['Job_Name']

					tmp_dict['job_state_full'] = qstat_flags[tmp_dict['job_state']]
					qstat_data["jobs"][job_id] = tmp_dict

					if qstat_flags[qstat_data["jobs"][job_id]['job_state']] not in qstat_data["status"]:
						qstat_data["status"][qstat_flags[qstat_data["jobs"][job_id]['job_state']]] = []

					qstat_data["status"][qstat_flags[qstat_data["jobs"][job_id]['job_state']]].append(job_id)

			except Exception as e:
				exc_type, exc_value, exc_tb = sys.exc_info()
				error_traceback = traceback.format_exception(exc_type, exc_value, exc_tb)

				qstat_data['error'] = {"type":error_traceback,"job_data":job_data}

		qstat_data["job_count"] = len(qstat_data["jobs"])

		return qstat_data

	def job_status(self,job_identifier,dataset_type):
		data_path = os.path.join(self.paths['data_path'],"webservices",dataset_type)
		results_file = os.path.join(data_path,job_identifier + ".results.json")
		options_file = os.path.join(data_path,job_identifier + ".options.json")
		
		qstat_data_response = {}
		qstat_data_response['job_state'] = "U"
		qstat_data_response['job_state_full'] = "Unknown Failure"
		qstat_data_response["job_queue_rank"] = 0
		
		if os.path.exists(results_file):
			qstat_data_response['job_state'] = "C"
			qstat_data_response['job_state_full'] = "Job is completed"
			qstat_data_response["job_queue_rank"] = 0

			if os.path.exists(results_file):
				with open(results_file) as json_file:
					data = json.load(json_file)
					qstat_data_response['data'] = data

			return qstat_data_response

		qstat_data = self.show_jobs()

		if job_identifier in qstat_data["jobs"]:
			qstat_data_response['job_identifier'] = job_identifier
			qstat_data_response['job_state_full'] = qstat_data["jobs"][job_identifier]['job_state_full']
			qstat_data_response['job_state'] = qstat_data["jobs"][job_identifier]['job_state']
			qstat_data_response["job_count"] = qstat_data["job_count"]
			qstat_data_response["job_queue_rank"] = qstat_data["status"][qstat_data_response['job_state_full']].index(qstat_data_response['job_identifier'])+1
		elif os.path.exists(options_file):
			qstat_data_response['job_identifier'] = job_identifier
			qstat_data_response["job_count"] = qstat_data["job_count"]
			qstat_data_response['job_state'] = "F"
			qstat_data_response['job_state_full'] = "Job started but failed to complete"
			qstat_data_response["job_queue_rank"] = 0
		else:
			qstat_data_response['job_identifier'] = job_identifier
			qstat_data_response["job_count"] = qstat_data["job_count"]
			qstat_data_response['job_state'] = "U"
			qstat_data_response['job_state_full'] = "Job Identifier cannot be found"
			qstat_data_response["job_queue_rank"] = 0
		
		del qstat_data
		return qstat_data_response

	def queue_director(self):
		try:
			status = "Success"
			if 'queue_tasks' in self.options:
				tasks = self.options['queue_tasks']
			else:
				return {
				"status":"Error",
				"data":"tasks option must be one of [delete_job|submit_jobs|show_nodes|show_jobs]"
				}

			qstat_data = {}

			if "clear_cache" in tasks:
				if self.options['logged_in']:
					if "dataset_type" in self.options:
						
						import glob
						delete_files = glob.glob(os.path.join(self.paths['data_path'],"webservices",self.options['dataset_type'], '*'))
						for delete_file in delete_files:
							os.remove(delete_file)
							
						qstat_data["delete_files"] = delete_files
					else:
						return {
						"status":"Error",
						"data":"delete_job tasks option must be accompanied by a dataset_type option"
						}
				else:
					return {
						"status":"Error",
						"data":"You don't have permission to do that"
						}
						
			elif "delete_job" in tasks:
				if "job_id" in self.options:
					qstat_data = {"jobs":{}}

					qstat_str = os.popen('qdel ' + self.options["job_id"]).read()
					qstat_data["jobs"] = 'qdel ' + self.options["job_id"]

					import glob
					delete_files = glob.glob(os.path.join(self.paths['data_path'],"webservices","*",self.options["job_id"] + '*'))
					for delete_file in delete_files:
						os.remove(delete_file)

					qstat_data["delete_files"] = delete_files
				else:
					return {
					"status":"Error",
					"data":"delete_job tasks option must be accompanied by a job_id option"
					}

			elif "clean_unfinished_jobs" in tasks:
			
				import glob
				file_counter = {}
				qstat_data["delete_files"] = []

				job_files = glob.glob(os.path.join(self.paths['data_path'],"webservices","*",'*'))
				for job_file in job_files:
					job_id = job_file.split(".")[0]

					if job_id not in file_counter :
						file_counter[job_id] = 0

					file_counter[job_id] += 1


				for job_id in file_counter:
					if file_counter[job_id] != 4:
						delete_files = glob.glob(os.path.join(self.paths['data_path'],"webservices","*",job_id + '*'))
						for delete_file in delete_files:
							os.remove(delete_file)
						
						qstat_data["delete_files"] += delete_files
		
			elif "submit_jobs" in tasks:

				job_options = {}
				job_options["tool"] = self.options["dataset_type"]
				job_options['cmd'] = "python3.7 /home/scripts/data_management/queryManager.py"
				job_options['cmd_options'] = {}
				job_options['verbose'] = False
				
				for cmd_option in self.options:
					if cmd_option not in self.removed_cmd_options:
						job_options['cmd_options'][cmd_option] = self.options[cmd_option]

				hash = self.params_to_hash(job_options)
				job_options['job_id'] = hash
				
				job_check = self.job_status(job_options['job_id'],self.options["dataset_type"])
				results_file = os.path.join(self.paths['data_path'],"webservices",job_options["tool"],job_options['job_id'] + ".results.json")
				
				#if not self.options["cached"]:
				#	import glob
				#	delete_files = glob.glob(os.path.join(self.paths['data_path'],"webservices",self.options['dataset_type'], job_options['job_id'] + '*'))
				#	for delete_file in delete_files:
				#		os.remove(delete_file)
				
				if os.path.exists(results_file):
					try:
						status = "Completed"
						with open(results_file) as json_file:
							data = json.load(json_file)
							qstat_data = data
					except:
						status = "Submitted"
						job_options['cmd_options']['outfile'] = results_file

						qstat_data['job_id'] = job_options['job_id']
						qstat_data["job_submission_info"] = self.create_jobfile(job_options)
						
						#qstat_data['job_queue'] = self.show_jobs()
						qstat_data = self.job_status(job_options["job_id"],self.options["dataset_type"])

				elif job_check['job_state'] != "U":
					qstat_data = job_check
					#qstat_data['job_queue'] = self.show_jobs()
					status = qstat_data['job_state_full']
				else:
					status = "Submitted"
					job_options['cmd_options']['outfile'] = results_file

					qstat_data['job_id'] = job_options['job_id']
					qstat_data["job_submission_info"] = self.create_jobfile(job_options)
					
					#qstat_data['job_queue'] = self.show_jobs()
					qstat_data = self.job_status(job_options["job_id"],self.options["dataset_type"])
					#qstat_data = job_check

					if self.options["remote_addr"] not in self.user_jobs:
						self.user_jobs[self.options["remote_addr"]] = {}

					self.user_jobs[self.options["remote_addr"]][job_options['job_id']] = self.options

			elif "show_nodes" in tasks:
				qstat_str = os.popen("pbsnodes -a").read()
				qstat_data = {"nodes":[]}
				for line in qstat_str.strip().split("\n"):
					for line_bits in line.split(","):
						qstat_data["nodes"].append(line_bits.strip())

			elif "job_id" in self.options and "dataset_type" in self.options and "job_status" in tasks:
				qstat_data = self.job_status(self.options["job_id"],self.options["dataset_type"])
			
			elif "show_jobs" in tasks:
				qstat_data = self.show_jobs()
		
			elif "show_users" in tasks:
				qstat_data = {"users":self.user_jobs}

			else:
				return {
				"status":"Error",
				"error_type":"It seems you didn't select a task"
				}

			return {
			"status":status,
			"options":self.options,
			"data":qstat_data
			}

		except Exception as e:
			exc_type, exc_value, exc_tb = sys.exc_info()
			error_traceback = traceback.format_exception(exc_type, exc_value, exc_tb)

			return {
			"status":"Error",
			"data":error_traceback,
			"options":self.options,
			"paths":self.paths
			}

	
