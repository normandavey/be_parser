import os,re,json,inspect,sys,pprint

from xml.dom import minidom
from xml.etree import cElementTree as elementtree

file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"../"))
sys.path.append(os.path.join(file_path,"../utilities/"))
import utilities_downloader
import utilities_basic
import config_reader

#-----
import logging
logger = logging.getLogger(__name__)
#-----

#
# https://github.com/ProteinsWebTeam/interpro7-api/blob/master/docs/modifiers.md
# https://www.ebi.ac.uk/interpro/api/entry/all/pfam/protein/reviewed/Q14289/
# https://www.ebi.ac.uk/interpro/api/set/pfam/Cl0632/
# https://www.ebi.ac.uk/interpro/api/entry/all/pfam/PF00373

class pfamDownloader():

	def __init__(self):
		self.options = {}
		self.options["wait_time"] = 0.01
		self.options["data_path"] = ""
		self.options["remake_age"] = 1800

		config_options = config_reader.load_configeration_options(sections=["general","structure_reader"])

		for config_option in config_options:
			self.options[config_option] = config_options[config_option]

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def check_directory(self):
		if not os.path.exists(os.path.join(self.options["data_path"], "pfam")):
			os.mkdir(os.path.join(self.options["data_path"], "pfam"))

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def grabPfamInterPro(self,pfam_id):
		self.check_directory()

		url = 'https://www.ebi.ac.uk/interpro/api/entry/all/pfam?search=' + pfam_id + '&extra_fields=short_name,description'
		out_path = os.path.join(self.options["data_path"], "pfam",pfam_id + ".pfam.json")
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=False)
		return status

	def grabPfamClanMappingInterPro(self,pfam_id):
		self.check_directory()
		url = 'https://www.ebi.ac.uk/interpro/api/entry/all/set/pfam/?search=' + pfam_id + '&extra_fields=short_name,description'
		out_path = os.path.join(self.options["data_path"], "pfam",pfam_id + ".clan.pfam.json")
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=False)
		return status

	def grabPfamClanInterPro(self,clan_id):
		self.check_directory()

		url = 'https://www.ebi.ac.uk/interpro/api/set/pfam/' + clan_id
		out_path = os.path.join(self.options["data_path"], "pfam",clan_id + ".clan.json")
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=False)
		return status

	def parsePfamClanName(self,clan_id):
		self.grabPfamClanInterPro(clan_id)

		json_path = os.path.join(self.options["data_path"], "pfam",clan_id + ".clan.json")
		if os.path.exists(json_path):
			with open(json_path) as outfile:
				pfam_json_content = json.load(outfile)

			return pfam_json_content["metadata"]["name"]["name"]
		else:
			return ""

	def parsePfam(self,pfam_id=None):
		if pfam_id == None:
			if 'pfam_id' in self.options:pfam_id = self.options['pfam_id']
			else: return {}

		status = self.grabPfamInterPro(pfam_id)

		if status['status'] == "Error":
			return status

		json_path  = os.path.join(self.options["data_path"], "pfam",pfam_id + ".pfam.json")

		try:
			if os.path.exists(json_path):
				try:
					try:
						with open(json_path) as outfile:
							pfam_json_content = json.load(outfile)
					except:
						logger.error(json_path + " not correctly formatted")
						pfam_json_content = {}

					try:
						description = utilities_basic.remove_html_tags(pfam_json_content['results'][0]['extra_fields']['description'][0])
					except:
						description = ""

					short_name = pfam_json_content['results'][0]['extra_fields']['short_name']
					domain_name = pfam_json_content['results'][0]['metadata']['name']

					try:
						status = self.grabPfamClanMappingInterPro(pfam_id)
						json_clan_mapping_path  = os.path.join(self.options["data_path"], "pfam",pfam_id + ".clan.pfam.json")
						with open(json_clan_mapping_path) as outfile:
							pfam_clan_mapping_json_content = json.load(outfile)

						clan_id = pfam_clan_mapping_json_content['results'][0]['set_subset'][0]['accession']
						clan_name = self.parsePfamClanName(clan_id)
					except:
						clan_id = ''
						clan_name = ''

					pfam_data = {"pfam_id":pfam_id,"pfam_acc":short_name,"domain_name":domain_name ,"domain_description":description,"clan_name":clan_name,"clan_id":clan_id}

					try:
						if 'integrated' in pfam_json_content['results'][0]['metadata']:
							pfam_data['interpro_id'] = pfam_json_content['results'][0]['metadata']['integrated']
					except:
						pass

					return {"status":"Added","data":pfam_data}
				except:
					return {"status":"Error","error_type":"File parsing failed"}
			else:
				return {"status":"Error","error_type":"File not found"}

		except Exception as e:
			raise
			return {"status":"Error","error_type":str(e)}
