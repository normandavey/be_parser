import sys,os,inspect,pprint,json

file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../"))
sys.path.append(os.path.abspath(os.path.join(file_path,"../utilities")))
sys.path.append(os.path.abspath(os.path.join(file_path,"../data_management")))

import queryRunner

import utilities_downloader
import utilities_basic

import config_reader
import option_reader

#-----
import logging
logger = logging.getLogger(__name__)
#-----

class databaseQueryRestManager():
	
	#####------------------------------------------------------------------#####
	
	def __init__(self):
		self.options = {
			"username":"",
			"password":""
		}

		self.options.update(config_reader.load_configeration_options(sections=['general']))
		self.options.update(option_reader.load_commandline_options(self.options,{},purge_commandline=True))
		self.options.update({
			'params': {}
		})
				
	############################################################################
	############################################################################
	############################################################################
	
	def get_architecture(self):
		self.options['architecture_data_path'] = os.path.join(os.path.join(self.options['data_path'],"uniprot","architecture"))

		if not os.path.exists(self.options['architecture_data_path']):
			logging.error("Making: " + self.options['architecture_data_path'])
			os.mkdir(self.options['architecture_data_path'])

		out_path = os.path.join(self.options['architecture_data_path'], str(self.options['accession']) + ".architecture.jsonc")
		
		if not os.path.exists(out_path):
			url =  "https://slim.icr.ac.uk/restapi/database/protein/" + str(self.options['accession']) + "/architecture"
			
			sessionDownloaderObj = utilities_downloader.sessionDownloader()
			status = sessionDownloaderObj.download_file(url,out_path,method="GET",zipped_JSON=True,replace_empty=True,verify=False)
			logger.debug(status)
		
		if os.path.exists(out_path):
			architecture = utilities_basic.read_from_json(out_path,zipped=True)
			return architecture['result']
		else:
			return {"status":"error","error_type":out_path + " - file does not exist"}

	#####------------------------------------------------------------------#####
 
	def get_features(self):
		self.options['features_data_path'] = os.path.join(os.path.join(self.options['data_path'],"uniprot","features"))

		if not os.path.exists(self.options['features_data_path']):
			logging.error("Making: " + self.options['features_data_path'])
			os.mkdir(self.options['features_data_path'])

		out_path = os.path.join(self.options['features_data_path'], str(self.options['accession']) + ".features.jsonc")
		
		if not os.path.exists(out_path):
			url =  "https://slim.icr.ac.uk/restapi/database/protein/" + str(self.options['accession']) + "/features"
			
			sessionDownloaderObj = utilities_downloader.sessionDownloader()
			status = sessionDownloaderObj.download_file(url,out_path,method="GET",zipped_JSON=True,replace_empty=True,verify=False)
			if status['status'] == "error":
				logger.debug(status)
		
		if os.path.exists(out_path):
			features = utilities_basic.read_from_json(out_path,zipped=True)
			if 'residue_centric' in self.options:
				if self.options['residue_centric']:
					return self.collapse_features_residue_centric(features['result'])

			return features['result']
		else:
			return {"status":"error","error_type":out_path + " - file does not exist"}

	def collapse_features_residue_centric(self,features,skip_features=['chain']):
		features_residue_centric = {}
		for feature in features:
			for i in range(int(feature['start']),int(feature['stop'])+1):
				try:
					if feature['type'] == 'region': 
						if feature['description']['description'].split(';')[0] == "chain": 
							continue
				except:
					pass

				if i not in features_residue_centric: features_residue_centric[i] = {}
				if feature['type'] not in features_residue_centric: features_residue_centric[i][feature['type']] = []
				
				features_residue_centric[i][feature['type']].append(feature)
		
		return features_residue_centric
	
	#####------------------------------------------------------------------#####
 
	def get_attributes(self):
		self.options['attributes_data_path'] = os.path.join(os.path.join(self.options['data_path'],"uniprot","attributes"))

		if not os.path.exists(self.options['attributes_data_path']):
			logging.error("Making: " + self.options['attributes_data_path'])
			os.mkdir(self.options['attributes_data_path'])

		out_path = os.path.join(self.options['attributes_data_path'], str(self.options['accession']) + ".attributes.jsonc")
		
		if not os.path.exists(out_path):
			url =  "https://slim.icr.ac.uk/restapi/database/protein/" + str(self.options['accession']) + "/attributes"
		
			sessionDownloaderObj = utilities_downloader.sessionDownloader()
			status = sessionDownloaderObj.download_file(url,out_path,method="GET",zipped_JSON=True,replace_empty=True,verify=False)
			logger.debug(status)
		
		if os.path.exists(out_path):
			features = utilities_basic.read_from_json(out_path,zipped=True)
			return features
		else:
			return {"status":"error","error_type":out_path + " - file does not exist"}

	#####------------------------------------------------------------------#####


	