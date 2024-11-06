import os,re,time,string,json,pprint, hashlib, sys, copy

from xml.dom import minidom
from xml.etree import cElementTree as elementtree

import inspect
file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"../"))
sys.path.append(os.path.join(file_path,"../utilities/"))
import config_reader

import pfamDownloader
import utilities_downloader
import utilities_error

#-----
import logging
logger = logging.getLogger(__name__)
#-----

class uniprotDownloader():

	##------------------------------------------------------------------##

	def __init__(self):
		self.options = {}
		self.options["wait_time"] = 0.01
		self.options["remake_age"] = 180
		self.options["use_features"] = []
		self.options["format"] = "list"

		config_options = config_reader.load_configeration_options(sections=["general"])
		self.options.update(config_options)

		self.settings = {}
		self.settings["disease_keywords"] = ["KW-0014","KW-0015","KW-0020","KW-0023","KW-0026","KW-0038","KW-0043","KW-0065","KW-0069","KW-0070","KW-0087","KW-0122","KW-0161","KW-0172","KW-0182","KW-0192","KW-0199","KW-0209","KW-0214","KW-0218","KW-0219","KW-0225","KW-0241","KW-0242","KW-0248","KW-0263","KW-0307","KW-0316","KW-0322","KW-0331","KW-0335","KW-0355","KW-0360","KW-0361","KW-0362","KW-0367","KW-0370","KW-0380","KW-0435","KW-0451","KW-0454","KW-0461","KW-0466","KW-0510","KW-0523","KW-0550","KW-0553","KW-0586","KW-0622","KW-0656","KW-0657","KW-0682","KW-0685","KW-0705","KW-0751","KW-0757","KW-0772","KW-0792","KW-0821","KW-0852","KW-0855","KW-0856","KW-0857","KW-0887","KW-0898","KW-0900","KW-0901","KW-0905","KW-0907","KW-0908","KW-0910","KW-0912","KW-0913","KW-0923","KW-0940","KW-0942","KW-0947","KW-0948","KW-0951","KW-0954","KW-0955","KW-0956","KW-0958","KW-0976","KW-0977","KW-0982","KW-0983","KW-0984","KW-0985","KW-0986","KW-0987","KW-0988","KW-0989","KW-0991","KW-0992","KW-0993","KW-1004","KW-1007","KW-1008","KW-1011","KW-1013","KW-1014","KW-1016","KW-1020","KW-1021","KW-1022","KW-1023","KW-1024","KW-1026","KW-1054","KW-1056","KW-1057","KW-1058","KW-1059","KW-1060","KW-1062","KW-1063","KW-1065","KW-1066","KW-1067","KW-1068","KW-1186","KW-1211","KW-1212","KW-1215","KW-1268","KW-1274"]

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def check_directory(self):
		if not os.path.exists(os.path.join(self.options["data_path"], "uniprot")):
			os.mkdir(os.path.join(self.options["data_path"], "uniprot"))

		for directory in ["pdb","genetree","uniref","mobidb",'taxonomy','xml','bulk','fasta','pfam','attributes',"obsolete",'keywords','history']: #"features",
			if not os.path.exists(os.path.join(self.options["data_path"], "uniprot",directory)):
				print(("Making",os.path.join(self.options["data_path"], "uniprot",directory)))
				os.mkdir(os.path.join(self.options["data_path"], "uniprot",directory))

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def check_protein_accession_files(self):
		from pathlib import Path

		for path in Path(os.path.join(self.options["data_path"], "uniprot")).rglob(self.options['accession'][0] + '*'):
			print(path)

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def check_accession(self,accession):
		accession = accession.split("-")[0]
		if re.match("\A([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\Z", accession):
			return True
		elif re.match("\A([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\-{0,1}[0-9]*\Z", accession):
			logger.error(str(accession) + " - isoform")
			return True
		else:
			logger.error(str(accession) + " - not a valid accession")
			return False

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def check_isoform(self,accession):
		return accession.split("-")[0]

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def grab_taxon_identifier_information_from_uniprot(self):
		self.check_directory()

		rerun = True

		while rerun:
			#url = "https://www.uniprot.org/taxonomy/?sort=score&desc=&query=id:" + str(self.options['taxon_id'])  + '&force=no&format=tab&columns=id'
			url = "https://rest.uniprot.org/taxonomy/search?query=id:" + str(self.options['taxon_id'])  + "&force=no&format=tsv&columns=id"
			out_path = os.path.join(self.options["data_path"], "uniprot","taxonomy", str(self.options['taxon_id']) + ".uniprot_info.tdt")

			sessionDownloaderObj = utilities_downloader.sessionDownloader()
			status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True)

			f  = open(out_path,'r')
			content = f.read()
			f.close()
			if 'Taxon' in content:
				rerun = False
			else:
				os.remove(out_path)


	def parse_taxon_identifier_information_from_uniprot(self):
		self.grab_taxon_identifier_information_from_uniprot()
		out_path = os.path.join(self.options["data_path"], "uniprot","taxonomy", str(self.options['taxon_id']) + ".uniprot_info.tdt")

		f = open(out_path,'r')
		content = f.read().split("\n")
		f.close()

		header = content[0].split("\t")

		line = content[1].split("\t")
		data = {}
		for idx, col in enumerate(header):
			data[col] = line[idx]
		data.update(self.parse_taxonomy_lineage())
		return data


	def grab_taxonomy_lineage_by_string(self):
		self.check_directory()

		url = "https://www.ebi.ac.uk/proteins/api/taxonomy/name/" + str(self.options['taxonomy']) + "?pageNumber=1&pageSize=1&searchType=EQUALSTO&fieldName=SCIENTIFICNAME"

		out_path = os.path.join(self.options["data_path"], "uniprot","taxonomy", str(self.options['taxonomy']) + ".taxonomy.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True)

	def parse_taxonomy_lineage_by_string(self):
		self.grab_taxonomy_lineage_by_string()
		json_path = os.path.join(self.options["data_path"], "uniprot","taxonomy", str(self.options['taxonomy']) + ".taxonomy.json")

		taxonomyId = ""
		if os.path.exists(json_path):
			with open(json_path) as outfile:
				json_data = json.load(outfile)

			if 'taxonomies' in json_data:
				taxonomyId = json_data['taxonomies'][0]['taxonomyId']

		return  taxonomyId

	##------------------------------------------------------------------##
	##------------------------------------------------------------------##

	def grab_taxonomy_lineage(self):
		self.check_directory()

		url = "https://www.ebi.ac.uk/proteins/api/taxonomy/lineage/" + str(self.options['taxon_id'])

		out_path = os.path.join(self.options["data_path"], "uniprot","taxonomy", str(self.options['taxon_id']) + ".lineage.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True)

	##------------------------------------------------------------------##

	def grab_taxonomy_information(self):
		self.check_directory()

		url = "https://www.ebi.ac.uk/proteins/api/taxonomy/id/" + str(self.options['taxon_id']) + ""

		out_path = os.path.join(self.options["data_path"], "uniprot","taxonomy", str(self.options['taxon_id']) + ".info.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True)


	##------------------------------------------------------------------##
	def parse_taxonomy_lineage(self,detailed=False):

		self.grab_taxonomy_lineage()
		json_path = os.path.join(self.options["data_path"], "uniprot","taxonomy", str(self.options['taxon_id']) + ".lineage.json")

		data = {
			'lineage':[],
			'taxonomies':[],
			'ranks':{}
			}

		if os.path.exists(json_path):
			with open(json_path) as outfile:
				json_data = json.load(outfile)

				if 'taxonomies' in json_data:
					data['taxonomies'] = json_data['taxonomies']

			for rank in data['taxonomies']:
				if rank['rank'] != 'no rank':
					data['ranks'][rank['rank']] = rank['scientificName']
				if rank['rank'] in ['species',
										'genus',
										'family',
										'order',
										'class',
										'subphylum',
										'phylum',
										'kingdom',
										'superkingdom']:
					data['lineage'].append(rank['scientificName'])

		if detailed == False:
			del data['taxonomies']

		data['lineage'].reverse()

		return data

	def parse_taxonomy_details(self):

		self.grab_taxonomy_information()

		json_path = os.path.join(self.options["data_path"], "uniprot","taxonomy", str(self.options['taxon_id']) + ".info.json")

		data = {
			'mnemonic':"",
			'scientificName':"",
			'commonName':"",
			'synonym':"",
			'rank':""
			}

		if os.path.exists(json_path):

			with open(json_path) as outfile:
				json_data = json.load(outfile)

			if 'mnemonic' in json_data: data['mnemonic'] = json_data['mnemonic']
			if 'scientificName' in json_data: data['scientificName'] = json_data['scientificName']
			if 'commonName' in json_data: data['commonName'] = json_data['commonName']
			if 'synonym' in json_data: data['synonym'] = json_data['synonym']
			if 'rank' in json_data: data['rank'] = json_data['rank']

		return data

	def parse_taxon_identifier_information(self):
		data = self.parse_taxonomy_details()
		data.update(self.parse_taxonomy_lineage())

		return data

	def parse_taxonomy_information(self):
		taxonomyId = self.parse_taxonomy_lineage_by_string()

		if taxonomyId != "":
			self.options['taxon_id'] = taxonomyId
			data = self.parse_taxon_identifier_information()
		else:
			return {"status":"Error","error_type":"Taxonomy not found"}
		return data

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def grabPTMvocabulary(self):
		self.check_directory()

		# url = "https://www.uniprot.org/docs/ptmlist.txt"
		url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist"
		out_path = os.path.join(self.options['data_path'], "uniprot", "ptmlist.txt")


		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,method="GET",remake_age=self.options['remake_age'])


	def parsePTMvocabulary(self):
		self.grabPTMvocabulary()
		out_path = os.path.join(self.options['data_path'], "uniprot", "ptmlist.txt")

		ptms = []
		f = open(out_path,'r')
		content = f.read().split("//")
		f.close()

		key_mapped = {
			"ID": "id",
			"AC": "accession",
			"FT": "feature",
			"TG": "target",
			"PA": "mod_position_aa",
			"PP": "mod_position_polypeptide",
			"CF": "correction_formula",
			"MM": "monoisotopic_mass_difference",
			"MA": "avg_mass_difference",
			"LC": "localisation",
			"TR": "taxonomic_range",
			"KW": "keywords",
			"DR": "xrefs"
		}

		for elem in content:
			items = elem.split("\n")
			item_start = [x for x in items if x.startswith("ID")]
			if not item_start:continue
			item_start_idx = items.index(item_start[0])
			items = items[item_start_idx:]
			# print("items", items)
			tmp = {}
			for item in items:
				item = item.strip()
				item_start = item[:2]
				if item_start not in key_mapped: continue
				key_name = key_mapped[item_start]
				if item_start in ['DR','KW','TR']:
					if key_name not in tmp: tmp[key_name] = []

					val = item[2:].strip()
					if val.endswith("."): val = val[:-1]

					if item_start == 'DR':
						source = val.split(";")[0].strip()
						id = val.split(";")[-1].strip()

						val = {"source": source, "id": id}
						tmp[key_name].append(val)
					elif item_start == 'KW':
						val = val.split(";")
						val = [x.replace(".","").strip() for x in val]
						tmp[key_name].extend(val)
						continue
					elif item_start == 'TR':
						if 'taxons' not in tmp: tmp['taxons'] = []
						items = val.split(";")
						clade = items[0].strip()
						tmp[key_name].append(clade)
						taxons = items[1].split(",")
						for taxon in taxons:
							taxon_id =re.findall("taxId:\d+",taxon)
							if taxon_id:
								taxon_id = taxon_id[0].split("taxId:")[-1].strip()
								name = re.findall("\((.*?)\)",taxon)[0]
								tmp['taxons'].append({"taxon_id": taxon_id, "name": name, "top": clade})


				else:
					tmp[key_name] = item[2:].strip()
					if tmp[key_name].endswith("."): tmp[key_name] = tmp[key_name][:-1]
			# if tmp['feature'] == 'CROSSLNK': continue
			ptms.append(tmp)


		# print("ptms",len(ptms), set([x['feature'] for x in ptms]))
		# print("\ncontent", content[:5])
		return ptms

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def get_keyword_proteins(self):
		self.check_directory()

		url = 'https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28keyword%3A' + self.options['keyword_id'] +'%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29'
		out_path = os.path.join(self.options['data_path'], "uniprot", "keywords",self.options['keyword_id'] + '.tdt')
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,method="GET",remake_age=self.options['remake_age'])
		return open(out_path).read().strip().split("\n")

	def parse_keywords_proteins(self):
		self.grab_keywords_proteins()
		out_path = os.path.join(self.options['data_path'], "uniprot", "kwlist.txt")

		kw = {}
		f = open(out_path,'r')
		content = f.read().split("\n")
		f.close()

		columns = content[0].split("\t")

		for elem in content[1:]:
			line = elem.split("\t")
			if len(line) < 2: continue
			id =line[columns.index('Keyword ID')]
			kw[id] = {}
			for idx, col in enumerate(columns):
				kw[id][col] = line[idx]

		return kw

	def grabKeywords(self):

		self.check_directory()

		url = "https://rest.uniprot.org/keywords/stream?fields=id%2Cname%2Ccategory%2Cgene_ontologies&format=tsv&query=%2A"
		out_path = os.path.join(self.options['data_path'], "uniprot", "kwlist.txt")


		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,method="GET",remake_age=self.options['remake_age'])

	def parseKeywordsDetails(self):
		self.grabKeywords()
		out_path = os.path.join(self.options['data_path'], "uniprot", "kwlist.txt")

		kw = {}
		f = open(out_path,'r')
		content = f.read().split("\n")
		f.close()

		columns = content[0].split("\t")

		for elem in content[1:]:
			line = elem.split("\t")
			if len(line) < 2: continue
			id =line[columns.index('Keyword ID')]
			kw[id] = {}
			for idx, col in enumerate(columns):
				kw[id][col] = line[idx]

		return kw

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def grabHumsavar(self):
		self.check_directory()

		# url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
		url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
		out_path = os.path.join(self.options["data_path"], "uniprot","humsavar.txt")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		# status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,method="FTP",remake_age=self.options['remake_age'])
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,method="GET",remake_age=self.options['remake_age'])

	def parseHumsavar(self,accession="all"):
		self.check_directory()

		self.grabHumsavar()
		out_path = os.path.join(self.options["data_path"], "uniprot","humsavar.txt")

		mutations = {}

		with open(out_path) as fileobject:
			for line in fileobject:
				line_bits = line.split()
				if len(line_bits) > 5:
					if line_bits[1] == accession or accession == "all":

						if line_bits[1] not in mutations:
							mutations[line_bits[1]] = {"mutations":{},"positions":{},"phenotype":{}}

						wildtype = line_bits[3][2:5]
						mutant = line_bits[3][-3:]
						position = line_bits[3][5:-3]

						mutations[line_bits[1]]["mutations"][line_bits[2]] = {
						"position":position,
						"wildtype":wildtype,
						"mutant":mutant,
						"FTId": line_bits[2],
						"dbSNP": line_bits[5],
						"variant_type": line_bits[4],
						"variant": line_bits[3],
						"description":" ".join(line_bits[6:])
						}

						if line_bits[4] not in mutations[line_bits[1]]["phenotype"]:
							mutations[line_bits[1]]["phenotype"][line_bits[4]] = {}

						if position not in mutations[line_bits[1]]["phenotype"][line_bits[4]]:
							mutations[line_bits[1]]["phenotype"][line_bits[4]][position] = []

						mutations[line_bits[1]]["phenotype"][line_bits[4]][position].append(line_bits[2])

						if position not in mutations[line_bits[1]]["positions"]:
							mutations[line_bits[1]]["positions"][position] = []

						mutations[line_bits[1]]["positions"][position].append(line_bits[2])

		if accession == "all":
			return mutations
		else:
			if accession in mutations:
				return mutations[accession]
			else:
				return {}

	"""
	def grabFeatures(self,accession):
		url = "http://slim.ucd.ie/rest/features/?uniprot_acc=" + accession + "&type=motif&type=modification"

		out_path = os.path.join(self.options["data_path"], "uniprot","features",accession + ".features.json")
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True)


	def parseFeatures(self,accession):
		try:
			self.grabFeatures(accession)

			out_path = os.path.join(self.options["data_path"], "uniprot","features",accession + ".features.json")
			with open(out_path) as outfile:
				return {"status":"Added","data":json.load(outfile)}
		except Exception,e:
			return {"status":"Error","error_type":str(e)}
	"""

	def grabAttributes(self,accession):
		self.check_directory()

		# url = "http://slim.icr.ac.uk/peptools_webservices/attributes/?uniprot_acc=" + accession + "&type=Conservation&type=IUPred&type=Anchor&type=PsiPred"
		url = "http://slim.icr.ac.uk/restapi/database/protein/" + accession + "/attributes"
		out_path = os.path.join(self.options["data_path"], "uniprot",'attributes',accession + ".attributes.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True)

	def parseAttributes(self,accession):
		try:
			self.grabAttributes(accession)

			out_path = os.path.join(self.options["data_path"], "uniprot",'attributes',accession + ".attributes.json")
			with open(out_path) as outfile:
				data = json.load(outfile)
				return {"status":"Added","data":data}
		except Exception as e:
			return {"status":"Error","error_type":str(e)}

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##
	def grabUniProtGenetree(self,accession):
		self.check_directory()

		# url = "https://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=" + accession + "&fil=&limit=1&force=no&preview=true&format=tab&columns=id,database(GeneTree)"
		url = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cxref_genetree&format=tsv&query=%28accession%3A" + accession + "%29"
		out_path = os.path.join(self.options["data_path"] , "uniprot","genetree",accession + ".genetree")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method='GET',replace_empty=True)

	def getUniProtGenetree(self,accession):
		self.grabUniProtGenetree(accession)

		out_path = os.path.join(self.options["data_path"] , "uniprot","genetree",accession + ".genetree")

		if os.path.exists(out_path):
			line_bits = open(out_path).read().split("\n")
			if len(line_bits) > 1:
				return line_bits[1].split("\t")[1].strip(";")
			else:
				return ""

	def getUniProtGenetreeRecursive(self,accession):

		uniProtUnirefMembers = self.getUniProtUnirefMembers(accession.split("-")[0])
		homologue_accessions = uniProtUnirefMembers['data'][0]['Cluster members'].split(";")
		genetree_accession = ""
		for homologue_accession in homologue_accessions:
			if homologue_accession[0:3] != "UPI":
				genetree_accession = self.getUniProtGenetree(homologue_accession)

			if genetree_accession != "":
					break

		return genetree_accession
	
	def getClosestReviewedProteinUniRef(self,accession):
		for identity in [1.0,0.9,0.5]:
			protein_cluster_data = self.getUniProtUnirefMembers(accession,identity)
			protein_data = self.parseBasic(protein_cluster_data['data'][0]['Cluster ID'].split("_")[1] )
			
			if protein_data['data']['reviewed']:
				return {"accession":protein_data['data']['accession'],"cluster_identity":identity}
			
		return {}
	
	##------------------------------------------------------------------##
	#"https://rest.uniprot.org/uniref/search?compressed=true&fields=id%2Cname%2Ctypes%2Ccount%2Corganism%2Clength%2Cidentity&format=tsv&query=%28uniprot_id%3A" + + "%29&size=500"

	def grabUniProtUnirefMembers(self,accession,identity=0.5):
		url = "https://rest.uniprot.org/uniref/search?fields=id%2Cname%2Ctypes%2Ccount%2Cmembers%2Corganism%2Clength%2Cidentity&format=tsv&query=%28uniprot_id%3A" + accession + "%29%20AND%20%28identity%3A" + str(identity) + "%29&size=500"

		if identity == 'all':
			url = "https://rest.uniprot.org/uniref/search?fields=id%2Cname%2Ctypes%2Ccount%2Cmembers%2Corganism%2Clength%2Cidentity&format=tsv&query=%28uniprot_id%3A" + accession + "%29&size=500"


		out_path = os.path.join(self.options["data_path"] , "uniprot","uniref", accession + "." + str(identity) + ".uniref")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,remake_age=self.options['remake_age'])

	def getUniProtUnirefMembers(self,accession,identity=0.5):
		self.grabUniProtUnirefMembers(accession,identity)

		out_path = os.path.join(self.options["data_path"] , "uniprot","uniref", accession + "." + str(identity) + ".uniref")

		if os.path.exists(out_path):
			f = open(out_path,'r')
			content = f.read().split("\n")
			f.close()

			if len(content) > 1:
				data = []
				header = content[0].split("\t")
				try:
					idx= header.index("")
					header[idx] ="Status"
				except:
					pass

				for line in content[1:]:
					data_tmp = {}
					line_splitted = line.split("\t")
					if len(line_splitted) > 1:
						for i in range(0, len(header)):
							data_tmp[header[i]] = line_splitted[i]

						data.append(data_tmp)

				return {
				"status": "Success",
				"data":data
				}
			else:
				return {"status": "Error","error_type":"No data returned"}
		else:
			return {"status": "Error", "error_type": "Couldn't download a file."}

	def getUniProtUnirefClusterName(self,accession,identity=0.5):
		self.grabUniProtUnirefMembers(accession,identity)

		out_path = os.path.join(self.options["data_path"] , "uniprot","uniref", accession + "." + str(identity) + ".uniref")

		if os.path.exists(out_path):
			line_bits = open(out_path).read().split("\n")

			if len(line_bits) > 1:
				return {
				"status": "Success",
				"data":line_bits[1].split("\t")[0]
				}
			else:
				return {"status": "Error","error_type":"No data returned"}

	#------------------------------------
	#------------------------------------
	"""

	CURRENTLY NOT WORKING
	def grabUniProtUnirefByTaxon(self, range="0.99",reviewed=False):
		url = "http://www.uniprot.org/uniref/?sort=score&desc=&compress=no"
		url += "&query=taxonomy:" + self.options["taxon_id"]
		url += "%20AND%20reviewed:yes"
		url += "&fil=identity:" + range
		url += "&force=no&preview=true&format=tab"
		url += "&columns=id,reviewed,name,count,members,organisms,length,identity"

		if reviewed:
			url += "&fil=reviewed%3Ayes"

		out_path = os.path.join(self.paths["uniprot_dir"] , self.options["taxon_id"]  + ".uniref")
		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True)

	#-----------------------------------------------
	#-----------------------------------------------

	def getUniProtUnirefByTaxon(self, range="0.99",reviewed=False):
		self.grabUniProtUniref(range,reviewed)

		out_path = os.path.join(self.paths["uniprot_dir"] , self.options["taxon_id"]  + ".uniref")

		self.mapping["uniprot_uniref"] = {}
		self.mapping["uniref_uniprot"] = {}

		if os.path.exists(out_path):
			for line in open(out_path).read().strip().split("\n"):
				line_bits = line.split("\t")

				print len(self.mapping["uniprot_uniref"]),len(self.mapping["uniref_uniprot"]),len(line_bits[4].split(";"))

				for accession in line_bits[4].split(";"):
					if accession.strip() in self.mapping["uniprot"]:
						self.mapping["uniprot_uniref"][accession.strip()] = {"cluster":line_bits[0],"cluster_name":line_bits[2],}
						self.mapping["uniref_uniprot"][line_bits[0]] = [x.strip() for x in line_bits[4].split(";")]
	"""
	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def grabMobiDB(self,accession,force=False):
		accession = self.check_isoform(accession)

		# self.check_directory()

		self.check_directory()
		url = "http://mobidb.bio.unipd.it/ws/entries/" + accession + "/disorder"
		out_path = os.path.join(self.options["data_path"], "uniprot","mobidb",accession + ".mobidb.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,remake_age=self.options['remake_age'])

	def parseMobiDB(self,accession):
		try:
			self.grabMobiDB(accession)

			out_path = os.path.join(self.options["data_path"], "uniprot","mobidb",accession + ".mobidb.json")
			with open(out_path) as outfile:
				return {"status":"Added","data":json.load(outfile)}

		except Exception as e:
			return {"status":"Error","error_type":str(e)}

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def get_uniprot_accession_taxa(self,taxon_id,format="list",reviewed=False,structures=False,nofragment=True, complete_proteome=False, include_version = False, protein_existence = False):

		self.check_directory()

		if structures == False:
			url = "https://rest.uniprot.org/uniprotkb/stream?query=(taxonomy_id=" + str(taxon_id) + ")"
		else:
			url = "https://rest.uniprot.org/uniprotkb/stream?query=(taxonomy_id=" + str(taxon_id) + ")" + " AND (keyword=KW-0002)"

		file_tag = ""

		if complete_proteome:
			url += " AND (keyword=KW-1185)"
			file_tag += ".complete_proteome"

		if nofragment:
			url += " AND (fragment=false)"
			file_tag += ".nofragment"
		if reviewed:
			url += " AND (reviewed=true)"
			file_tag += ".reviewed"

		if structures:
			file_tag += ".structures"

		if protein_existence:
			file_tag += ".protein_existence"
			url += " NOT(existence:4) NOT(existence:5)"


		if include_version:
			url += "&fields=accession%2Cid%2Cgene_primary%2Cversion%2Cdate_modified%2Csequence_version"
			format = "tsv"
			file_tag += '.versions'


		url += "&format=" + format

		out_path = os.path.join(self.options["data_path"], "uniprot",'taxonomy',str(taxon_id) + file_tag + "." +  format)

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method='GET',replace_empty=True,remake_age=30)

	##------------------------------------------------------------------##

	def parse_uniprot_accession_taxa(self,taxon_id,format="list",reviewed=False,structures=False,nofragment=True, include_version = False, protein_existence = False):
		file_tag = ""

		if nofragment:
			file_tag += ".nofragment"
		if reviewed:
			file_tag += ".reviewed"

		if structures:
			file_tag += ".structures"

		if protein_existence:
			file_tag += ".protein_existence"

		if include_version:
			format = "tsv"
			file_tag += '.versions'

		if format=="path":
			return {"path":os.path.join(self.options["data_path"], "uniprot",'taxonomy',str(taxon_id) + file_tag + "." + format )}
		else:
			self.options['taxonomy_path'] = os.path.join(self.options["data_path"], "uniprot",'taxonomy',str(taxon_id) + file_tag + "." + format )

		self.get_uniprot_accession_taxa(taxon_id,format,reviewed,structures,nofragment, include_version = include_version, protein_existence = protein_existence)

		if format=="list":
			print(self.options['taxonomy_path'])
			if os.path.exists(self.options['taxonomy_path'] ):
				return open(self.options['taxonomy_path'] ).read().split()
			else:
				return []

		elif format=="fasta":
			if os.path.exists(self.options['taxonomy_path'] ):
				proteins = {}
				for protein in open(self.options['taxonomy_path'] ).read().split(">")[1:-1]:
					try:
						proteins[protein.split("\n")[0].split("|")[1]] = "".join(protein.split("\n")[1:])
					except:
						pass

				return proteins
			else:
				return []

		elif format == "tsv":
			if os.path.exists(self.options['taxonomy_path'] ):
				proteins = {}
				f = open(self.options['taxonomy_path'] ,'r')
				content = f.read().split("\n")
				f.close()
				header = ["accession", "id", "gene_name", "version","date_mod","sequence_version"]
				# header = content[0].split("\t")
				for line in content[1:]:
					tmp = {}
					item = line.split("\t")
					if len(item) < 2: continue
					for idx, col in enumerate(header):
						tmp[col] = item[idx]
					proteins[tmp['accession']] = tmp
				return proteins
			else:
				return {}
		else:
			return self.options['taxonomy_path'] 


		##------------------------------------------------------------------##

	def parse_region_annotation(self):
		protein = self.options["accession"][0]
		start = self.options["region_start"]
		end = self.options["region_end"]
		return self.parseUniProtRegion(protein,start,end)

	##------------------------------------------------------------------##

	def parse_region_features(self,regions,flank_length=0,pad_peptide=True,format="JSON"):

		regions_features = {}
		for region in regions:
			region_bits = region.split(";")

			if region_bits[0] not in regions_features:
				regions_features[region_bits[0]] = {'regions':{}}

			regions_features[region_bits[0]]['regions'][region_bits[1]] = {
				"Start":int(region_bits[1]),
				"End":int(region_bits[2])
			}

		rows = []
		for protein in regions_features:
			for region in regions_features[protein]['regions']:
				start = regions_features[protein]['regions'][region]['Start']
				end = regions_features[protein]['regions'][region]['End']

				regions_features[protein] = self.parseUniProtRegion(protein,start,end)


		return regions_features
	##------------------------------------------------------------------##

	def parse_region_sequence(self,regions,flank_length=0,pad_peptide=True,format="JSON"):

		regions_sequences = {}
		for region in regions:
			region_bits = region.split(";")

			if region_bits[0] not in regions_sequences:
				regions_sequences[region_bits[0]] = {'regions':{}}

			regions_sequences[region_bits[0]]['regions'][region_bits[1]] = {
				"Start":int(region_bits[1]),
				"End":int(region_bits[2])
			}

		proteins = self.get_uniprot_mapping_details(list(regions_sequences.keys()), mapping_from='UniProtKB_AC-ID', mapping_to="UniProtKB", columns='accession,id,reviewed,protein_name,gene_names,organism_name,length,sequence')

		rows = []
		for protein in proteins:
			for region in regions_sequences[protein]['regions']:

				sequence = proteins[protein]['Sequence']
				start = regions_sequences[protein]['regions'][region]['Start']
				end = regions_sequences[protein]['regions'][region]['End']

				region_sequence = sequence[max(0,start-1-flank_length):min(end+flank_length,int(proteins[protein]['Length']))]

				if pad_peptide and start-1-flank_length < 0:
					pad_length = -(start-1-flank_length)
					region_sequence = "-"*pad_length + region_sequence

				if pad_peptide and end+flank_length > int(proteins[protein]['Length']):
					pad_length = (end+flank_length) - int(proteins[protein]['Length'])
					region_sequence = region_sequence + "-"*pad_length
					region_sequence = region_sequence[:flank_length*2 + 1]

				regions_sequences[protein]['Entry'] = proteins[protein]['Entry']
				regions_sequences[protein]['Entry_Name'] = proteins[protein]['Entry Name']
				regions_sequences[protein]['Gene_Name'] = ""
				if  proteins[protein]['Gene Names'] !="":
					regions_sequences[protein]['Gene_Name'] = proteins[protein]['Gene Names'].split()[0]

				regions_sequences[protein]['Protein_Name'] = proteins[protein]['Protein names']
				regions_sequences[protein]['regions'][region]['Region_Sequence'] = region_sequence

				row = [region_sequence,str(start),str(end)]
				for header in ['Entry','Entry_Name','Gene_Name','Protein_Name']:
					row.append(str(regions_sequences[protein][header]))

				rows.append(row)

		del proteins

		if format == "list":
			return rows
		elif format == "tdt":
			return "\n".join(["\t".join(row) for row in rows])
		else:
			return regions_sequences

	##------------------------------------------------------------------##

	def get_uniprot_mapping(self,accessions,mapping_from='UniProtKB_AC-ID',mapping_to="UniProtKB"):
		return self.get_uniprot_mapping_details(accessions,mapping_from=mapping_from,mapping_to=mapping_to,columns='accession')

	##------------------------------------------------------------------##

	def get_genename_to_uniprot_mapping_details(self,accessions,out_path, format="tsv"):

		##
		## FILTER WEIRD GENE NAME - They raise an error at uniprot
		##

		while("" in accessions) :
			accessions.remove("")

		#add_accessions = []
		remove_accessions = []
		required_format = "^[A-Za-z0-9-]{2,10}$"
#		
		"""
		for accession in accessions:
			if accession.count("-") > 0:
				remove_accessions.append(accession)
				add_accessions += accession.split("-")

		accessions = accessions + add_accessions
		"""

		for accession in accessions:
			if not re.match(required_format, accession):
				remove_accessions.append(accession)
				logger.error(accession + " does not match gene name format")

		for remove_accession in remove_accessions:
			if remove_accession in accessions:
				accessions.remove(remove_accession)

		prefix = ".".join(out_path.split(".")[:-1])
		suffix = out_path.split(".")[-1]

		mappings = {}
		chunk_size = 50
		for i in range(0,len(accessions),chunk_size):
			out_path = prefix +  "." + str(i) + "-" + str(i+chunk_size) + "." + suffix + "." + format
			if not os.path.exists(out_path):
				logging.debug("Making: " + accessions[i:i+chunk_size][0] + " -> " + accessions[i:i+chunk_size][-1] + " " + str(len(accessions[i:i+chunk_size])) + " of " + str(len(accessions)))
				logging.debug(out_path)
				# gene_list = 'gene_exact:' +  " OR gene_exact:".join(accessions[i:i+100])

				gene_list = "((gene:" + ") OR (gene:".join(accessions[i:i+chunk_size]) + "))" + " AND organism_id=9606"
				url = "https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession,id,reviewed,protein_name,gene_names,organism_name,organism_id&size=500&format=" + format + "&query=" + gene_list

				sessionDownloaderObj = utilities_downloader.sessionDownloader()
				status = sessionDownloaderObj.download_file(url,out_path,method='GET',replace_empty=True, remake_age=30)

			if os.path.exists(out_path):
				logging.debug("Reading: " + accessions[i:i+chunk_size][0] + " -> " + accessions[i:i+chunk_size][-1] + " " + str(len(accessions[i:i+chunk_size])) + " of " + str(len(accessions)))
				logging.debug(out_path)
		
				if format == "json":
					f = open(out_path,'r')
					content = json.loads(f.read())
					f.close()

					column_mapper =  {
						"primaryAccession": "Entry",
						"uniProtkbId": "Entry Name"
					}

					for item in content['results']:
						gene_names = []

						for gene in item['genes']:
							if 'geneName' in gene:
								gene_names.append(gene['geneName']['value'])
							if 'synonyms' in gene:
								for other_name in gene['synonyms']:
									gene_names.append(other_name['value'])

						for gene_name in gene_names:
							print(gene_name)
							if gene_name in accessions:
								if gene_name not in mappings:
									mappings[gene_name] = {"Gene Names": " ".join(gene_names)}
									for k, v in column_mapper.items():
										mappings[gene_name][v] = item[k]
									mappings[gene_name]['Reviewed'] = "reviewed" if 'Swiss-Prot' in item['entryType'] else "unreviewed"
									mappings[gene_name]['Organism'] = item['organism']['scientificName']
									if 'commonName' in item['organism']:
										if item['organism']['commonName']:
											mappings[gene_name]['Organism'] += " (" + item['organism']['commonName'] + ")"
									mappings[gene_name]['Organism (ID)'] = item['organism']['taxonId']
									if 'recommendedName' not in item['proteinDescription']:
										other_names = []
										for submitted_name in item['proteinDescription']['submissionNames']:
											other_names.append(submitted_name['fullName']['value'])
										mappings[gene_name]['Protein names'] = " ".join(other_names)
									else:
										mappings[gene_name]['Protein names'] = item['proteinDescription']['recommendedName']['fullName']['value']

									if 'alternativeNames' in item['proteinDescription']:
										for other_name in item['proteinDescription']['alternativeNames']:
											mappings[gene_name]['Protein names'] += " (" + other_name['fullName']['value'] + ")"
				else:
					mappings_data = open(out_path).read()

					lines = mappings_data.strip().split("\n")
					headers = lines[0].split("\t")
					for line in lines[1:]:

						mappings_tmp = {}
						for i in range(0,len(headers)):
							line_bits = line.split("\t")

							mappings_tmp[headers[i]] = line_bits[i]

						gene_names = mappings_tmp['Gene Names']

						for input_identifier in gene_names.split():
							if input_identifier in accessions:
								if input_identifier not in mappings:
									mappings[input_identifier] = mappings_tmp
			else:
				print(accessions)
				return {"status":"Error"}

		missing = set(accessions).difference(list(mappings.keys()))

		if len(missing) > 0:
			logger.error("Missing gene names")
			logger.error(missing)

		#sys.exit()

		return mappings


	def get_uniprot_mapping_details(self,accessions,mapping_from='UniProtKB_AC-ID',mapping_to="UniProtKB",format='tsv',columns='accession,id,reviewed,protein_name,gene_names,organism_name,length', best_hit = True, remake_age = 60):

		accessions = list(set(accessions))

		self.check_directory()

		url = 'https://www.uniprot.org/uploadlists/'

		accessions.sort()

		hash_accession = ";".join(accessions) + mapping_to + mapping_from + format + columns

		hash = str(hashlib.md5(hash_accession.encode()).hexdigest())


		out_path = os.path.join(self.options["data_path"], "uniprot",'bulk',hash + ".list." + format)

		if mapping_from == "GENENAME":
			return self.get_genename_to_uniprot_mapping_details(accessions,out_path)

		if os.path.exists(out_path):
			if (time.time() - os.path.getmtime(out_path))/60/60/24 > remake_age:
				os.system("rm " + out_path)

		if not os.path.exists(out_path) or self.options['remake']:
			# url = 'https://www.uniprot.org/uploadlists/'
			url = "https://rest.uniprot.org/idmapping/run"


			old_remap = {
				"ACC+ID": "UniProtKB_AC-ID",
				"ACC": "UniProtKB"
			}

			if mapping_from in old_remap:
				mapping_from = old_remap[mapping_from]
			if mapping_to in old_remap:
				mapping_to = old_remap[mapping_to]

			if format == "tab": format="tsv"

			### SUBMIT
			params = {
				'from':mapping_from,
				'to':mapping_to,
				'ids':" ".join(accessions)
			}

			logger.debug(params)

			sessionDownloaderObj = utilities_downloader.sessionDownloader()

			resp =sessionDownloaderObj.session.post(url = url, data = params)
			content = resp.json()
			job_id = content['jobId']

			## CHECK STATUS
			url_status = "https://rest.uniprot.org/idmapping/status/" + job_id

			is_running = True
			while is_running:
				resp = sessionDownloaderObj.session.get(url_status)
				content = resp.json()
				logger.debug(content)
				if "jobStatus" in content:
					if content['jobStatus'].lower() == 'running':
						time.sleep(2)
					else:
						is_running = False
				else:
					is_running = False


			url = "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/" + job_id + "?fields=" + columns + "&format=" + format

			## GRAB RESULTS
			sessionDownloaderObj = utilities_downloader.sessionDownloader()
			status = sessionDownloaderObj.download_file(url,out_path,method='GET',replace_empty=True,remake_age=self.options['remake_age'])
			# print(status)

		if os.path.exists(out_path):

			mappings_data = open(out_path).read()
			mappings = {}

			lines = mappings_data.strip().split("\n")
			headers = lines[0].split("\t")

			# print("headers",headers)
			for line in lines[1:]:
				if line.split("\t") == 2:
					input_identifier = line.split("\t")[0]
					if input_identifier not in mappings: mappings[input_identifier] =[]
					mappings[input_identifier].append(line.split("\t")[1])
				else:
					input_identifier = "-"
					line_bits = line.split("\t")
					if 'From' in headers:
						input_idx = headers.index("From")
						input_identifier = line_bits[input_idx]

					if input_identifier not in mappings:
						mappings[input_identifier]= []
					mappings_tmp = {}
					for i in range(0,len(headers)):
						if len(line_bits) > i:

							if headers[i].split(":")[0] == "yourlist":
								input_identifier = line_bits[i]
							elif headers[i].split(":")[0] == "isomap":
								pass
							else:
								if columns == "accession":
									mappings_tmp = line_bits[i]
								else:
									mappings_tmp[headers[i]] = line_bits[i]

					mappings[input_identifier].append(mappings_tmp)

			if best_hit:
				for k, v in mappings.items():
					mappings[k] = v[0]
			return mappings
		else:
			return {"status":"Error"}

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def annotate_list_uniprot(self,identifiers,mapping_from=None):

		mapping_to="UniProtKB"
		format='tsv'

		if mapping_from != None:
			if mapping_from in ["ID","UniProtKB_AC-ID","ACC+ID"]:
				mapping_from = 'UniProtKB_AC-ID'
				mappings = self.get_uniprot_accession_list(identifiers,mapping_from='UniProtKB_AC-ID')

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def grabUniProtByPDB(self,pdb_id):
		self.check_directory()

		url = "https://rest.uniprot.org/uniprotkb/stream?&format=tsv&query=%28xref%3Apdb-" + pdb_id + "%29" + "&fields=accession,id,reviewed,protein_name,gene_names,organism_name"
		out_path = os.path.join(self.options["data_path"],"uniprot","pdb",pdb_id + ".uniprot.tdt")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method='GET',replace_empty=True)

	"""
	def parseUniProtByPDBs(self,pdb_ids,detailed=False):
		if detailed == False:
			protein = []
			for pdb_id in pdb_ids:
				protein += self.parseUniProtByPDB(pdb_id).keys()

			return list(set(protein))
		else:
			return self.parseUniProtByPDB(pdb_id)


	def parseUniProtByPDB(self,pdb_id,detailed=False):
		self.grabUniProtByPDB(pdb_id)

		tdt_path = os.path.join(self.options["data_path"],"uniprot","pdb",pdb_id + ".uniprot.tdt")
		proteins = []

		if os.path.exists(tdt_path):
			if detailed:
				uniprotDownloaderObj = uniprotDownloader.uniprotDownloader()
				uniprotDownloaderObj.options["data_path"] = self.options["data_path"]

			tdt_data = open(tdt_path).read().strip().split("\n")

			headers = tdt_data[0].split("\t")


			for protein in tdt_data[1:]:
				accession = protein.split("\t")[0]
				proteins.append(accession)

			if len(proteins) > 20: return {"status":"Failed - Too many proteins","data":{}}

			if detailed:
				annotation_json_data = {}
				for accession in proteins:

					protein_data = uniprotDownloaderObj.parseBasic(accession)

					if protein_data['status'] == "Success":
						annotation_json_data[accession] = protein_data['data']

				return {"status":"Success","data":annotation_json_data}
			else:
				return proteins
		else:
			print(("error",pdb_id))
			return proteins
	"""
	##------------------------------------------------------------------##

	def grabUniProtPfamInterpro(self,accession,force=False):
		self.check_directory()

		url = "https://www.ebi.ac.uk/interpro/api/entry/all/pfam/protein/uniprot/" + accession + "/"
		out_path = os.path.join(self.options["data_path"], "uniprot","pfam",accession + ".pfam.json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",JSON=True,replace_empty=True)
		

	def grabUniProtPfam(self,accession,force=False):
		self.check_directory()
		url = "http://pfam.xfam.org/protein/" + accession + "?output=xml"
		out_path = os.path.join(self.options["data_path"], "uniprot","pfam",accession + ".pfam.xml")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True)

	##------------------------------------------------------------------##

	def parseUniProtPfam(self,accession,region_start=None,region_end=None):
		if region_start == None:
			if 'region_start' in self.options:region_start = self.options['region_start']
		if region_end == None:
			if 'region_end' in self.options:region_end = self.options['region_end']
		
		self.grabUniProtPfamInterpro(accession)

		pfamDownloaderObj = pfamDownloader.pfamDownloader()
		pfamDownloaderObj.options["data_path"] = self.options["data_path"]

		json_path = os.path.join(self.options["data_path"], "uniprot","pfam",accession + ".pfam.json")

		domains_data = {}

		try:
			if os.path.exists(json_path):
				try:
					with open(json_path) as outfile:
						pfam_json_content = json.load(outfile)
				except:
					pfam_json_content = {}

				"""
				if len(pfam_json_content) == 0:
					os.remove(json_path)
					self.grabUniProtPfamInterpro(accession)
					
					with open(json_path) as outfile:
						try:
							pfam_json_content = json.load(outfile)
						except:
							pfam_json_content = {}

					if len(pfam_json_content) == 0:
						return  {"status":"Error","error_type":"no data returned"}
				"""
				protein_sequence = self.parseSequence(accession)

				for domain_type in pfam_json_content['results']:
					
					if not domain_type['proteins'][0]['entry_protein_locations']:
						continue			
					pfam_accession = domain_type['metadata']['accession']
					pfam_id = domain_type['metadata']['accession']
					pfam_type = domain_type['metadata']['type']

					add_domain = True
					for domain_instance in domain_type['proteins'][0]['entry_protein_locations']:
						for fragment in domain_instance['fragments']:
							start = fragment['start']
							end = fragment['end']
							evalue = domain_instance['score']

							if region_start == None and region_end == None:
								if pfam_accession not in domains_data:
									domains_data[pfam_accession] = {"matches":[]}

								domains_data[pfam_accession]["matches"].append({"start":start,"end":end,"evalue":evalue,'sequence':protein_sequence[int(start):int(end)]})
							else:
								try:
									domain_length = int(end) - int(start)
									region_length = int(region_end) - int(region_start)

									overlap =  set(range(int(start),int(end))).intersection(list(range(int(region_start),int(region_end))))

									if region_length == 0 or domain_length == 0:
										add_domain = False
									else:
										if float(len(overlap))/domain_length > 0.75 or float(len(overlap))/region_length > 0.75:
											if pfam_accession not in domains_data:
												domains_data[pfam_accession] = {"matches":[]}

											domains_data[pfam_accession]["matches"].append({"start":start,"end":end,"evalue":evalue,'sequence':protein_sequence[int(start):int(end)],"region_overlap":float(len(overlap))/domain_length,"region_ratio":float(region_length)/domain_length})
											add_domain = True
										else:
											add_domain = False
								except:
									add_domain = False


						if add_domain:
							if "pfam_id" not in domains_data[pfam_accession]:

								pfam_details = pfamDownloaderObj.parsePfam(pfam_accession)
								
								domains_data[pfam_accession]["pfam_accession"] = pfam_accession
								domains_data[pfam_accession]["pfam_type"] = pfam_type
								
								if 'data' in pfam_details:
									for tag in ['clan_id','clan_name','pfam_acc','domain_name','domain_description']:
										domains_data[pfam_accession][tag] = pfam_details['data'][tag]

								if 'pfam_acc' in domains_data[pfam_accession]:
									domains_data[pfam_accession]["pfam_id"] = domains_data[pfam_accession]['pfam_acc']
								else:
									domains_data[pfam_accession]["pfam_id"] = pfam_id
									
				return {"status":"Added","data":domains_data}
			else:
				return {"status":"Error","error_type":"File not found","data":{}}

		except Exception as e:
			print("Error - parseUniProtPfam @", e,accession,json_path)
			return {"status":"Error","error_type":str(e),"data":{}}


	def parseUniProtPfamInstances(self,accession):
		domain_data = self.parseUniProtPfam(accession)
		domain_instances_data = []
		if 'data' in domain_data:
			for pfam_id in domain_data['data']:
				for match in domain_data['data'][pfam_id]['matches']:
					try:
						domain_instances_data.append({
							'domain_start':match['start'],
							'domain_end':match['end'],
							'domain_domain_name':domain_data['data'][pfam_id]['domain_name'],
							'domain_pfam_accession':domain_data['data'][pfam_id]['pfam_accession'],
							'domain_pfam_id':domain_data['data'][pfam_id]['pfam_id']
						})
					except:
						pass

		return domain_instances_data

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def grabUniProtFasta(self,accession,force=False):
		self.check_directory()

		url = "https://rest.uniprot.org/uniprotkb/" + accession + ".fasta"
		out_path = os.path.join(self.options["data_path"] , "uniprot","fasta",accession + ".fasta")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,force=force,remake_age=self.options['remake_age'])

		return status

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def parseUniProtFasta(self,accession):
		out_path = os.path.join(self.options["data_path"] , "uniprot","fasta",accession + ".fasta")
		protein_data = {}
		gnRe = re.compile("GN=\S+")

		if not os.path.exists(out_path):
			self.grabUniProtFasta(accession)

		if os.path.exists(out_path):
			try:
				line_bits = open(out_path).read().split("\n")

				sequence = "".join(line_bits[1:])

				accessions = line_bits[0].split()[0].split("|")


				if len(accessions) > 1:
					accession = accessions[1]
					identifier = accessions[2]

					fullname = " ".join(line_bits[0].split("OS=")[0].split()[1:])
					gene = ""
					geneBits = gnRe.findall(line_bits[0])
					if len(geneBits) > 0:
						gene = geneBits[0][3:]
				else:
					accession = line_bits[0][1:].split()[0]
					identifier = accession
					gene = accession
					fullname = accession

				protein_data = {
					"identifier":identifier,
					"gene":gene,
					"fullname":fullname,
					"sequence":sequence,
					"header":line_bits[0][1:]
					}
			except Exception as e:
				protein_data = {"status":"Error","error_type":str(e)}

		return protein_data


	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def get_protein_topology(self,accession):

		protein_topology = []

		topological_domains = [
			"Cytoplasmic",
			"Nuclear",
			"Extracellular",
			"Lumenal",
			"Mitochondrial intermembrane",
			"Perinuclear space",
			"Signal Peptide"
		]

		try:
			features_response = self.parseFeatures(accession)

			if 'initiator methionine' in features_response['data']['features']:
				for region in features_response['data']['features']['initiator methionine' ]:
					if region['type'] == 'initiator methionine':
						protein_topology.append(region)

			if 'transmembrane region' in features_response['data']['features']:
				for region in features_response['data']['features']['transmembrane region']:
					if region['type'] == 'transmembrane region':
						protein_topology.append(region)

			if 'intramembrane region' in features_response['data']['features']:
				for region in features_response['data']['features']['intramembrane region']:
					if region['type'] == 'intramembrane region':
						region['description'] = 'intramembrane region'
						protein_topology.append(region)

			if 'signal peptide' in features_response['data']['features']:
				for region in features_response['data']['features']['signal peptide']:
					if region['type'] == 'signal peptide':
						protein_topology.append(region)

			if 'propeptide' in features_response['data']['features']:
				for region in features_response['data']['features']['propeptide']:
					if region['type'] == 'propeptide':
						protein_topology.append(region)

			if 'topological domain' in features_response['data']['features']:
				for region in features_response['data']['features']['topological domain']:
					if region['type'] == 'topological domain':
						#print(region)
						for topological_domain in topological_domains:
							if region['description'].count(topological_domain) > 0:
								region['description'] = topological_domain
								protein_topology.append(region)
		except:
			logger.error(accession)

		return protein_topology

	def get_protein_topology_string(self,accession):

		protein_sequence = self.parseBasic(accession)['data']['sequence']
		protein_topology = ['-']*len(protein_sequence)

		topological_domains = {
			"Cytoplasmic":"-",
			"Nuclear":"-",
			"Extracellular":"E",
			"Lumenal":"L",
			"Mitochondrial intermembrane":"I",
			"Perinuclear space":"P",
			"signal peptide":"S",
			"propeptide":"s",
			"transmembrane region":"M",
			"intramembrane region":"I",
			"initiator methionine":"m"
		}

		try:
			protein_topology_response = self.get_protein_topology(accession)
			for region in protein_topology_response:
				try:
					for i in range(int(region['start'])-1,int(region['end'])):

						if region['type'] in topological_domains:
							protein_topology[i] = topological_domains[region['type']]
						elif region['description'] in topological_domains:
							protein_topology[i] = topological_domains[region['description']]
						else:
							logger.error(region,topological_domains.keys())
				except:
					logger.debug(region)
		except:
			logger.error("Can't make protein topology string")

		return "".join(protein_topology)

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def check_uniprot_xml_file_exists(self,accession):
		xml_path = os.path.join(self.options["data_path"] , "uniprot", "xml",accession + ".xml")

		if os.path.exists(xml_path):
			return True
		else:
			return False

	def check_uniprot_fasta_file_exists(self,accession):
		fasta_path = os.path.join(self.options["data_path"] , "uniprot", "fasta",accession + ".fasta")

		if os.path.exists(fasta_path):
			return True
		else:
			return False

	def check_uniprot_versions(self,accession,accession_version_details,force=False):
		accession_status = {}

		uniprot_xml_file_exists = self.check_uniprot_xml_file_exists(accession)
		if not uniprot_xml_file_exists:
			logger.info(accession + " xml file not found")
			accession_status["xml_file_status"] = "not found"
			return accession_status
		else:
			accession_status["xml_file_status"] = "found"

		uniprot_fasta_file_exists = self.check_uniprot_fasta_file_exists(accession)
		if not uniprot_fasta_file_exists:
			logger.info(accession + " fasta file not found")
			accession_status["fasta_file_status"] = "not found"
			return accession_status
		else:
			accession_status["fasta_file_status"] = "found"

		uniprot_info = self.parseBasic(accession)
		uniprot_version_local = uniprot_info['data']['version']
		uniprot_sequence_version_local = uniprot_info['data']['sequence_version']

		accession_status["uniprot_version_local"] = uniprot_version_local
		accession_status["uniprot_sequence_version_local"] = uniprot_sequence_version_local
		accession_status["uniprot_version_remote"] = accession_version_details['version']
		accession_status["uniprot_sequence_version_remote"] = accession_version_details['sequence_version']
		accession_status["uniprot_version_uptodate"] = accession_version_details['version'] == uniprot_version_local
		accession_status["uniprot_sequence_version_uptodate"] = accession_version_details['sequence_version'] == uniprot_sequence_version_local

		if not (accession_status["uniprot_version_uptodate"]):
			if self.options["update"]:
				logger.debug("Version mismatch - Updating xml for " + accession )
				self.grabUniProt(accession,force=True)
				accession_status["xml_file_status"] = "updated"
				uniprot_info = self.parseBasic(accession)
				fasta_data = self.parseUniProtFasta(accession)
				if fasta_data['sequence'] != uniprot_info['data']['sequence']:
					logger.debug("Updating fasta for " + accession )
					self.grabUniProtFasta(accession,force=True)
					accession_status["fasta_file_status"] = "updated"
		else:
			if self.options["update"]:
				fasta_data = self.parseUniProtFasta(accession)
				if fasta_data['sequence'] != uniprot_info['data']['sequence']:
					logger.debug("Sequence mismatch - Updating fasta for " + accession )
					self.grabUniProtFasta(accession,force=True)
					accession_status["fasta_file_status"] = "updated"

		logger.info("\t".join([
			accession,
			accession_status["uniprot_version_local"],
			accession_status["uniprot_version_remote"],
			accession_status["uniprot_sequence_version_local"],
			accession_status["uniprot_sequence_version_remote"],
			str(accession_status["uniprot_version_uptodate"]),
			str(accession_status["uniprot_sequence_version_uptodate"])
			]))

		return accession_status

	def check_uniprot_versions_taxa(self,taxon_id):

		accession_status = {}
		accession_version_details = self.parse_uniprot_accession_taxa(taxon_id,reviewed=self.options['reviewed'],structures=False,nofragment=True,include_version=True)

		self.options['remake'] = False
		self.options["remake_age"] = 100000

		logger.debug(str(len(accession_version_details)) + " uniprot accessions")

		for accession in accession_version_details:
			accession_status[accession] = self.check_uniprot_versions(accession,accession_version_details[accession])

		return accession_status

	def checkUniProtEntry(self,accession,force=False):
		try:
			self.check_directory()
			
			if self.options['remake'] == False:
				out_xml_path = os.path.join(self.options["data_path"] , "uniprot","xml",accession + ".xml")
				if os.path.exists(out_xml_path):
					title_pattern = re.compile("<title>(.*?)</title>")
					title = title_pattern.findall(open(out_xml_path).read())
					
					if len(title) > 0:
						if title[0] == 'Error':
							return {"status": "Error", "error_type": "Uniprot website"}
				else:
					return{"status": "Error", "error_type": "No file"}

				if os.path.getsize(out_xml_path) == 0:
					return{"status": "Success", "error_type": "None"}


			out_path = os.path.join(self.options["data_path"] , "uniprot","obsolete",accession + ".json")

			if os.path.exists(out_path):
				if force or self.options['remake'] :
					logger.debug("Forced - Deleting: " + out_path)
					os.remove(out_path)

			if os.path.exists(out_path):
				if (time.time() - os.path.getctime(out_path))/60/60/24 > self.options["remake_age"] :
					if self.options["verbose"]:print(("Outdated - Deleting",out_path,(time.time() - os.path.getctime(out_path))/60/60/24))

					logger.debug("Outdated - Deleting: " + out_path + " " + str((time.time() - os.path.getctime(out_path))/60/60/24))
					os.remove(out_path)

			if not os.path.exists(out_path):
				logger.debug("Grabbing: " + out_path)

				try:
					url = "https://rest.uniprot.org/uniprotkb/search?&query=" + accession + " AND active:false"

					sessionDownloaderObj = utilities_downloader.sessionDownloader()
					status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True,JSON=True)
				except Exception as e:
					logger.error(str(e))

			try:
				try:
					with open(out_path) as outfile:
						status = json.load(outfile)
				except:
					logging.error(out_path + " format error")
					if os.path.exists(out_path):
						os.remove(out_path)

				if len(status['results']) == 0:
					return {"status":"Success", "error_type": "None"}

				if 'inactiveReason' in status['results'][0]:
					if status['results'][0]['inactiveReason']['inactiveReasonType'] == "DELETED":
						return {"status":"Error","error_type":"obsolete","data":status['results'][0]}
					if 'mergeDemergeTo' in status['results'][0]['inactiveReason']:
						response = status['results'][0]['inactiveReason']
						updated_accession = response['mergeDemergeTo'][0]
						response.update({"status":"Error","error_type":"renamed", "updated_accession":updated_accession})

						return response
				else:
					return {"status":"Success", "error_type": "None","data":status['results'][0]}
			except Exception as e:
			

				return  {"status":"Error","error_type":str(e)}
		except Exception as e:
			return  {"status":"Error","error_type":str(e)}

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def grabProteinHistory(self,accession):
		self.check_directory()
		url = "https://rest.uniprot.org/unisave/" + accession + "?format=json"
		out_path = os.path.join(self.options["data_path"] , "uniprot","history",accession + ".json")

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,method="GET",replace_empty=True,JSON=True,remake_age=1)


	def parseProteinHistory(self, accession):
		self.grabProteinHistory(accession = accession)
		out_path = os.path.join(self.options["data_path"] , "uniprot","history",accession + ".json")

		f = open(out_path, 'r')
		data = json.loads(f.read())
		f.close()
		return data



	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def grabUniProt(self,accession,force=False):
		self.check_directory()

		url = "https://www.uniprot.org/uniprot/" + accession + ".xml"
		out_path = os.path.join(self.options["data_path"] , "uniprot","xml",accession + ".xml")

		if os.path.exists(out_path):
			if force:
				if self.options["verbose"]: print("Forced - Deleting",out_path)
				os.remove(out_path)

			elif (time.time() - os.path.getctime(out_path))/60/60/24 > self.options["remake_age"]:
				if self.options["verbose"]:print("Outdated - Deleting",out_path,(time.time() - os.path.getctime(out_path))/60/60/24)
				os.remove(out_path)

			# elif os.path.getsize(out_path) == 0:
			else:
				entry_check = self.checkUniProtEntry(accession)

				if entry_check['error_type'] == "obsolete":
					if self.options["verbose"]: print("Obsolete accession:",accession,entry_check)
					return entry_check
				elif entry_check['error_type'] == "not found":
					if self.options["verbose"]: print("Not found accession:",accession)
					return entry_check
				elif entry_check['error_type'] == "renamed":
					if self.options["verbose"]: print("Renamed accession:",accession,"->",entry_check['updated_accession'])
					return entry_check
				elif entry_check['error_type'] == 'Uniprot website':
					if self.options['verbose']: print("UniProt website is not available. Try to download again every 10 secs for 0.5 h, then simply give up.")
					is_downloading = True
					start_time = time.time()
					while is_downloading:
						if os.path.exists(out_path): os.remove(out_path)
						time.sleep(10)
						if self.options['verbose']: print("Redownloading protein file:", accession)
						sessionDownloaderObj = utilities_downloader.sessionDownloader()
						status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,remake_age=self.options['remake_age'])
						entry_recheck = self.checkUniProtEntry(accession)
						if entry_recheck['error_type'] != 'Uniprot website':
							is_downloading = False
							if self.options['verbose']: print("Redownloading protein file - successful.")
						if time.time() - start_time > 1800:
							if self.options['verbose']: print("Redownloading protein file - failed. Uniprot website still is not available.")
							is_downloading = False
							return entry_recheck
				else:
					if os.path.getsize(out_path) == 0:
						if self.options["verbose"]: print("Empty - Deleting",out_path,(time.time() - os.path.getctime(out_path))/60/60/24)
						os.remove(out_path)

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,remake_age=self.options['remake_age'])

		if status['status'] == 'Error':
			status = self.checkUniProtEntry(accession)
		else:
			if os.path.getsize(out_path) < 1000:
				logger.debug("File exists but is empty: " + str(os.path.getsize(out_path)))
				status = self.checkUniProtEntry(accession)

		return status

	##------------------------------------------------------------------##
	def parseGO(self,accession,parse_keywords=True,generic=False):
		data = self.parseUniProt(accession,parse_keywords=True,parse_generic=True)

		output = {"go":{}}

		if 'data' in data:
			if 'GO' in data['data']:
				output['go'] = data['data']['GO']

		return output

	def parseGODetailed(self,accession):
		if not os.path.exists(os.path.join(self.options["data_path"], "uniprot","go")):
			os.mkdir(os.path.join(self.options["data_path"], "uniprot","go"))

		out_path = os.path.join(self.options["data_path"] , "uniprot","go",accession + ".json")
		url = "http://golr-aux.geneontology.io/solr/select?fq=document_category:%22annotation%22&q=*:*&fq=bioentity:%22UniProtKB:" + accession + "%22&wt=json&rows=10000"

		sessionDownloaderObj = utilities_downloader.sessionDownloader()
		status = sessionDownloaderObj.download_file(url,out_path,replace_empty=True,remake_age=self.options['remake_age'])

		with open(out_path) as outfile:
			json_data = json.load(outfile)

		go_data = {
			"terms":{},
			"hierarchy_terms":{}
		}

		for go_term in json_data['response']['docs']:
			if go_term['annotation_class'] in go_data["terms"]: continue

			go_data["terms"][go_term['annotation_class']] = {"name":go_term['annotation_class_label'],"class":go_term["aspect"],"is_a":{}}

			go_data['hierarchy_terms'][go_term['annotation_class']] = {"name":go_term['annotation_class_label'],"type":"annotated","descendent_name":[],"descendent_id":[],"class":go_term["aspect"]}


			for i in range(0,len(go_term["isa_partof_closure"])):
				go_data["terms"][go_term['annotation_class']]["is_a"][go_term["isa_partof_closure"][i]] = {"name":go_term["isa_partof_closure_label"][i]}

				if go_term["isa_partof_closure"][i] not in go_data['hierarchy_terms']:
					go_data['hierarchy_terms'][go_term["isa_partof_closure"][i]] = {"name":go_term["isa_partof_closure_label"][i],"type":"ancestor","descendent_name":[go_term['annotation_class_label']],"descendent_id":[go_term['annotation_class']],"class":go_term["aspect"]}
				else:
					go_data['hierarchy_terms'][go_term["isa_partof_closure"][i]]["descendent_name"].append(go_term['annotation_class_label'])
					go_data['hierarchy_terms'][go_term["isa_partof_closure"][i]]["descendent_id"].append(go_term['annotation_class'])

		return go_data


	##------------------------------------------------------------------##

	def parseLocalisation(self,accession,parse_keywords=True,generic=False):
		data = self.parseUniProt(accession,parse_keywords=True,parse_generic=True)

		output = {}

		if 'data' in data:
			localisation = {}
			if 'localisation' in data['data']:
				localisation["localisation_keyword"] = data['data']['localisation']

				if data['data']['localisation_isoform'] != {}:
					localisation["localisation_isoform"] = data['data']['localisation_isoform']


				localisation["basic_localisation"] = []
				
				for subcellular_compartment_basic in ['Nucleus','Cytoplasm',"Membrane","Secreted","Mitochondrion"]:
					for subcellular_compartment in localisation["localisation_keyword"]:
						if subcellular_compartment.lower().count(subcellular_compartment_basic.lower()):
							localisation["basic_localisation"].append(subcellular_compartment_basic)
				
				localisation["basic_localisation"] = list(set(localisation["basic_localisation"]))

			localisation["Transmembrane"] = False
			if 'keywords' in data['data']:
				if 'Transmembrane helix' in data['data']['keywords'] or 'Transmembrane' in data['data']['keywords'] or 'Transmembrane' in data['data']['localisation']:
					localisation["basic_localisation"].append('Transmembrane')
					localisation["Transmembrane"] = True


			for header in ['id','protein_name','family','gene_name']:
				if header in data['data']:
					output[header] = data['data'][header]

			localisation["localisation_go"] = {}
			if 'GO' in data['data']:
				if 'C' in data['data']['GO']:
					for compartment in data['data']['GO']['C']:
						localisation["localisation_go"][compartment['term']] = compartment

			output = {"localisation":localisation}

		return output

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##


	def parseDiseases(self,accession,generic=False):
		data = self.parseUniProt(accession,parse_keywords=True,parse_generic=generic)

		try:
			if 'diseases' in data['data']:
				return data['data']['diseases']
		except:
			pass

		return {}

	def parseKeywords(self,accession,generic=False):
		return self.parseUniProt(accession,parse_keywords=True,parse_generic=generic)

	def parseIsoforms(self,accession,generic=False):
		return self.parseUniProt(accession,parse_isoforms=True,parse_generic=generic)

	def parseSNPs(self,accession,generic=False):
		return self.parseUniProt(accession,parse_snps=True,parse_generic=generic)

	def parseSequence(self,accession):
		return self.parseUniProt(accession,parse_sequence=True)

	def parseFeatures(self,accession,generic=False):
		return self.parseUniProt(accession,parse_features=True,parse_generic=generic)

	def parseSecondaryStructure(self,accession):
		self.options['use_features'] = ['helix','strand','turn']
		data = self.parseUniProt(accession,parse_features=True,parse_generic=False)
		return data['data']['features']

	def parseMutagenesis(self,accession):
		self.options['use_features'] = ['mutagenesis site']
		data = self.parseUniProt(accession,parse_features=True,parse_generic=False)
		return data['data']['features']

	def parsePTMs(self,accession):
		self.options['use_features'] = ['modified residue']
		data = self.parseUniProt(accession,parse_features=True,parse_generic=False)
		return data['data']['features']

	def parseRegionsOfInterest(self,accession):
		self.options['use_features'] = ['region of interest','short sequence motif','metal ion-binding site','site','DNA-binding region']
		data = self.parseUniProt(accession,parse_features=True,parse_generic=False)
		return data['data']['features']

	def parseDomains(self,accession):
		data = self.parseUniProt(accession)

		if 'data' in data:
			if 'Pfam' in data['data']:
				return {"status":"Success",'data':data['data']['Pfam']}
			else: return {"status":"Success",'data':{}}
		else: return {"status":"Error",'data':{}}


	def parsePDBs(self,accession):
		if not self.check_accession(accession):
			return {"status":"Error",'data':"Not a UniProt Accession:" + accession}

		data = self.parseUniProt(accession,parse_generic=True)

		if 'data' in data:
			if 'PDB' in data['data']:
				return {"status":"Success",'data':data['data']['PDB']}
			else: return {"status":"Success",'data':{}}
		else: return data

	def parseBasic(self,accession,basic_data_types = ["accession","description","family",'gene_name','gene_names','id','protein_name','protein_names','sequence','species_common','species_scientific','taxon_id','taxonomy','version','modified','reviewed','sequence_version']):
		if not self.check_accession(accession):
			return {"status":"Error",'data':"Not a UniProt Accession" + accession}

		data = self.parseUniProt(accession,parse_generic=True)
		basic_data = {}
		if 'data' in data:
			for data_type in basic_data_types:
				if data_type in data['data']:
					basic_data[data_type] =  data['data'][data_type]

			return {"status":"Success",'data':basic_data}
		else: return {"status":"Error",'data':data}

	def parseOrganismHost(self,accession):
		if not self.check_accession(accession):
			return {"status":"Error",'data':"Not a UniProt Accession" + accession}

		data = self.parseUniProt(accession,parse_generic=True)
		basic_data = {}
		if 'data' in data:
			organism_host = []
			if "organism_host" in data['data']:
				organism_host = data['data']['organism_host']

			return {"status":"Success",'data':organism_host}
		else: return {"status":"Error",'data':[]}


	def parseName(self,accession):
		return self.parseBasic(accession,basic_data_types = ["accession",'gene_name','protein_name'])

	def parseDbXref(self,accession,db_names=None):
		protein_data = self.parseUniProt(accession,parse_db_reference=True)
		if db_names == None:
			db_xrefs = protein_data['data']['db_xref']
		else:
			db_xrefs = {}
			for db_name in db_names:
				if db_name in protein_data['data']['db_xref']:
					db_xrefs[db_name] = protein_data['data']['db_xref'][db_name]

		return db_xrefs

	def parsePathways(self,accession):
		return self.parseDbXref(accession,db_names=['Reactome'])

	def parseDrugbank(self,accession):
		return self.parseDbXref(accession,db_names=['DrugBank'])

	def parseSecondaryAccessions(self,accession):
		data = self.parseUniProt(accession,parse_generic=True)

		if 'data' in data:
			if "secondary_accessions" in data['data']:
				return {"status":"Success",'data':{"secondary_accessions":data['data']["secondary_accessions"]}}
			else:
				return {"status":"Error",'data':{}}
		else: return {"status":"Error",'data':{}}

	def parseUniProtRegion(self,accession,start,end,generic=False):
		start = int(start)
		end = int(end)

		data = self.parseUniProt(accession,parse_generic=True)

		region_data = {}
		region_data['sequence'] = data['data']['sequence'][start:end]
		protein_data = self.parseUniProt(accession,parse_features=True)
		domain_data = self.parseUniProtPfam(accession,start,end)

		region_data['PDB'] = {}
		region_data['Pfam'] = {}
		region_data['features'] = {}

		if 'data' in protein_data:
			if 'PDB' in protein_data['data']:
				for pdb in protein_data['data']['PDB']:
					tmp_data = {}
					for chain in protein_data['data']['PDB'][pdb]:
						region_start = int(protein_data['data']['PDB'][pdb][chain]['start'])
						region_end = int(protein_data['data']['PDB'][pdb][chain]['end'])

						overlap =  set(range(int(start),int(end))).intersection(list(range(int(region_start),int(region_end))))

						if len(overlap) > 0:
							tmp_data[chain] = protein_data['data']['PDB'][pdb][chain]

					if len(tmp_data) > 0:
						region_data['PDB'][pdb] = tmp_data

			if 'features' in protein_data['data']:
				for feature_type in protein_data['data']['features']:
					tmp_data = []
					for feature in protein_data['data']['features'][feature_type]:

						region_start = int(feature['start'])
						region_end = int(feature['end'])
						overlap =  set(range(int(start),int(end))).intersection(list(range(int(region_start),int(region_end))))

						if len(overlap) > 0:
							tmp_data.append(feature)

					if len(tmp_data) > 0:
						region_data['features'][feature_type] = tmp_data

		if 'data' in domain_data:
			for domain in domain_data['data']:
				tmp_data = []
				for match in domain_data['data'][domain]['matches']:

					region_start = int(match['start'])
					region_end = int(match['end'])

					overlap =  set(range(int(start),int(end))).intersection(list(range(int(region_start),int(region_end))))

					if len(overlap) > 0:
						tmp_data.append(match)

				if len(tmp_data) > 0:
					region_data['Pfam'][domain] = copy.deepcopy(domain_data['data'][domain])
					region_data['Pfam'][domain]['matches'] = tmp_data

		return {"status":"Success",'data':region_data}

	def parseUniProtBulk(self,accessions):
		if not isinstance(accessions,(list)):
			accessions = accessions.split(",")

		protein_data = {}
		for accession in accessions:
			try:
				protein_data[accession] = self.parseUniProt(accession)
			except Exception as e:
				protein_data[accession]  = {"status":"Error",'error_type':str(e)}

		return protein_data


	def parseUniProtByPDB(self,pdb_id):
		import pdbDownloader
		dataDownloaderObj = pdbDownloader.pdbDownloader()
		dataDownloaderObj.options['pdbID'] = pdb_id
		accessions = dataDownloaderObj.getPDBProteins()
		return self.parseUniProtBulk(accessions)

	def parseUniProt(self,accession,parse_features=False,parse_names=False,parse_keywords=False,parse_snps=False,parse_attributes=False,parse_disorder=False,parse_isoforms=False,parse_sequence=False,parse_generic=True,parse_db_reference=False, parse_caution=True,force=False):

		if not self.check_accession(accession):
			return {"status":"Error",'data':"Not a UniProt Accession" + accession}

		if 'parse_features' in self.options:
			parse_features = self.options['parse_features']
		if 'parse_names' in self.options:
			parse_names = self.options['parse_names']
		if 'parse_keywords' in self.options:
			parse_keywords = self.options['parse_keywords']
		if 'parse_snps' in self.options:
			parse_snps = self.options['parse_snps']
		if 'parse_attributes' in self.options:
			parse_attributes = self.options['parse_attributes']
		if 'parse_disorder' in self.options:
			parse_disorder = self.options['parse_disorder']
		if 'parse_isoforms' in self.options:
			parse_isoforms = self.options['parse_isoforms']
		if 'parse_sequence' in self.options:
			parse_sequence = self.options['parse_sequence']
		if 'parse_generic' in self.options:
			parse_generic = self.options['parse_generic']
		if 'parse_caution' in self.options:
			parse_caution = self.options['parse_caution']

		#accession = self.check_isoform(accession)

		download_status = self.grabUniProt(accession,force)

		if download_status["status"] == "Error":
			if "updated_accession" in download_status:

				logging.error("Accession updated: " + accession + " " + download_status["updated_accession"]) 
				logging.error(download_status)
				
				accession = download_status["updated_accession"]

				if not self.check_accession(accession):
					return {"status":"Error",'data':"Not a UniProt Accession" + accession}
			else:
				return download_status

		proteins_data = {
		"accession":accession,
		"isoform":False
		}

		if len(accession.split("-")) > 1:
			proteins_data["isoform"] = True

		xml_path = os.path.join(self.options["data_path"] , "uniprot", "xml",accession + ".xml")

		error_pattern = re.compile("<error>.+</error>")

		try:
			if os.path.exists(xml_path):
				if len(open(xml_path).read()) == 0:
					return  {"status":"Error","error_type":"no data returned"}

				tree = elementtree.parse(xml_path)
				root = tree.getroot()

				for entry in tree.iter('{http://uniprot.org/uniprot}entry'):
					#----

					references = {}
					for reference in tree.iter('{http://uniprot.org/uniprot}evidence'):
						for dbReference in reference.iter('{http://uniprot.org/uniprot}dbReference'):
							if dbReference.attrib['type'] == 'PubMed':
								references[reference.attrib['key']] = dbReference.attrib['id']

					#----

					if parse_generic or parse_names:
						#try:
						#	proteins_data["accession"] = entry.find('{http://uniprot.org/uniprot}accession').text#version
						#except:
						#	proteins_data["accession"] = ""

						try:
							proteins_data["modified"] = entry.attrib['modified'] #modified date (timestamp for database)
						except:
							proteins_data["modified"] = ""

						try:
							proteins_data["reviewed"] = True if entry.attrib['dataset']=='Swiss-Prot' else False #reviewed
						except:
							proteins_data["reviewed"] = ""


						try:
							proteins_data["id"] = entry.find('{http://uniprot.org/uniprot}name').text#version
						except:
							proteins_data["id"] = ""

						try:
							proteins_data["gene_name"] = entry.find('.//{http://uniprot.org/uniprot}gene/{http://uniprot.org/uniprot}name').text
						except:
							proteins_data["gene_name"] = ""

						try:
							proteins_data["protein_name"] = entry.find('.//{http://uniprot.org/uniprot}recommendedName/{http://uniprot.org/uniprot}fullName').text
						except:
							proteins_data["protein_name"] = ""

						if proteins_data["protein_name"] == "":
							try:
								proteins_data["protein_name"] = entry.find('.//{http://uniprot.org/uniprot}submittedName/{http://uniprot.org/uniprot}fullName').text
							except:
								proteins_data["protein_name"] = ""

						try:
							proteins_data["species_common"] = entry.find("{http://uniprot.org/uniprot}organism/{http://uniprot.org/uniprot}name[@type='common']").text#http://uniprot.org/uniprot}organism/{http://uniprot.org/uniprot}name").text
						except:
							proteins_data["species_common"] = ""

						try:
							proteins_data["species_scientific"] = entry.find("{http://uniprot.org/uniprot}organism/{http://uniprot.org/uniprot}name[@type='scientific']").text.split("(")[0]#http://uniprot.org/uniprot}organism/{http://uniprot.org/uniprot}name").text
						except:
							proteins_data["species_scientific"] = ""


						#<<<<<<<<<<<<<<<<<<<<<<
						if parse_names:
							return proteins_data
						#<<<<<<<<<<<<<<<<<<<<<<

						try:
							proteins_data["secondary_accessions"] = []
							for accession in tree.iter('{http://uniprot.org/uniprot}accession'):
								proteins_data["secondary_accessions"].append(accession.text)
						except:
								proteins_data["secondary_accessions"] = []

						try:
							proteins_data["fragment"] = False

							if "fragment" in entry.find('{http://uniprot.org/uniprot}sequence').attrib:
								proteins_data["fragment"] = True
						except:
							proteins_data["fragment"] = False

						try:
							proteins_data["version"] = entry.attrib['version'] #version
						except:
							proteins_data["version"] = ""

						try:
							proteins_data["sequence_version"] =entry.find('{http://uniprot.org/uniprot}sequence').attrib['version']
						except:
							proteins_data["sequence_version"] = ""

						try:
							proteins_data["protein_existence"] = entry.find('{http://uniprot.org/uniprot}proteinExistence').attrib['type']#version
						except:
							proteins_data["protein_existence"] = ""


						try:
							proteins_data["sequence"] = entry.find('{http://uniprot.org/uniprot}sequence').text.replace("\n","")
							if parse_sequence and not proteins_data["isoform"]:
								return proteins_data["sequence"]
						except:
							proteins_data["sequence"] = ""


						try:
							proteins_data["protein_names"] = []
							for protein_name in entry.findall('.//{http://uniprot.org/uniprot}alternativeName/{http://uniprot.org/uniprot}fullName'):
								proteins_data["protein_names"].append(protein_name.text)
						except:
							proteins_data["protein_names"] = []


						try:
							proteins_data["gene_names"] = []
							for gene_name in entry.findall('.//{http://uniprot.org/uniprot}gene/{http://uniprot.org/uniprot}name'):
								proteins_data["gene_names"].append(gene_name.text)
						except:
							proteins_data["gene_names"] = []



						try:
							proteins_data["taxon_id"] = entry.find('{http://uniprot.org/uniprot}organism/{http://uniprot.org/uniprot}dbReference').attrib['id']
						except:
							proteins_data["taxon_id"] = ""

						try:
							proteins_data["description"] = entry.find("{http://uniprot.org/uniprot}comment[@type='function']/{http://uniprot.org/uniprot}text").text
						except:
							proteins_data["description"] = ""

						try:
							proteins_data["caution"] = entry.find("{http://uniprot.org/uniprot}comment[@type='caution']/{http://uniprot.org/uniprot}text").text
						except:
							proteins_data["caution"] = ""

						try:
							proteins_data["sequence caution"] = entry.find("{http://uniprot.org/uniprot}comment[@type='sequence caution']/{http://uniprot.org/uniprot}text").text
						except:
							proteins_data["sequence caution"] = ""

						try:
							proteins_data["family"] = ""
							for comment in entry.findall("{http://uniprot.org/uniprot}comment[@type='similarity']/{http://uniprot.org/uniprot}text"):
								if comment.text.replace('.','').split()[-1].replace('.','') == "family":
									proteins_data["family"] = comment.text.replace('.','').replace('Belongs to the ','')
						except:
							proteins_data["family"] = ""
						
						try:
							proteins_data["cofactor"] = entry.find("{http://uniprot.org/uniprot}comment[@type='cofactor']/{http://uniprot.org/uniprot}cofactor/{http://uniprot.org/uniprot}name").text
						except:
							proteins_data["cofactor"] = ""

						try:
							proteins_data["catalytic_activity"] = entry.find("{http://uniprot.org/uniprot}comment[@type='catalytic activity']/{http://uniprot.org/uniprot}reaction/{http://uniprot.org/uniprot}text").text
						except:
							proteins_data["catalytic_activity"] = ""

						proteins_data["taxonomy"] = []

						for elem in tree.iter('{http://uniprot.org/uniprot}taxon'):
							proteins_data["taxonomy"].append(elem.text)

						for elem in tree.iter('{http://uniprot.org/uniprot}organismHost'):
							if 'organism_host' not in proteins_data:
								proteins_data['organism_host'] = []

							try:
								scientific_name = elem.find("{http://uniprot.org/uniprot}name[@type='scientific']").text.split("(")[0]
							except:
								scientific_name = ""
							try:
								common_name = elem.find("{http://uniprot.org/uniprot}name[@type='common']").text
							except:
								common_name = ""
							try:
								taxon_id = elem.find("{http://uniprot.org/uniprot}dbReference").attrib['id']
							except:
								taxon_id = ""

							host = {
								"scientific_name":scientific_name,
								"common_name":common_name,
								"taxon_id":taxon_id
							}

							proteins_data["organism_host"].append(host)

						proteins_data["taxonomy"] = "|".join(proteins_data["taxonomy"])

						######################################

						if proteins_data["isoform"]:
							self.grabUniProtFasta(proteins_data["accession"])
							fasta_data = self.parseUniProtFasta(proteins_data["accession"])

							if "sequence" in fasta_data:
								proteins_data["sequence_primary_isoform"] = proteins_data["sequence"]
								proteins_data["sequence"] = fasta_data["sequence"]

							else:
								proteins_data["sequence"] = ""
								proteins_data["isofrom_error"] = fasta_data

							if parse_sequence:
								return proteins_data["sequence"]

							isoform_names = []

							for isoform in tree.iter('{http://uniprot.org/uniprot}isoform'):

								if isoform.find('{http://uniprot.org/uniprot}id').text == proteins_data["accession"]:
									# isoform.find('{http://uniprot.org/uniprot}sequence').attrib["type"]
									# isoform.find('{http://uniprot.org/uniprot}sequence').attrib["ref"]

									for name in isoform.iter('{http://uniprot.org/uniprot}name'):
										isoform_names.append(name.text)

									proteins_data["isoform_names"] = isoform_names

									if len(isoform_names) > 1:
										proteins_data["protein_name"] = proteins_data["protein_name"] + " - Isoform " + "/".join(isoform_names[1:])
									elif len(isoform_names) == 1:
										proteins_data["protein_name"] = proteins_data["protein_name"] + " - Isoform " + isoform_names[0]
									else:
										proteins_data["protein_name"] = fasta_data["fullname"]

						if not proteins_data["isoform"]:
							proteins_data["Pfam"] = []
							proteins_data["GO"] = {}
							proteins_data["PDB"] = {}

							for dbReference in tree.iter('{http://uniprot.org/uniprot}dbReference'):
								if dbReference.attrib["type"] == "PDB":
									pdb_method, pdb_resolution, pdb_chains = None, None, []
									for dbReferenceProperty in dbReference.iter('{http://uniprot.org/uniprot}property'):

										if dbReferenceProperty.attrib["type"] == "chains":
											proteins_data["PDB"][dbReference.attrib["id"]] = {}

											for chain in dbReferenceProperty.attrib['value'].split(","):
												try:
													chainBits = chain.split("=")
													chainIds = chainBits[0].strip()

													for chainId in chainIds.split("/"):
														start = int(chainBits[1].split("-")[0])
														end = int(chainBits[1].split("-")[1].strip(" ."))
														proteins_data["PDB"][dbReference.attrib["id"]][chainId] = {"start":start,"end":end}
														pdb_chains.append(chainId)
												except:
													pass
										if dbReferenceProperty.attrib["type"] == "method":
											pdb_method = dbReferenceProperty.attrib['value']
										if dbReferenceProperty.attrib["type"] == "resolution":
											pdb_resolution = dbReferenceProperty.attrib['value']

									for chainId in pdb_chains:
										if pdb_method != None:
											proteins_data["PDB"][dbReference.attrib["id"]][chainId]['method'] = pdb_method
										if pdb_resolution != None:
											proteins_data["PDB"][dbReference.attrib["id"]][chainId]['resolution'] = pdb_resolution

								elif dbReference.attrib["type"] == "Pfam":
									proteins_data["Pfam"].append(dbReference.attrib["id"])

								elif dbReference.attrib["type"] == "GeneTree":
									proteins_data["GeneTree"] = dbReference.attrib["id"]
								elif dbReference.attrib["type"] == "GO":
									tmp = {"id":dbReference.attrib["id"]}
									go_class = ""
									for dbReferenceProperty in dbReference.iter('{http://uniprot.org/uniprot}property'):
										if dbReferenceProperty.attrib["type"] == "term":
											go_class =  dbReferenceProperty.attrib["value"].split(':')[0]
											go_term_name =  dbReferenceProperty.attrib["value"].split(':')[1]

											tmp['term'] = go_term_name
											tmp['class'] = go_class
										else:
											tmp[dbReferenceProperty.attrib["type"]] = dbReferenceProperty.attrib["value"]

									if go_class not in proteins_data["GO"]:
										proteins_data["GO"][go_class] = []

									proteins_data["GO"][go_class].append(tmp)
								else:
									if parse_db_reference:
										if 'db_xref' not in proteins_data:
											proteins_data['db_xref'] = {}

										if dbReference.attrib["type"] not in proteins_data['db_xref']:
											proteins_data['db_xref'][dbReference.attrib["type"]] = {}
											
										dbReferenceData = {}
										for dbReferenceProperty in dbReference.iter('{http://uniprot.org/uniprot}property'):
											dbReferenceData[dbReferenceProperty.attrib["type"]] = dbReferenceProperty.attrib['value']

										for dbReferenceProperty in dbReference.iter('{http://uniprot.org/uniprot}molecule'):
											dbReferenceData['molecule'] = dbReferenceProperty.attrib['id']

										proteins_data['db_xref'][dbReference.attrib["type"]][dbReference.attrib["id"]] = dbReferenceData

						#####
						#####
						#####

					if parse_isoforms or proteins_data["isoform"]:
						isoform_pattern = re.compile("[iI]soform [0-9]+")
						proteins_data["isoforms"] = {
						"alternative_exons":[],
						"isoform_details":{},
						"isoform_comments":{},
						"isoform_name_map":{}
						}

						for isoform in tree.iter('{http://uniprot.org/uniprot}isoform'):
							isoform_id = isoform.find('{http://uniprot.org/uniprot}id').text

							proteins_data["isoforms"]["isoform_details"][isoform_id] = {
							"id":isoform_id,
							"names":[],
							"sequence_type":"",
							"sequence_ref":"",
							"description":"",
							"comments":{}
							}

							for name in isoform.iter('{http://uniprot.org/uniprot}name'):

								proteins_data["isoforms"]["isoform_details"][isoform_id]["names"].append(name.text)
								proteins_data["isoforms"]['isoform_name_map'][name.text] = isoform_id

								isoform_sequence = isoform.find('{http://uniprot.org/uniprot}sequence')
								if isoform_sequence != None:
									if 'type' in isoform_sequence.attrib: proteins_data["isoforms"]["isoform_details"][isoform_id]["sequence_type"] = isoform_sequence.attrib['type']
									if 'ref' in isoform_sequence.attrib: proteins_data["isoforms"]["isoform_details"][isoform_id]["sequence_ref"] = isoform_sequence.attrib['ref']

								isoform_text = isoform.find('{http://uniprot.org/uniprot}text')

								if isoform_text != None:
									proteins_data["isoforms"]["isoform_details"][isoform_id]["description"] = isoform_text.text

						for variant in tree.iter("{http://uniprot.org/uniprot}feature"):
							if variant.attrib['type'] == 'splice variant':
								description = variant.attrib['description']

								location_begin_tag = variant.find('{http://uniprot.org/uniprot}location/{http://uniprot.org/uniprot}begin')
								location_end_tag = variant.find('{http://uniprot.org/uniprot}location/{http://uniprot.org/uniprot}end')
								location_position_tag = variant.find('{http://uniprot.org/uniprot}location/{http://uniprot.org/uniprot}position')
								original_tag = variant.find('{http://uniprot.org/uniprot}original')
								variation_tag = variant.find('{http://uniprot.org/uniprot}variation')

								start = None
								end = None
								original = None
								variation = None

								if location_begin_tag != None: start = location_begin_tag.attrib['position']
								if location_end_tag != None: end = location_end_tag.attrib['position']
								if location_position_tag != None: position = location_position_tag.attrib['position']
								if original_tag != None: original = original_tag.text
								if variation_tag != None: variation = variation_tag.text

								if location_position_tag != None:
									start = location_position_tag.attrib['position']
									end = location_position_tag.attrib['position']

								isoform_ids = []
								for isoform_match in isoform_pattern.findall(description):
									if isoform_match.split()[1] in proteins_data["isoforms"]['isoform_name_map']:
										isoform_id = proteins_data["isoforms"]['isoform_name_map'][isoform_match.split()[1]]
										isoform_ids.append(isoform_id)

								alternative_exon = {
								"id": variant.attrib['id'],
								"start": start ,
								"end": end,
								"original": original,
								"variation": variation,
								"removed": original == None and variation == None,
								"description":description,
								"isoform_ids":isoform_ids
								}

								proteins_data["isoforms"]["alternative_exons"].append(alternative_exon)

								for isoform_id in isoform_ids:
									if isoform_id in proteins_data["isoforms"]["isoform_details"]:
										if 'variation' not in proteins_data["isoforms"]["isoform_details"][isoform_id ]:
											proteins_data["isoforms"]["isoform_details"][isoform_id]['variation' ] = []

										proteins_data["isoforms"]["isoform_details"][isoform_id]['variation'].append(alternative_exon)


						for comment in tree.iter("{http://uniprot.org/uniprot}comment"):
							molecule_tag = comment.find('{http://uniprot.org/uniprot}molecule')
							if molecule_tag != None:
								isoform_id = ""
								for isoform_match in isoform_pattern.findall(molecule_tag.text):
									if isoform_match.split()[1] in proteins_data["isoforms"]['isoform_name_map']:
										isoform_id = proteins_data["isoforms"]['isoform_name_map'][isoform_match.split()[1]]

								text_tag = comment.find('{http://uniprot.org/uniprot}text')

								if text_tag != None:
									if isoform_id in proteins_data["isoforms"]["isoform_details"]:
										if comment.attrib["type"] not in proteins_data["isoforms"]["isoform_details"][isoform_id]['comments']:
											 proteins_data["isoforms"]["isoform_details"][isoform_id]['comments'][comment.attrib["type"]] = []

										proteins_data["isoforms"]["isoform_details"][isoform_id]['comments'][comment.attrib["type"]].append(text_tag.text)

						if proteins_data["isoform"] and proteins_data["accession"] in proteins_data["isoforms"]["isoform_details"]:
							proteins_data["isoform_details"] = proteins_data["isoforms"]["isoform_details"][proteins_data["accession"]]
							if not parse_isoforms: del proteins_data["isoforms"]
						else:
							proteins_data["isoform_details"] = {}

					#####
					#####
					#####

					if parse_keywords:
						proteins_data["keywords"] = {}
						proteins_data["disease_keywords"] = {}
						proteins_data["diseases"] = []
						proteins_data["localisation"] = []
						proteins_data["localisation_isoform"] = {}

						for keyword in tree.iter('{http://uniprot.org/uniprot}keyword'):
							proteins_data["keywords"][keyword.text] = keyword.attrib['id']
							if keyword.attrib['id'] in self.settings["disease_keywords"]:
								proteins_data["disease_keywords"][keyword.text] = keyword.attrib['id']

						#for subcellularLocation in tree.iter('{http://uniprot.org/uniprot}subcellularLocation'):


						for comment in tree.iter("{http://uniprot.org/uniprot}comment"):
							if comment.attrib["type"] == "disease":
								for disease in comment.findall('{http://uniprot.org/uniprot}disease'):
									disease_id = disease.attrib["id"]
									disease_name = disease.find('{http://uniprot.org/uniprot}name').text
									disease_acronym = disease.find('{http://uniprot.org/uniprot}acronym').text
									disease_description = disease.find('{http://uniprot.org/uniprot}description').text
									disease_db_reference = {"source":disease.find('{http://uniprot.org/uniprot}dbReference').attrib["type"],"id":disease.find('{http://uniprot.org/uniprot}dbReference').attrib["id"]}

									proteins_data["diseases"].append({
										"disease_id":disease_id,
										"disease_name":disease_name,
										"disease_acronym":disease_acronym,
										"disease_description":disease_description,
										"disease_db_reference":disease_db_reference
									})

							if comment.attrib["type"] == "subcellular location":
								for subcellularLocation in comment.findall('{http://uniprot.org/uniprot}subcellularLocation'):
									location_tag = subcellularLocation.find('{http://uniprot.org/uniprot}location')
									topology_tag = subcellularLocation.find('{http://uniprot.org/uniprot}topology')
									orientation_tag = subcellularLocation.find('{http://uniprot.org/uniprot}orientation')
									molecule_tag = comment.find('{http://uniprot.org/uniprot}molecule')

									localisation = []
									location = ""
									topology = ""
									orientation = ""

									if location_tag != None:
										location = location_tag.text
										localisation.append(location)

									if topology_tag != None:
										topology = topology_tag.text
										localisation.append(topology)

									if orientation_tag != None:
										orientation = orientation_tag.text
										localisation.append(orientation)

									if len(localisation) > 1:
										if localisation[0] not in proteins_data["localisation"]:
											proteins_data["localisation"].append(localisation[0])

										if "|".join(localisation) not in proteins_data["localisation"]:
											proteins_data["localisation"].append("|".join(localisation))
									else:
										if "|".join(localisation) not in proteins_data["localisation"]:
											proteins_data["localisation"].append("|".join(localisation))

									if molecule_tag != None:
										molecule = molecule_tag.text
										if molecule not in proteins_data["localisation_isoform"]:
											proteins_data["localisation_isoform"][molecule] = []

										localisation = "|".join(localisation)
										if localisation not in proteins_data["localisation"]:
											proteins_data["localisation_isoform"][molecule].append(localisation)

						if parse_isoforms:
							pass

					#####
					#####
					#####

					if parse_features:
						proteins_data["features"] = {}

						for feature in tree.iter('{http://uniprot.org/uniprot}feature'):

							try:
								if feature.attrib["type"] in self.options['use_features'] or len(self.options['use_features']) == 0:
									pmid = []
									if "evidence" in feature.attrib:
										for evidence in feature.attrib["evidence"].split():
											if evidence in references:
												pmid.append(references[evidence])

										del feature.attrib["evidence"]

									if feature.attrib["type"] not in proteins_data["features"]:
										proteins_data["features"][feature.attrib["type"]] = []

									original = feature.find('{http://uniprot.org/uniprot}original')
									variation = feature.find('{http://uniprot.org/uniprot}variation')
									start = feature.find('{http://uniprot.org/uniprot}location/{http://uniprot.org/uniprot}begin')
									end = feature.find('{http://uniprot.org/uniprot}location/{http://uniprot.org/uniprot}end')

									if original != None:
										feature.attrib['original'] = original.text
									if variation != None:
										feature.attrib['variation'] = variation.text

									if start != None:
										feature.attrib['start'] = start.attrib['position']
									else:
										feature.attrib['start'] = feature.find('{http://uniprot.org/uniprot}location/{http://uniprot.org/uniprot}position').attrib['position']

									if end != None:
										feature.attrib['end'] = end.attrib['position']
									else:
										feature.attrib['end'] = feature.find('{http://uniprot.org/uniprot}location/{http://uniprot.org/uniprot}position').attrib['position']

									feature.attrib['pmid'] = pmid
									proteins_data["features"][feature.attrib["type"]].append(feature.attrib)
							except:
								if self.options['verbose']: print(("Can't parse " + feature.attrib["type"]))
								pass
				#####
				#####
				#####

				if parse_snps:
					proteins_data["mutations"] = self.parseHumsavar(proteins_data["accession"])
					proteins_data['disease_mutations'] = None

					if 'phenotype' in proteins_data["mutations"]:
						if 'Disease' in proteins_data["mutations"]['phenotype']:
							proteins_data['disease_mutations'] = proteins_data["mutations"]['phenotype']['Disease']

				if parse_attributes:
					response = {}# self.parseAttributes(proteins_data["accession"])
					proteins_data['disorder'] = None
					proteins_data['conservation'] = None
					proteins_data['anchor'] = None

					if 'data' in response:
						if 'IUPred' in  response['data']:
							proteins_data['disorder'] = response['data']['IUPred']

						if 'Conservation' in  response['data']:
							proteins_data['conservation'] = response['data']['Conservation']

						if 'Anchor' in  response['data']:
							proteins_data['anchor'] = response['data']['Anchor']

				if parse_disorder:
					response= self.parseMobiDB(proteins_data["accession"])
					proteins_data['disprot_consensus'] = None
					proteins_data['disprot'] = None

					if 'data' in response:
						if 'consensus' in response['data']:
							proteins_data['disorder_consensus'] = response['data']['consensus']['full']

							if 'disprot' in response['data']:
								proteins_data['disprot']  = response['data']['consensus']['disprot']


				return {"status":"Success","data":proteins_data}
			else:
				return {"status":"Error","error_type":"File not found"}

		except Exception as e:
			logger.error("Error:" + xml_path)
			logger.error(utilities_error.getError())

			raise
			return {"status":"Error","error_type":str(e)}

	##------------------------------------------------------------------##
	##
	##------------------------------------------------------------------##

	def getUniProtRest(self,rest_options):

		data = {}
		status = "Success"
		error = {"status":"Error","error_type":"undefined"}

		option_list = [
			"parse_uniprot",#parseUniProt
			"parse_mobidb",#parseMobiDB
			"parse_sequence",#parseSequence
			"parse_uniprot_pfam",#parseUniProtPfam,
			"parse_attributes",#parseAttributes,
			"parse_features",#parseFeatures
			"parse_basic",#parseBasic
			"parse_domains",#parseDomains
			"parse_go",#parseGO
			"parse_isoforms",#parseIsoforms
			"parse_keywords",#parseKeywords
			"parse_localisation",#parseLocalisation
			"parse_pdb",#parsePDBs
			"parse_secondary_accessions",#parseSecondaryAccessions
			"parse_snps"#parseSNPs
		]

		required = ['accession']

		defaults = {
		}

		defaults.update(rest_options)
		rest_options = defaults

		try:
			###----###
			###----###

			if rest_options["task"] == "help":
				data = option_list

			if "task" not in rest_options:
				status = "Error"
				error["error_type"] = "No task selected: " + ",".join(option_list)
			else:
				if rest_options["task"] not in option_list:
					status = "Error"
					error["error_type"] = "Task not in task list: " + ",".join(option_list)

			for rest_option in required:
				if rest_option not in rest_options:
					status = "Error"
					error["error_type"] = "No " + rest_option + " supplied"

			###----###
			###----###

			if status != "Error":
				if rest_options["task"] == "parse_uniprot":
					data = self.parseUniProt(rest_options['accession'])
				if rest_options["task"] == "parse_mobidb":
					data = self.parseMobiDB(rest_options['accession'])
				if rest_options["task"] == "parse_sequence":
					data = self.parseSequence(rest_options['accession'])
				if rest_options["task"] == "parse_uniprot_pfam":
					data = self.parseUniProtPfam(rest_options['accession'])
				if rest_options["task"] == "parse_attributes":
					data = self.parseAttributes(rest_options['accession'])
				if rest_options["task"] == "parse_features":
					data = self.parseFeatures(rest_options['accession'])
				if rest_options["task"] == "parse_basic":
					data = self.parseBasic(rest_options['accession'])
				if rest_options["task"] == "parse_domains":
					data = self.parseDomains(rest_options['accession'])
				if rest_options["task"] == "parse_go":
					data = self.parseGO(rest_options['accession'])
				if rest_options["task"] == "parse_isoforms":
					data = self.parseIsoforms(rest_options['accession'])
				if rest_options["task"] == "parse_keywords":
					data = self.parseKeywords(rest_options['accession'])
				if rest_options["task"] == "parse_localisation":
					data = self.parseLocalisation(rest_options['accession'])
				if rest_options["task"] == "parse_pdb":
					data = self.parsePDBs(rest_options['accession'])
				if rest_options["task"] == "parse_secondary_accessions":
					data = self.parseSecondaryAccessions(rest_options['accession'])
				if rest_options["task"] == "parse_snps":
					data = self.parseSNPs(rest_options['accession'])

			###----###
			###----###

			if 'data' in data:
				data = data['data']

			if 'status' in data:
				if data['status'] == "Error":
					status = "Error"
					error["error_type"] = "?"

			if status == "Error":
				return error
			else:
				return {
				"status": status,
				"options":rest_options,
				"data":data
				}
		except Exception as e:
			raise
			return {
			"status":"Error",
			"error_type":str(e)
			}

