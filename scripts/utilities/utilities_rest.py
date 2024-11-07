import urllib.request, urllib.error, urllib.parse, json, time, os, pprint, sys, optparse, requests

sys.path.append(os.path.join(os.path.dirname(__file__), "../"))
import config_reader

#-----
import logging
logger = logging.getLogger(__name__)
#-----

class RestSubmission():
	def __init__(self):
		pass
		self.headers = {}
		self.token = ''
		
		self.options = config_reader.load_configeration_options(sections=["general"])
        
		self.options.update({
			'verbose':False,
			'waiting_reponses':["In progress","Submitted",'Running', "Submitted to Kafka"]
		})


	def login(self):
		logger.debug('logging in - ' + self.login_url)
		login_response = urllib.request.urlopen(self.login_url).read()
		login_response = json.loads(login_response)

		self.headers = login_response["headers"]

	#--------------------------#
	#--------------------------#

	def loginToken(self):

		params = {
			"username": self.username,
			"password": self.password
		}

		return self.restRequestPOST(os.path.join(self.options["server_url"], "restapi/users/login/access-token-and-user"), params)

	#--------------------------#
	#--------------------------#

	def checkToken(self):
		return self.restRequestPOST(os.path.join(self.options["server_url"], "restapi/users/login/check-token/" + self.token),{})

	#--------------------------#
	#--------------------------#

	def whoAmI(self):
		return self.restRequestPOST(os.path.join(self.options["server_url"], "restapi/users/me"),{})

	#--------------------------#
	#--------------------------#

	def waiter(self,job_status):
		isRunning = True
		logger.info("Still running - " + job_status)
		if job_status in self.options['waiting_reponses']:
			time.sleep(2)
		else:
			isRunning = False
		return isRunning

	#--------------------------#
	#--------------------------#

	def restRequestGET(self,urlRequest,verify=True):
		opener = urllib.request.build_opener()

		for header in self.headers:
			opener.addheaders.append(('Cookie', str(header[1])))

		logger.debug("GET:" + urlRequest)
		f = opener.open(urlRequest).read()
		return json.loads(f)

	#--------------------------#
	#--------------------------#

	def restRequestPOST(self,urlRequest,params,verify=True):
		opener = urllib.request.build_opener()

		self.session = requests.Session()

		logger.debug("POST:" + urlRequest)

		if len(self.token) > 0:
			self.headers['Content-Type'] = 'application/json'
			self.headers['Authorization'] = 'Bearer {}'.format(self.token)
			
		resp = self.session.post(urlRequest, data=params,headers=self.headers,verify=verify)
		del self.session

		if not resp.ok:
			return {"status":"Error","error_type":resp.status_code,"response":resp.text}

		return resp.json()

	#--------------------------#
	#--------------------------#

	def run_query(self,  path,params={},submission="POST",verify=True):
		url = urllib.parse.urljoin(self.server_url,path)
	
		logger.debug("Submitting: " + submission + ' ' + url)
		logger.debug(params)
		try:
			if submission == "POST":
				response = self.restRequestPOST(url,params,verify=verify)
			else:
				response = self.restRequestGET(url,verify=verify)

			return response

		except Exception as e:
			return {'status': 'Error',"error_type":str(e)}

	#--------------------------#
	#--------------------------#

if __name__ == "__main__":

	op = optparse.OptionParser()
	options = {
	"populate_from_file":False,
	"populator_file":"",
	"get_reference":False,
	"get_structure":False,
	"get_elm_instance":False,
	"get_elm_class":False,
	"get_protein":False,
	"identifier":"",
	"identifier_type":"",
	"format":"",
	"taxon_id":"9606",
	"flanks":5,
	"production":False
	}

	op.add_option("--identifier",
		action="store",
		dest="identifier",
		default=False,
		help="")
	op.add_option("--identifier_type",
		action="store",
		dest="identifier_type",
		default=False,
		help="")

	op.add_option("--populate_from_file",
		action="store",
		dest="populate_from_file",
		default=False,
		help="Add ids from file.")
	op.add_option("--populator_file",
		action="store",
		dest="populator_file",
		default="",
		help="Add ids from file path.")

	op.add_option("--get_reference",
		action="store",
		dest="get_reference",
		default=False,
		help="PMIDS")
	op.add_option("--get_structure",
		action="store",
		dest="get_structure",
		default=False,
		help="Structures")
	op.add_option("--get_elm_instance",
		action="store",
		dest="get_elm_instance",
		default=False,
		help="ELM")
	op.add_option("--get_elm_class",
		action="store",
		dest="get_elm_class",
		default=False,
		help="ELM")
	op.add_option("--get_protein",
		action="store",
		dest="get_protein",
		default=False,
		help="get_protein")

	op.add_option("--format",
		action="store",
		dest="format",
		default=False,
		help="format")

	op.add_option("--production",
		action="store",
		dest="production",
		default=False,
		help="Populated the production server.")


	opts, args = op.parse_args()

	for option in options:
		try:
			options[option] = getattr(opts, option)
		except:
			pass

	if options["populate_from_file"]:
		if not os.path.exists(options['populator_file']):
			print(("[populator_file] option:", options['populator_file'], "does not exist"))
			sys.exit()

	#--------------------------#

	restSubmitter = RestSubmission()

	#--------------------------#

	if options["production"]:
		restSubmitter.server_url = "http://slim.ucd.ie/webservices"
	else:
		restSubmitter.server_url = "http://localhost:6543"

	restSubmitter.login_url = restSubmitter.server_url + "/login/?login=slim&password=enterslim"
	restSubmitter.login()

	#--------------------------#

	if options["get_reference"]:
		if options["populate_from_file"]:
			pmids = open(options['populator_file']).read().strip().split("\n")
		else:
			pmids = options["identifier"].split(",")

		for pmid in pmids:
			try:
				params = {"pmid":pmid}
				response = restSubmitter.run_query('/reference/get/',params, "POST")
				print(response)
			except:
				raise

	#--------------------------#

	if options["get_structure"]:
		if options["populate_from_file"]:
			pdbids = open(options['populator_file']).read().strip().split("\n")
		else:
			pdbids = options["identifier"].split(",")

		for pdbid in pdbids:
			try:
				params = {"pdbid":pdbid}
				params = {"pdbid":pdbid,"tasks":["get_contacts_detailed"]} #"get_pockets","get_motifs", "get_details", get_proteins, get_contacts
				response = restSubmitter.run_query('/structure/get/',params, "POST")

				pprint.pprint(response)
			except:
				raise

	#--------------------------#

	if options["get_protein"]:
		if options["populate_from_file"]:
			proteinids = open(options['populator_file']).read().strip().split("\n")
		else:
			proteinids = options["identifier"].split(",")

		for proteinid in proteinids:
			try:
				params = {"accession":proteinid,"tasks":options['identifier_type'].split(",")}
				response = restSubmitter.run_query('/protein/get/',params, "POST")
				pprint.pprint(response)
			except:
				raise

	#--------------------------#

	if options["get_elm_class"]:
		if options["populate_from_file"]:
			pmids = open(options['populator_file']).read().strip().split("\n")
		else:
			elmids = options["identifier"].split(",")

		for elmid in elmids:
			try:
				params = {"identifier":elmid,"identifier_type":options['identifier_type'],"tasks":"get_elm_classes"}
				response = restSubmitter.run_query('/elm/get/',params, "POST")

				pprint.pprint(response)
			except:
				raise

	#--------------------------#

	if options["get_elm_instance"]:
		if options["populate_from_file"]:
			pmids = open(options['populator_file']).read().strip().split("\n")
		else:
			elmids = options["identifier"].split(",")

		for elmid in elmids:
			try:
				params = {"identifier":elmid,"identifier_type":options['identifier_type'],"tasks":"get_elm_instances","flanks":options['flanks'],"format":options['format'],"aligned_type":"peptides_padded_pssm_aligned","max_gap_proportion":0.25}
				response = restSubmitter.run_query('/elm/get/',params, "POST")

				pprint.pprint(response)
			except:
				raise
