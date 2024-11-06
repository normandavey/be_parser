import sys, os, inspect, time, requests, string

###################################################################
### download file class with sessions
###################################################################

import utilities_error
import utilities_basic

#-----
import logging
logger = logging.getLogger(__name__)
#-----

class sessionDownloader():
	#=============================================================================================

	def __init__(self):
		self.session = requests.Session()

	def download_file(self,url,out_path,method="POST",params={},json_params={},data={},force=False,JSON=False,zipped_JSON=False,remake_age=365,replace_empty=True,wait_time=0.5,zipped=False,save_http_errors=False,headers=None,print_content=False,verify=True):
		if os.path.exists(out_path):
			if force:
				logger.debug("Forced - Deleting: " + out_path)
				os.remove(out_path)

			elif (time.time() - os.path.getmtime(out_path))/60/60/24 > remake_age:
				logger.debug("Outdated - Deleting" + out_path + " - age:" + "%1.2f"%((time.time() - os.path.getmtime(out_path))/60/60/24))
				os.remove(out_path)

			elif (replace_empty and os.path.getsize(out_path) == 0):
				logger.debug("Empty - Deleting: " + out_path)
				os.remove(out_path)
			else:
				logger.debug("File exists: " + out_path)
				return {"status":"Success","status_type":"File exists: " + out_path}

		if not os.path.exists(out_path):
			try:
				if zipped:
					unzipped_out_path = out_path
					out_path = out_path + ".gz"

					logger.debug("Downloading zipped file to" + out_path)
					
				if method == "FTP":
					import ftplib

					download_server = url.split("/")[2]
					download_path = "/".join(url.split("/")[3:-1])
					download_file = url.split("/")[-1]
					ftp_usr = ""
					ftp_pass = ""

					logger.debug("FTP: " + download_server)
					logger.debug("FTP path: " + download_path)
					logger.debug("FTP file: " + download_file)

					ftp = ftplib.FTP(download_server)

					ftp.login(ftp_usr, ftp_pass)
					ftp.cwd(download_path)
					ftp.retrbinary("RETR " + download_file ,open(out_path, 'wb').write)
					ftp.quit()

					logger.debug("FTP loaded")

					return {"status":"Success"}

				if len(params) != 0:
					if headers == None:
						resp = self.session.post(url, params=params)
					else:
						resp = self.session.post(url, params=params,headers=headers) 
				elif len(json_params) != 0:
					if headers == None:
						resp = self.session.post(url, json=json_params)
					else:
						resp = self.session.post(url, json=json_params,headers=headers) 
			
				elif len(data) != 0:
					if headers == None:
						resp = self.session.post(url, data=data)
					else:
						resp = self.session.post(url, data=data,headers=headers) 
				
				else:
					if method == "POST":
						logger.debug("POST: " + url)
						if headers == None:
							resp = self.session.post(url,verify=verify)
						else:
							resp = self.session.post(url,headers=headers,verify=verify) 
					else:
						logger.debug("GET: " + url)
						if headers == None:
							resp = self.session.get(url,verify=verify)
						else:
							resp = self.session.get(url,headers=headers,verify=verify) 
				
				
				status = resp.status_code

				if status not in [200,204]:
					logger.debug("GET: " + url + ' - ' + str(resp.status_code))
					import pprint
					if save_http_errors == False:
						try:
							status_text = resp.text
						except:
							status_text = ""

						logger.debug("GET: " + url + ' - ' + str(status_text))
						return {
							"status":"Error",
							"status_code":resp.status_code,
							"status_text":status_text,
							"details":str(resp),
							"error_type":"Error downloading file: " + os.path.basename(out_path) + " from " + url
							}

				if print_content:
					print(resp.text)
					
				if JSON:
					import json
					content = json.loads(resp.text)
					
					if 'data' in content:
						content = content['data']
					
					status = utilities_basic.write_to_json(out_path,content,zipped=zipped_JSON,normalise_json=False)
					if status['status'] == "Error":logger.error(status)
				else:
					content = resp.text

					logger.debug("Writing: " + out_path)
					try:
						open(out_path,"w").write(content)
					except:
						content = "".join([x for x in content if x in string.printable])
						open(out_path,"w").write(content)
				
				if zipped:
					logger.debug("Unzipping zipped file to" + unzipped_out_path)
					import gzip
					import shutil
					with gzip.open(out_path, 'rb') as f_in:
						with open(unzipped_out_path, 'wb') as f_out:
							shutil.copyfileobj(f_in, f_out)

				time.sleep(wait_time)

				return {"status":"Success","content":content}
			except Exception:
				return {"status":"Error","error_type":utilities_error.getError()}

if __name__ == "__main__":
	test = "test2"

	if test == "test1":
		pdb_id = "2AST"
		pdb_file_path = os.path.dirname(inspect.stack()[0][1])

		url = "http://files.rcsb.org/download/" + pdb_id + ".pdb"
		out_path = os.path.abspath(os.path.join(pdb_file_path, pdb_id + ".pdb"))

		sessionDownloaderObj = sessionDownloader()
		response = sessionDownloaderObj.download_file(url,out_path)

		print((open(out_path).read()))
		print(response)

	if test == "test2":
		url = "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt"
		url = "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/README"
		out_path = "./README.txt"
		sessionDownloaderObj = sessionDownloader()
		response = sessionDownloaderObj.download_file(url,out_path,method="FTP")

		print(response)
