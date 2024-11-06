import math,sys,traceback,pprint,os,hashlib,logging
import zlib,json, base64
import warnings
warnings.filterwarnings("ignore", message="the sets module is deprecated")

###################################################################
### Errors
###################################################################

def writeError(e):
	print(e)
	error_type = sys.exc_info()[0]
	error_value = sys.exc_info()[1]
	error_traceback = traceback.extract_tb(sys.exc_info()[2])
	print(error_traceback)

	print("\n")
	print('Error in routine\n<br>')
	print(('Error Type       : ' + str(error_type)))
	print(('Error Value      : ' + str(error_value)))
	print(('File             : ' + str(error_traceback[-1][0])))
	print(('Method           : ' + str(error_traceback[-1][2])))
	print(('Line             : ' + str(error_traceback[-1][1])))
	print(('Error            : ' + str(error_traceback[-1][3])))



###################################################################
### remove_html_tags
###################################################################
def remove_html_tags(text):
	"""Remove html tags from a string"""
	import re
	clean = re.compile('<.*?>')
	return re.sub(clean, '', text)
	
	
###################################################################
### Checks
###################################################################

def isInt(x):
	try:
		int(x)
		return True
	except:
		return False

###################################################################
### files
###################################################################

def gunzip_shutil(source_filepath, dest_filepath, block_size=65536):
	import gzip
	import shutil
	with gzip.open(source_filepath, 'rb') as s_file, open(dest_filepath, 'wb') as d_file:
		shutil.copyfileobj(s_file, d_file, block_size)

def fileAgeDays(filepath):
	from datetime import datetime
	if os.path.exists(filepath):
		stat = os.stat(filepath)
		fileage = datetime.fromtimestamp(stat.st_mtime)
		now = datetime.now()
		delta = now - fileage

		return delta.days
	else:
		return -1
		
def fileChecker(filepath,fileDesc,exit=False):
	if os.path.exists(filepath):
		return True
	else:
		print(("Error @ file check : " + fileDesc + " - " + filepath + " does not exist"))

		if exit:
			sys.exit()

		return False

def params_to_hash(params):
	param_str = ""
	param_keys = list(params.keys())
	param_keys.sort()

	for param in param_keys:
		param_str += param+"="+str(params[param])+"\n"

	hash = str(hashlib.md5(param_str.encode()).hexdigest())
	return hash

def list_to_hash(list_for_hashing):
	list_for_hashing.sort()
	hash = str(hashlib.md5(str(list_for_hashing).encode()).hexdigest())
	return hash

import json
import numpy as np

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

def write_to_json(out_file_json,json_data,zipped=False,normalise_json=False):
	if normalise_json:
		json_data = json.loads(json.dumps(json_data, cls=NpEncoder))
	
	if zipped:
		json_data = json_zip(json_data)

	try:
		with open(out_file_json ,"w") as outfile:
			json.dump(json_data, outfile)
		return {"status":"Success","data":out_file_json}
	except Exception as e:
		print(str(e))
		return {"status":"Error","error_type":"File writing failed","error_details":str(e)}

def read_from_json(out_file_json,zipped=False):
	
	if os.path.exists(out_file_json):
		try:
			if zipped:
				return json_unzip(json.loads(open(out_file_json).read()))
			else:
				with open(out_file_json) as data_file:
					return json.load(data_file)
		except:
			logging.error("Error reading: " + out_file_json)
			raise
	else:
		return {"status":"Error","error_type":"File does not exist","error":str(sys.exc_info()[0])}

###################################################################
### JSON ZIP
###################################################################


ZIPJSON_KEY = 'base64(zip(o))'

def json_zip(j):
    j = {
        ZIPJSON_KEY: base64.b64encode(
            zlib.compress(
                json.dumps(j).encode('utf-8')
            )
        ).decode('ascii')
    }

    return j


def json_unzip(j, insist=False):
    try:
        assert (j[ZIPJSON_KEY])
        assert (set(j.keys()) == {ZIPJSON_KEY})
    except:
        if insist:
            raise RuntimeError("JSON not in the expected format {" + str(ZIPJSON_KEY) + ": zipstring}")
        else:
            return j

    try:
        j = zlib.decompress(base64.b64decode(j[ZIPJSON_KEY]))
    except:
        raise RuntimeError("Could not decode/unzip the contents")

    try:
        j = json.loads(j)
    except:
        raise RuntimeError("Could interpret the unzipped contents")

    return j


###################################################################
### Strings
###################################################################
def alignMotif(str1,str2,returnScore=False,similarity=0.5):
	origStr1 = str1
	origStr2 = str2

	strLength = len(str2)
	offset = -len(str2)

	maxOffset = [0,0]

	str1 = "-"*len(str2) + str1  + "-"*len(str2)

	while len(str2) <= len(str1) :

		matches = [str1[i] == str2[i] for i in range(0,len(str2))].count(True)

		#print offset,matches,str2.count("."),strLength - str2.count("."),strLength

		if matches > maxOffset[1]:
			maxOffset = [offset,matches]

		offset += 1
		str2 = "*" + str2

	if float(maxOffset[1])/(strLength - str2.count(".")) < similarity:
		#print float(maxOffset[1])/strLength
		if returnScore:
			return ["False",float(maxOffset[1])/(strLength - str2.count("."))]
		else:
			return "False"
	else:
		"""print
		if maxOffset[0] >= 0:
			print origStr1
			print "-"*maxOffset[0] + origStr2
		else:
			print "-"*abs(maxOffset[0]) + origStr1
			print origStr2#"""

		if returnScore:
			return [maxOffset[0],float(maxOffset[1])/(strLength - str2.count("."))]
		else:
			return maxOffset[0]


def alignStrings(str1,str2,returnScore=False):
	origStr1 = str1
	origStr2 = str2

	strLength = len(str2)
	offset = -len(str2)

	maxOffset = [0,0]

	str1 = "-"*len(str2) + str1 + "-"*len(str2)

	while len(str2) <= len(str1):
		matches = [str1[i] == str2[i] for i in range(0,len(str2))].count(True)

		if matches > maxOffset[1]:
			maxOffset = [offset,matches]

		offset += 1
		str2 = "*" + str2

	if strLength > 0:
		if float(maxOffset[1])/strLength < 0.5:
			if returnScore:
				return ["False",float(maxOffset[1])/strLength]
			else:
				return "False"
		else:
			if returnScore:
				return [maxOffset[0],float(maxOffset[1])/strLength]
			else:
				return maxOffset[0]

###################################################################
###  Lists
###################################################################

def diffLists(list1,list2):
	for item in list1:
		if item in list2:
	  		list2.remove(item)

	return list2

def removeRedundency(list): # Not order preserving
	tmp = {}
	for val in list:
		tmp[val] = 1
	return list(tmp.keys())

def removeRedundencyOrdered(list):
	set = {}
	return [set.setdefault(e,e) for e in list if e not in set]

def binList(list,logValue=10,log=False,normaliser=1):
	binDict = {}

	for val in list:
		val = float(val)
		if log:
			try:
				bin = abs(int(math.log(val,logValue)))
				if bin in binDict:
					binDict[bin] += 1
				else:
					binDict[bin] = 1
			except:
				pass
		else:
			#print val,int(val*normaliser)
			try:
				if int(val*normaliser) in binDict:
					binDict[int(val*normaliser)] += 1
				else:
					binDict[int(val*normaliser)] = 1
			except:
				pass

	return binDict

def plotList(plot,adjuster=100):
	for val in range(0,len(plot)):
		print(((val + 1),"\t","%1.3f"%plot[val],"\t","*"*int(plot[val]*adjuster)))

def inFrame(start,stop,rangeStart,rangeStop):
	start = int(start)
	stop = int(stop)
	rangeStart = int(rangeStart)
	rangeStop = int(rangeStop)

	if start <= rangeStop and start >= rangeStart or stop <= rangeStop and stop >= rangeStart or stop >= rangeStop and start <= rangeStart :
		return True
	else:
		return False

def featureDistance(start,stop,rangeStart,rangeStop):
	start = int(start)
	stop = int(stop)
	rangeStart = int(rangeStart)
	rangeStop = int(rangeStop)

	minDistance = 10000
	for x in [start - rangeStart, stop - rangeStart, start - rangeStop, stop - rangeStop]:
		if abs(x) < abs(minDistance):
			minDistance = x

	return minDistance

def overlap(start,stop,rangeStart,rangeStop):
	return set(list(range(int(start),int(stop)+1))).intersection(set(list(range(int(rangeStart),int(rangeStop)+1))))

def within(start,stop,offset):
	return offset >= start and offset < stop



###################################################################
### Dicts
###################################################################

def printDict(tmpDict):
	pp = pprint.PrettyPrinter(indent=4)
	pp.pprint(tmpDict)


def plotDist(dist):
	sorter = list(dist.keys())
	sorter.sort()
	sorter.reverse()

	samples = sum(dist.values())


	outStr = ""
	summer = 0
	for v in sorter:
		summer += dist[v]
		outStr += str(10**-v) + "\t"
		outStr += str(dist[v]) + "\t"
		outStr += str(summer) + "\t"
		outStr += "%1.3f"%(float(dist[v])/samples) + "\t"
		outStr += "%1.3f"%(float(summer)/samples) +"\t"
		outStr += str(int((float(dist[v])/samples)*100)*"*")
		outStr += "\n"

	return outStr

###################################################################
###  Probability related
###################################################################

def getProb(offsets,scores):
	prob = 1
	for offset in offsets:
		prob *=  rlcProb(scores[offset],0,1)

	return prob

def upd(u,n):
	return (float(pow(-1,n-1))/factorial(n - 1))*pow(math.log(u),(n-1))


def std_dev(list_temp):
	try:
		mean = sum(list_temp)/len(list_temp)

		sd = 0

		for v in list_temp:
			sd += (v**2 - mean**2)

		sd = sd/len(list_temp)
		return math.sqrt(sd)
	except:
		return -10

def gauss(x):
	u = 0
	s = 1
	return (1/(s*math.sqrt(2*math.pi)))*math.exp(-((x - u)**2)/2*(s**2))


def informationContent(column,IC=""):
	AAs = "ACDEFGHIKLMNPQRSTVWY"
	default_equiv = 'AGS,FYW,FYH,ILV,ILMV,ILMVF,KR,KRH,DE,ST'
	default_equiv_set = list(AAs) + default_equiv.split(",")

	pattern = ""

	countAA = {}
	best = ("",0)

	for aa in set(column):
		countAA[aa] = float(column.count(aa))/len(column)

	for equiv in default_equiv_set:
		sum = 0
		for val in list(equiv):
			if val in countAA:
				sum += countAA[val]

		if best[1] < sum:
			best = (equiv, sum)

	if best[1] > 0.90:
		if len(best[0]) > 1:
			pattern = "[" + best[0] + "]"
		else:
			pattern = best[0]
	else:
		pattern = "."

	pattern = pattern.replace("-","")

	return pattern
