import sys,copy,pprint

###################################################################
### Errors
###################################################################

basepair = {"A":"T","T":"A","G":"C","C":"G","N":"N","X":"X","-":"-"}

def check_oligo_sequence_features(oligo,avoided_sequences={},gc_content_cutoff={"min":0,"max":1,"strict_min":0,"strict_max":1},max_complementarity_oligo=10):

	oligo_GC_content =  calculate_GC_content(oligo)

	oligo_GC_content_check = oligo_GC_content >= gc_content_cutoff['min'] and oligo_GC_content <= gc_content_cutoff['max']	
	oligo_GC_content_strict_check = oligo_GC_content >= gc_content_cutoff['strict_min'] and oligo_GC_content <= gc_content_cutoff['strict_max']	
	
	oligo_details = copy.deepcopy(avoided_sequences)
	oligo_details['oligo'] = oligo 
	oligo_details['GC_content'] = oligo_GC_content 
	oligo_details['GC_content_check'] = oligo_GC_content_check 
	oligo_details['GC_content_strict_check'] = oligo_GC_content_strict_check
	oligo_details['max_complementarity_library_check'] = oligo_GC_content_check 
	oligo_details['max_complementarity_oligo_check'] = oligo_GC_content_check 
	
	max_complementarity_library_check = True #check_max_complementarity_library(oligo)
	max_complementarity_oligo_check = check_max_complementarity_oligo(oligo,max_complementarity_oligo)

	oligo_details['check_status'] = {"all":[],"strict":[],"soft":[],"failed_check_type":[]}
	
	oligo_details['all'] = {
		"check":{
			"strict":[
				oligo_GC_content_strict_check
			],
			"soft":[
				oligo_GC_content_check,
				max_complementarity_library_check,
				max_complementarity_oligo_check
			]
		}
	}

	for avoided_oligo_sequence in avoided_sequences:
		oligo_details[avoided_oligo_sequence]['check'] = oligo.count(avoided_oligo_sequence) <= oligo_details[avoided_oligo_sequence]['count']

		if oligo_details[avoided_oligo_sequence]['check'] == False:
			oligo_details['check_status']["failed_check_type"].append(oligo_details[avoided_oligo_sequence]['name'])

		oligo_details[avoided_oligo_sequence]['matches'] = oligo.count(avoided_oligo_sequence)
		
		if oligo_details[avoided_oligo_sequence]['strict']:
			oligo_details['all']['check']["strict"].append(oligo.count(avoided_oligo_sequence) <= oligo_details[avoided_oligo_sequence]['count'])
		else:
			oligo_details['all']['check']["soft"].append(oligo.count(avoided_oligo_sequence) <= oligo_details[avoided_oligo_sequence]['count'])

	oligo_details['check_status']['soft'] = oligo_details['all']['check']["soft"].count(True) == len(oligo_details['all']['check']["soft"])
	oligo_details['check_status']['strict'] = oligo_details['all']['check']["strict"].count(True) == len(oligo_details['all']['check']["strict"])
	oligo_details['check_status']['all'] = oligo_details['check_status']['strict'] and oligo_details['check_status']['soft']
			
	return oligo_details

def complementry(oligo):
	complementry_oligo = ""
	for nucleotide in oligo:
		complementry_oligo += basepair[nucleotide]
	
	return complementry_oligo

def calculate_GC_content(oligo):
	return float(oligo.count('G') +  oligo.count('C'))/len(oligo)

def check_max_complementarity_oligo(oligo,max_complementarity_oligo):

	# PMID:15335214 -- longer loops are unstable. 4-5 bases are optimal, what happens if you limit the distance to 2-8 or 2-10 bases?

	complementry_nucleotides = complementry(oligo)
	complementry_nucleotides = complementry_nucleotides[::-1]

	complementarity_oligo_check = True

	for i in range(0,len(oligo)- max_complementarity_oligo + 1):
		tmp_seq = oligo[i:i+max_complementarity_oligo]

		if complementry_nucleotides.count(tmp_seq) > 0:
			tmp_seq_complementry = complementry(tmp_seq)[::-1]

			split_complementry_nucleotides = complementry_nucleotides.split(tmp_seq_complementry) 
			if len(split_complementry_nucleotides) == 2 and len(complementry_nucleotides.split(tmp_seq)) == 2:

				split_bits = split_complementry_nucleotides[0].split(tmp_seq) + split_complementry_nucleotides[1].split(tmp_seq)
		
				if len(split_bits) > 2:
					if len(split_bits[1]) > 2 and len(split_bits[1]) <= 10:
						#print("Warning Hairpin - Spaced Oligo Complementary - distance:" + str(len(split_bits[1])) + " " + oligo + " " + tmp_seq + " " + tmp_seq_complementry)
						complementarity_oligo_check = False
						#logger.debug("Warning Oligo Overlapping Complementary: " + oligo + " " + tmp_seq + " " + str(len(split_bits[1])))
					
	return complementarity_oligo_check


if __name__ == "__main__":
	avoided_sequences = {
		"AAAAA":{"count":0,"strict":False,"name":"A low complexity","final_strict":False},
		"GGGGG":{"count":0,"strict":False,"name":"G low complexity","final_strict":False},
		"CCCCC":{"count":0,"strict":False,"name":"C low complexity","final_strict":False},
		"TTTTT":{"count":0,"strict":True,"name":"termination signal for RNA pol III","final_strict":False},
		"CGTCTC":{"count":0,"strict":True,"name":"Esp3I/BsmBI site","final_strict":False},
		"GAGACG":{"count":0,"strict":True,"name":"Esp3I/BsmBI site reverse complement","final_strict":False}
	}
	pprint.pprint(check_oligo_sequence_features("TACCAGAAGGTTTTTTCTATGGGAC",avoided_sequences=avoided_sequences))
	