import sys, traceback, pprint, re, os

###################################################################
###    readAlignment readAlignment readAlignment readAlignment  ### 
###################################################################

#-----
import logging
logger = logging.getLogger(__name__)
#-----

sys.path.append(os.path.join(os.path.dirname(__file__), "../data_management"))
import queryRunner

###################################################################

def read_blosum_matrix(matrix_name="blosum62"):
	return read_blast_similarity_matrix(os.path.join(os.path.dirname(__file__), "../matrices/" + matrix_name + ".bla"))

###################################################################

def read_similarity_matrix(similarity_matrix_path):
	similarity_matrix = {}

	file_content = open(similarity_matrix_path).read().split("\n")

	header = file_content[0].strip().split('\t')
	idx_to_aa = dict(list(zip(list(range(0,len(header))), header)))

	max_score = 0
	for line in file_content[1:]:
		fields = line.strip().split('\t')
		from_aa = fields[0]
		similarity_matrix[from_aa] = {}

		for idx, score in enumerate(fields):
			if idx == 0:
				continue

			to_aa = idx_to_aa[idx]
			similarity_matrix[from_aa][to_aa] = float(score)

			if to_aa not in similarity_matrix:
				similarity_matrix[to_aa] = {}

			similarity_matrix[to_aa][from_aa] = float(score)

			if float(score) > max_score:
				max_score = float(score)

		similarity_matrix[from_aa]["-"] = -4

	similarity_matrix["-"] = {}

	for aa in similarity_matrix:
		if aa == "-":
			similarity_matrix["-"][aa] = -1
		else:
			similarity_matrix["-"][aa] = -4
	
	return similarity_matrix

##################################################################################################

def read_blast_similarity_matrix(similarity_matrix_path,delimiter=" "):
	similarity_matrix = {}

	file_content = open(similarity_matrix_path).read().strip().split("\n")

	header = None
	idx_to_aa = None
	max_score = 0
	row_counter = 0
	for line in file_content:
		
		if line[0] == "#": continue

		if header == None:
			if delimiter == " ":
				header = line.strip().split()
			else:
				header = line.strip().split(delimiter)

			idx_to_aa = dict(list(zip(list(range(0,len(header))), header)))
		else:
			if delimiter == " ":
				fields = line.strip().split()
			else:
				fields = line.strip().split(delimiter)


			from_aa = idx_to_aa[row_counter]
			row_counter += 1
			similarity_matrix[from_aa] = {}

			for idx, score in enumerate(fields):
				to_aa = idx_to_aa[idx]
				similarity_matrix[from_aa][to_aa] = float(score)

				if to_aa not in similarity_matrix:
					similarity_matrix[to_aa] = {}

				similarity_matrix[to_aa][from_aa] = float(score)

				if float(score) > max_score:
					max_score = float(score)

			similarity_matrix[from_aa]["-"] = -4

	similarity_matrix["-"] = {}

	for aa in similarity_matrix:
		if aa == "-":
			similarity_matrix["-"][aa] = -1
		else:
			similarity_matrix["-"][aa] = -4
	
	return similarity_matrix

##################################################################################################

def sequenceIdentity(seqA,seqB):
	return "".join(["+" if (seqA[x] == seqB[x] and seqA[x] != "-") else " " for x in range(0,len(seqA))])

def sequenceSimilarity(seqA,seqB,similarity_matrix_path=None):
	seqA = seqA.replace('U',"X")
	seqB = seqB.replace('U',"X")
	
	similarity_matrix = read_blast_similarity_matrix(os.path.join(os.path.dirname(__file__), "../matrices/blosum62.bla")) #read_similarity_matrix(similarity_matrix_path)

	identity = [True if (seqA[x] == seqB[x] and seqA[x] != "-") else False for x in range(0,len(seqA))]
	similarity = [True if (similarity_matrix[seqA[x]][seqB[x]] > 0) else False for x in range(0,len(seqA))]

	similarity_string  = ""
	for i in range(0,len(seqA)):
		if identity[i]:
			similarity_string += "+"
		elif similarity[i]:
			similarity_string += "."
		else:
			similarity_string += " "

	return similarity_string


def sequenceSimilarityScores(seqA,seqB):
	similarity_string = sequenceSimilarity(seqA,seqB)
	return {
		"identity":similarity_string.count("+")/len(similarity_string),
		"similarity":(similarity_string.count("+")+similarity_string.count("."))/len(similarity_string)
	}

##################################################################################################

def sequencesToFasta(sequences):
	fasta_data = ""
	for sequence in sequences:
		if "start" in sequence and "end" in sequence:
			fasta_data += ">" + "\n".join([sequence['identifier']+"/"+str(sequence['start'])+"-"+str(sequence['end']),sequence['sequence']]) + "\n"
		else:
			fasta_data += ">" + "\n".join([sequence['identifier'],sequence['sequence']]) + "\n"

	return fasta_data

def alignSequences(alignment_input_path, alignment_output_path, alignment_algorithm,alignment_algorithm_path,remake=False,installed=False):
	try:
		if not os.path.exists(alignment_algorithm_path) or remake:
			if not os.path.exists(alignment_algorithm_path) and installed == False:
				return  {"status":"Error","type":"Alignment algorithm does not appear to exist"}

			if alignment_algorithm == "muscle":
				cmd = alignment_algorithm_path + " -in " + alignment_input_path+ " -out " + alignment_output_path  + ""

			if alignment_algorithm == "clustalo":
				cmd = alignment_algorithm_path + " --force --infile=" + alignment_input_path + " --outfile=" + alignment_output_path  + ""

			if alignment_algorithm == "clustalw":
				cmd = alignment_algorithm_path + " -ALIGN -INFILE=" + alignment_input_path + " -OUTFILE=" + alignment_output_path  + " -OUTPUT=FASTA -OUTORDER=INPUT"

			if alignment_algorithm == "fsa":
				cmd = alignment_algorithm_path + " " + alignment_input_path + " > " + alignment_output_path + ""

			if alignment_algorithm == "mafft":
				#os.system('MAFFT_BINARIES="/Users/normandavey/Documents/Work/Tools/mafft/libexec"; export MAFFT_BINARIES;')
				cmd = alignment_algorithm_path + " --quiet --anysymbol --auto " + alignment_input_path + " > " + alignment_output_path  + ""

			if alignment_algorithm == "probcons":
				cmd = alignment_algorithm_path +  " " + alignment_input_path + " > " + alignment_output_path  + ""

			if len(cmd) > 0:
				logger.debug(cmd)
				os.popen(cmd).read()
			else:
				print("Can't run",alignment_input_path)

			return {"status":"Completed","data":alignment_output_path}
		else:
			return {"status":"Exists"}
	except:
		return {"status":"Error","error_type":utilities_error.printError()}

##################################################################################################


def printAlignment(alignment_path):
	logger.debug(alignment_path)
	print(open(alignment_path).read().strip())

def readAlignment(alignment_path):
	logger.debug(alignment_path)
	alignment_content = open(alignment_path).read().strip()
	alignment_format = checkAlignmentFormat(alignment_content)

	if alignment_format == "FASTA":
		return  readFastaAlignment(alignment_content)
	elif alignment_format == "STOCKHOLM":
		return readStockholmAlignment(alignment_content)
	else:
		return {}

def readStockholmAlignment(alignment_content):
	alignmentBits = alignment_content.split("\n")

	alignment_data = {}

	for alignmentBit in alignmentBits:
		try:
			if len(alignmentBit.strip()) == 0: continue
			if alignmentBit[0] != "#" and alignmentBit[0] != "/":
				proteinBits = alignmentBit.split()
				name = proteinBits[0]
				sequence = proteinBits[1]
				if name not in alignment_data:
					alignment_data[name] = []

				alignment_data[name].append(sequence)
		except:
			raise
	
	for name in alignment_data:		
		alignment_data[name] = "".join(alignment_data[name]).replace(".","-").upper()
		
	return alignment_data

def readFastaAlignment(alignment_content):
	alignmentBits = alignment_content[1:].split("\n>")

	alignment_data = {}

	for alignmentBit in alignmentBits:
		try:
			proteinBits = alignmentBit.split("\n")
			name = proteinBits[0]
			sequence = "".join(proteinBits[1:])
			alignment_data[name] = sequence
		except:
			pass
	
	return alignment_data

##################################################################################################

def checkAlignmentFormat(alignment_content):
	try:
		if alignment_content[0] == ">":
			return "FASTA"
		elif alignment_content[0] == "#":
			return "STOCKHOLM"
		else:
			return "UNKNOWN"
	except:
		logger.error("Error reading: " + alignment_content)

def mapAlignmentColumnsToOffets(alignment_path):
	logger.debug(alignment_path)
	
	alignment_data_mapping = {
		"mapping":{}
	}

	alignmentData = readAlignment(alignment_path)

	for alignmentSequenceIdentifier in alignmentData:
		try:
			sequence = alignmentData[alignmentSequenceIdentifier]
			mapping = {}
			query_protein_offet = 0
			for pos in range(0,len(sequence)):
				if sequence[pos] == "-":
					pass
				else:
					mapping[query_protein_offet] = pos
					query_protein_offet += 1

			alignment_data_mapping["mapping"][alignmentSequenceIdentifier] = mapping
		except:
			raise
	
	return alignment_data_mapping

##################################################################################################

def parseAlignmentColumnsByOffset(alignment_path):
	logger.debug(alignment_path)

	alignment_data_mapping = {
		"column":{},
		"order":[]
	}
	
	alignmentData = readAlignment(alignment_path)

	for alignmentSequenceIdentifier in alignmentData:
		try:
			sequence = alignmentData[alignmentSequenceIdentifier]

			alignment_data_mapping['order'].append(alignmentSequenceIdentifier)
			
			for pos in range(0,len(sequence)):
				if pos not in alignment_data_mapping["column"]:
					alignment_data_mapping["column"][pos] = []
				alignment_data_mapping["column"][pos].append(sequence[pos])

		except:
			raise

	return alignment_data_mapping

##################################################################################################


def parseAlignmentByProtein(alignment_path):
	logger.debug(alignment_path)
	
	alignment_data = {
	"proteins":{},
	"order":[]
	}

	alignmentData = readAlignment(alignment_path)

	for alignmentSequenceIdentifier in alignmentData:
		try:
			sequence = alignmentData[alignmentSequenceIdentifier].upper().replace(".","-")
			
			uniprot_accession = find_uniprot_accession(alignmentSequenceIdentifier)
			taxon_id = find_taxon_id(alignmentSequenceIdentifier)
			species = find_species(alignmentSequenceIdentifier)
			genename = find_genename(alignmentSequenceIdentifier)

			try:
				
				region_boundaries = find_region_boundaries(alignmentSequenceIdentifier)
				region_boundary_start = region_boundaries[0]
				region_boundary_end = region_boundaries[1]

				alignment_data['order'].append(uniprot_accession)
				if uniprot_accession not in alignment_data["proteins"]:
					alignment_data["proteins"][uniprot_accession] = []

				alignment_data["proteins"][uniprot_accession].append({
				"header":alignmentSequenceIdentifier,
				"uniprot_accession":uniprot_accession,
				"taxon_id":taxon_id,
				"species":species,
				"genename":genename,
				"region_boundary_start":region_boundary_start,
				"region_boundary_end":region_boundary_end,
				"sequence":sequence,
				"sequence_length":len(sequence.replace("-","").replace(".",""))
				})
			except:
				alignment_data['order'].append(uniprot_accession)
				if uniprot_accession not in alignment_data["proteins"]:
					alignment_data["proteins"][uniprot_accession] = []

				alignment_data["proteins"][uniprot_accession].append({
				"header":alignmentSequenceIdentifier,
				"uniprot_accession":uniprot_accession,
				"taxon_id":taxon_id,
				"species":species,
				"genename":genename,
				"sequence":sequence,
				"sequence_length":len(sequence.replace("-","").replace(".",""))
				})

		except:
			raise

	return alignment_data


##################################################################################################

def degapSequences(seq_a,seq_b,type="either"):
	seq_a_degapped = ''
	seq_b_degapped = ''
	for i in range(0,len(seq_a)):
		if type == "either":
			if seq_a[i] != "-" or seq_b[i] != "-":
				seq_a_degapped += seq_a[i]
				seq_b_degapped += seq_b[i]
		if type == "both":
			if seq_a[i] != "-" and seq_b[i] != "-":
				seq_a_degapped += seq_a[i]
				seq_b_degapped += seq_b[i]
		if type == "relative":
			if seq_a[i] != "-":
				seq_a_degapped += seq_a[i]
				seq_b_degapped += seq_b[i]
			
	return [seq_a_degapped,seq_b_degapped]

def degapAlignment(alignment_path,query_uniprot_accession=False):
	logger.debug(alignment_path)

	query_protein = True

	alignment_data = {
	"query":{
		"accession":"",
		"mapping":{},
		"query_gaps":[],
		"query_aas":[]
		},
	"proteins":{},
	"order":[]
	}

	alignmentData = readAlignment(alignment_path)

	for alignmentSequenceIdentifier in alignmentData:
		try:
			sequence = alignmentData[alignmentSequenceIdentifier].upper().replace(".","-")
			uniprot_accession = find_uniprot_accession(alignmentSequenceIdentifier)
			
			if (query_protein and query_uniprot_accession == False) or query_uniprot_accession == uniprot_accession:
				query_protein = False
				alignment_data["query"]['accession'] = uniprot_accession

				query_protein_offet = 0
				for pos in range(0,len(sequence)):
					if sequence[pos] == "-":
						alignment_data["query"]["query_gaps"].append(pos)
					else:
						alignment_data["query"]["query_aas"].append(pos)
						alignment_data["query"]["mapping"][query_protein_offet] = pos
						query_protein_offet += 1

			alignment_data['order'].append(uniprot_accession)

			alignment_data["proteins"][uniprot_accession] = {
			"sequence":sequence,
			"sequence_length":len(sequence.replace("-","").replace(".",""))
			}

		except:
			raise

	for protein in alignment_data['order']:
		degapped_sequence = ''
		for pos in alignment_data["query"]["query_aas"]:
			degapped_sequence += alignment_data["proteins"][protein]['sequence'][pos]

		alignment_data["proteins"][protein]["degapped_sequence"] = degapped_sequence
		alignment_data["proteins"][protein]["unaligned_residues_query"] = degapped_sequence.count('-')
		alignment_data["proteins"][protein]["aligned_residues"] = len(degapped_sequence) - degapped_sequence.count('-')
		alignment_data["proteins"][protein]["unaligned_residues"] = alignment_data["proteins"][protein]["sequence_length"] - alignment_data["proteins"][protein]["aligned_residues"]

	return alignment_data

##################################################################################################

def find_uniprot_accession(header):
	try:
		if re.findall("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})", header):
			return re.findall("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})", header)[0][0]
		elif re.match("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\-{0,1}[0-9]*", header):
			return re.match("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\-{0,1}[0-9]*", header)[0][0]
		else:
			return header.split("/")[0]
	except:	
		logger.error("Can't parse:" + header)
		return ""

def find_taxon_id(header):
	try:
		if re.findall("OX=[0-9]+", header):
			return re.findall("OX=[0-9]+", header)[0].split("=")[-1]
		else:
			return ""
	except:	
		logger.error("Can't parse:" + header)
		return ""

def find_genename(header):
	try:
		if re.findall("GN=[a-zA-Z0-9]+", header):
			return re.findall("GN=[a-zA-Z0-9]+", header)[0].split("=")[-1]
		else:
			return ""		
	except:	
		logger.error("Can't parse:" + header)
		return ""

def find_species(header):
	try:
		if re.findall("OS=[^\=]+=", header):
			return" ".join(re.findall("OS=[^\=]+=", header)[0].split("=")[1:-1])[:-3]
		else:
			return ""
	except:
		logger.error("Can't parse:" + header)
		return ""

##################################################################################################

def find_region_boundaries(header):
	if re.findall("\/", header):
		region_boundary = re.findall("\/[0-9]+\-[0-9]+",header)[0]
		region_boundary_start = region_boundary[1:].split("-")[0]
		region_boundary_end = region_boundary[1:].split("-")[1]
		return [region_boundary_start,region_boundary_end]
	else:
		return [-1,-1]

##################################################################################################

def make_alignment_input(accessions,alignment_path=None):
	input_fasta = ""
	for accession in accessions:
		logger.debug("Adding " + accession + " to alignment")
		data = queryRunner.queryRunner("uniprot","parse_uniprot_fasta",accession=accession).run()
		input_fasta  += ">" + data['data']['header'] + "\n"
		input_fasta  += data['data']['sequence'] + "\n"

	if alignment_path != None:
		open(alignment_path,"w").write(input_fasta)

	return input_fasta

##################################################################################################

def replaceNonCommonResidues(peptide):
	processed_peptide  = list(peptide)
	common_aas = list("ACDEFGHIKLMNPQRSTVWY")
	
	for i in range(0,len(peptide)):
		if processed_peptide[i] not in common_aas:
			processed_peptide[i] = "X"

	return "".join(processed_peptide)

def alignPeptides(peptide1,peptide2,gap_open=-10,gap_extend=-1,matrix_type="BLOSUM62"):
	import utilities_error
	
	try:
		from Bio import Align
		from Bio.Align import substitution_matrices

		matrix = substitution_matrices.load(matrix_type)
		aligner = Align.PairwiseAligner()
		aligner.internal_open_gap_score = gap_open
		aligner.internal_extend_gap_score = gap_extend
		aligner.substitution_matrix = matrix

		peptide1 = replaceNonCommonResidues(peptide1)
		peptide2 = replaceNonCommonResidues(peptide2)
		alignments = aligner.align(peptide1, peptide2)
		
		peptide1_aligned = [peptide1[i] if i != -1 else "-" for i in alignments[0].indices[0] ]
		peptide2_aligned = [peptide2[i] if i != -1 else "-" for i in alignments[0].indices[1] ]

		return [[peptide1_aligned,peptide2_aligned]]
	except:
		utilities_error.printError()
		logging.error("Bio.Align import error - trying pairwise2")	

	try:	
		from Bio import pairwise2
		from Bio.SubsMat import MatrixInfo as matlist
		matrix = matlist.blosum62
		alignments = pairwise2.align.globalds( peptide1, "" + peptide2 + "" , matrix, gap_open, gap_extend)
		return alignments
	except:
		logging.error("Bio.SubsMat")
		logging.error([peptide1,peptide2])

	return []

##################################################################################################

def make_alignment_dict(alignment_path):
	import networkx as nx
	G = nx.Graph()

	alignment = parseAlignmentByProtein(alignment_path)
	
	alignment_dict = {}
	for accession in alignment['proteins']:
		aa_offset = 0
		alignment_offset = 0
		
		for position in range(0,len(alignment['proteins'][accession][0]['sequence'])):
			if alignment['proteins'][accession][0]['sequence'][position] != "-":
				aa_offset += 1
			
			if alignment_offset not in alignment_dict:
				alignment_dict[alignment_offset] = {}

			alignment_dict[alignment_offset][accession] = {
				"aa_offset":aa_offset+int(alignment['proteins'][accession][0]['region_boundary_start']),
				"aa":alignment['proteins'][accession][0]['sequence'][alignment_offset]
			}

			alignment_offset += 1
		
	return alignment_dict

def make_alignment_graph(alignment_input_path, alignment_output_path, alignment_algorithm,alignment_algorithm_path,alignment_data):

	open(alignment_input_path,"w").write(sequencesToFasta(alignment_data))
	alignSequences(alignment_input_path, alignment_output_path,alignment_algorithm,alignment_algorithm_path,remake=True,installed=True)
	
	alignment_dict = make_alignment_dict(alignment_output_path)

	import networkx as nx
	import matplotlib.pyplot as plt

	G = nx.Graph()
	
	accessions = list(alignment_dict[0].keys())
	alignment_offsets = list(alignment_dict.keys())

	for alignment_offset in alignment_offsets:
		for accession_iter in range(0,len(accessions)):
			accession = accessions[accession_iter]
			try:
				aa = alignment_dict[alignment_offset][accession]['aa']
				aa_offset = alignment_dict[alignment_offset][accession]['aa_offset']

				if aa == "-":
					G.add_node(accession + "_" + str(alignment_offset), accession=accession, aa=aa, aa_offset=False ,alignment_offset=alignment_offset)#label=accession + ":" + str(aa_offset) + ":" + aa,
				else:
					G.add_node(accession + "_" + str(alignment_offset), accession=accession, aa=aa, aa_offset=aa_offset,alignment_offset=alignment_offset)#label=accession + ":" + str(aa_offset) + ":" + aa,
				
				if accession_iter+1 < len(accessions):
					accession_edge = accessions[accession_iter+1]
					G.add_edge(accession + "_" + str(alignment_offset), accession_edge + "_" + str(alignment_offset), type="homologue")
			except:
				pass

	return G

def find_homologous_residues(accession_search,aa_offset_search_list,alignment_path=None,G=None,group_by="aa_offset"):
	import networkx as nx

	if G == None:
		G = make_alignment_graph(alignment_path)

	nodeDict = dict(G.nodes(data=True))

	mapped_homologous_residues = {}

	filtered_protein_nodes = [x for x,y in G.nodes(data=True) if y['aa_offset'] in aa_offset_search_list and y['accession']==accession_search]	

	filtered_protein_connected_nodes = {}

	for node_id in filtered_protein_nodes:
		filtered_protein_connected_nodes[node_id] = list(nx.node_connected_component(G,node_id))
		
	for node_id in filtered_protein_connected_nodes:
		for filtered_homologue_node_id in filtered_protein_connected_nodes[node_id]:
			filtered_homologue_node = [nodeDict[filtered_homologue_node_id],nodeDict[node_id]]

			if filtered_homologue_node[0]['accession'] not in mapped_homologous_residues:
				mapped_homologous_residues[filtered_homologue_node[0]['accession']] = {}

			mapped_homologous_residues[filtered_homologue_node[0]['accession']][filtered_homologue_node[0][group_by]] = {
				"aa":filtered_homologue_node[0]['aa'],
				"offset":filtered_homologue_node[0]['aa_offset'],
				"accession":filtered_homologue_node[0]['accession'],
				"alignment_offset":filtered_homologue_node[0]['alignment_offset'],
				"query_aa_offset":filtered_homologue_node[1]['aa_offset'],
				"query_aa":filtered_homologue_node[1]['aa'],
				"query_accession":filtered_homologue_node[1]['accession'],
				"alignment_offset":filtered_homologue_node[1]['alignment_offset'],
			}
		
	del G
	del filtered_protein_nodes

	return mapped_homologous_residues

def draw_network(G):
	import networkx as nx
	import matplotlib.pyplot as plt

	pos = nx.circular_layout(G)
	labels = nx.get_node_attributes(G, "name")
	edge_labels = nx.get_edge_attributes(G, "type")

	nx.draw_networkx_labels(G, pos, labels=labels, font_size=10, font_color="white")
	nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10)

	nx.draw(G)
	plt.show()
	

if __name__ == "__main__":
	alignment_path = sys.argv[-1]
	find_homologous_residues(alignment_path,accession_search="P12004",aa_offset_search_list=[100,201,110,200,134])
	sys.exit()

	accessions = ["Q24141","Q5FBB7","Q562F6","F4J3S1","Q0WTB8","Q6FMT2","Q18412","Q4KLP8","Q08490","Q9P7A0","O13734"]
	accessions = list(set(accessions))

	searchdb_options = {
			"dataset_type": "evolution",
			"orthdb_taxon_id": "qfo",
			"verbose": False,
			"logfile": False,
			"debug":False,
			"show_output": False
		}
	
	data = queryRunner.queryRunner("evolution","get_searchdb_taxons",options=searchdb_options).run()
	taxon_ids = data['data']
	taxon_ids = ["9606","10239","7227"]
	for taxon_id in taxon_ids:	
		data = queryRunner.queryRunner("hmm","run_hmm_search_accessions",options={"accessions":accessions,"taxon_id":taxon_id}).run()
		pprint.pprint(data['data'])

