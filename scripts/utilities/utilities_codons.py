import  random, os, pprint, sys

sys.path.append(os.path.join(os.path.dirname(__file__), "../"))

import utilities_basicReader

basepair = {"A":"T","T":"A","G":"C","C":"G","N":"N","-":"-"}

codon_mapping_path = os.path.join(os.path.dirname(__file__), "./codon_table.tdt")

#-----
import logging
logger = logging.getLogger(__name__)
#-----


def get_codon_table():
	codon_table = utilities_basicReader.readTableFile(codon_mapping_path,key="Codon",relationship="binary",byColumn=False,hasHeader=True)
	return codon_table

def translate_oligo(oligo):
	translated_peptide = ''
	codon_table = get_codon_table()
	
	#if len(oligo)%3 != 0:
	#	logger.error("oligo not divisible by 3")
	
	for i in range(0,len(oligo),3):
		if len(oligo[i:i+3]) == 3:
			try:
				translated_peptide += codon_table[oligo[i:i+3]]['AA']
			except:
				translated_peptide += ""
		else:
			translated_peptide += "-"
			
	return translated_peptide

def find_oligo_frame(oligo,peptide):
	oligo_translated_peptides = translate_oligo_frames(oligo)
	oligo_frame_peptide = ""
	
	for oligo_translated_peptide in oligo_translated_peptides:
		if oligo_translated_peptide.count(peptide) > 0:
			oligo_frame_peptide = oligo_translated_peptide.strip("-")
	
	oligo_frame_peptide_bits = oligo_frame_peptide.split(peptide)
	oligo_frame_peptide_response = oligo_frame_peptide_bits[0].lower() + peptide + oligo_frame_peptide_bits[1].lower()

	return oligo_frame_peptide_response

def translate_oligo_frames(oligo,add_details=False,strand=None, frame=None):
	translated_peptides = {}

	for i in range(0,3):
		frame_oligo = oligo[i:]
		frame_oligo_translated = translate_oligo(frame_oligo)

		if strand != None and frame == None:
			if strand == "forward":
				translated_peptides[frame_oligo_translated] = {"strand":"forward","frame":i+1}
		elif strand == None and frame != None:
			if frame == i+1:
				translated_peptides[frame_oligo_translated] = {"strand":"forward","frame":i+1}
		elif strand != None and frame != None:
			if strand == "forward" and frame == i+1:
				translated_peptides[frame_oligo_translated] = {"strand":"forward","frame":i+1}
		else:
			translated_peptides[frame_oligo_translated] = {"strand":"forward","frame":i+1}

	complementry_oligo = get_reverse_complementry_oligo(oligo)

	for i in range(0,3):
		frame_oligo = complementry_oligo[i:]
		frame_oligo_translated = translate_oligo(frame_oligo)	

		if strand != None and frame == None:
			if strand == "reverse":
				translated_peptides[frame_oligo_translated] = {"strand":"reverse","frame":i+1}
		elif strand == None and frame != None:
			if frame == i+1:
				translated_peptides[frame_oligo_translated] = {"strand":"reverse","frame":i+1}
		elif strand != None and frame != None:
			if strand == "reverse" and frame == i+1:
				translated_peptides[frame_oligo_translated] = {"strand":"reverse","frame":i+1}
		else:
			translated_peptides[frame_oligo_translated] = {"strand":"reverse","frame":i+1}


	if add_details:
		return translated_peptides
	else:
		return list(translated_peptides.keys())

def translate_codon(codon):
	codon_table = get_codon_table()

	if len(codon) != 3:
		logger.error("codon not divisible by 3")
		return "*"
	else:
		return codon_table[codon]['AA']
	
def get_complementry_oligo(oligo):
	complementry_oligo = ""
	for nucleotide in oligo:
		complementry_oligo += basepair[nucleotide]
		
	return complementry_oligo

def get_reverse_complementry_oligo(oligo):
	return get_complementry_oligo(oligo)[::-1]

def edit_guide(guide,editor,edit_window):
	editable_nucleotides = {
		"ABE":{"C":"c","A":"G","T":"t","G":"g"},
		"CBE":{"C":"T","A":"a","T":"t","G":"g"}
	}

	edited_guide = list(guide.lower())
	
	for edit_window_position in edit_window:
		edit_window_position_list = edit_window_position - 1 
		wt_nucleotide = guide[edit_window_position_list]

		if wt_nucleotide not in editable_nucleotides[editor]: continue

		mut_nucleotide =  editable_nucleotides[editor][wt_nucleotide]
		edited_guide[edit_window_position_list] = mut_nucleotide

	return edited_guide

def find_edited_mutations(guide,protein_sequence,editor,edit_window):

	frame = ""
	strand = "" 
	oligo_translated_peptides = translate_oligo_frames(guide,add_details=True)
	
	translated_sequence = None
	for oligo_translated_peptide in oligo_translated_peptides:
		if protein_sequence.count(oligo_translated_peptide.strip("-")) and len(oligo_translated_peptide.strip("-")) > 0:
			translated_sequence = oligo_translated_peptides[oligo_translated_peptide]
			translated_sequence['peptide'] = oligo_translated_peptide.strip("-")
			translated_sequence['subpeptide'] = False

	if translated_sequence == None:
		for i in range(1,3):
			if translated_sequence != None: 
				break

			for direction in ['forward','reverse']:
				for oligo_translated_peptide in oligo_translated_peptides:
					if direction == 'forward':
						if protein_sequence.count(oligo_translated_peptide.strip("-")[:-i]):
							translated_sequence = oligo_translated_peptides[oligo_translated_peptide]
							translated_sequence['peptide'] = oligo_translated_peptide.strip("-")
							translated_sequence['subpeptide'] = oligo_translated_peptide.strip("-")[:-i]
					if direction == 'reverse':
						if protein_sequence.count(oligo_translated_peptide.strip("-")[i:]):
							translated_sequence = oligo_translated_peptides[oligo_translated_peptide]
							translated_sequence['peptide'] = oligo_translated_peptide.strip("-")
							translated_sequence['subpeptide'] = oligo_translated_peptide.strip("-")[i:]

	frame = translated_sequence['frame']
	strand = translated_sequence['strand']
	
	edited_guide = edit_guide(guide,editor,edit_window)

	wt_peptide = translate_oligo_frames("".join(guide).upper(),strand=strand,frame=frame)[0].strip("-")
	mut_peptide = translate_oligo_frames("".join(edited_guide).upper(),strand=strand,frame=frame)[0].strip("-")
	
	print(editor,wt_peptide,mut_peptide,edited_guide)

	edit_change = []
	changes_dict = {}
	for i in range(0,len(wt_peptide)):
		if wt_peptide[i] != mut_peptide[i]:
			edit_change.append(wt_peptide[i] + str(i+1) + "->" + mut_peptide[i])
			changes_dict[i+1] = {"wt":wt_peptide[i],"mut":mut_peptide[i]}

	return {"wt":wt_peptide,"mutant":mut_peptide,"changes":edit_change,"changes_dict":changes_dict,"frame":frame,"strand":strand,'subpeptide':translated_sequence['subpeptide']}
