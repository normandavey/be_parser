import os, inspect, sys, pprint, bisect,logging

import numpy as np
file_path = os.path.dirname(inspect.stack()[0][1])

sys.path.append(os.path.join(file_path,"../"))
import config_reader

sys.path.append(os.path.join(file_path,"../utilities"))
import utilities_error

sys.path.append(os.path.join(file_path,"../structure"))
from structure.dssp.dssp import Dssp

sys.path.append(os.path.join(file_path,"../data_management/"))
import queryRunner

sys.path.append(os.path.join(file_path,"../pocket_discovery/"))
from pdb_parser import read_clean_pdb


from Bio.PDB import PDBParser, Selection

try:
	from Bio.PDB.Polypeptide import three_to_one
except:
	from Bio.PDB.Polypeptide import three_to_index
	from Bio.PDB.Polypeptide import index_to_one

from Bio import Align

try:
	from Bio.Align import substitution_matrices
except:
	logging.error("Error from 'from Bio.Align import substitution_matrices'")

from scipy.spatial.distance import euclidean
from scipy.spatial import Delaunay

#-----
import logging
logger = logging.getLogger(__name__)
#-----

#source https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html
substitutions = {
    '2AS':'ASP', '3AH':'HIS', '5HP':'GLU', 'ACL':'ARG', 'AGM':'ARG', 'AIB':'ALA', 'ALM':'ALA', 'ALO':'THR', 'ALY':'LYS', 'ARM':'ARG',
    'ASA':'ASP', 'ASB':'ASP', 'ASK':'ASP', 'ASL':'ASP', 'ASQ':'ASP', 'AYA':'ALA', 'BCS':'CYS', 'BHD':'ASP', 'BMT':'THR', 'BNN':'ALA',
    'BUC':'CYS', 'BUG':'LEU', 'C5C':'CYS', 'C6C':'CYS', 'CAS':'CYS', 'CCS':'CYS', 'CEA':'CYS', 'CGU':'GLU', 'CHG':'ALA', 'CLE':'LEU', 'CME':'CYS',
    'CSD':'ALA', 'CSO':'CYS', 'CSP':'CYS', 'CSS':'CYS', 'CSW':'CYS', 'CSX':'CYS', 'CXM':'MET', 'CY1':'CYS', 'CY3':'CYS', 'CYG':'CYS',
    'CYM':'CYS', 'CYQ':'CYS', 'DAH':'PHE', 'DAL':'ALA', 'DAR':'ARG', 'DAS':'ASP', 'DCY':'CYS', 'DGL':'GLU', 'DGN':'GLN', 'DHA':'ALA',
    'DHI':'HIS', 'DIL':'ILE', 'DIV':'VAL', 'DLE':'LEU', 'DLY':'LYS', 'DNP':'ALA', 'DPN':'PHE', 'DPR':'PRO', 'DSN':'SER', 'DSP':'ASP',
    'DTH':'THR', 'DTR':'TRP', 'DTY':'TYR', 'DVA':'VAL', 'EFC':'CYS', 'FLA':'ALA', 'FME':'MET', 'GGL':'GLU', 'GL3':'GLY', 'GLZ':'GLY',
    'GMA':'GLU', 'GSC':'GLY', 'HAC':'ALA', 'HAR':'ARG', 'HIC':'HIS', 'HIP':'HIS', 'HMR':'ARG', 'HPQ':'PHE', 'HTR':'TRP', 'HYP':'PRO',
    'IAS':'ASP', 'IIL':'ILE', 'IYR':'TYR', 'KCX':'LYS', 'LLP':'LYS', 'LLY':'LYS', 'LTR':'TRP', 'LYM':'LYS', 'LYZ':'LYS', 'MAA':'ALA', 'MEN':'ASN',
    'MHS':'HIS', 'MIS':'SER', 'MLE':'LEU', 'MPQ':'GLY', 'MSA':'GLY', 'MSE':'MET', 'MVA':'VAL', 'NEM':'HIS', 'NEP':'HIS', 'NLE':'LEU',
    'NLN':'LEU', 'NLP':'LEU', 'NMC':'GLY', 'OAS':'SER', 'OCS':'CYS', 'OMT':'MET', 'PAQ':'TYR', 'PCA':'GLU', 'PEC':'CYS', 'PHI':'PHE',
    'PHL':'PHE', 'PR3':'CYS', 'PRR':'ALA', 'PTR':'TYR', 'PYX':'CYS', 'SAC':'SER', 'SAR':'GLY', 'SCH':'CYS', 'SCS':'CYS', 'SCY':'CYS',
    'SEL':'SER', 'SEP':'SER', 'SET':'SER', 'SHC':'CYS', 'SHR':'LYS', 'SMC':'CYS', 'SOC':'CYS', 'STY':'TYR', 'SVA':'SER', 'TIH':'ALA',
    'TPL':'TRP', 'TPO':'THR', 'TPQ':'ALA', 'TRG':'LYS', 'TRO':'TRP', 'TYB':'TYR', 'TYI':'TYR', 'TYQ':'TYR', 'TYS':'TYR', 'TYY':'TYR'
}

#source patchfinder codebase
#Detection of Functionally Important Regions in “Hypothetical Proteins” of Known Structure
#http://patchfinder.tau.ac.il/
#https://www.dropbox.com/s/6cej7i7tm2u9keq/patchfinder_del.zip?dl=0
rad_siz = {'ALAN':1.65, 'ALACA':1.87, 'ALAC':1.76, 'ALAO':1.40, 'ALACB':1.87, 'ARGN':1.65, 'ARGCA':1.87, 'ARGC':1.76, 'ARGO':1.40, 'ARGCB':1.87, 'ARGCG':1.87, 'ARGCD':1.87, 'ARGNE':1.65, 'ARGCZ':1.76, 'ARGNH1': 1.65, 'ARGNH2': 1.65, 'ASPN':1.65, 'ASPCA':1.87, 'ASPC':1.76, 'ASPO':1.40, 'ASPCB':1.87, 'ASPCG':1.76, 'ASPOD1': 1.40, 'ASPOD2': 1.40, 'ASNN':1.65, 'ASNCA':1.87, 'ASNC':1.76, 'ASNO':1.40, 'ASNCB':1.87, 'ASNCG':1.76, 'ASNOD1': 1.40, 'ASNND2': 1.65, 'CYSN':1.65, 'CYSCA':1.87, 'CYSC':1.76, 'CYSO':1.40, 'CYSCB':1.87, 'CYSSG':1.85, 'GLUN':1.65, 'GLUCA':1.87, 'GLUC':1.76, 'GLUO':1.40, 'GLUCB':1.87, 'GLUCG':1.87, 'GLUCD':1.76, 'GLUOE1': 1.40, 'GLUOE2': 1.40, 'GLNN':1.65, 'GLNCA':1.87, 'GLNC':1.76, 'GLNO':1.40, 'GLNCB':1.87, 'GLNCG':1.87, 'GLNCD':1.76, 'GLNOE1': 1.40, 'GLNNE2': 1.65, 'GLYN':1.65, 'GLYCA':1.87, 'GLYC':1.76, 'GLYO':1.40, 'HISN':1.65, 'HISCA':1.87, 'HISC':1.76, 'HISO':1.40, 'HISCB':1.87, 'HISCG':1.76, 'HISND1': 1.65, 'HISCD2': 1.76, 'HISCE1': 1.76, 'HISNE2': 1.65, 'ILEN':1.65, 'ILECA':1.87, 'ILEC':1.76, 'ILEO':1.40, 'ILECB':1.87, 'ILECG1': 1.87, 'ILECG2': 1.87, 'ILECD1': 1.87, 'LEUN':1.65, 'LEUCA':1.87, 'LEUC':1.76, 'LEUO':1.40, 'LEUCB':1.87, 'LEUCG':1.87, 'LEUCD1': 1.87, 'LEUCD2': 1.87, 'LYSN':1.65, 'LYSCA':1.87, 'LYSC':1.76, 'LYSO':1.40, 'LYSCB':1.87, 'LYSCG':1.87, 'LYSCD':1.87, 'LYSCE':1.87, 'LYSNZ':1.50, 'METN':1.65, 'METCA':1.87, 'METC':1.76, 'METO':1.40, 'METCB':1.87, 'METCG':1.87, 'METSD':1.85, 'METCE':1.87, 'PHEN':1.65, 'PHECA':1.87, 'PHEC':1.76, 'PHEO':1.40, 'PHECB':1.87, 'PHECG':1.76, 'PHECD1': 1.76, 'PHECD2': 1.76, 'PHECE1': 1.76, 'PHECE2': 1.76, 'PHECZ':1.76, 'PRON':1.65, 'PROCA':1.87, 'PROC':1.76, 'PROO':1.40, 'PROCB':1.87, 'PROCG':1.87, 'PROCD':1.87, 'SERN':1.65, 'SERCA':1.87, 'SERC':1.76, 'SERO':1.40, 'SERCB':1.87, 'SEROG':1.40, 'THRN':1.65, 'THRCA':1.87, 'THRC':1.76, 'THRO':1.40, 'THRCB':1.87, 'THROG1': 1.40, 'THRCG2': 1.87, 'TRPN':1.65, 'TRPCA':1.87, 'TRPC':1.76, 'TRPO':1.40, 'TRPCB':1.87, 'TRPCG':1.76, 'TRPCD1': 1.76, 'TRPCD2': 1.76, 'TRPNE1': 1.65, 'TRPCE2': 1.76, 'TRPCE3': 1.76, 'TRPCZ2': 1.76, 'TRPCZ3': 1.76, 'TRPCH2': 1.76, 'TYRN':1.65, 'TYRCA':1.87, 'TYRC':1.76, 'TYRO':1.40, 'TYRCB':1.87, 'TYRCG':1.76, 'TYRCD1': 1.76, 'TYRCD2': 1.76, 'TYRCE1': 1.76, 'TYRCE2': 1.76, 'TYRCZ':1.76, 'TYROH':1.40, 'VALN':1.65, 'VALCA':1.87, 'VALC':1.76, 'VALO':1.40, 'VALCB':1.87, 'VALCG1': 1.87, 'VALCG2': 1.87}

class simpleMappingAccessibilityUtilities():
	def __init__(self):
		self.options = {}

		config_options = config_reader.load_configeration_options(sections=["general","accessibility"])

		for config_option in config_options:
			self.options[config_option] = config_options[config_option]

	def get_pdb_id(self):
		if isinstance(self.options["pdb_id"],(list)):
			return self.options["pdb_id"][0]
		else:
			return self.options["pdb_id"]	

	def alignPeptides(self, peptide1, peptide2, gap_open=-10, gap_extend=-1):
		aligner = Align.PairwiseAligner()
		aligner.internal_open_gap_score = gap_open
		aligner.internal_extend_gap_score = gap_extend
		aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
		alignments = aligner.align(peptide1, peptide2)
		
		# @Hazem  Found a bug here so updated code
		# alignments_str = str(alignments[0]).split('\n')
		
		pdb_sequence_aligned = [peptide1[i] if i != -1 else "-" for i in alignments[0].indices[0] ]
		uniprot_sequence_aligned = [peptide2[i] if i != -1 else "-" for i in alignments[0].indices[1] ]
	
		
		return {
	 		'num_alignments': len(alignments),
			'pdb_sequence_aligned':pdb_sequence_aligned,
			'uniprot_sequence_aligned':uniprot_sequence_aligned
		}

	def map_pdb_offsets_to_uniprot(self,pdb_id,chain_id,offset_keys,offset_keys_no_icode,offset_sequence,offset_standard_aa,return_alignment=False):
		logger.debug("map_pdb_offsets_to_uniprot " + pdb_id + " " + chain_id)
	
		pdb_id = self.get_pdb_id()
		uniprot_mapping_structure_data = queryRunner.queryRunner("pdb","get_sequence_structure",{"pdb_id":pdb_id,"chain_id":chain_id.split(".")[-1]}).run()

		residue_mapping = {}

		try:
			if 'data' not in uniprot_mapping_structure_data:
				return residue_mapping
				
			offset_keys_no_icode_sorted = np.sort(offset_keys_no_icode)
			offset_keys_no_icode_arg_sorted = np.argsort(offset_keys_no_icode)

			for uniprot_accession in uniprot_mapping_structure_data['data']:
				for region in uniprot_mapping_structure_data['data'][uniprot_accession]:
					uniprot_sequence = region['uniprot_sequence']
					uniprot_start = region['unp_start']
					pdb_sequence_start = str(region['start']['author_residue_number'])+str(region['start']['author_insertion_code'])
					pdb_sequence_end = str(region['end']['author_residue_number'])+region['end']['author_insertion_code']
					
					try:
						try:
							pdb_start = offset_keys.index(pdb_sequence_start)
						except:
							# Handle the case when the start residue is missing in the pdb file
							#if region['start']['author_insertion_code']:
								#print("Potential bug (minor) in map_pdb_offsets_to_uniprot in finding region start", pdb_id,chain_id)
							
							pdb_start = bisect.bisect_left(offset_keys_no_icode_sorted, region['start']['author_residue_number'])
							
							if pdb_start >= len(offset_keys):
								#happens due to problems in pdb files and/or the response from SIFTS api
								continue
							pdb_start = offset_keys_no_icode_arg_sorted[pdb_start]


						try:
							pdb_end = offset_keys.index(pdb_sequence_end)
						except:
							# Handle the case when the end residue is missing in the pdb file
							#if region['end']['author_insertion_code']:
								#print("Potential bug (minor) in map_pdb_offsets_to_uniprot in finding region end", pdb_id,chain_id)
								
							pdb_end = bisect.bisect_right(offset_keys_no_icode_sorted, region['end']['author_residue_number'])-1
							
							if pdb_end < 0:
								#happens due to problems in pdb files and/or the response from SIFTS api
								continue
							pdb_end = offset_keys_no_icode_arg_sorted[pdb_end]
					except:
						pdb_start = region['pdb_start'] - 1
						pdb_end = region['pdb_end']
					
					if pdb_end < pdb_start:
						#happens due to problems in pdb files and/or the response from SIFTS api
						continue
						
					pdb_sequence = offset_sequence[pdb_start:pdb_end+1]
					
					sequence_alignment = self.alignPeptides(pdb_sequence,uniprot_sequence)
					pdb_sequence_aligned = sequence_alignment['pdb_sequence_aligned']
					uniprot_sequence_aligned = sequence_alignment['uniprot_sequence_aligned']
					
					# This problem can be solved by inserting missing residues before performing the alignment.
					#if sequence_alignment['num_alignments'] > 3:
						#print("Skipped a region due to multiple pairwise alignments", chain_id, sequence_alignment)
						#print(pdb_sequence,uniprot_sequence)
						
					pdb_offset = pdb_start
					uniprot_offset = uniprot_start
					
					for i in range(0,len(pdb_sequence_aligned)):
						try:
							pdb_sequence_aligned_aa = pdb_sequence_aligned[i]
							uniprot_sequence_aligned_aa = uniprot_sequence_aligned[i]
							
							if pdb_sequence_aligned_aa != '-':
								pdb_key = offset_keys[pdb_offset]
								pdb_aa = offset_sequence[pdb_offset]
								is_standard = offset_standard_aa[pdb_offset]
								
								if pdb_aa != pdb_sequence_aligned_aa:
									logger.debug(["This is bug in map_pdb_offsets_to_uniprot and the code shouldn't enter here",str(i),chain_id])
								
							aligned = pdb_sequence_aligned_aa == uniprot_sequence_aligned_aa
							
							#if aligned:
								#if pdb_sequence_aligned_aa == 'X':
									#print("This is bug in map_pdb_offsets_to_uniprot due to bio pairwise alignment: Found X aligned", pdb_id,chain_id)
									#print(pdb_sequence, uniprot_sequence)
									#continue
							if pdb_sequence_aligned_aa != "-":	
								if pdb_key not in residue_mapping:
									residue_mapping[pdb_key] = {}

								residue_mapping[pdb_key][uniprot_accession] = {
									"uniprot_offset":uniprot_offset,
									"uniprot_sequence_aligned_aa":uniprot_sequence_aligned_aa,
									"pdb_standard_aa":is_standard,
									"aligned":aligned
								}
							
							if pdb_sequence_aligned_aa != "-": pdb_offset += 1
							if uniprot_sequence_aligned_aa != "-": uniprot_offset += 1
						except:
							logger.error(["This is bug in map_pdb_offsets_to_uniprot and the code shouldn't enter here",i,pdb_id,chain_id])

			logger.debug("Complete")
		except:
			utilities_error.printError()
			logger.error(chain_id)

		if return_alignment:
			return {
				"residue_mapping":residue_mapping,
				"region_mapping":uniprot_mapping_structure_data['data']
			}
		else:
			return { 
				"residue_mapping": residue_mapping
			} 
	
	#get_surface_accessibility_data in dssp_utilities returns different accessibility scores from monomers (single chains), joined (chains in complex), and the difference between them.
	#This service was built to get normalised accessibility scores based on monomers (the only required score in pocket discovery) from dssp after cleaning/fixing some PDB issues (use_cleaner_fixer=True).
	#Now, pocket discovery is using get_tessellation_accessibility, but this function will remain for other services which might become relying on it.
	
	def get_surface_accessibility_normalised(self):
			
		pdb_id = self.get_pdb_id()
		pdb_structure_data = queryRunner.queryRunner("pdb","download_pdb_structure",{"pdb_id":pdb_id}).run()
		
		pdb_file = pdb_structure_data['data']
		dssp_runner = Dssp(self.options['dssp_path'])

		try:
			dssp_splitted, _ = dssp_runner.run(pdb_file, pdb_id, splitted=True, use_cleaner_fixer=True)
		except:
			logger.error("ERROR: There was a problem running DSSP for pdb: {}".format(pdb_file))
			raise

		if 'chain_id' in self.options:
			pdb_chains = [pdb_id + "." + chain_id for chain_id in self.options['chain_id']]
		else:
			pdb_chains = list(dssp_splitted.to_dict().keys())
		

		self.data = {
		"chain_details":{}
		}

		for pdb_chain in pdb_chains:
		
			sa_normalised_splitted = dssp_splitted.chain(pdb_chain).get_summary().surface_accessibility_normalised if pdb_chain in dssp_splitted else {}

			chain_details = {}
			chain_details['surface_accessibility_normalised'] = sa_normalised_splitted
			
			self.data['chain_details'][pdb_chain] = chain_details
			
		return self.data
		
	#Built to replace get_pdb_uniprot_mapping in dssp_utilities.
	#Get the pdb sequence by reading the pdb file residues. The previous service was using dssp offsets and polymer_entities to get the pdb sequence.
	#Solved some issues related to mapping multiple mapping regions returned from www.ebi.ac.uk APIs.
	#Used in pocket discovery.
	def get_pdb_uniprot_simple_mapping(self):
			
		pdb_id = self.get_pdb_id()
		pdb_structure_data = queryRunner.queryRunner("pdb","download_pdb_structure",{"pdb_id":pdb_id}).run()
		
		pdb_file = pdb_structure_data['data']
		
		p = PDBParser()
		structure = p.get_structure(pdb_id, pdb_file)
		model = structure[0]
		pdb_chains = [c.id for c in Selection.unfold_entities(model, "C")]

		if 'chain_id' in self.options:
			pdb_chains = [chain_id for chain_id in self.options['chain_id']]
		

		self.data = {
		"chain_details":{}
		}

		for pdb_chain in pdb_chains:
		
			chain = model[pdb_chain]
			
			residues = Selection.unfold_entities(chain, "R")
			offset_keys = []
			offset_keys_no_icode = []
			offset_sequence = ""
			offset_standard_aa = []
			for r in residues:

				resseq = r.get_id()[1]
				icode = r.get_id()[2].strip()
				key = str(resseq)+icode

				standard = True if r.get_resname() not in substitutions else False
				resname = substitutions[r.get_resname()] if not standard else r.get_resname()

				try:
					try:
						one_letter = three_to_one(resname)	
					except:
						one_letter = index_to_one(three_to_index(resname))
				except:
					all_atoms = [atom.get_name().strip() for atom in r.get_atoms()]
					# skip residues lacking a backbone atom
					if "C" not in all_atoms or "CA" not in all_atoms or "N" not in all_atoms:
						continue
					one_letter = "X"
					standard = False

				offset_keys.append(key)
				offset_keys_no_icode.append(resseq)
				offset_sequence += one_letter
				offset_standard_aa.append(standard)

			simple_uniprot_mapping = self.map_pdb_offsets_to_uniprot(pdb_id, pdb_id + "." + pdb_chain, offset_keys, offset_keys_no_icode, offset_sequence, offset_standard_aa, return_alignment=False)
			chain_details = simple_uniprot_mapping
			
			
			self.data['chain_details'][pdb_id + "." + pdb_chain] = chain_details
		
		return self.data

	def get_uniprot_pdb_simple_mapping(self):
		pdb_uniprot_simple_mapping = self.get_pdb_uniprot_simple_mapping()
		uniprot_pdb_simple_mapping = {}
		for chain_id in pdb_uniprot_simple_mapping['chain_details']:
			for pdb_offset in pdb_uniprot_simple_mapping['chain_details'][chain_id]['residue_mapping'].keys():
				for uniprot_accession in pdb_uniprot_simple_mapping['chain_details'][chain_id]['residue_mapping'][pdb_offset].keys():
					mapping = pdb_uniprot_simple_mapping['chain_details'][chain_id]['residue_mapping'][pdb_offset][uniprot_accession]
					mapping['pdb_offset'] = pdb_offset
					uniprot_offset = mapping['uniprot_offset'] 
					del  mapping['uniprot_offset']
					
					if uniprot_accession not in uniprot_pdb_simple_mapping:
						uniprot_pdb_simple_mapping[uniprot_accession] = {}

					if uniprot_offset not in uniprot_pdb_simple_mapping[uniprot_accession]:
						uniprot_pdb_simple_mapping[uniprot_accession][uniprot_offset] = {}
					
					uniprot_pdb_simple_mapping[uniprot_accession][uniprot_offset][chain_id.split(".")[1]] = mapping

		return uniprot_pdb_simple_mapping
		
	def get_tessellation_accessibility(self, pdb_id=None, pdb_file=None, alphafold=False, domain_residues=None,filter_inaccessibility=True,add_atom_details=True,skip_atoms=[]):	#the 4 parameters are used in the case of alphafold, the 4th one is used when called from pocket_discovery
		peptide_offset_shift_start = 0
		peptide_offset_shift_end = 0

		if alphafold:
			self.options['chain_id'] = ["A"]
		else:
			pdb_id = self.get_pdb_id()
			pdb_structure_data = queryRunner.queryRunner("pdb","download_pdb_structure",{"pdb_id":pdb_id}).run()
			
			pdb_file = pdb_structure_data['data']
		
		if 'chain_id' in self.options:
			pdb_chains = self.options['chain_id']
		else:
			parser = PDBParser()
			structure = parser.get_structure(pdb_id, pdb_file)
			model = structure[0]
			pdb_chains = [c.id for c in Selection.unfold_entities(model, "C")]
		
		self.data = {
		"chain_details":{}
		}

		if alphafold:
			## Reads the Uniprot strat stop positions of the PDB file
			db_ref = [line for line in open(pdb_file).read().split("\n") if line[0:5] == "DBREF"]
			
			if len(db_ref) > 0:
				logging.debug(pdb_file + " " + db_ref[0])
				peptide_offset_shift_start = int(db_ref[0].split()[-2]) - 1
				peptide_offset_shift_end = int(db_ref[0].split()[-1]) - 1
				
		for pdb_chain in pdb_chains:
			model = read_clean_pdb(pdb_id, pdb_file, pdb_chain, domain_residues)
			
			atm_keys = []
			coords = []

			for r in model.get_residues():
				res_id = r.get_full_id()[3]
				res_key = str(res_id[1]+peptide_offset_shift_start)+res_id[2].strip()
				for atom in r.get_atoms():
					if atom.get_name().strip() != "H" and atom.get_name().strip() not in skip_atoms:
						atm_keys.append(res_key+'_'+r.get_resname().strip()+'_'+atom.get_name().strip())
						coords.append(atom.get_coord())
			
			tessellationUtilitiesObj = tessellationUtilities(np.array(coords), np.array(atm_keys),filter_inaccessibility=filter_inaccessibility,add_atom_details=add_atom_details)

			accessible_residues, direct_neighbors, expanded_neighbors = tessellationUtilitiesObj.get_tessellation_accessible_residues_and_their_neighbors()
					
			chain_details = {
				'accessible_residues': accessible_residues,
				'direct_neighbors': direct_neighbors,
				'expanded_neighbors': expanded_neighbors,
				"peptide_offset_shift_centre":(peptide_offset_shift_start + peptide_offset_shift_end - 1)/2,
				"peptide_offset_shift_start":peptide_offset_shift_start,
				"peptide_offset_shift_end":peptide_offset_shift_end
			}
			
			self.data['chain_details'][pdb_chain] = chain_details
			
		return self.data

#referece: patchfinder
class tessellationUtilities():
	def __init__(self, coords, atm_keys,filter_inaccessibility=True,add_atom_details=True,distance_cutoff = 2.8):
		self.faces = {}

		self.distance_cutoff = distance_cutoff

		self.filter_inaccessibility = filter_inaccessibility
		self.add_atom_details = add_atom_details
		self.faces_tetrahedrons1 = {}
		self.faces_tetrahedrons2 = {}
		self.removed_faces = {}
		
		self.coords = coords
		self.atm_keys = atm_keys
	
	def get_tessellation_accessible_residues_and_their_neighbors(self):
		
		accessible_residues = {}
		direct_neighbors = {}
		expanded_neighbors = {}
		tri = Delaunay(self.coords)
		updated_cavity_set = True
		
		for indx, tetrahedron in enumerate(tri.simplices):
			tetrahedron.sort()
			for i in range(4):
				mask = np.ones(4, np.bool)
				mask[i] = 0
				self.fill_data_dicts(" ".join(np.char.mod('%d', tetrahedron[mask])), indx)
		
		while updated_cavity_set:
			updated_cavity_set = False
			for face, num_tetrahedrons in self.faces.items():
				if num_tetrahedrons == 1:
					a, b, c = face.split()
					if self.VDW_distance(a, b) >= self.distance_cutoff or self.VDW_distance(a, c) >= self.distance_cutoff or self.VDW_distance(c, b) >= self.distance_cutoff:
						tetrahedron = tri.simplices[self.faces_tetrahedrons1[face]]
						tetrahedron_indx = self.faces_tetrahedrons1[face]
						for i in range(4):
							mask = np.ones(4, np.bool)
							mask[i] = 0
							self.update_data_dicts(" ".join(np.char.mod('%d', tetrahedron[mask])), tetrahedron_indx)
						updated_cavity_set = True
		
		for face, num_tetrahedrons in self.faces.items():
			if num_tetrahedrons == 1 or self.removed_faces[face] == 1:
				a, b, c = face.split()
				self.add_accessible_atm(a, accessible_residues)
				self.add_accessible_atm(b, accessible_residues)
				self.add_accessible_atm(c, accessible_residues)
		
		#get tri direct neighbors
		for face, num_tetrahedrons in self.faces.items():
			if self.filter_inaccessibility == False:
				a, b, c = face.split()
				if self.VDW_distance(a, b) <= self.distance_cutoff:
					self.connect_neighbor_atms(a, b, direct_neighbors)
				if self.VDW_distance(a, c) <= self.distance_cutoff:
					self.connect_neighbor_atms(a, c, direct_neighbors)
				if self.VDW_distance(b, c) <= self.distance_cutoff:
					self.connect_neighbor_atms(b, c, direct_neighbors)
			else:
				if num_tetrahedrons == 1:
					a, b, c = face.split()
					self.connect_neighbor_atms(a, b, direct_neighbors)
					self.connect_neighbor_atms(a, c, direct_neighbors)
					self.connect_neighbor_atms(b, c, direct_neighbors)
				
		#combine tri expanded neighbors
		for res_key1 in direct_neighbors:
			expanded_neighbors[res_key1] = set()
			expanded_neighbors[res_key1].update(direct_neighbors[res_key1])
			for res_key2 in direct_neighbors[res_key1]:
				expanded_neighbors[res_key1].update(direct_neighbors[res_key2])
			expanded_neighbors[res_key1].discard(res_key1)
				
		#sets to lists to make them compatible with json format
		for res_key in accessible_residues:
			accessible_residues[res_key]['accessible_atms'] = list(accessible_residues[res_key]['accessible_atms'])

		if self.add_atom_details == False:
			for res_key in direct_neighbors:
				direct_neighbors[res_key] = list(direct_neighbors[res_key])
		
		for res_key in expanded_neighbors:
			expanded_neighbors[res_key] = list(expanded_neighbors[res_key])
		
		#add non-accessible residues to accessible_residues with zero scores
		for atm_key in self.atm_keys:
			res_key, res_name, atm = atm_key.split('_')

			try:
				one_letter = three_to_one(res_name)	
			except:
				one_letter = index_to_one(three_to_index(res_name))
				
			if res_key not in accessible_residues:
				accessible_residues[res_key] = {'accessible_atms': [], 'side_chain_score': 0, 'any_atm_score': 0,"aa":one_letter}
			else:
				accessible_residues[res_key]["aa"] = one_letter
		
		return accessible_residues, direct_neighbors, expanded_neighbors

	def fill_data_dicts(self, face, indx):
		if face in self.faces:
			self.faces[face]+=1
			self.faces_tetrahedrons2[face] = indx
		else:
			self.faces[face]=1
			self.faces_tetrahedrons1[face] = indx
			self.faces_tetrahedrons2[face] = -1
			self.removed_faces[face] = 0

	def VDW_rad(self, indx):
		res_key, res_name, atm = self.atm_keys[int(indx)].split('_')
		rad_key = res_name+atm
		return rad_siz[rad_key] if rad_key in rad_siz else 1.7
		
	def VDW_distance(self, atm1_indx, atm2_indx):
		return euclidean(self.coords[int(atm1_indx)], self.coords[int(atm2_indx)])-self.VDW_rad(atm1_indx)-self.VDW_rad(atm2_indx)

	def update_data_dicts(self, face, indx):
		self.faces[face]-=1
		if self.faces[face] == 0:
			self.removed_faces[face] = 1
		if self.faces_tetrahedrons1[face] == indx:
			self.faces_tetrahedrons1[face] = self.faces_tetrahedrons2[face]

	def add_accessible_atm(self, atm_indx, accessible_residues):
		atm_key = self.atm_keys[int(atm_indx)]
		res_key, res_name, atm = atm_key.split('_')
		if res_key not in accessible_residues:
			accessible_residues[res_key] = {}
			accessible_residues[res_key]['accessible_atms'] = set()
			accessible_residues[res_key]['any_atm_score'] = 0
			accessible_residues[res_key]['side_chain_score'] = 0
		
		if atm not in accessible_residues[res_key]['accessible_atms']:
			accessible_residues[res_key]['any_atm_score'] += 1
			if atm not in {'CA', 'C', 'N', 'O'}:
				accessible_residues[res_key]['side_chain_score'] += 1

		accessible_residues[res_key]['accessible_atms'].add(atm)

	def pair_type(self,res1_name,res2_name,atm1,atm2):
		atm1_type = "sidechain"
		atm2_type = "sidechain"

		if res1_name == "GLY":
			if atm1 in ['C', 'N', 'O']:atm1_type="backbone"
		else:
			if atm1 in ['CA', 'C', 'N', 'O']:atm1_type="backbone"
			
		if res2_name == "GLY":
			if atm2 in ['C', 'N', 'O']:atm2_type="backbone"
		else:
			if atm2 in ['CA', 'C', 'N', 'O']:atm2_type="backbone"

		return [atm1_type,atm2_type]
	
	def connect_neighbor_atms(self, atm1_indx, atm2_indx, direct_neighbors):
		atm1_key = self.atm_keys[int(atm1_indx)]
		res1_key, res1_name, atm1 = atm1_key.split('_')
		atm2_key = self.atm_keys[int(atm2_indx)]
		res2_key, res2_name, atm2 = atm2_key.split('_')
		
		if self.add_atom_details:
			if res1_key != res2_key:
				if res1_key not in direct_neighbors: direct_neighbors[res1_key] = {}
				if res2_key not in direct_neighbors[res1_key]: direct_neighbors[res1_key][res2_key] = {}
				if 	atm1+'-'+atm2 not in direct_neighbors[res1_key][res2_key]:
					
					direct_neighbors[res1_key][res2_key][atm1+'-'+atm2] = "-".join(self.pair_type(res1_name,res2_name,atm1,atm2)) #{"atm1":atm1,"atm2":atm2,"pair_type":"-".join(self.pair_type(atm1,atm2))}
	
				if res2_key not in direct_neighbors: direct_neighbors[res2_key] = {}
				if res1_key not in direct_neighbors[res2_key]: direct_neighbors[res2_key][res1_key] = {}
				if 	atm1+'-'+atm2 not in direct_neighbors[res2_key][res1_key]:
					direct_neighbors[res2_key][res1_key][atm2+'-'+atm1] = "-".join(self.pair_type(res1_name,res2_name,atm1,atm2)) #{"atm1":atm1,"atm2":atm2,"pair_type":"-".join(self.pair_type(atm1,atm2))}

				"""
				if 	atm1+'-'+atm2 not in direct_neighbors[res1_key][res2_key]:
					direct_neighbors[res1_key][res2_key][atm1+'-'+atm2] = {
						"type":"-".join(self.pair_type(res1_name,res2_name,atm1,atm2)),
						"distance":self.VDW_distance(atm1_indx,atm2_indx)
					 }

				if res2_key not in direct_neighbors: direct_neighbors[res2_key] = {}
				if res1_key not in direct_neighbors[res2_key]: direct_neighbors[res2_key][res1_key] = {}
				if 	atm1+'-'+atm2 not in direct_neighbors[res2_key][res1_key]:
					direct_neighbors[res2_key][res1_key][atm2+'-'+atm1] = {
						"type":"-".join(self.pair_type(res1_name,res2_name,atm1,atm2)),
						"distance":self.VDW_distance(atm1_indx,atm2_indx)
					 }
				"""
		else:
			if res1_key != res2_key:
				if res1_key not in direct_neighbors:
					direct_neighbors[res1_key] = set()
				direct_neighbors[res1_key].add(res2_key)
				if res2_key not in direct_neighbors:
					direct_neighbors[res2_key] = set()
				direct_neighbors[res2_key].add(res1_key)

if __name__ == "__main__":
	simpleMappingAccessibilityUtilitiesOj = simpleMappingAccessibilityUtilities()

	pprint.pprint("Nothing to do...")
