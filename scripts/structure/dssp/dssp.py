from os.path import dirname
from os.path import join
from os.path import splitext
from os.path import exists
from os import remove, system as _run_cmd
from Bio.PDB import is_aa
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBIO import Select
from Bio.PDB.PDBParser import PDBParser
from subprocess import check_output
from subprocess import PIPE
from subprocess import CalledProcessError
from structure.dssp.parser import Parser

#-----
import logging
logger = logging.getLogger(__name__)
#-----

pdbfixer_loaded = True

"""
"""

'''
try:
	from structure.dssp.pdbfixer.pdbfixer import PDBFixer
	from simtk.openmm.app import PDBFile

	pdbfixer_loaded = True
except:
	logger.error("Issue loading pdbfixer")
'''


class Dssp(object):
	def __init__(self, path_to_executable="dssp"):
		"""
		Creates a new Dssp object, the path to dssp executable can be supplied. If not, it is
		assumed that the executable name is 'dssp' and is accessible in current enviroment.
			:param self:
			:param path_to_executable="dssp":
		"""
		self.dssp = path_to_executable

	def _save_dssp(self, save_to_file, dssp_results):
		saved_files = {c: "{}.dssp".format(splitext(pdb)[0])
			for c,(pdb, _) in dssp_results.items()} if save_to_file else {}
		for c, dssp_file in saved_files.items():
			with open(dssp_file, 'w') as fh:
				_, dssp_output = dssp_results[c]

				try:
					logger.debug("Writing:" + dssp_file)

					if len(dssp_output) > 0:
						fh.write(dssp_output.decode("utf-8"))
					else:
						fh.write("".decode("utf-8"))
				except:
					logger.error("Error writing" + dssp_file + " - " + dssp_output)

		return saved_files

	def _parse_dssp_results(self, dssp_results, pdb_id, splitted):
		parsed = {c: Parser().dssp_parse(results, [], pdb_id, True, True)
			for c,(pdb, results) in dssp_results.items()}

		if splitted:
			chain_ids = list(parsed.keys())
			result = parsed[chain_ids[0]]
			for c in chain_ids[1:]:
				if parsed[c].chains():
					chain_id = list(parsed[c].chains())[0]
					result.add_chain(chain_id, parsed[c].chain(chain_id))
			return result
		else:
			return parsed['all']

	def _rm_tmp_pdb_files(self, plain_results, rm_tmp, splitted):
		if rm_tmp and splitted:
			files_to_remove = [f for _, (f, _) in plain_results.items()]
			for f in files_to_remove:
				if exists(f):
					remove(f)

	def run(self, pdb_file, pdb_id, splitted=True,  chain_ids=[], save_to_file=True, parse=True, rm_tmp=True, use_cleaner_fixer=False):
		"""
		Runs dssp on a PDB file. PDB can be used as a whole or can be splitted in each different
		chain and computed separetally. Results from PDB can be optionally parsed into a
		DsspDataPdb object.
			:param self:
			:param pdb_file: The path of the inut pdb file.
			:param pdb_id: The pdb id of the input file.
			:param splitted=True: If true, dssp is run for each chain separetally.
			:param save_to_file=True: Save dssp results to file (or many files if splitted).
			:param parse=True: Parse dssp output into a DsspDataPdb object.
			:param rm_tmp=True: Remove intermideate PDB files generated when splitted option is True.

		Returns a tuple containg the dssp results and a dict with the path of saved files.

		Dssp results are plain text or a DsspDataPdb object, if output is plain text, the result is
		a dict where chain ids are keys when option splitted is True or 'all' is the key when
		splitted is False.
		"""
		pdbs_to_run = self._split_pdb(pdb_file, chain_ids, use_cleaner_fixer) if splitted else {'all': pdb_file}
		plain_results = {c:(pdb, self._run(pdb)) for c, pdb in pdbs_to_run.items()}

		saved_files = self._save_dssp(save_to_file, plain_results)
		result = self._parse_dssp_results(plain_results, pdb_id, splitted) if parse else {c: r for c,(_, r) in plain_results.items()}
		self._rm_tmp_pdb_files(plain_results, rm_tmp, splitted)

		return result, saved_files

	def _run(self, pdb_file):
		#
		# Updated to fail by running other version of dssp in case a different dssp is installed. 
		# This is not a clean way to do this
		#
		try:
			# -i is invalid option from dssp v.4.4
			# result = check_output([self.dssp, '-i', pdb_file])
			result = check_output([self.dssp, pdb_file, '--output-format', 'dssp'])
		except CalledProcessError:
			try:
				logger.error("Error running DSSP")
				result = check_output([self.dssp, '-i', pdb_file])
			except CalledProcessError:
				result = ""

		return result

	def _adapt_pdb_file(self, original_file_name, chain):
		m = splitext(original_file_name)
		case = "Upper" if chain.isupper() else "Lower"
		return "{}.{}_{}.pdb".format(m[0], chain.upper(), case)

	def _split_pdb(self, pdb_file, chain_ids, use_cleaner_fixer):
		pdb_header = "HEADER    DNA BINDING PROTEIN                     20-JUL-17   XXXX"
		parser = PDBParser()
		structure = parser.get_structure("XXXX", pdb_file)
		first_model = structure[0]
		result = {}
		for c in first_model.get_chains():
			if len(chain_ids) == 0 or c.id in chain_ids:
				io = PDBIO()
				io.set_structure(structure)
				new_pdb_filename = self._adapt_pdb_file(pdb_file, c.id)
				if not use_cleaner_fixer or not pdbfixer_loaded:
					io.save(new_pdb_filename, ChainSelect(c.id))
				else:
					io.save(new_pdb_filename, CleanerSelect(c.id))
					#fixer = PDBFixer(filename=new_pdb_filename)
					#fixer.missingResidues = {}
					#fixer.findMissingAtoms()
					#fixer.addMissingAtoms()
					#PDBFile.writeFile(fixer.topology, fixer.positions, open(new_pdb_filename, 'w'), keepIds=True)

				## add header so it looks like PDB format
				cmd = "sed -i '1i " + pdb_header + "' " + new_pdb_filename
				out = _run_cmd(cmd)

				result.update({c.id: new_pdb_filename})
		return result


class ChainSelect(Select):
	"""
	Custom Bio.PDB.PDBIO.Select class to select a single chain
		:param Select: parent Bio.PDB.PDBIO.Select
	"""
	def __init__(self, chain):
		"""
		Creates a new ChainSelect class, that selects a single givenc chain
			:param self:
			:param chain: a character string representing a PDB chain.
		"""
		self.accepted = chain

	def accept_model(self, model):
		"""
		Accepts only the first model of a PDB
			:param self:
			:param model:
		"""
		return 1 if model.id == 0 else 0

	def accept_chain(self, chain):
		"""
		Accepts a single chain
			:param self:
			:param chain:
		"""
		return 1 if chain.id == self.accepted else 0


class CleanerSelect(Select):
	"""
	Custom Bio.PDB.PDBIO.Select class to select a single cleaned chain
		:param Select: parent Bio.PDB.PDBIO.Select
	"""
	def __init__(self, chain):
		"""
		Creates a new CleanerSelect class, that selects a single given chain
			:param self:
			:param chain: a character string representing a PDB chain.
		"""
		self.accepted = chain

	def accept_model(self, model):
		"""
		Accepts only the first model of a PDB
			:param self:
			:param model:
		"""
		return 1 if model.id == 0 else 0

	def flag_residue_for_deletion(self, residue):
		for atom in residue.child_list:
			if atom.is_disordered():
				atom = atom.disordered_get()
				atom.altloc = "@"
			else:
				print('Error:', chain.id, key)

	def accept_chain(self, chain):
		"""
		Accepts a single chain, and flag disordered residues for deletion (to solve issues in 1AW8, 2HAL, 4AON, 4Z0Y, 6RXH)
			:param self:
			:param chain:
		"""

		disordered_dict = {}
		for res in chain.get_residues():
			if res.is_disordered():
				resseq = str(res.get_id()[1])
				icode = str(res.get_id()[2].strip())
				key = resseq+icode
				if key in disordered_dict:
					if disordered_dict[key].get_resname() != res.get_resname():
						if not is_aa(disordered_dict[key].get_resname(), standard=True) and is_aa(res.get_resname(), standard=True):
							self.flag_residue_for_deletion(disordered_dict[key])
						else:
							self.flag_residue_for_deletion(res)
				else:
					disordered_dict[key] = res

		return 1 if chain.id == self.accepted else 0

	def accept_residue(self, residue):
		"""
		Accepts all residues, but flag the disordered atoms which should be kept (to solve dropping atoms ex: 1F7A, 3FMA) (Not rely on pdbfixer in doing so because it is failed in 3FMA)
			:param self:
			:param residue:
		"""
		for atom in residue.child_list:
			if atom.is_disordered():
				atom = atom.disordered_get()
				if atom.altloc != "@":
					atom.altloc = " "

		'''
		# CONVERT SELENOMETHIONINES TO METHIONINES
		if residue.resname == 'MSE' or residue.resname == 'MET':

			residue.resname = 'MET'

			for atom in residue.child_list:
                    		# IF DISORDERED ATOM GET THE ONE WHICH WILL BE KEPT
				if atom.is_disordered():
					atom = atom.disordered_get()


				if atom.name == 'SE' and atom.element == 'SE':
					atom.name = 'SD'
					atom.element = 'S'
		'''
		return 1

	def accept_atom(self, atom):
		"""
		Reject disordered atoms which are not flagged
			:param self:
			:param atom:
		"""
		if atom.is_disordered() and not atom.altloc.isspace():
			return 0
		return 1
