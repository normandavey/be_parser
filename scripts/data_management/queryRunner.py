import sys

if sys.version_info[0] != 3 or sys.version_info[1] < 5:
    print("This script requires Python version 3.5 or greater")
    sys.exit(1)

import os
import inspect
file_path = os.path.dirname(inspect.stack()[0][1])

from queryManager import queryManager

class queryRunner:

	##------------------------------------------------------------------##

	def __init__(self,dataset_type, task, options={}, accession=None, pdb_id =None, pmid = None, pfam_id = None):
		self.queryManagerObj = queryManager()
		self.queryManagerObj.options['dataset_type'] = dataset_type
		self.queryManagerObj.options['task'] = task
		self.queryManagerObj.options['is_superuser'] = False
		self.queryManagerObj.options.update(options)

		if accession != None:	
			self.queryManagerObj.options['accession'] = accession

		if pdb_id != None:
			self.queryManagerObj.options['pdb_id'] = pdb_id
			
		if pmid != None: 
			self.queryManagerObj.options['pmid'] = pmid
			
		if pfam_id != None:
			self.queryManagerObj.options['pfam_id'] = pfam_id
	
	##------------------------------------------------------------------##

	def run(self):
		self.data = self.queryManagerObj.main()
		return self.data
	
	##------------------------------------------------------------------##
	
