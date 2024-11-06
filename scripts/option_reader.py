import sys,os,inspect

from optparse import (OptionParser,BadOptionError,AmbiguousOptionError)

file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.abspath(os.path.join(file_path,"./utilities")))

import utilities_error

#-----
import logging
logger = logging.getLogger(__name__)
#-----

def isInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def isFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class MyOptionParser(OptionParser):
	def error(self, msg):
		pass#logger.debug("option error:" + str(msg))
		
def load_commandline_options(options,options_help,purge_commandline=False,show_errors=False,verbose=False,debug=False):

	## TODO: This still has some issues with dropping options when a unknown option is passed

	commandline_options = {}

	op = MyOptionParser()
	
	if show_errors:
		op.verbose = show_errors

	options_sorted = list(options.keys())
	options_sorted.sort()

	for option in options_sorted:
		help_text = "--Not set--"

		if option in options_help:
			help_text = str(options_help[option])

		op.add_option("--" + option,
			action="store",
			dest=option,
			default=options[option],
			help=help_text)
	
	opts, args = op.parse_args()
	
	try:
		for option in vars(opts):
			if debug:
				print(option,getattr(opts, option),args)

			try:
				if getattr(opts, option) == None:
					commandline_options[option] = None
				elif getattr(opts, option) == "":
					commandline_options[option] = ""
				elif isinstance(getattr(opts, option), (dict)):
					commandline_options[option] = getattr(opts, option)
				elif isinstance(getattr(opts, option), (list)):
					commandline_options[option] = getattr(opts, option)
				elif (isFloat(getattr(opts, option)) or isInt(getattr(opts, option))) and not isinstance(getattr(opts, option), bool):
					try:
						if str(getattr(opts, option)).count(".") > 0 or str(getattr(opts, option)).count("e") > 0:
							commandline_options[option] = float(getattr(opts, option))
						else:
							commandline_options[option] = int(getattr(opts, option))
					except:
						commandline_options[option] = getattr(opts, option)
				elif getattr(opts, option) == "True" or getattr(opts, option) == True:
					commandline_options[option] = True
				elif getattr(opts, option) == "False" or getattr(opts, option) == False:
					commandline_options[option] = False
				else:
					commandline_options[option] = getattr(opts, option)
			except:
				logger.debug((option,getattr(opts, option)))
				
	except:
		if show_errors:
			utilities_error.printError()
	
	if purge_commandline:
		sys.argv = sys.argv[:1]

	return commandline_options
