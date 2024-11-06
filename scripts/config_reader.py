import os,inspect

try:
    import configparser
except:
    import ConfigParser as configparser

#-----
import logging
logger = logging.getLogger(__name__)
#-----

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

def load_configeration_options(sections=[],verbose=False):
	options = {}
	config = configparser.SafeConfigParser()
	
	if os.path.exists('/home/scripts/config.ini'):
		config_path = '/home/scripts/config.ini'
	elif os.path.exists(os.path.join(os.path.dirname(inspect.stack()[0][1]),'./config_local.ini')):
		config_path = os.path.abspath(os.path.join(os.path.dirname(inspect.stack()[0][1]),'./config_local.ini'))
	else:
		config_path = os.path.abspath(os.path.join(os.path.dirname(inspect.stack()[0][1]),'./config.ini'))
	
	loading_issues = False
	if os.path.exists(config_path):
		config.read(config_path)
		
		for each_section in config.sections():
			if len(sections) == 0 or each_section in sections:
				for (each_key, each_val) in config.items(each_section):
					each_val = each_val.strip()
					if each_val == "None" or each_val == "none":
						options[each_key] = None
					elif each_val == "True" or each_val == "true":
						options[each_key] = True
					elif each_val == "False" or each_val == "false":
						options[each_key] = False
					elif isfloat(each_val):
						if str(each_val).count(".") > 0 or str(each_val).count("e") > 0:
							try:
								options[each_key] =  float(each_val)
							except:
								options[each_key] =  each_val
						else:
							options[each_key] = int(each_val)
					else:
						options[each_key] = each_val

		for option in options:
			if option[-4:] == 'path':
				if not os.path.exists(options[option]):
					if verbose:
						logger.debug("In config.ini - " + option + "=" + options[option] + "- path issue - path does not exist.")
						loading_issues = True
	else:
		logger.debug("Config File Path not found:" + config_path)

	if loading_issues:
		logger.debug("Config loading issues" + config_path)
	
	return options
