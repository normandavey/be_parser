import sys,traceback,pprint

###################################################################
### Errors
###################################################################

def getError(description=None, email_format=False):
	try:
		exc_type, exc_value, exc_tb = sys.exc_info()
		exc_file = "-"
		exc_line = "-"

		format_exc = traceback.format_exc().strip().split("\n")

		if len(format_exc) > 0:
			exc_file = format_exc[-3].strip()

		if len(format_exc) > 1:
			exc_line = format_exc[-2].strip()


		if email_format:
			msg = []
			msg.append("Some error occured. Details below.")

			error_msg = {
				"Description": description,
				"Value":str(exc_value),
				"Exc_file":exc_file,
				"Exc_line":exc_line,
				"Format_exc":traceback.format_exc()
			}

			for k, v in error_msg.items():
				msg.append("<span style='padding-left: 2px;'><b>"+k+"</b>" + ": " + str(v) + "</span>")

			return "<br>".join(msg)

		return {
			"type":str(exc_type),
			"value":str(exc_value),
			"exc_file":exc_file,
			"exc_line":exc_line,
			"format_exc":traceback.format_exc(),
			"description":description
			}
	except:
		return {"type":"Unknown error"}

def printError(description=None):
	pprint.pprint(getError(description=description))
