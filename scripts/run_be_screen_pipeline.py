import os,sys,inspect,logging

#-----

file_path = os.path.dirname(inspect.stack()[0][1])
sys.path.append(os.path.join(file_path,"./base_editing_screen_tools/"))
import baseEditingScreenParser

#-----

if __name__ == "__main__":
	baseEditingParserObj = baseEditingScreenParser.baseEditingParser()
	baseEditingParserObj.process_screen()