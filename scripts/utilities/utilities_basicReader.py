import os,sys,csv
import utilities_basic as basic

def readTableHeader(tablePath,delimiter="\t",hasHeader=True,key="",byColumn=True,relationship="binary",stripQuotes=False):
	if os.path.exists(tablePath):
		rows = open(tablePath).read().replace("\n","\r").strip("\n\r").split("\r")

		if hasHeader:
			start_row = 1
			headers = rows[0].split(delimiter)
		else:
			start_row = 0
			headers = list(range(0,len(rows[0].split(delimiter))))

	return headers

def readTableFile(tablePath,delimiter="\t",hasHeader=True,key="",byColumn=True,relationship="binary",stripQuotes=False,verbose=False):
	if os.path.exists(tablePath):
		rows = open(tablePath).read()

		if rows[0] == "<":
			return {}

		if stripQuotes:
			rows = rows.replace('"',"")

		#rows = rows.replace("\r","\n").strip("\n").strip("\r")
		rows = rows.splitlines()

		start_row = 0
		while (rows[start_row][0:1] == "#"):
			start_row += 1

		if hasHeader:
			headers = [header.strip() for header in rows[start_row].split(delimiter)]
			start_row += 1
		else:
			headers = list(range(0,len(rows[0].split(delimiter))))

		tableData = {}

		if key == "":
			if byColumn == True:
				for header in headers:
					tableData[header] = []
		else:
			if key not in headers:

				print(("Key not in header",tablePath))
				return {}

			keyIndex = headers.index(key)

		########################################################################
		with open(tablePath) as csvfile:
			csv_reader = csv.reader(csvfile, delimiter=delimiter, quotechar='"')
			line_counter = 0 
			for bits in csv_reader:
				if line_counter >= start_row:
					try:
						if len(bits) > 0:
							#bits = row.split(delimiter)
							if verbose: print((len(bits),len(headers),bits))

							if byColumn == True:
								if key == "":
									for i in range(0,len(headers)):
										tableData[headers[i]].append(bits[i].strip())
								else:
									tableData[bits[keyIndex].strip()] = {}
									for i in range(0,len(headers)):
										if key != headers[i]:
											tableData[bits[keyIndex].strip()][headers[i]] = bits[i].strip()
							else:
								if key == "":
									tableData[len(tableData)] = {}
									for i in range(0,len(headers)):
										tableData[len(tableData) -1][headers[i]] = bits[i].strip()
								else:
									if relationship == "binary":
										tableData[bits[keyIndex].strip()] = {}
									else:
										if bits[keyIndex].strip() not in tableData:
											tableData[bits[keyIndex].strip()] = []

										tableDataTmp = {}

									for i in range(0,len(headers)):
										if key != headers[i]:
											if relationship == "binary":
												tableData[bits[keyIndex].strip()][headers[i]] = bits[i].strip()
											else:
												tableDataTmp[headers[i]] = bits[i].strip()

									if relationship != "binary":
										tableData[bits[keyIndex].strip()].append(tableDataTmp)
						
							#if bits[keyIndex].strip() == "P04637":
							#	import pprint
							#	pprint.pprint(tableData[bits[keyIndex].strip()])
							#	pprint.pprint(bits)
							#	sys.exit

					except Exception as e:
						print(e)
						print(("Skipping",row))
						return {}

				line_counter += 1

		return tableData
	else:
		print(("File not found:" + tablePath))


if __name__ == "__main__":
	table = sys.argv[1]

	readTableFile(table,key='elminstanceid')

	tableData = readTableFile(table)

	print((len(basic.removeRedundency(tableData['uniprotid']))))
	print((basic.removeRedundency(tableData['uniprotid'])))
