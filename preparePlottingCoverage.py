#!/usr/bin/env python3

import sys

def printRanges (chrName, start, end, dist, step, svN, chrs, backwards = 0, inside = 0):
	start, end, dist, step = int(start), int(end), int(dist), int(step)
	if inside:
		backwards = 0
		dist = end - start # This needs to be done before changing end
		end = start

	if backwards:
		for i in range(dist, 0, -step):
			if start - i < 0:
				continue
			yield f"{chrName}\t{start - i}\t{start - i + step}\tSV_{svN}-{'{:,}'.format(i)}" #svN is set outside function
	else:
		for i in range(0, dist, step):
			if end + i + step > chrs[chrName]:
				break
			if inside and start + i + step > start + dist: #if is going to be out of range, we get a smaller interval, and finish
				yield f"{chrName}\t{end + i}\t{start + dist}\tSV_{svN}_in_{'{:,}'.format(i + step)}"
				break
			yield f"{chrName}\t{end + i}\t{end + i + step}\t" + (f"SV_{svN}+{'{:,}'.format(i + step)}" if not inside else \
				f"SV_{svN}_in_{'{:,}'.format(i + step)}")

def writeIndexFile (refGenome, inputBed, idxFile):
	chrs = {}
	
	with open(refGenome) as refGenome:
		for line in refGenome.readlines():
			chrName, start, end = line.rstrip().split("\t")[0:3]
			chrs[chrName] = int(end) #Start is always 0 so no point in storing it

	idxContent = []
	with open(inputBed) as f, open(idxFile, "w") as idx:
		svN = 0

		for line in f.readlines():
			line = line.rstrip()
			lineFields = line.split("\t")

			chrName, svStart, svEnd = lineFields[0:3]
			if chrName not in chrs.keys():
				continue
			svType = ""
			if (len(lineFields) >= 4):
				svType = lineFields[3]

			svCode = "" # Will be similar to "SV_1_DEL"
			if not svType:
				svCode = f"SV_{svN}"
			else:
				svCode  = f"SV_{str(svN) + '_' + svType}"

			# For index file
			svFullRange = f"{chrName}\t{int(svStart) - 1000000}\t{int(svEnd) + 1000000}\t{svCode}\n"
			idxContent.append(svFullRange)
			svN += 1
		idx.write("".join(idxContent))

		return idxContent

def writeMosdepthInput (refGenome, inputBed, workingDir, sample, infix):
	"""
	Creates a bed file as side effect and yields a list of str
	refGenome: genome.bed
	"""
	chrs = {} # To check we dont go out of range
	with open(refGenome) as refGenome:
		for line in refGenome.readlines():
			chrName, start, end = line.rstrip().split("\t")[0:3]
			chrs[chrName] = int(end) #Start is always 0 so no point in storing it

	idxContent = []
	with open(inputBed) as f:
		svN = 0

		for line in f.readlines():
			line = line.rstrip()
			lineFields = line.split("\t")

			chrName, svStart, svEnd = lineFields[0:3]
			if chrName not in chrs.keys():
				continue
			svType = ""
			if (len(lineFields) >= 4):
				svType = lineFields[3]

			svCode = "" # Will be similar to "SV_1_DEL"
			if not svType:
				svCode = f"SV_{svN}"
			else:
				svCode  = f"SV_{str(svN) + '_' + svType}"

			# For index file
			svFullRange = f"{chrName}\t{int(svStart) - 1000000}\t{int(svEnd) + 1000000}\t{svCode}\n"
			idxContent.append(svFullRange)
			svN += 1

	def linesGenerator(): # So the function can have the side effect of witing the index file
		with open(inputBed) as f:
			svN = 0

			for chrom, end in chrs.items():
					yield f"{chrom}\t0\t{end}\t{chrom}\n"
			
			for line in f.readlines():
				line = line.rstrip()
				lineFields = line.split("\t")

				chrName, svStart, svEnd = lineFields[0:3]
				if chrName not in chrs.keys():
					continue
				svType = ""
				if (len(lineFields) >= 4):
					svType = lineFields[3]

				svCode = "" # Will be similar to "SV_1_DEL"
				if not svType:
					svCode = f"SV_{svN}"
				else:
					svCode  = f"SV_{str(svN) + '_' + svType}"

				# For index file
				svFullRange = f"{chrName}\t{int(svStart) - 1000000}\t{int(svEnd) + 1000000}\t{svCode}\n"
				idxContent.append(svFullRange)

				# For mosdepth

				# For the t-test
				
				if (not svType): 
					yield f"{line}\t{svCode}\n" # Avg. coverage of SV
				else:
					yield '\t'.join([chrName, svStart, svEnd]) + f"\t{svCode}\n" # Avg. coverage of SV

				# for x in printRanges(chrName, svStart, svEnd, 1e6, 1, svN, chrs, backwards = 1): # Before the SV
				# 	yield x + "\n"
				# for x in printRanges(chrName, svStart, svEnd, 0, 1, svN, chrs, inside = 1): # Inside the SV, dist is ignored, thus the 0
				# 	yield x + "\n"
				# for x in printRanges(chrName, svStart, svEnd, 1e6, 1, svN, chrs): # After the SV
				# 	yield x + "\n"

				# For plotting
				for x in printRanges(chrName, svStart, svEnd, 1e6, 1e4, svN, chrs, backwards = 1): # Before the SV
					yield x + "_plot\n"
				for x in printRanges(chrName, svStart, svEnd, 0, 1e4, svN, chrs, inside = 1): # Inside the SV, dist is ignored, thus the 0
					yield x + "_plot\n"
				for x in printRanges(chrName, svStart, svEnd, 1e6, 1e4, svN, chrs): # After the SV
					yield x + "_plot\n"

				svN += 1

	return linesGenerator()