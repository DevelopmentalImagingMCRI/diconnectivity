#!/usr/bin/python3

import numpy
import numpy.linalg

import os
import re
import struct

# converts the MIF datatypes
# Bit	bitwise data
# Int8	signed 8-bit (char) integer
# UInt8	unsigned 8-bit (char) integer
# Int16	signed 16-bit (short) integer
# UInt16	unsigned 16-bit (short) integer
# Int16LE	signed 16-bit (short) integer (little-endian)
# UInt16LE	unsigned 16-bit (short) integer (little-endian)
# Int16BE	signed 16-bit (short) integer (big-endian)
# UInt16BE	unsigned 16-bit (short) integer (big-endian)
# Int32	signed 32-bit int
# UInt32	unsigned 32-bit int
# Int32LE	signed 32-bit int (little-endian)
# UInt32LE	unsigned 32-bit int (little-endian)
# Int32BE	signed 32-bit int (big-endian)
# UInt32BE	unsigned 32-bit int (big-endian)
# Float32	32-bit floating-point
# Float32LE	32-bit floating-point (little-endian)
# Float32BE	32-bit floating-point (big-endian)
# Float64	64-bit (double) floating-point
# Float64LE	64-bit (double) floating-point (little-endian)
# Float64BE	64-bit (double) floating-point (big-endian)
#
# into fields suitable for struct.pack and struct.unpack
# returns a tuple,
#	the first element is the endianness,
#		'>' for explicit big-endian
#		'<' for explicit big-endian
#		'' for native
#	the second element is the data type
#		c 	char 	string of length 1 	1 	 
#		b 	signed char 	integer 	1 	(3)
#		B 	unsigned char 	integer 	1 	(3)
#		h 	short 	integer 	2 	(3)
#		H 	unsigned short 	integer 	2 	(3)
#		i 	int 	integer 	4 	(3)
#		I 	unsigned int 	integer 	4 	(3)
#		l 	long 	integer 	4 	(3)
#		L 	unsigned long 	integer 	4 	(3)
#		q 	long long 	integer 	8 	(2), (3)
#		Q 	unsigned long long 	integer 	8 	(2), (3)
#		f 	float 	float 	4 	(4)
#		d 	double

def dataTypeGetStructStrings(datatype):
	pat = re.compile('(U)?((Int|Float)(\d+))(LE|BE)?')
	mat = pat.match(datatype)
	
	if mat == None:
		return None
	else:
		T = {'Int8': 'b', 'Int16': 'h', 'Int32': 'i', 'Float32': 'f', 'Float64': 'd'}
		
		if mat.group(2) in T:
			dataletter = T[mat.group(2)]
			if mat.group(1) != None:
				dataletter = dataletter.upper()
		size = int(int(mat.group(4)) / 8)
		if mat.group(5) == None:
			endian = ''
		else:
			E = {'LE': '<', 'BE': '>'}
			endian = E[mat.group(5)]
		#print datatype + ": " + dataletter + " " + endian + " " + str(size)
		return (dataletter, endian, size)

def structStringsGetDataType(letter, endian):

	T = {'b': 'Int8', 'h': 'Int16', 'i': 'Int32', 'f': 'Float32', 'd': 'Float64'}
	
	if letter == letter.upper():
		signedPart = 'U'
	else:
		signedPart = ''
	
	datatypePart = T[letter.lower()]
	
	if len(endian) == 1:
		E = {'<': 'LE', '>': 'BE'}
		endianPart = E[endian]
	else:
		endianPart = ''
	#print signedPart + datatypePart + endianPart
	return signedPart + datatypePart + endianPart

#letter, endian, size = dataTypeGetStructStrings('Int8'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('UInt8'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Int16'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('UInt16'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Int16LE'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('UInt16LE'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Int16BE'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('UInt16BE'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Int32'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('UInt32'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Int32LE'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('UInt32LE'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Int32BE'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('UInt32BE'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Float32'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Float32LE'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Float32BE'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Float64'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Float64LE'); structStringsGetDataType(letter, endian);
#letter, endian, size = dataTypeGetStructStrings('Float64BE'); structStringsGetDataType(letter, endian);
#
def loadTrackFile(fileName, origVectorMode = False):

# reads an mrtrix .tck file, returns the header and tracks from the the file "fileName"
#	(dict): with keys
# 	"init_threshold": the lowest FA to initialise tracks
# 	"lmax": the spherical harmonic order
# 	"max_dist": maximum length of any tracts
# 	"max_num_attempts": maximum number of tracts sent out
# 	"max_num_tracks": maximum number of accepted tracts
# 	"max_trials": number of probabilistic trials at each point
# 	"method": tractography method
# 	"min_curv": ???
# 	"min_dist": minimum track distance
# 	"no_mask_interp": whether we DONT interpolate the masks
# 	"sh_precomputed": Legendre polynomials precomputed?
# 	"source": tensor or CSD image used for tracking
# 	"step_size": pixel step size
# 	"stop_when_included": whether we stopped when we went into an include zone
# 	"threshold": FA threshold
# 	"unidirectional": whether we sent out in one direction
# 	"roi":
#		"type (list of strings)": seed or mask
#		"file (list of strings)": filenames
#	"datatype (dict)":
#		"letter": string that has the pack/unpack letter for the datatype 
#		"endian": '<' for little or '>' for big or '' for native
#		"size": the size, in bytes, of each element
# 	"count": number of accepted tracks
# 	"total_count": number of tracks sent out
#	"tracks (list of numpy arrays)": the track data
	
	#print(fileName)
	if not os.path.isfile(fileName):
		print("warning: could not find track file: " + fileName)
		return None
	else:
		FID = open(fileName, 'rb');
		# first line should be "mrtrix image"
		firstLine = FID.readline().decode().rstrip()
		#print(firstLine)

		if firstLine != 'mrtrix tracks':
			print("file " + fileName + " not mrtrix track file")
			FID.close()
			return None
		else:
			del firstLine
			

			trackStruct = dict()

			trackStruct['init_threshold'] = None
			trackStruct['lmax'] = None
			trackStruct['max_dist'] = None
			trackStruct['max_num_attempts'] = None
			trackStruct['max_num_tracks'] = None
			trackStruct['max_trials'] = None
			trackStruct['method'] = None
			trackStruct['min_curv'] = None
			trackStruct['min_dist'] = None
			trackStruct['no_mask_interp'] = None
			trackStruct['sh_precomputed'] = None
			trackStruct['source'] = None
			trackStruct['step_size'] = None
			trackStruct['stop_when_included'] = None
			trackStruct['threshold'] = None
			trackStruct['unidirectional'] = None
			trackStruct['roi'] = dict()
			trackStruct['roi']['type'] = list();
			trackStruct['roi']['file'] = list();
			trackStruct['datatype'] = None
			trackStruct['count'] = None
			trackStruct['total_count'] = None
			trackStruct['tracks'] = list()
			pat = re.compile('^([a-z_]+): (.*)$')

			while True:
				curLine = FID.readline().decode().rstrip()
				if curLine == "END":
					break
				else:
					mat = pat.match(curLine)
					if mat != None:
						
						curKeyword = mat.group(1)
						curValue = mat.group(2)

						if curKeyword in ['init_threshold', 'max_dist', 'min_curv', 'min_dist', 'step_size', 'threshold']:
							trackStruct[curKeyword] = float(curValue)
						elif curKeyword in ['lmax', 'max_num_attempts', 'max_num_tracks', 'max_trials', 'no_mask_interp', 'sh_precomputed', 'threshold', 'unidirectional', 'count', 'total_count']:
							trackStruct[curKeyword] = int(curValue)
						elif curKeyword in ['method', 'source']:
							trackStruct[curKeyword] = curValue[:]
						elif curKeyword == 'datatype':
							trackStruct['datatype'] = dict()
							trackStruct['datatype']['letter'], trackStruct['datatype']['endian'], trackStruct['datatype']['size'] = dataTypeGetStructStrings(curValue)
						elif curKeyword == 'roi':
							roimat = re.match('(\S+)\s+(\S+)', curValue)
							if roimat != None:
								trackStruct['roi']['type'].append(roimat.group(1))
								trackStruct['roi']['file'].append(roimat.group(2))
							del roimat
						elif curKeyword == 'file':
							trackStruct['file'] = curValue[:]
							filemat = re.match('\.\s+(\d+)', curValue)
							if filemat != None:
								trackStruct['dataoffset'] = int(filemat.group(1))
							del filemat
		#if firstLine != 'mrtrix tracks':
		#print(trackStruct)
		if trackStruct['datatype'] != None:
			FID.seek(0, 2) # 2 means the end of the file
			fileSize = FID.tell()
			
			if trackStruct['dataoffset'] >= 0 and trackStruct['dataoffset'] < fileSize:
				FID.seek(trackStruct['dataoffset'], 0) # 0 means absolute position
				
				numElements = (fileSize - trackStruct['dataoffset']) / trackStruct['datatype']['size']

				if numElements % 3 != 0:
					print("The data section of " + fileName + " does not have a multiple of 3 elements")
				else:
					trackStruct['tracks'] = list()
					
					## original fully vector reading of the file, read ALL data at once
					if origVectorMode:
						T = FID.read()
						
						FID.close()
						
						T = struct.unpack(trackStruct['datatype']['endian'] + trackStruct['datatype']['letter'] * numElements, T)
						T = numpy.reshape(numpy.array(T), (3, numElements / 3), order='F')

						separatorIDX = numpy.where(numpy.any(numpy.logical_not(numpy.isfinite(T)), axis = 0))[0]
						#print(separatorIDX)
						leftIDX = 0
						
						for z in range(separatorIDX.size - 1):
							trackStruct['tracks'].append(numpy.array(T[:, leftIDX:separatorIDX[z]]))
							leftIDX = separatorIDX[z] + 1
						## end original fully vector mode
					else:
						# block reading mode
						numTriplesInFile = int(numElements / 3)

						# how many triples to read in each time
						tripleBlockSize = int(5000)
						
						numTriplesRead = int(0)
						
						# curTrack is used to store tracks that straddle block boundaries
						curTrack = list()
						#allSeparatorIDX = list()
						while numTriplesRead < numTriplesInFile:

							B = FID.read(tripleBlockSize * 3 * trackStruct['datatype']['size'])
							curNumTriplesRead = int(len(B) / 3 / trackStruct['datatype']['size'])
							#print "curNumTriplesRead: " + str(curNumTriplesRead)
							
							T = struct.unpack(trackStruct['datatype']['endian'] + trackStruct['datatype']['letter'] * (curNumTriplesRead * 3), B)
							del B
							T = numpy.reshape(numpy.array(T), (3, curNumTriplesRead), order='F')
							
							#print T
							separatorIDX = numpy.where(numpy.any(numpy.logical_not(numpy.isfinite(T)), axis = 0))[0]
							#print separatorIDX	
							#allSeparatorIDX.append(separatorIDX + numTriplesRead)
							numTriplesRead += curNumTriplesRead
							
							if numpy.size(separatorIDX) == 0: # no separators in this set, just append the whole thing to the current track
								curTrack.append(numpy.array(T))
							else:
								if separatorIDX[0] == 0 and len(curTrack) > 0: # the separator is at the start, push the current track onto the main list
									trackStruct['tracks'].append(numpy.concatenate(curTrack, axis = 1))
									#print "separatoridx[0] == 0, adding track"
								elif separatorIDX[0] > 0: # otherwise append and then add the track
									curTrack.append(T[:, :separatorIDX[0]])
									trackStruct['tracks'].append(numpy.concatenate(curTrack, axis = 1))
								
								curTrack = list()
								#print "separatoridx[0] > 0, concatenating and adding track"
								
								if numpy.size(separatorIDX) > 1:
									# add all the intermediate tracks
									for z in range(numpy.size(separatorIDX) - 1):
										if separatorIDX[z + 1] > separatorIDX[z] + 1:
											#print "adding intermediate track: " + str(separatorIDX[z]) + ":" + str(separatorIDX[z + 1])
											trackStruct['tracks'].append(T[:, (separatorIDX[z] + 1):separatorIDX[z + 1]])
								
								if separatorIDX[-1] < T.shape[1] - 1:
									curTrack.append(T[:, (separatorIDX[-1] + 1):])
									#print "adding last bit to curTrack"

							del T
							#print str(len(trackStruct['tracks']))
						#while numTriplesRead < numTriplesInFile:
						#allSeparatorIDX = numpy.concatenate(allSeparatorIDX)
						#print allSeparatorIDX
						#trackStruct['tracks']

						FID.close()
					numTracks = len(trackStruct['tracks'])
					#print "numTracks = " + str(numTracks)
					if numTracks != trackStruct['count']:
						print("Warning: the number of tracks in the file " + fileName + ", which was " + str(numTracks) + ", does not match the number specified by the count header field, which was " + str(trackStruct['count']))
					#print trackStruct['tracks']
					#print T[:, 0:10]
					#print T.shape
				#print numElements
		return trackStruct
	#if not os.path.isfile(fileName):
#def loadTrackFile(fileName):

def saveTrackFile(trackStruct, fileName):
	
# DESCRIPTION
#	Saves a mrtrix track file, the track format used by mrtrix.
#
# PARAMETERS
#	trackStruct (dict)
#	MANDATORY FIELDS:
#	"datatype" (dict):
#		"letter": string that has the pack/unpack letter for the datatype 
#		"endian": '<' for little or '>' for big or '' for native
#		"size": the size, in bytes, of each element
#	"count": number of accepted tracks
#	"tracks": the track data
#	OPTIONAL FIELDS:
# 	"init_threshold": the lowest FA to initialise tracks
# 	"lmax": the spherical harmonic order
# 	"max_dist": maximum length of any tracts
# 	"max_num_attempts": maximum number of tracts sent out
# 	"max_num_tracks": maximum number of accepted tracts
# 	"max_trials": number of probabilistic trials at each point
# 	"method": tractography method
# 	"min_curv": ???
# 	"min_dist": minimum track distance
# 	"no_mask_interp": whether we DONT interpolate the masks
# 	"sh_precomputed": Legendre polynomials precomputed?
# 	"source": tensor or CSD image used for tracking
# 	"step_size": pixel step size
# 	"stop_when_included": whether we stopped when we went into an include zone
# 	"threshold": FA threshold
# 	"unidirectional": whether we sent out in one direction
# 	"total_count": number of tracks sent out

# check the must have fields in trackStruct
	mustHaveFields = ['datatype', 'tracks']

	for z in range(len(mustHaveFields)):
		if not mustHaveFields[z] in list(trackStruct.keys()):
			print("The field: " + mustHaveFields[z] + " was not in the structure")
			return
	
	if not (trackStruct['datatype']['letter'] == 'f' or trackStruct['datatype']['letter'] == 'd'):
		print("The tracks must have a floating point datatype")
		return

	try:
		FID = open(fileName, 'wb')
	except e:
		print("Could not open " + fileName + " for writing")
		return
	
	FID.write("mrtrix tracks\n".encode())
	
	for curField in sorted(trackStruct.keys()):
		if curField == 'roi':
			for z in range(len(trackStruct[curField]['type'])):
				FID.write(('roi: %s %s\n' % (trackStruct[curField]['type'][z], trackStruct[curField]['file'][z])).encode())
		elif curField == 'datatype':
			FID.write(('datatype: %s\n' % structStringsGetDataType(trackStruct[curField]['letter'], trackStruct[curField]['endian'])).encode())
		elif curField != 'tracks':
			FID.write(("%s: %s\n" % (curField, str(trackStruct[curField]))).encode())
	
	if not 'count' in trackStruct:
		FID.write(("count: %d\n" % (len(trackStruct['tracks']))).encode())
	FID.write("file: ".encode())
	offset = FID.tell()
	offset = offset + 14
	FID.write((". %d\nEND\n" % offset).encode())
	curPos = FID.tell()
	
	bytesToWrite = offset - curPos

	if bytesToWrite > 0:
		FID.write(('\x00' * bytesToWrite).encode())
	
	if trackStruct['datatype']['letter'] == 'f':
		curPadding = numpy.single(numpy.array([numpy.nan, numpy.nan, numpy.nan]))
	else:
		curPadding = numpy.double(numpy.array([numpy.nan, numpy.nan, numpy.nan]))
	
	curPadding.tolist()
	for z in range(len(trackStruct['tracks'])):
		
		if trackStruct['datatype']['letter'] == 'f':
			curTrack = numpy.single(trackStruct['tracks'][z]).flatten(order='F')
		else:
			curTrack = numpy.double(trackStruct['tracks'][z]).flatten(order='F')
		
		# serial version
		#for k in curTrack.flat:
		#	FID.write(struct.pack(trackStruct['datatype']['endian'] + trackStruct['datatype']['letter'], k))
		#for k in curPadding.flat:
		#	FID.write(struct.pack(trackStruct['datatype']['endian'] + trackStruct['datatype']['letter'], k))

		# "vectorised" version
		curTrack = curTrack.tolist()
		#print(trackStruct['datatype']['endian'] + trackStruct['datatype']['letter'] * len(curTrack))
		#print(trackStruct['datatype']['endian'] + trackStruct['datatype']['letter'], curTrack[0])
		FID.write(struct.pack(trackStruct['datatype']['endian'] + trackStruct['datatype']['letter'] * len(curTrack), *curTrack))
		FID.write(struct.pack(trackStruct['datatype']['endian'] + trackStruct['datatype']['letter'] * len(curPadding), *curPadding))
		del curTrack

	#FID.write(struct.pack(trackStruct['datatype']['endian'] + trackStruct['datatype']['letter'] * len(curPadding), *curPadding))

	#for k in curPadding.flat:
	#	FID.write(struct.pack(trackStruct['datatype']['endian'] + trackStruct['datatype']['letter'], k))
			
	FID.close()

def getMrtrixTransformFromNIBabelHeader(NIIHeader):
	
	QForm, QFormCode = NIIHeader.get_qform(coded = True)
	SForm, SFormCode = NIIHeader.get_sform(coded = True)
	
	#print QForm
	#print QFormCode
	#print SForm
	#print SFormCode
	pixDims = NIIHeader.get_header()['pixdim'][1:4]

	if QFormCode <= 0 or SFormCode > 0:
		transform = numpy.array(SForm)
		voxelSizes = numpy.sqrt(numpy.sum(transform[0:3, 0:3] * transform[0:3, 0:3], axis = 0))
		transform[0:3, 0:3] = transform[0:3, 0:3] / numpy.atleast_2d(voxelSizes)
	else:
		# use the qcode, just divide by the voxel sizes
		transform = numpy.array(QForm)
		transform[0:3, 0:3] = transform[0:3, 0:3] / numpy.atleast_2d(pixDims)
		pass
	
	#transform[0, 0] = -transform[0, 0]
	#transform[0:2, :] = numpy.concatenate((numpy.atleast_2d(transform[1, :]), numpy.atleast_2d(transform[0, :])), axis = 0)
	
	permutation = numpy.argmax(numpy.abs(transform[0:3, 0:3] * transform[0:3, 0:3]), axis = 1)
	
	flip = transform[(numpy.arange(3), permutation)]
	flip = (flip < 0)
	#print permutation
	
	if not numpy.array_equal(permutation, numpy.array([0, 1, 2])) or numpy.any(flip):
		origTransform = numpy.array(transform)
		transform[0:3, 0:3] = transform.take(permutation, axis = 1).take([0, 1, 2], axis = 0)
		for curColumn in range(3):
			dimLength = pixDims[curColumn] * (NIIHeader.shape[curColumn] - 1)

			if flip[curColumn] == True:
				transform[0:3, curColumn] = -transform[0:3, curColumn]
				for curRow in range(3):
					transform[curRow, 3] = transform[curRow, 3] + origTransform[curRow, permutation[curColumn]] * dimLength
	#if not numpy.array_equal(permutation, numpy.array([0, 1, 2])) or numpy.any(flip):
	return transform

# NIIHeader *MUST* be direct output from nibabel.load

def tracksWorldToIMG(tracksWorldSpace, NIIHeader):
	tracksIMGSpace = list()
	transform = getMrtrixTransformFromNIBabelHeader(NIIHeader)
	
	print(transform)
	invTransform = numpy.matrix(numpy.linalg.inv(transform))

	invTransform = invTransform[0:3, :]
	
	tracksSZ = numpy.zeros((len(tracksWorldSpace)), dtype = numpy.int64)
	for z in range(len(tracksWorldSpace)):
		tracksSZ[z] = tracksWorldSpace[z].shape[1]
	
	tracksT = numpy.concatenate(tracksWorldSpace, axis = 1)

	tracksT = numpy.array(invTransform * numpy.matrix(numpy.concatenate((tracksT, numpy.ones((1, tracksT.shape[1]))), axis = 0)))
	tracksT = tracksT / numpy.atleast_2d(NIIHeader.get_header()['pixdim'][1:4]).T
	tracksT[1] = (NIIHeader.shape[1] - 1) - tracksT[1]
		
	leftIDX = 0
	tracksIMGSpace = list()
	for z in range(len(tracksWorldSpace)):
		tracksIMGSpace.append(tracksT[:, leftIDX:(leftIDX + tracksSZ[z])])
		leftIDX += tracksSZ[z]

#	for curTrack in tracksWorldSpace:
#		
#		# apply the transformation
#		T = invTransform * numpy.matrix(numpy.concatenate((curTrack, numpy.ones((1, curTrack.shape[1]))), axis = 0))
#
#		# remove the voxel spacing
#		T = T / numpy.atleast_2d(NIIHeader.get_header()['pixdim'][1:4]).T
#		
#		T[1, :] = (NIIHeader.shape[1] - 1) - T[1, :]
#		if numpy.any(numpy.logical_not(numpy.isfinite(T))):
#			print T
#		#print T
#		#quit()
#		tracksIMGSpace.append(numpy.array(T))
#		del T
	return tracksIMGSpace

if __name__ == "__main__":
	trackFile = 'bigtracks_culled'
	#MIF = loadTrackFile(trackFile + '.tck')
	MIF = loadTrackFile(trackFile + '.tck', origVectorMode = False)
	#print(len(MIF['tracks']))
	saveTrackFile(MIF, trackFile + "_saved.tck")
