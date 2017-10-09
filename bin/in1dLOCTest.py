#!/usr/bin/python

import numpy

# like in1d but returns the index in B of each element of A rather than a binary vector

def in1dLOC(A, B):
	
	# concatenate the vectors together
	ABConcat = numpy.concatenate((B, A))
	ABConcatSortedIDX = numpy.argsort(ABConcat)
	
	print ABConcat
	print ABConcatSortedIDX
	ABConcatSorted = ABConcat[ABConcatSortedIDX]
	print ABConcatSorted

	# the locations of the elements of B are those whose indices are less than or equal to B.size
	BIDXInSorted = numpy.where(ABConcatSortedIDX < B.size)[0]
	
	print BIDXInSorted
	
	# make a mask that is 1 for the locations of the B elements
	BIDXInSortedMask = numpy.zeros_like(ABConcatSortedIDX)
	BIDXInSortedMask[BIDXInSorted] = 1

	# perform a cumulative sum so that BIDXInSortedCS[I] is the index we are looking for that element in B
	
	BIDXInSortedCS = numpy.cumsum(BIDXInSortedMask) - 1
	print BIDXInSortedMask
	print BIDXInSortedCS
	InvalidIDX = BIDXInSortedCS < 0

	BIDXInSortedCS[InvalidIDX] = 0
	
	# now check for equality between ABConcat BIDXInSortedCS
	AIDXInSortedMask = (ABConcatSortedIDX >= B.size)
	print AIDXInSortedMask	
	AIDXEqualToB = (ABConcatSorted == B[BIDXInSortedCS])
	#print ABConcatSorted[AIDXInSortedMask]
	#print B[BIDXInSortedCS]
	#print AIDXEqualToB
	

A = numpy.array([1, 4, 1, 2, 3, 4])
B = numpy.array([2, 3])

C = in1dLOC(A, B)
