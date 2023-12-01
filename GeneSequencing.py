#!/usr/bin/python3
import math

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		if banded:
			scor, align1, align2 = self.gene_banded(seq1,seq2)
			align1 = align1[:100]
			align2 = align2[:100]
		else:
			scor, align1, align2 = self.gene_algo(seq1,seq2)
			align1 = align1[:100]
			align2 = align2[:100]
			x = 0

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		score = scor
		alignment1 = align1.format(len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = align2.format(len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}


	def gene_algo(self,seq1,seq2):
		x = (min(len(seq1),self.MaxCharactersToAlign))+1
		y = (min(len(seq2),self.MaxCharactersToAlign))+1
		# initialize matrices of scores and back pointers
		E = [[0] * (y) for _ in range(x)]
		prev = [[0]*(y)for _ in range(x)]
		returnAl1 = ""
		returnAl2 = ""
		# initialize top row
		for i in range(x):
			E[i][0] = INDEL*i
			prev[i][0] = 'U'
		# 	initialize left column
		for j in range(y):
			E[0][j] = INDEL*j
			prev[0][j] = 'L'
		# fill in scores and back pointers until bottom corner reached
		for i in range(1,x):
			for j in range(1,y):
				if seq1[i-1] == seq2[j-1]:
					match = E[i-1][j-1] + MATCH
				else:
					match = E[i-1][j-1] + SUB
				insert = E[i][j-1] + INDEL
				delete = E[i-1][j] + INDEL
				point = min(match,insert,delete)
				E[i][j] = point
				if point == insert:
					prev[i][j] = 'L'
				elif point == delete:
					prev[i][j] = 'U'
				else:
					prev[i][j] = 'D'
		score = E[x-1][y-1]
		# recurse backwards until you've reached the top left corner
		i = x-1
		j = y-1
		while not(i == 0 and j == 0):
			next = prev[i][j]
			if next == 'L':
				returnAl1 = '-' + returnAl1
				returnAl2 = seq2[j-1] + returnAl2
				j-=1
			elif next == 'U':
				returnAl1 = seq1[i-1] + returnAl1
				returnAl2 = '-' + returnAl2
				i-=1
			else:
				returnAl1 = seq1[i - 1] + returnAl1
				returnAl2 = seq2[j - 1] + returnAl2
				j-=1
				i-=1
		# return score and the two alignments
		return score, returnAl1, returnAl2


	def gene_banded(self, seq1, seq2):
		length = len(seq1)-len(seq2)
		absolute = math.fabs(length)
		if absolute > 1000:
			return math.inf, "No Alignment Possible", "No Alignment Possible"

		n = max(len(seq1),len(seq2))
		n = min(n,self.MaxCharactersToAlign)+1
		k = 7
		d = 3
		returnAl1 = ""
		returnAl2 = ""
		score = math.inf
		#initialize k*n matrices
		E = [[math.inf] * (k) for _ in range(n)]
		prev = [['X'] * (k) for _ in range(n)]
		test  = 0
		# initialize top and left
		for i in range(d+1):
			E[i][0] = i * INDEL
			prev[i][0] = 'U'
		for j in range(0,d+1):
			E[0][j] = j * INDEL
			prev[0][j] = 'L'
		# for the first d arrays, use this formula to fill
		for i in range(1,4):
			for j in range(1,4+i):
				if seq1[i - 1] == seq2[j - 1]:
					match = E[i - 1][j - 1] + MATCH
				else:
					match = E[i - 1][j - 1] + SUB
				insert = E[i][j - 1] + INDEL
				delete = E[i - 1][j] + INDEL
				point = min(match, insert, delete)
				E[i][j] = point
				if point == insert:
					prev[i][j] = 'L'
				elif point == delete:
					prev[i][j] = 'U'
				else:
					prev[i][j] = 'D'
		# for the middle part from d to almost the end use this formula to populate
		i = 4
		gap = 1
		endReached = False
		while((7+gap) <= n) and not endReached:
			for j in range(7):
				if j == 0:
					#if [i == len(seq1) or (j+gap == len(seq2))]
					if seq1[i - 1] == seq2[j+(gap) - 1]:
						match = E[i-1][j] + MATCH
					else:
						match = E[i-1][j] + SUB
					insert = E[i-1][j+1] + INDEL
					point = min(match, insert)
					E[i][j] = point
					if point == insert:
						prev[i][j] = 'L'
					else:
						prev[i][j] = 'D'
				elif j == 6:
					if seq1[i - 1] == seq2[j+(gap) - 1]:
						match = E[i-1][j] + MATCH
					else:
						match = E[i-1][j] + SUB
					delete = E[i][j-1] + INDEL
					point = min(match, delete)
					E[i][j] = point
					if point == delete:
						prev[i][j] = 'U'
					else:
						prev[i][j] = 'D'
				else:
					if seq1[i - 1] == seq2[j+(gap) - 1]:
						match = E[i-1][j] + MATCH
					else:
						match = E[i-1][j] + SUB
					insert = E[i-1][j + 1] + INDEL
					delete = E[i][j-1] + INDEL
					point = min(match, insert, delete)
					E[i][j] = point
					if point == insert:
						prev[i][j] = 'L'
					elif point == delete:
						prev[i][j] = 'U'
					else:
						prev[i][j] = 'D'
			if (i) == min(n-1,len(seq1)):
				endReached = True
				score = E[i][6]
			i += 1
			gap += 1
		x = 0
		# populate the end of the matrix, the last few arrays
		while (x < 6) and not endReached:
			gap  = n - (6)
			i = n-(3-x)
			for j in range((x+1),7):
				if seq1[i - 1] == seq2[j+gap - 2]:
					match = E[i - 1][j - 1] + MATCH
				else:
					match = E[i - 1][j - 1] + SUB
				insert = E[i][j - 1] + INDEL
				delete = E[i - 1][j] + INDEL
				point = min(match, insert, delete)
				E[i][j] = point
				if point == insert:
					prev[i][j] = 'L'
				elif point == delete:
					prev[i][j] = 'U'
				else:
					prev[i][j] = 'D'
			if (i) == min(n-1,len(seq1)):
				endReached = True
			score = E[i][6]
			x +=1
		test = 3

		# begin recursing backwards from bottom corner to top in sections
		# recursing section 1
		i = (min(n-1,len(seq1)))
		j = k - 1
		while E[i][0] == math.inf:
			next = prev[i][j]
			if next == 'L':
				returnAl1 = '-' + returnAl1
				returnAl2 = seq2[j+gap -2] + returnAl2
				j -= 1
			elif next == 'U':
				returnAl1 = seq1[i - 1] + returnAl1
				returnAl2 = '-' + returnAl2
				i -= 1
			elif next == 'D':
				returnAl1 = seq1[i - 1] + returnAl1
				returnAl2 = seq2[j+gap -2] + returnAl2
				j -= 1
				i -= 1
		# 		recursing the majority of the matrix
		while E[i-1][6] != math.inf:
			next = prev[i][j]
			if next == 'U':
				returnAl1 = '-' + returnAl1
				returnAl2 = seq2[j+gap -2] + returnAl2
				j -= 1
			elif next == 'L':
				returnAl1 = seq1[i - 1] + returnAl1
				returnAl2 = '-' + returnAl2
				j+=1
				gap -=1
				i -=1
			elif next == 'D':
				returnAl1 = seq1[i - 1] + returnAl1
				returnAl2 = seq2[j+gap -2] + returnAl2
				i -= 1
				gap-=1
		# 		recurse the end to the top left corner
		while not(i == 0 and j == 0):
			next = prev[i][j]
			if next == 'L':
				returnAl1 = '-' + returnAl1
				returnAl2 = seq2[j-1] + returnAl2
				j -= 1
			elif next == 'U':
				returnAl1 = seq1[i - 1] + returnAl1
				returnAl2 = '-' + returnAl2
				i -= 1
			elif next == 'D':
				returnAl1 = seq1[i - 1] + returnAl1
				returnAl2 = seq2[j-1] + returnAl2
				j -= 1
				i -= 1
		# return score and optimal alignments
		return score, returnAl1,returnAl2