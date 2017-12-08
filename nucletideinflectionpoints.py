"""
Script to run mean AT conent inflection point directionality alignment

Wren Saylor
December 2017

Other: Analog and discrete signals, Fourier series, Spectral analysis, Fourier transform, Discrete Fourier Transform, Nyquist-Shannon sampling theorem, aliasing

Copyright 2017 Harvard University, Wu Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

import argparse
import pandas as pd
from collections import defaultdict
import pybedtools as pbt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import math
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from scipy.stats import chisquare
from scipy.interpolate import splrep, splev
import scipy.stats as ss
from scipy.stats import mstats
import time


# set command line arguments
def get_args():
	# file lists
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=str,help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("-r","--randomfile",type=argparse.FileType('rU'),help="A file containing a list of paths to the random regions equable with your elements to plot in contrast") # option to not use random at all, or generate randoms?

	# columns in element file - all 0 based
	parser.add_argument("-lc", "--labelcolumn",type=int,help='column in the element file where the label (exonic, intergenic, intronic) is')
	parser.add_argument("-ic", "--idcolumn",type=int,help='column in the element file where the id is, if not provided will be generated for the sliding window collection')

	# genome files
	parser.add_argument("-g","--genome",type=str, default="hg19.genome")
	parser.add_argument("-f","--fasta",type=str,default="hg19.fa")

	# integer parameters
	parser.add_argument("-t","--total",type=int,default="600",help='total size of region to look at (region + flanks), should be an even number, suggested to be at least triple your element')
	parser.add_argument("-e","--element",type=int,help='size of your element (region without flanks), should be an even number, if not provided will use the smallest size of your input data')
	parser.add_argument("-i","--inset",type=int,default="50",help='size into your element from the boundaries, should be an even number, will also be used for number bp out side of you element - to get a boundary for finding inflection points')
	parser.add_argument("-w","--window",type=int,default="11",help='size of sliding window, should be an odd number, previous studies have used 11')

	# plot filename addition
	parser.add_argument('-s',"--stringname",type=str,help='string to add to the outfile name')

	# directionality parameters
	parser.add_argument('-c',"--reversecomplement",action='store_true',help='if you want the reverse complement to be plotted')
	parser.add_argument("-n", "--numberrandomassignments",type=int,help='the number of times to to randomly assign direction, will only be used when "--reversecomplement" is "random"')

	# Add additional descriptive file name information
	return parser.parse_args()

# set all args that will be used throughout the script
def set_global_variables(args):
	# Integer parameters
	global num
	global elementsize
	global inuce
	global window
	global halfwindow
	global fillX
	num = args.total
	elementsize = args.element
	inuce = args.inset
	window = args.window
	halfwindow = ((window/2)+1)
	fillX = range(0,(num-window))

	# Element, random regions and methylation files
	global eFiles
	global rFiles
	eFiles = args.efile
	if args.randomfile:
		rFiles = [line.strip() for line in args.randomfile]
	else:
		rFiles = None
	
	# Column labels in element file
	global labelcolumn
	global idcolumn
	labelcolumn = args.labelcolumn
	idcolumn = args.idcolumn

	# Genome files from UCSC
	global sizeGenome
	global faGenome
	sizeGenome = args.genome
	faGenome = args.fasta

	# Lists with the types and directions to use
	global nucList
	nucList = ['A','T']
	
	# A string to add to the out file name in case you want to set up runs and let be
	global stringName
	stringName = args.stringname
	
	global reverseComplement
	global randomassignments
	reverseComplement = args.reversecomplement
	randomassignments = args.numberrandomassignments
	
	print 'collected global parameters'

def set_ploting_parameters():
	# Locations for plotting with sliding window
	global plotLineLocationOne # upstream element boundary
	global plotLineLocationTwo # downstream element boundary
	global plotLineLocationThree # upstream element inset
	global plotLineLocationFour # downstream element inset
	plotLineLocationOne = (((num-uce)/2)+(inuce-halfwindow))
	plotLineLocationTwo = (((num-uce)/2)+(uce-inuce-halfwindow))
	plotLineLocationThree = (((num-uce)/2)-halfwindow)
	plotLineLocationFour = (((num-uce)/2)+uce-halfwindow)

	global centerelementpoint
	centerelementpoint = ((num-uce)/2)+(uce/2)-halfwindow
	print 'center point', centerelementpoint

	print 'set plotting parameters'

# get bt features
def get_bedtools_features(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

# get the correct range for fang evaluation
def collect_coordinates_for_element_boundaries(btFeatures):
	midFeatures = pd.read_table(btFeatures.fn, header=None)
	midFeatures['middle'] = midFeatures.loc[:,1:2].mean(axis=1).astype(int)
	midFeatures['size'] = midFeatures.loc[:,2].astype(int)-midFeatures.loc[:,1].astype(int)
	midFeatures.columns.values[0]='chr'
	midFeatures.columns.values[1]='start'
	midFeatures.columns.values[2]='end'
	global uce
	if elementsize:
		uce = elementsize
	else:
		getmin = midFeatures['size'].min()
		if (getmin % 2) == 0:
			uce = getmin
		else:
			uce = getmin + 1
	flankSize = (num - uce)/2
	inregion = uce-(inuce*2)
	midFeatures['sInset'] = midFeatures['start'].astype(int) + inuce
	midFeatures['eInset'] = midFeatures['end'].astype(int) - inuce
	midFeatures['sOutset'] = midFeatures['start'].astype(int) - inuce
	midFeatures['eOutset'] = midFeatures['end'].astype(int) + inuce
	midFeatures=midFeatures.drop(['middle'],axis=1)
	if labelcolumn:
		midFeatures['type'] = midFeatures.loc[:,labelcolumn]
	if idcolumn:
		midFeatures['id'] = midFeatures.loc[:,idcolumn]
	else:
		midFeatures.insert(len(midFeatures.columns),'id',range(0,0+len(midFeatures)))
	print 'collected the coordinates to create the sequence strings'
	return midFeatures

# check that the coordinates with the surrounding regions don't fall off the end of genome, if they do, print and skip
# future might just fill in with N's
def check_coords_beyond_genome(rangeFeatures):
	genomeBedtool = get_bedtools_features(sizeGenome)
	genomeFeatures = pd.read_table(genomeBedtool.fn, header=None)
	genomeFeatures['chr'] = genomeFeatures.loc[:,0]
	genomeFeatures['end'] = genomeFeatures.loc[:,1]
	genomeFeatures['start'] = 0
	initiallength=len(rangeFeatures.index)
	chrList = rangeFeatures['chr'].unique()
	bychromosome = []
	for chr in chrList:
		end = genomeFeatures.loc[(genomeFeatures['chr'] == chr),'end'].values[0]
		belowthresh = rangeFeatures[(rangeFeatures['chr'] == chr) & (rangeFeatures['eOutset'] <= end)]
		bychromosome.append(belowthresh)
		# can check if any starts are negative
	catFeatures = pd.concat(bychromosome,axis=0)
	checklength = initiallength - len(catFeatures.index)
	print "there were {0} out of {1} total beyond the end of the genome".format(checklength, initiallength)
	return catFeatures

# get the strings for sliding window regions
def get_fasta_for_element_boundary_coordinates(rangeFeatures):
	rangeFeatures = check_coords_beyond_genome(rangeFeatures)
	rangeFeatures['sSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','sOutset','sInset']].values.astype(str).tolist()))
	rangeFeatures['sSeq'] = rangeFeatures['sSeq'].str.upper()
	rangeFeatures['eSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','eInset','eOutset']].values.astype(str).tolist()))
	rangeFeatures['eSeq'] = rangeFeatures['eSeq'].str.upper()
	print 'collected the sequence strings for running the sliding window'
	return rangeFeatures

# get the reverse complement
def reverse_complement_dictionary(sequence):
	seqDict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
	return "".join([seqDict[base] for base in reversed(sequence)])

# used in get_fasta_for_element_coordinates to extract just the fasta strings
def get_just_fasta_sequence_for_feature(inFeature):
	seqFeature = inFeature.sequence(fi=faGenome)
	outFeature = pd.read_table(seqFeature.seqfn)
	outSequence = outFeature[::2]
	outSequence = outSequence.reset_index(drop=True)
	return outSequence

# get coordinates with flanking regions
def collect_element_coordinates(fileName):
	btFeatures = get_bedtools_features(fileName)
	subsetFeatures = collect_coordinates_for_element_boundaries(btFeatures)
	rangeFeatures = get_fasta_for_element_boundary_coordinates(subsetFeatures)
	return rangeFeatures

# run the sliding window for a string in a column
def run_sliding_window(feature):
	outCollect = []
	n = inset*2
	s = 1
	start, end = 0, window
	while end < n:
		current = feature[start:end]
		percentage = float(100*current.count(key)/len(current))
		outCollect.append(percentage)
		start, end = start + s, end + s
	return outCollect

# sliding window as list in column
def collect_sliding_window_values(rangeFeatures):
	rangeFeatures['sWindow'] = rangeFeatures.apply(lambda row: (run_sliding_window(row['sSeq'])),axis=1)
	rangeFeatures['eWindow'] = rangeFeatures.apply(lambda row: (run_sliding_window(row['eSeq'])),axis=1)
	return rangeFeatures

# def boundary_inflection_point_locations(rangeFeatures):


# With the results from compare_boundaries_size_n per each element, evaluate directionality into new column
def assign_directionality_from_arg_or_boundary(rangeFeatures,fileName):
	rangeFeatures = collect_sliding_window_values(rangeFeatures)
	print rangeFeatures

	return rangeFeatures



def main():
	starttime = time.time()
	
	# Get and set parameters
	args = get_args()
	set_global_variables(args)
	paramlabels = '{0}_{1}_{2}_{3}_{4}_{5}'.format(elementsize,inuce,num,window,eFiles,stringName)
	
	# Get coords and strings for elements
	rangeFeatures = collect_element_coordinates(eFiles)
	directionFeatures = assign_directionality_from_arg_or_boundary(rangeFeatures,eFiles)
	
	# If label column is supplied to split df by type
	if labelcolumn:
		typeList = directionFeatures['type'].unique()
		for type in typeList:
			typeBool,typeWindow,typeNames = separate_dataframe_by_group(type,directionFeatures,'type',eFiles)
			percentage_at_for_element(typeBool['combineString'],'{0}, {1}'.format(eFiles,type))
			probOptionstype = make_probabilites_for_direction(typeBool,'directionality')
			spreadRandomtype,spreadRandomtypeRC,denseRandomtype,denseRandomtypeRC,lengthrandom=[],[],[],[],[]
			
			# If random files are supplied
			if rFiles:
				lengthrandom.append(len(rFiles))
				lengthrandom.append('randomfiles')
				for randomFile in rFiles:
					randomFeatures = collect_element_coordinates(randomFile)
					randirFeatures= assign_directionality_from_arg_or_boundary(randomFeatures,randomFile)
					rantypeBool,rantypeWindow,rantypeNames = separate_dataframe_by_group(type,randirFeatures,'type',randomFile)
					spreadRandomtype.append(rantypeWindow)
					if reverseComplement:
						typeWindowRCrandom,typeNamesRCrandom = sort_elements_by_directionality(rantypeBool,'directionality')
						spreadRandomtypeRC.append(typeWindowRCrandom)
			
			# If random assignments of elements used
			if randomassignments:
				lengthrandom.append(randomassignments)
				lengthrandom.append('randomassingments')
				for i in range(randomassignments):
					typeBool['randomDirectiontype'] = np.random.choice(dirOptions,len(typeBool.index),p=probOptionstype)
					typedirWindow, typedirNames = sort_elements_by_directionality(typeBool,'randomDirectiontype')
					spreadRandomtype.append(typedirWindow)
					spreadRandomtypeRC.append(typedirWindow)
			
			# Collect average of either type of randoms
			if any([rFiles,randomassignments]):
				denseRandomtype = sliding_window_df_to_collect_all_random(spreadRandomtype,typeNames)
			if reverseComplement:
				typeWindowRC,typeNamesRC = sort_elements_by_directionality(typeBool,'directionality')
				denseRandomtypeRC = sliding_window_df_to_collect_all_random(spreadRandomtypeRC,typeNames)
				graph_element_line_means_with_rc_sorted(typeWindow,typeNames,typeWindowRC,'{0}_rc_{1}_{2}'.format(type,paramlabels,lengthrandom),spreadRandomtype,spreadRandomtypeRC,denseRandomtype,denseRandomtypeRC)
			else:
				graph_element_line_means(typeWindow,typeNames,'{0}_{1}_{2}'.format(type,paramlabels,lengthrandom),spreadRandomtype,denseRandomtype)
	
	# If all elements are to run together
	else:
		allWindow,allNames = sliding_window_wrapper(directionFeatures['combineString'],directionFeatures['id'])
		spreadRandom,spreadRandomRC,denseRandom,denseRandomRC,lengthrandom=[],[],[],[],[]
		
		# If random files are supplied
		if rFiles:
			lengthrandom.append(len(rFiles))
			lengthrandom.append('randomfiles')
			for randomFile in rFiles:
				randomFeatures = collect_element_coordinates(randomFile)
				randirFeatures= assign_directionality_from_arg_or_boundary(randomFeatures,randomFile)
				sWindow,sNames = sliding_window_wrapper(randirFeatures['combineString'],randirFeatures['id'])
				spreadRandom.append(sWindow)
				if reverseComplement:
					ranrevWindow, ranrevNames = sort_elements_by_directionality(randirFeatures,'directionality')
					spreadRandomRC.append(ranrevWindow)
		
		# If random assignments of elements used
		if randomassignments:
			lengthrandom.append(randomassignments)
			lengthrandom.append('randomassingments')
			for i in range(randomassignments):
				directionFeatures['randomDirection'] = np.random.choice(dirOptions,len(directionFeatures.index),p=probOptions)
				randirWindow, randirNames = sort_elements_by_directionality(directionFeatures,'randomDirection')
				spreadRandom.append(randirWindow)
				spreadRandomRC.append(randirWindow)
		
		# Collect average of either type of randoms
		if any([rFiles,randomassignments]):
			denseRandom = sliding_window_df_to_collect_all_random(spreadRandom,allNames)
		if reverseComplement:
			revWindow,revNames = sort_elements_by_directionality(directionFeatures,'directionality')
			denseRandomRC = sliding_window_df_to_collect_all_random(spreadRandomRC,allNames)
			graph_element_line_means_with_rc_sorted(allWindow,allNames,revWindow,'all_rc_{0}_{1}'.format(paramlabels,lengthrandom),spreadRandom,spreadRandomRC,denseRandom,denseRandomRC)
		else:
			graph_element_line_means(allWindow,allNames,'all_{0}_{1}'.format(paramlabels,lengthrandom),spreadRandom,denseRandom)
	
	endtime = time.time()
	print 'total time elapsed is {0}'.format(endtime-starttime)

if __name__ == "__main__":
	main()