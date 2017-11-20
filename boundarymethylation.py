"""
Script to run methylation analysis

Wren Saylor
September 19 2017

To Do:
unittest
check frequency calculation
change directionality arg column
streamline and remove excess columns
.index len to break up large datasets

May be graphing incorrect value; add frequency and add cpg count before dividing

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

# set command line arguments
def get_args():
	# File lists
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=str,help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("-r","--randomfile",type=argparse.FileType('rU'),help="A file containing a list of paths to the random regions equable with your elements to plot in contrast") # option to not use random at all, or generate randoms?
	parser.add_argument("-m","--methylationfile",type=argparse.FileType('rU'),help="A file containing a list of paths to the methlylation bedfiles")

	# Columns in element file - all 0 based
	parser.add_argument("-lc", "--labelcolumn",type=int,help='column in the element file where the label (exonic, intergenic, intronic) is')
	parser.add_argument("-dc", "--directionalitycolumn",type=int,help='column in the element file where directionality is, if not supplied will infer by AT content')
	parser.add_argument("-ic", "--idcolumn",type=int,help='column in the element file where the id is')

	# Genome Files
	parser.add_argument("-g","--genome",type=str, default="hg19.genome")
	parser.add_argument("-f","--fasta",type=str,default="hg19.fa")

	# Integer Parameters for element
	parser.add_argument("-p","--periphery",type=int,default="10",help='number bp from your boundary you want to include in the analysis')
	parser.add_argument("-b","--bin",type=int,default="100",help='size of bins used to compare element ends and determine directionality')
	parser.add_argument("-e","--element",type=int,help='size of your element (region without flanks), should be an even number, if not provided will use the smallest size of your input data')#,default="200"
	
	# Integer Parameters for methlation
	parser.add_argument("-mp", "--methylationthresholdpercentage", type=int, default="10", help='size to threshold % methylation data')
	parser.add_argument("-mc", "--methylationthresholdcoverage", type=int, default="10", help='size to threshold uncapped coverage of methylation data to send to % methylation, field often uses 10')

	# Plot parameters
	parser.add_argument('-s',"--stringname",type=str,help='string to add to the outfile name')
	parser.add_argument('-c',"--reversecomplement",action='store_true',help='if you want the reverse complement to be plotted')
	parser.add_argument("-n", "--numberrandomassignments",type=int,default="1",help='the number of times to to randomly assign direction, will only be used when "--reversecomplement" is "random"')
	parser.add_argument('-t',"--twographs",action='store_true',help='if you want to see the upstream and downstream boundaries separately')
	parser.add_argument('-d',"--motifdirectionality",type=str,help='if directionality should be assigned by motif presence')

	# Add additional descriptive file name information
	return parser.parse_args()

# set all args that will be used throughout the script
def set_global_variables(args):
	# Integer parameters
	global periphery
	global binDir
	global elementsize
	periphery = args.periphery
	binDir = args.bin
	elementsize = args.element

	# Methylation percentage and coverage thresholds
	global methPerThresh
	global methCovThresh
	methPerThresh = args.methylationthresholdpercentage
	methCovThresh = args.methylationthresholdcoverage

	# Element, random regions and methylation files
	global eFiles
	global rFiles
	global mFiles
	eFiles = args.efile
	mFiles = [line.strip() for line in args.methylationfile]
	if args.randomfile:
		rFiles = [line.strip() for line in args.randomfile]
	else:
		rFiles = None
	
	# Column labels in element file
	global labelcolumn
	global directionalitycolumn
	global idcolumn
	labelcolumn = args.labelcolumn
	directionalitycolumn = args.directionalitycolumn
	idcolumn = args.idcolumn

	# Genome files from UCSC
	global sizeGenome
	global faGenome
	sizeGenome = args.genome
	faGenome = args.fasta

	# A string to add to the out file name in case you want to set up runs and let be
	global stringName
	global splitgraphs
	stringName = args.stringname
	splitgraphs = args.twographs
	
	global reverseComplement
	global randomassignments
	global motifdirectionality
	reverseComplement = args.reversecomplement
	randomassignments = args.numberrandomassignments
	motifdirectionality = args.motifdirectionality

# get bt features
def get_bedtools_features(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

# get the correct range for fang evaluation
def collect_coordinates_for_element_positions(btFeatures):
	midFeatures = pd.read_table(btFeatures.fn, header=None)
	midFeatures['chr'] = midFeatures.loc[:,0]
	midFeatures['start'] = midFeatures.loc[:,1]
	midFeatures['end'] = midFeatures.loc[:,2]
	midFeatures['size'] = midFeatures.loc[:,2].astype(int)-midFeatures.loc[:,1].astype(int)
	global uce
	if elementsize:
		uce = elementsize
	else:
		getmin = midFeatures['size'].min()
		if (getmin % 2) == 0:
			uce = getmin
		else:
			uce = getmin + 1
	midFeatures['startup'] = midFeatures.loc[:,1] - periphery
	midFeatures['startdown'] = midFeatures.loc[:,1] + periphery
	midFeatures['endup'] = midFeatures.loc[:,2] - periphery
	midFeatures['enddown'] = midFeatures.loc[:,2] + periphery
	if idcolumn:
		midFeatures['id'] = midFeatures.loc[:,idcolumn]
	else:
		midFeatures.insert(0,'id',range(0,0+len(midFeatures)))
	if labelcolumn:
		midFeatures['type'] = midFeatures.loc[:,labelcolumn]
	if directionalitycolumn:
		midFeatures['directionality'] = midFeatures.loc[:,directionalitycolumn]
	return midFeatures

# get the strings for sliding window regions
def get_fasta_for_element_coordinates(rangeFeatures):#rangeFeatures,faGenome
	rangeFeatures = check_coords_beyond_genome(rangeFeatures)
	rangeFeatures['upstreamsequence'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','startup','startdown']].values.astype(str).tolist()))
	rangeFeatures['upstreamsequence'] = rangeFeatures['upstreamsequence'].str.upper()
	rangeFeatures['downstreamsequence'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','endup','enddown']].values.astype(str).tolist()))
	rangeFeatures['downstreamsequence'] = rangeFeatures['downstreamsequence'].str.upper()
	return rangeFeatures

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
		belowthresh = rangeFeatures[(rangeFeatures['chr'] == chr) & (rangeFeatures['enddown'] <= end)]
		bychromosome.append(belowthresh)
		# can check if any starts are negative
	catFeatures = pd.concat(bychromosome,axis=0)
	checklength = initiallength - len(catFeatures.index)
	print "there were {0} out of {1} total beyond the end of the genome".format(checklength,initiallength)
	return catFeatures

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
	subsetFeatures = collect_coordinates_for_element_positions(btFeatures)
	rangeFeatures = get_fasta_for_element_coordinates(subsetFeatures)
	return rangeFeatures

# convert bedtool to panda
# If there is nothing in the btobject, will it read the data from the previous itteration!?
def convert_bedtool_to_panda(btobject):
	save_bedtool_as_bedfile(btobject,'temp.bed')
	pdObject = pd.read_table(btobject.fn, header=None)
	return pdObject

# convert panda to bedtool
def convert_panda_to_bed_format(panda):
	arArFeatures = panda.values.tolist()
	btoutFeatures = get_bedtools_features(arArFeatures)
	return btoutFeatures

# save file from bedtool
def save_bedtool_as_bedfile(btObject,strFilename):
	btObject.saveas(strFilename)

# Do the comparison between boundaries to get + - or =
def calculate_nucleotides_at(element,size):
	start = element[:size]
	end = element[-size:]
	perSize = []
	perSize.append(eval('100*float(start.count("A") + start.count("a") + start.count("T") + start.count("t"))/len(start)'))
	perSize.append(eval('100*float(end.count("A") + end.count("a") + end.count("T") + end.count("t"))/len(end)'))
	return perSize

# Get the number of times a motif appears in each boundary
def locate_boundary_with_motif(element,size,motif):
	rcmotif = reverse_complement_dictionary(motif)
	start = element[:size]
	end = element[-size:]
	perSize = []
	perSize.append(eval('100*float(start.count(motif) + start.count(rcmotif))/len(start)'))
	perSize.append(eval('100*float(end.count(motif) + end.count(rcmotif))/len(end)'))
	return perSize

# Directionality, as inferred by comparing first and last n base pairs from input parameters
def compare_boundaries_size_n(element,size):
	if motifdirectionality:
		perSize = locate_boundary_with_motif(element,size,motifdirectionality)
	else:
		perSize = calculate_nucleotides_at(element,size)
	# give + - = depending on which side has larger AT content
	if perSize[0] > perSize[1]: outList = '+'
	if perSize[1] > perSize[0]: outList = '-'
	if perSize[1] == perSize[0]: outList = '='
	return outList

# With the results from compare_boundaries_size_n per each element, evaluate directionality into new column
def assign_directionality_from_arg_or_boundary(rangeFeatures,fileName):
	if not directionalitycolumn:
		rangeFeatures['feature'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','start','end']].values.astype(str).tolist()))
		rangeFeatures['feature'] = rangeFeatures['feature'].str.upper()
		rangeFeatures['directionality'] = rangeFeatures.apply(lambda row: (compare_boundaries_size_n(row['feature'],binDir)),axis=1)
		rangeFeatures=rangeFeatures.drop(['feature'],axis=1)
		numval = rangeFeatures['directionality'].value_counts()
	print '{0} in {1}'.format(numval,fileName)
	return rangeFeatures

# get the reverse complement
def reverse_complement_dictionary(sequence):
	seqDict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
	return "".join([seqDict[base] for base in reversed(sequence)])

# Concat the negative and positive groups together and get the updated frequency
def concat_positive_and_negative_directionality_and_get_frequency(posmethylation,newnegmethylation):
	# Concat pos and revised neg meth dfs
	frames = [posmethylation,newnegmethylation]
	methlationconcat = pd.concat(frames)
	submethylation = methlationconcat[['chr','id','methylationlocation','methylationpercentage','cytosine','tissue']]
	submethylation['methylationcount'] =  submethylation.groupby(['tissue'])['tissue'].transform('count')
	return submethylation

# Correct the Sequence of elements assigned negative directionality
def negative_directionality_corrected_features(negmethylation):
	# Zip reversed range to make a dictionary for replacing the location of the neg methylation
	originalRange = range(0,periphery*2)
	reverseRange = originalRange[::-1]
	rangeDict = dict(zip(originalRange,reverseRange))
	# Zip reverse complement sequence for replacing the nucleotides for neg methylation
	seqDict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
	negmethylation['newmethylationlocation'] = negmethylation.methylationlocation.map(rangeDict)
	negmethylation['newcytosine'] = negmethylation.cytosine.map(seqDict)
	newnegmethylation = negmethylation[['chr','id','newmethylationlocation','methylationpercentage','newcytosine','tissue']]
	newnegmethylation.columns = ['chr','id','methylationlocation','methylationpercentage','cytosine','tissue']
	return newnegmethylation

# Calculate sums for cpg and cg
def get_c_and_g_sequency_info(rangeFeatures):
	rangeFeatures['upcpgsequencecount'] = rangeFeatures.apply(lambda row: float(row['upstreamsequence'].count("CG")),axis=1)
	upcpgsequencecountsum = rangeFeatures['upcpgsequencecount'].sum()
	rangeFeatures['downcpgsequencecount'] = rangeFeatures.apply(lambda row: float(row['downstreamsequence'].count("CG")),axis=1)
	downcpgsequencecountsum = rangeFeatures['downcpgsequencecount'].sum()
# 	rangeFeatures['upcandgsequencecount'] = rangeFeatures.apply(lambda row: (row['upstreamsequence'].count("G")+row['upstreamsequence'].count("C")),axis=1)
# 	upcandgsequencecountsum = rangeFeatures['upcandgsequencecount'].sum()
# 	rangeFeatures['downcandgsequencecount'] = rangeFeatures.apply(lambda row: (row['downstreamsequence'].count("G")+row['downstreamsequence'].count("C")),axis=1)
# 	downcandgsequencecountsum = rangeFeatures['downcandgsequencecount'].sum()
	return upcpgsequencecountsum,downcpgsequencecountsum

# Add columns for c and g content
def columns_for_nucleotide_and_methylation_content(rangeFeatures,upstreamcverify,downstreamcverify):
	upcpgsequencecountsum,downcpgsequencecountsum = get_c_and_g_sequency_info(rangeFeatures)
	upstreamcverify['cpgsequencecountsum'] = upcpgsequencecountsum
	upstreamcverify['cpgsequencecountsumcombineboundary'] = upcpgsequencecountsum + downcpgsequencecountsum
	downstreamcverify['cpgsequencecountsum'] = downcpgsequencecountsum
	downstreamcverify['cpgsequencecountsumcombineboundary'] = upcpgsequencecountsum + downcpgsequencecountsum
# 	downstreamcverify['methylationcountcombineboundary'] = len(upstreamcverify.index) + len(downstreamcverify.index)
# 	downstreamcverify['percentageunmethylatedcpg'] = downstreamcverify['methylationcount']/downstreamcverify['cpgsequencecountsum']
# 	downstreamcverify['percentageunmethylatedcpgcombineboundary'] = downstreamcverify['methylationcountcombineboundary']/downstreamcverify['cpgsequencecountsumcombineboundary']
# 	upstreamcverify['methylationcountcombineboundary'] = len(upstreamcverify.index) + len(downstreamcverify.index)
# 	upstreamcverify['percentageunmethylatedcpg'] = upstreamcverify['methylationcount']/upstreamcverify['cpgsequencecountsum']
# 	upstreamcverify['percentageunmethylatedcpgcombineboundary'] = upstreamcverify['methylationcountcombineboundary']/upstreamcverify['cpgsequencecountsumcombineboundary']
	return upstreamcverify,downstreamcverify

def make_methylation_data_frame_and_verify_cytosine(methdf,name,features,column):
	if methdf is not None:
		methdf['tissue'] = name.replace('.bed','')
		subfeatures = features[['id',column]]
		merge = pd.merge(methdf,subfeatures,how='left',on='id')
		merge['methLocBEnd'] = merge['methLoc'] - 1
		merge['methLocCEnd'] = merge['methLoc'] + 1
		merge['methLocEnd'] = merge['methLoc'] + 2
		merge['Cytosine'] = merge.apply(lambda row: row[column][row['methLoc']:row['methLocCEnd']],axis=1)
		merge['Context'] = merge.apply(lambda row: row[column][row['methLocCEnd']:row['methLocEnd']],axis=1)
		merge['BackContext'] = merge.apply(lambda row: row[column][row['methLocBEnd']:row['methLoc']],axis=1)
		merge['ContextCheck'] = merge.apply(lambda row: row[column][row['methLocBEnd']:row['methLocEnd']],axis=1)
		merge['Nuc'] = merge['Nuc'].str.upper()
		merge['sameSeq'] = merge['Nuc'] == merge['ContextCheck']
		# If the nucleotide in the cytosine column is 'G', make the context the other 
		# direction (reverse complement later, in graphing, in order to differentiate 
		# between strands)
		merge.loc[merge['Cytosine'] == 'G', 'Context'] = merge['BackContext']
		# sameSeq might be 'False' if 1) the c is at the end border for the downstream 
		# boundary, 2) the sequence bridges the sequence split for the upstream boundary
		falseMeth = (merge[merge['sameSeq'] == False])
		# Conditionally update contexts where False for matches between sequence and 
		# methylation nucleotides- c get context, and g gets backcontext
		merge.loc[merge['sameSeq'] == False,'Cytosine'] = merge['NucCytosine']
		merge.loc[(merge['sameSeq'] == False) & (merge['NucCytosine'] == 'C'),'Context'] = merge['NucContext']
		merge.loc[(merge['sameSeq'] == False) & (merge['NucCytosine'] == 'G'),'Context'] = merge['NucBackContext']
		merge['methylationcount'] = len(merge.index)
		subMeth = merge[['chr','id','methLoc','methPer','Cytosine','methylationcount','tissue']]
		subMeth.columns = ['chr','id','methylationlocation','methylationpercentage','cytosine','methylationcount','tissue']
	else:
		subMeth = pd.DataFrame(np.nan,index=[0],columns=['chr','id','methylationlocation','methylationpercentage','cytosine','methylationcount'])
		subMeth['tissue'] = name.replace('.bed','')
	return subMeth

# Threshold methylation data by coverage and percentage
def threshold_methylation_data(methFeatures,methName):
	pdmethFeatures = convert_bedtool_to_panda(methFeatures)
	pdmethThresh = (pdmethFeatures[(pdmethFeatures.loc[:,3] >= methCovThresh) & (pdmethFeatures.loc[:,4] >= methPerThresh)])
	btmethThresh = convert_panda_to_bed_format(pdmethThresh)
	print 'methylation coverage is being thresholded at {0} and percentage at {1}'.format(methCovThresh, methPerThresh)
	initiallength = len(pdmethFeatures.index)
	checklength = initiallength - len(pdmethThresh.index)
	print 'there were {0} out of {1} total removed by thresholding {2}'.format(checklength,initiallength,methName)
	return btmethThresh

# Intersect regions from the methylation data with element regions
def intersect_methylation_data_by_element(rangeFeatures,methFeature,methName):
	initiallength = len(methFeature)
	methylationupstreamboundary = methFeature.intersect(rangeFeatures[['chr','startup','startdown','id']].values.astype(str).tolist(),wb=True,wa=True)
	checklengthup = len(methylationupstreamboundary)
	print 'there were {0} intersections out of {1} in {2}'.format(checklengthup,initiallength,methName)
	if len(methylationupstreamboundary) != 0:
		pdmethup=convert_bedtool_to_panda(methylationupstreamboundary)
		pdmethup['int']=0
		pdmethup.columns = ['mchr','mstart','mstop','methCov','methPer','chr','sBoundary','sEdge','id','int']
		pdmethup['Context'] = pdmethup['mstop'] + 1
		pdmethup['BackContext'] = pdmethup['mstart'] -1
		pdmethup['Nuc'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethup[['mchr','BackContext','Context']].values.astype(str).tolist()))
		pdmethup['NucContext'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethup[['mchr','mstop','Context']].values.astype(str).tolist()))
		pdmethup['NucCytosine'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethup[['mchr','mstart','mstop']].values.astype(str).tolist()))
		pdmethup['NucBackContext'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethup[['mchr','BackContext','mstart']].values.astype(str).tolist()))
		pdmethup['methLoc'] = pdmethup['int'].astype(int)+(pdmethup['mstart'].astype(int)-pdmethup['sBoundary'].astype(int))
		outupstreammethylation = pdmethup[['chr','mstart','mstop','sBoundary','sEdge','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']]
		outupstreammethylation.columns = ['chr','methStart','methStop','eleStart','eleStop','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']
		sortoutupstreammethylation = outupstreammethylation.sort_values(['methLoc'],ascending=True)
	else:
		sortoutupstreammethylation = None
	methylationdownstreamboudnary = methFeature.intersect(rangeFeatures[['chr','endup','enddown','id']].values.tolist(),wb=True,wa=True)
	checklengthdown = len(methylationdownstreamboudnary)
	print 'there were {0} intersections out of {1} in {2}'.format(checklengthdown,initiallength,methName)
	if len(methylationdownstreamboudnary) != 0:
		pdmethdown=convert_bedtool_to_panda(methylationdownstreamboudnary)
		pdmethdown['int']=0
		pdmethdown.columns=['mchr','mstart','mstop','methCov','methPer','chr','eEdge','eBoundary','id','int']
		pdmethdown['Context']=pdmethdown['mstop'] + 1
		pdmethdown['BackContext']=pdmethdown['mstart'] -1
		pdmethdown['Nuc']=get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethdown[['mchr','BackContext','Context']].values.astype(str).tolist()))
		pdmethdown['NucContext']=get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethdown[['mchr','mstop','Context']].values.astype(str).tolist()))
		pdmethdown['NucCytosine']=get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethdown[['mchr','mstart','mstop']].values.astype(str).tolist()))
		pdmethdown['NucBackContext']=get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethdown[['mchr','BackContext','mstart']].values.astype(str).tolist()))
		pdmethdown['methLoc'] = pdmethdown['int'].astype(int)+(pdmethdown['mstart'].astype(int)-pdmethdown['eEdge'].astype(int))
		outdownstreammethylation=pdmethdown[['chr','mstart','mstop','eEdge','eBoundary','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']]
		outdownstreammethylation.columns =['chr','methStart','methStop','eleStart','eleStop','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']
		sortoutdownstreammethylation=outdownstreammethylation.sort_values(['methLoc'],ascending=True)
	else:
		sortoutdownstreammethylation= None
	return sortoutupstreammethylation,sortoutdownstreammethylation

# Run the analysis to extract percentage, frequency, coverage, location, context, and direction
def collect_methylation_data_by_element(rangeFeatures):
	upstreamcapture, downstreamcapture = [], []
	for methName in mFiles:
		methFeatures=get_bedtools_features(methName)
		pdmethThresh=threshold_methylation_data(methFeatures,methName)
		upstreamposition,downstreamposition=intersect_methylation_data_by_element(rangeFeatures,pdmethThresh,methName)
		upstreamcverify=make_methylation_data_frame_and_verify_cytosine(upstreamposition,methName,rangeFeatures,'upstreamsequence')
		downstreamcverify=make_methylation_data_frame_and_verify_cytosine(downstreamposition,methName,rangeFeatures,'downstreamsequence')
		upstreamnuc,downstreamnuc=columns_for_nucleotide_and_methylation_content(rangeFeatures,upstreamcverify,downstreamcverify)
		upstreamcapture.append(upstreamnuc)
		downstreamcapture.append(downstreamnuc)
	upstreamconcat = pd.concat(upstreamcapture)
	downstreamconcat = pd.concat(downstreamcapture)
	return upstreamconcat,downstreamconcat

# Methylation rcsorting
def sort_methylation_by_directionality(negStr,posStr):
	posmethylationupstream,posmethylationdownstream = collect_methylation_data_by_element(posStr)
	negmethylationupstream,negmethylationdownstream = collect_methylation_data_by_element(negStr)
	newnegmethylationupstream = negative_directionality_corrected_features(negmethylationupstream)
	newnegmethylationdownstream = negative_directionality_corrected_features(negmethylationdownstream)
	sortedmethylationupstream = concat_positive_and_negative_directionality_and_get_frequency(posmethylationupstream,newnegmethylationupstream)
	sortedmethylationdownstream = concat_positive_and_negative_directionality_and_get_frequency(posmethylationdownstream,newnegmethylationdownstream)
	posSeq = posStr[['chr','start','end','upstreamsequence','downstreamsequence']]
	negSeq = negStr[['chr','start','end','upstreamsequence','downstreamsequence']]
	negSeq['newupstreamsequence'] = negSeq.apply(lambda row: reverse_complement_dictionary(row['upstreamsequence']),axis=1)
	negSeq['newdownstreamsequence'] = negSeq.apply(lambda row: reverse_complement_dictionary(row['downstreamsequence']),axis=1)
	newnegSeq = negSeq[['chr','start','end','newupstreamsequence','newdownstreamsequence']]
	newnegSeq.columns = ['chr','start','end','upstreamsequence','downstreamsequence']
	catStr = pd.concat([posSeq,newnegSeq],axis=0,ignore_index=True)
	sortedmethylationupstream,sortedmethylationdownstream = columns_for_nucleotide_and_methylation_content(catStr,sortedmethylationupstream,sortedmethylationdownstream)
	return sortedmethylationupstream,sortedmethylationdownstream

# Separate on plus and minus orientation, rcsort and return methylation
def sort_elements_by_directionality(directionFeatures,columnCompare):
	negassingment = (directionFeatures[(directionFeatures[columnCompare] == '-')])
	posassingment = (directionFeatures[(directionFeatures[columnCompare] == '+')])
	upstream,downstream = sort_methylation_by_directionality(negassingment,posassingment)
	return upstream,downstream

# For type groups, separate the groups and run the analyses
def separate_dataframe_by_group(List,directionFeatures,typecolumn,fileName):
	bool = (directionFeatures[directionFeatures[typecolumn] == List])
	upstream,downstream = collect_methylation_data_by_element(bool)
	return bool,upstream,downstream

# If want to calculate frequency - might still need to remove duplicates
def calculate_methylation_frequency_remove_duplicates(methylationdf):
	methylationdf['methylationfrequency'] = methylationdf.groupby(['methylationcount','tissue','group'])['methylationcount'].transform('count')
	methylationsort = methylationdf.sort_values(by=['methylationcount','tissue','group'],ascending=True)
	methylationdup = methylationsort.drop_duplicates(['methylationfrequency','tissue','group'],keep='last')#methylationlocation,cytosine
	return methylationdup

# Make the color where hue is determined
def hue_by_group_color(df):
	methlationdf = df.dropna(axis=0)
	methlationdf['barcolor'] = np.where(methlationdf['group'] == 'element','Element','Random')
	return methlationdf

def get_percentage_cpg_methylated(methylationdf):
	methylationdf['percentcpgmethylated'] = methylationdf
	return methylationdf

# Make graphs for fangs
def graph_boundary_methylation(upstream,downstream,fileName):
	info = str(fileName)# + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"
	sns.set_style('ticks')
	pp = PdfPages('Methylation_{0}.pdf'.format(fileName))
	plt.figure(figsize=(14,7))
	plt.suptitle(info,fontsize=10)
	sns.set_palette("husl",n_colors=2)
	
	upstreamhue = hue_by_group_color(upstream)
	downstreamhue = hue_by_group_color(downstream)
	
	upstreamcalc = calculate_methylation_frequency_remove_duplicates(upstreamhue)
	downstreamcalc = calculate_methylation_frequency_remove_duplicates(downstreamhue) 
	
# 	firsttissue = upstreamcalc.iloc[0]['tissue']
# 	firsttissueup = (upstreamcalc[upstreamcalc['tissue']==firsttissue])
# 	firsttissuedown = (downstreamcalc[downstreamcalc['tissue']==firsttissue])
	
	surroundingfang = periphery*2

	if not splitgraphs:
		gs = gridspec.GridSpec(1,1,height_ratios=[1])
		gs.update(hspace=.8)
		frames = [upstreamhue,downstreamhue]
		catstreams = pd.concat(frames)
		catstreams['methylationfrequencyboundaries'] = catstreams.groupby(['tissue','group'])['methylationfrequency'].transform('sum')
		ax0 = plt.subplot(gs[0,:])
		sns.barplot(data=catstreams,x='tissue',y='methylationfrequencyboundaries',hue='barcolor',ax=ax0) #percentageunmethylatedcpgcombineboundary
		ax0.set_title('Fraction CpGs Methylated Across {0}bp Surrounding Fang'.format(surroundingfang),size=8)
		ax0.set_ylabel('Fraction',size=8)
# 		tissueframes = [firsttissueup,firsttissuedown]
# 		cattissues = pd.concat(tissueframes)
# 		for p in ax0.patches:
# 			ax0.annotate(str(p.get_height()), (p.get_x(), p.get_height()))
# 		ax1 = plt.subplot(gs[1,0])
# 		sns.barplot(data=cattissues,x='tissue',y='cpgsequencecountsum',hue='barcolor',ax=ax1)
# 		ax1.set_title('Total CpGs',size=8)
# 		ax1.set_ylabel('Counts',size=8)
# 		ax2 = plt.subplot(gs[1,1])
# 		sns.barplot(data=cattissues,x='tissue',y='candgsequencecountsum',hue='barcolor',ax=ax2)
# 		ax2.set_title('Total C and Gs',size=8)
# 		ax2.set_ylabel('Counts',size=8)
# 		ax1.legend_.remove()
# 		ax2.legend_.remove()

	else:
		gs = gridspec.GridSpec(2,1,height_ratios=[1,1],width_ratios=[1])
		gs.update(hspace=.8)
		ax0 = plt.subplot(gs[0,:])
		sns.barplot(data=upstreamcalc,x='tissue',y='methylationfrequency',hue='barcolor',ax=ax0) #percentageunmethylatedcpg
		ax0.set_title('Fraction CpGs Methylated Across {0}bp Surrounding Upstream Fang'.format(surroundingfang),size=8)
		ax0.set_ylabel('Fraction',size=8)
# 		ax1 = plt.subplot(gs[0,1])
# 		sns.barplot(data=firsttissueup,x='tissue',y='cpgsequencecountsum',hue='barcolor',ax=ax1)
# 		ax1.set_title('Total CpGs',size=8)
# 		ax1.set_ylabel('Counts',size=8)
# 		ax2 = plt.subplot(gs[0,0])
# 		sns.barplot(data=firsttissueup,x='tissue',y='candgsequencecountsum',hue='barcolor',ax=ax2)
# 		ax2.set_title('Total C and Gs',size=8)
# 		ax2.set_ylabel('Counts',size=8)
# 		ax1.legend_.remove()
# 		ax2.legend_.remove()
		ax3 = plt.subplot(gs[1,:])
		sns.barplot(data=downstreamcalc,x='tissue',y='methylationfrequency',hue='barcolor',ax=ax3) #percentageunmethylatedcpg
		ax3.set_title('Fraction CpGs Methylated Across {0}bp Surrounding Downstream Fang'.format(surroundingfang),size=8)
		ax3.set_ylabel('Fraction',size=8)
# 		ax4 = plt.subplot(gs[1,1])
# 		sns.barplot(data=firsttissuedown,x='tissue',y='cpgsequencecountsum',hue='barcolor',ax=ax4)
# 		ax4.set_title('Total CpGs',size=8)
# 		ax4.set_ylabel('Counts',size=8)
# 		ax5 = plt.subplot(gs[1,0])
# 		sns.barplot(data=firsttissuedown,x='tissue',y='candgsequencecountsum',hue='barcolor',ax=ax5)
# 		ax5.set_title('Total C and Gs',size=8)
# 		ax5.set_ylabel('Counts',size=8)
# 		ax4.legend_.remove()
# 		ax5.legend_.remove()

	sns.despine()
	pp.savefig()
	pp.close()

# Get the empirical probability for each direction classification
def make_probabilites_for_direction(directionFeatures,probabilitycolumn):
	lenAll = float(len(directionFeatures.index))
	numPlus = (directionFeatures[probabilitycolumn] == '+').sum()/lenAll
	numMinus = (directionFeatures[probabilitycolumn] == '-').sum()/lenAll
	numEqual = (directionFeatures[probabilitycolumn] == '=').sum()/lenAll
	probOptions = [numMinus,numPlus,numEqual]
	print 'made probabilities for df: {0} for +, {1} for - and {2} for ='.format(numPlus,numMinus,numEqual)
	return probOptions

def save_panda_data_frame(df,filename):
	df.to_csv(filename,sep="\t",header=True)

def main():
	# Get and set parameters
	args = get_args()
	set_global_variables(args)
	paramlabels = '{0}_{1}_{2}_{3}_{4}_{5}'.format(elementsize,periphery,binDir,eFiles,motifdirectionality,stringName)

	# Get coords and strings for elements
	rangeFeatures = collect_element_coordinates(eFiles)
	directionFeatures = assign_directionality_from_arg_or_boundary(rangeFeatures,eFiles)
	
	# Get the probability for each directional assignment, and use to randomly assign the correct number of random directions
	dirOptions = ['-','+','=']
	probOptions = make_probabilites_for_direction(directionFeatures,'directionality')

	collectupstream,collectdownstream = [],[]
	collectreversecomplementupstream,collectreversecomplementdownstream = [],[]
	
	numbertissues = []
	numbertissues.append(len(mFiles))
	numbertissues.append('numtissues')
	
	if labelcolumn:
		lenthrandom =[]
		typeList = directionFeatures['type'].unique()
		for type in typeList:
			typecollectupstream,typecollectdownstream = [],[]
			typecollectreversecomplementupstream,typecollectreversecomplementdownstream = [],[]
			typeBool,typeupstreammethylation,typedownstreammethylation = separate_dataframe_by_group(type,rangeFeatures,'type',eFiles)
			typeupstreammethylation['group'] = 'element'
			typedownstreammethylation['group'] = 'element'
			typecollectupstream.append(typeupstreammethylation)
			typecollectdownstream.append(typedownstreammethylation)
			if directionalitycolumn:
				probOptionstype = make_probabilites_for_direction(typeBool,'directionality')
			else:
				probOptionstype = make_probabilites_for_direction(typeBool,'compareBoundaries')
			if reverseComplement:
				if directionalitycolumn:
					revtypeupstreammethylation,revtypedownstreammethylation = sort_elements_by_directionality(typeBool,'directionality')
					revtypeupstreammethylation['group'] = 'element'
					revtypedownstreammethylation['group'] = 'element'
					typecollectreversecomplementupstream.append(revtypeupstreammethylation)
					typecollectreversecomplementdownstream.append(revtypedownstreammethylation)
				else:
					revtypeupstreammethylation,revtypedownstreammethylation = sort_elements_by_directionality(typeBool,'compareBoundaries')
					revtypeupstreammethylation['group'] = 'element'
					revtypedownstreammethylation['group'] = 'element'
					typecollectreversecomplementupstream.append(revtypeupstreammethylation)
					typecollectreversecomplementdownstream.append(revtypedownstreammethylation)
			if rFiles:
				lengthrandom.append(len(rFiles))
				lengthrandom.append('randomfiles')
				for randomFile in rFiles:
					randomFeatures = collect_element_coordinates(randomFile)
					randirFeatures= evaluate_boundaries_size_n(randomFeatures,randomFile)
					rantypeBool,randomtypeupstreammethylation,randomtypedownstreammethylation = separate_dataframe_by_group(type,randirFeatures,'type',randomFile)
					randomtypeupstreammethylation['group'] = 'random{0}'.format(randomFile)
					randomtypedownstreammethylation['group'] = 'random{0}'.format(randomFile)
					typecollectupstream.append(randomtypeupstreammethylation)
					typecollectdownstream.append(randomtypedownstreammethylation)
					if reverseComplement:
						randomrevtypeupstreammethylation,randomrevtypedownstreammethylation = sort_elements_by_directionality(rantypeBool,'compareBoundaries')
						randomrevtypeupstreammethylation['group'] = 'random{0}'.format(i)
						randomrevtypedownstreammethylation['group'] = 'random{0}'.format(i)
						typecollectreversecomplementupstream.append(randomrevtypeupstreammethylation)
						typecollectreversecomplementdownstream.append(randomrevtypedownstreammethylation)
			else:
				lengthrandom =[]
				lengthrandom.append(randomassignments)
				lengthrandom.append('randomassingments')
				for i in range(randomassignments):
					typeBool['randomDirectiontype'] = np.random.choice(dirOptions,len(typeBool.index),p=probOptionstype)
					typedirupstreammethylation,typedirdownstreammethylation = sort_elements_by_directionality(typeBool,'randomDirectiontype')
					typedirupstreammethylation['group'] = 'element'
					typedirdownstreammethylation['group'] = 'element'
					typecollectreversecomplementupstream.append(typedirupstreammethylation)
					typecollectreversecomplementdownstream.append(typedirdownstreammethylation)
			typeconcatupstream = pd.concat(typecollectupstream)
			typeconcatdownstream = pd.concat(typecollectdownstream)
			if reverseComplement:
				typeconcatupstreamreverse = pd.concat(typecollectreversecomplementupstream)
				typeconcatdownstreamreverse = pd.concat(typecollectreversecomplementdownstream)
				graph_boundary_methylation(typeconcatupstreamreverse,typeconcatdownstreamreverse,'{0}_rc_{1}_{2}_{3}'.format(type,paramlabels,lengthrandom,numbertissues))
			else:
				graph_boundary_methylation(typeconcatupstream,typeconcatdownstream,'{0}_{1}_{2}_{3}'.format(type,paramlabels,lengthrandom,numbertissues))

	else:
		lengthrandom =[]
		allmethylationupstream,allmethylationdownstream = collect_methylation_data_by_element(rangeFeatures)
		allmethylationupstream['group'] = 'element'
		allmethylationdownstream['group'] = 'element'
		collectupstream.append(allmethylationupstream)
		collectdownstream.append(allmethylationdownstream)
		if reverseComplement:
			if directionalitycolumn:
				revmethylationupstream,revmethylationdownstream = sort_elements_by_directionality(directionFeatures,'directionality')
				revmethylationupstream['group'] = 'element'
				revmethylationdownstream['group'] = 'element'
				collectreversecomplementupstream.append(revmethylationupstream)
				collectreversecomplementdownstream.append(revmethylationdownstream)
			else:
				revmethylationupstream,revmethylationdownstream = sort_elements_by_directionality(directionFeatures,'directionality')
				revmethylationupstream['group'] = 'element'
				revmethylationdownstream['group'] = 'element'
				collectreversecomplementupstream.append(revmethylationupstream)
				collectreversecomplementdownstream.append(revmethylationdownstream)
		if rFiles:
			lengthrandom.append(len(rFiles))
			lengthrandom.append('randomfiles')
			for randomFile in rFiles:
				randomFeatures = collect_element_coordinates(randomFile)
				randirFeatures= assign_directionality_from_arg_or_boundary(randomFeatures,randomFile)
				allrandommethylationupstream,allrandommethylationdownstream = collect_methylation_data_by_element(randirFeatures)
				allrandommethylationupstream['group'] = 'random{0}'.format(randomFile)
				allrandommethylationdownstream['group'] = 'random{0}'.format(randomFile)
				collectupstream.append(allrandommethylationupstream)
				collectdownstream.append(allrandommethylationdownstream)
				if reverseComplement:
					randomrevmethylationupstream,randomrevmethylationdownstream = sort_elements_by_directionality(randirFeatures,'directionality')
					randomrevmethylationupstream['group'] = 'random{0}'.format(i)
					randomrevmethylationdownstream['group'] = 'random{0}'.format(i)
					collectreversecomplementupstream.append(randomrevmethylationupstream)
					collectreversecomplementdownstream.append(randomrevmethylationdownstream)
		else:
			lengthrandom.append(randomassignments)
			lengthrandom.append('randomassingments')
			for i in range(randomassignments):
				directionFeatures['randomDirection'] = np.random.choice(dirOptions,len(directionFeatures.index),p=probOptions)
				randirrevmethylationupstream,randirrevmethylationdownstream = sort_elements_by_directionality(directionFeatures,'randomDirection')
				randirrevmethylationupstream['group'] = 'random{0}'.format(i)
				randirrevmethylationdownstream['group'] = 'random{0}'.format(i)
				collectupstream.append(randirrevmethylationupstream)
				collectdownstream.append(randirrevmethylationdownstream)
				collectreversecomplementupstream.append(randirrevmethylationupstream)
				collectreversecomplementdownstream.append(randirrevmethylationdownstream)
		concatupstream = pd.concat(collectupstream)
		concatdownstream = pd.concat(collectdownstream)
		if reverseComplement:
			concatupstreamreverse = pd.concat(collectreversecomplementupstream)
			concatdownstreamreverse = pd.concat(collectreversecomplementdownstream)
			graph_boundary_methylation(concatupstreamreverse,concatdownstreamreverse,'all_rc_{0}_{1}_{2}'.format(paramlabels,lengthrandom,numbertissues))
		else:
			graph_boundary_methylation(concatupstream,concatdownstream,'all_{0}_{1}_{2}'.format(paramlabels,lengthrandom,numbertissues))

if __name__ == "__main__":
	main()