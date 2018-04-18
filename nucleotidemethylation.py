"""
Script to run methylation analysis as heatmap

Wren Saylor
March 2018

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
from collections import Counter
from scipy import stats

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile",type=str,help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("-r","--randomfile",required=True,type=argparse.FileType('rU'),help="A file containing a list of paths to the random regions equable with your elements to plot in contrast")
	parser.add_argument("-m","--methylationfile",type=argparse.FileType('rU'),help="A file containing a list of paths to the methlylation bedfiles, where coverage is in the 4th column and percentage in the 5th. Files names in format; 'tissue_tissue-number.filetype'")
	parser.add_argument("-l", "--labelcolumn",type=int,help='column in the element file where the label (exonic, intergenic, intronic) is') # way to only read in certain entries, like only read in if 'intergenic'
	parser.add_argument("-d", "--directionalitycolumn",type=int,help='column in the element file where directionality is, if not supplied will infer by AT content')
	parser.add_argument("-i", "--idcolumn",type=int,help='column in the element file where the id is')
	parser.add_argument("-g","--genome",type=str, default="hg19.genome")
	parser.add_argument("-f","--fasta",type=str,default="hg19.fa")
	parser.add_argument("-t","--total",type=int,default="2200",help='total size of region to look at (region + flanks), should be an even number, suggested to be at least triple your element')
	parser.add_argument("-p","--inset",type=int,default="50",help='size into your element from the boundaries, should be an even number')
	parser.add_argument("-b","--bin",type=int,default="100",help='size of bins used to compare element ends and determine directionality')
	parser.add_argument("-e","--element",type=int,help='size of your element (region without flanks), should be an even number, if not provided will use the smallest size of your input data')
	parser.add_argument("-mu","--methylationthresholdpercentageupper",type=int, default="100",help='size to threshold percentage methylation data upper bound')
	parser.add_argument("-ml","--methylationthresholdpercentagelower",type=int, default="1",help='size to threshold percentage methylation data lower bound')
	parser.add_argument("-mc","--methylationthresholdcoverage",type=int,default="10",help='size to threshold uncapped coverage of methylation data to send to percentage methylation, field often uses 10')
	parser.add_argument("-ms","--combinemethylationsamples",action='store_true',help='whether to combine those samples with the same tissue/cell type or leave as seperate')
	parser.add_argument('-s',"--stringname",type=str,help='string to add to the outfile name')
	parser.add_argument('-c',"--reversecomplement",action='store_true',help='if you want the reverse complement to be plotted')
	return parser.parse_args()

# set all args that will be used throughout the script
def set_global_variables(args):
	global periphery
	global binDir
	global elementsize
	global methPerThreshupper
	global methPerThreshlower
	global methCovThresh
	global eFiles
	global rFiles
	global mFiles
	global labelcolumn
	global directionalitycolumn
	global idcolumn
	global sizeGenome
	global faGenome
	global stringName
	global reverseComplement
	global combinesamples
	global num
	num = args.total
	periphery = args.inset
	binDir = args.bin
	elementsize = args.element
	methPerThreshupper = args.methylationthresholdpercentageupper
	methPerThreshlower = args.methylationthresholdpercentagelower
	methCovThresh = args.methylationthresholdcoverage
	eFiles = args.efile
	mFiles = [line.strip() for line in args.methylationfile]
	if args.randomfile:
		rFiles = [line.strip() for line in args.randomfile]
	else:
		rFiles = None
	labelcolumn = args.labelcolumn
	directionalitycolumn = args.directionalitycolumn
	idcolumn = args.idcolumn
	sizeGenome = args.genome
	faGenome = args.fasta
	stringName = args.stringname
	reverseComplement = args.reversecomplement
	combinesamples = args.combinemethylationsamples

def set_ploting_parameters():
	# Locations for plotting with sliding window
	global plotLineLocationOne # upstream element boundary
	global plotLineLocationTwo # downstream element boundary
	global plotLineLocationThree # upstream element inset
	global plotLineLocationFour # downstream element inset
	plotLineLocationOne = ((num-uce)/2)+periphery
	plotLineLocationTwo = ((num-uce)/2)+(uce-periphery)
	plotLineLocationThree = (num-uce)/2
	plotLineLocationFour = ((num-uce)/2)+uce
	global centerelementpoint
	centerelementpoint = ((num-uce)/2)+(uce/2)
	print 'center point', centerelementpoint
	print 'set plotting parameters'

# get bt features
def get_bedtools_features(strFileName):
	return pbt.BedTool(strFileName)

# get the correct range for fang evaluation
def collect_coordinates_for_element_positions(btFeatures):
	midFeatures = pd.read_table(btFeatures.fn, header=None)
	midFeatures['middle'] = midFeatures.loc[:,1:2].mean(axis=1).astype(int)
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
	if idcolumn:
		midFeatures['id'] = midFeatures.loc[:,idcolumn]
	else:
		midFeatures.insert(0,'id',range(0,0+len(midFeatures)))
	if labelcolumn:
		midFeatures['type'] = midFeatures.loc[:,labelcolumn]
	if directionalitycolumn:
		midFeatures['directionality'] = midFeatures.loc[:,directionalitycolumn]
	global flankSize
	flankSize = (num - uce)/2
	global inregion
	inregion = uce-(periphery*2)
	midFeatures['sCenter'] = midFeatures['middle'].astype(int) - (inregion/2)
	midFeatures['eCenter'] = midFeatures['middle'].astype(int) + (inregion/2)
	midFeatures['sEdge'] = midFeatures['start'].astype(int) + periphery
	midFeatures['eEdge'] = midFeatures['end'].astype(int) - periphery
	midFeatures['sBoundary'] = midFeatures['start'].astype(int) - flankSize
	midFeatures['eBoundary'] = midFeatures['end'].astype(int) + flankSize
	return midFeatures

# get the strings for sliding window regions
def get_fasta_for_element_coordinates(rangeFeatures):#rangeFeatures,faGenome
	rangeFeatures['sBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','sBoundary','start']].values.astype(str).tolist()))
	rangeFeatures['sEdgeSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','start','sEdge']].values.astype(str).tolist()))
	rangeFeatures['MiddleSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','sCenter','eCenter']].values.astype(str).tolist()))
	rangeFeatures['eEdgeSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','eEdge','end']].values.astype(str).tolist()))
	rangeFeatures['eBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','end','eBoundary']].values.astype(str).tolist()))
	rangeFeatures['combineString'] = rangeFeatures['sBoundarySeq'].astype(str) + rangeFeatures['sEdgeSeq'].astype(str) + rangeFeatures['MiddleSeq'].astype(str) + rangeFeatures['eEdgeSeq'].astype(str) + rangeFeatures['eBoundarySeq'].astype(str)
	rangeFeatures['combineString'] = rangeFeatures['combineString'].str.upper()
	rangeFeatures=rangeFeatures.drop(['sBoundarySeq','sEdgeSeq','MiddleSeq','eEdgeSeq','eBoundarySeq'],axis=1)
	return rangeFeatures

# used in get_fasta_for_element_coordinates to extract just the fasta strings
def get_just_fasta_sequence_for_feature(inFeature):
	seqFeature = inFeature.sequence(fi=faGenome)
	outFeature = pd.read_table(seqFeature.seqfn)
	outSequence = outFeature[::2]
	return outSequence.reset_index(drop=True)

# get coordinates with flanking regions
def collect_element_coordinates(fileName):
	btFeatures = get_bedtools_features(fileName)
	subsetFeatures = collect_coordinates_for_element_positions(btFeatures)
	return get_fasta_for_element_coordinates(subsetFeatures)

# Do the comparison between boundaries to get + - or =
def calculate_nucleotides_at(element,size):
	start,end,perSize = element[:size],element[-size:],[]
	perSize.append(eval('100*float(start.count("A") + start.count("a") + start.count("T") + start.count("t"))/len(start)'))
	perSize.append(eval('100*float(end.count("A") + end.count("a") + end.count("T") + end.count("t"))/len(end)'))
	return perSize

# Directionality, as inferred by comparing first and last n base pairs from input parameters
def compare_boundaries_size_n(element,size):
	perSize = calculate_nucleotides_at(element,size)
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
		numval = rangeFeatures['directionality'].value_counts()
		rangeFeatures.drop(labels='feature',axis=1,inplace=True)
	print '{0} in {1}'.format(numval,fileName)
	return rangeFeatures

# save file from bedtool
def save_bedtool_as_bedfile(btObject,strFilename):
	btObject.saveas(strFilename)

# convert bedtool to panda
def convert_bedtool_to_panda(btobject):
	save_bedtool_as_bedfile(btobject,'temp.bed')
	return pd.read_table(btobject.fn, header=None)

# convert panda to bedtool
def convert_panda_to_bed_format(panda):
	arArFeatures = panda.values.tolist()
	return get_bedtools_features(arArFeatures)

# Threshold methylation data by coverage and percentage
def threshold_methylation_data(methFeatures,methName):
	pdmethFeatures = convert_bedtool_to_panda(methFeatures)
	pdmethThresh = (pdmethFeatures[(pdmethFeatures.loc[:,3] >= methCovThresh) & (pdmethFeatures.loc[:,4] <= methPerThreshupper) & (pdmethFeatures.loc[:,4] >= methPerThreshlower)])
	btmethThresh = convert_panda_to_bed_format(pdmethThresh)
	print 'methylation coverage is being thresholded at {0} and percentage between {1} and {2}'.format(methCovThresh,methPerThreshlower,methPerThreshupper)
	initiallength = len(pdmethFeatures.index)
	checklength = initiallength - len(pdmethThresh.index)
	print 'there were {0} out of {1} total removed by thresholding {2}'.format(checklength,initiallength,methName)
	return btmethThresh

# Intersect regions from the methylation data with element regions
def intersect_methylation_data(btFeatures,meFeature,meName,eStart,eEnd,addvalue):
	initiallength = len(meFeature)
	intersectboundary = meFeature.intersect(btFeatures[['chr',eStart,eEnd,'id']].values.astype(str).tolist(),wb=True,wa=True)
	checklength = len(intersectboundary)
	print 'there were {0} intersections out of {1} in {2}'.format(checklength,initiallength,meName)
	if len(intersectboundary) != 0:
		pdmeth=convert_bedtool_to_panda(intersectboundary)
		pdmeth.columns = ['mchr','mstart','mstop','coverage','percentage','chr','estart','estop','id']
		pdmeth['strand'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmeth[['mchr','mstart','mstop']].values.astype(str).tolist()))
		pdmeth['strand'] = pdmeth['strand'].str.upper()
		pdmeth['methlocation'] = pdmeth['estop'].astype(int)-pdmeth['mstart'].astype(int)+addvalue
		outmeth = pdmeth[['chr','mstart','mstop','estart','estop','id','percentage','methlocation','strand']]
		outmeth['Tissue'] = meName.split(".")[0]
		subfeatures = btFeatures[['id','directionality']]
		merge = pd.merge(outmeth,subfeatures,how='left',on='id')
		sortoutmeth = merge[['chr','id','directionality','methlocation','percentage','strand','Tissue']]
	else:
		sortoutmeth = pd.DataFrame(np.nan,index=[0],columns=['chr','id','directionality','methlocation','percentage','strand'])
		sortoutmeth['Tissue'] = meName.split(".")[0]
	return sortoutmeth

# Run the analysis to extract methylation properties
def collect_methylation_data_by_element(pdfeatures):
	capture = []
	for methname in mFiles:
		mfeatures=get_bedtools_features(methname)
		pdthresh=threshold_methylation_data(mfeatures,methname)
		upboundary=intersect_methylation_data(pdfeatures,pdthresh,methname,'sBoundary','sEdge',0)
		middle=intersect_methylation_data(pdfeatures,pdthresh,methname,'sCenter','eCenter',flankSize+periphery)
		downboundary=intersect_methylation_data(pdfeatures,pdthresh,methname,'eEdge','eBoundary',flankSize+(uce-periphery))
		capture.append(upboundary)
		capture.append(middle)
		capture.append(downboundary)
	totalconcat = pd.concat(capture)
	totalconcat.reset_index(drop=True,inplace=True)
	return totalconcat

# invert methlocation column
def invert_location(pdfeatures):
	originalRange = range(0,num)#-periphery,periphery+1
	reverseRange = originalRange[::-1]
	rangeDict = dict(zip(originalRange,reverseRange))
	pdfeatures['newmethlocation'] = pdfeatures.loc[:,'methlocation'].map(rangeDict)
	pdfeatures.drop(labels='methlocation',axis=1,inplace=True)
	pdfeatures.rename(columns={'newmethlocation':'methlocation'},inplace=True)
	return pdfeatures

# Reverse those column that need to be reversed in case of - directionality
def negative_directionality_columns_to_modify(negFeatures):
	invFeatures = invert_location(negFeatures)
	seqDict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
	invFeatures['newstrand'] = invFeatures.loc[:,'strand'].map(seqDict)
	newnegFeatures = invFeatures[['chr','id','directionality','methlocation','percentage','newstrand','Tissue','group']]
	newnegFeatures.columns = ['chr','id','directionality','methlocation','percentage','strand','Tissue','group']
	return newnegFeatures

# Separate on plus and minus orientation, rcsort and return methylation
def modify_negative_directionality_elements(btFeatures):
	if not btFeatures['directionality'].isnull().all():
		negFeatures = (btFeatures.loc[btFeatures['directionality']=='-'])
		nonFeatures = btFeatures.loc[btFeatures['directionality']!='-']
		newnegFeatures = negative_directionality_columns_to_modify(negFeatures)
		outFeatures = pd.concat([newnegFeatures,nonFeatures])
	else:
		outFeatures = btFeatures
		print 'There are no features with "-" directionality, will use the original unsorted set for the analysis'
	return outFeatures

# get the total cpg in a column of strings
def collect_total_avail_cpg_in_column(btFeatures,column):
	btFeatures['cpgcount'] = btFeatures.apply(lambda row: float(row[column].count("CG")),axis=1)
	btFeatures['cgcount'] = btFeatures.apply(lambda row: (row[column].count("G")+row[column].count("C")),axis=1)
	return btFeatures['cpgcount'].sum(),btFeatures['cgcount'].sum()

# get the reverse complement
def reverse_complement_dictionary(sequence):
	seqDict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
	return "".join([seqDict[base] for base in reversed(sequence)])

# collect the values for up/down/rev up/rev down stream total cpg counts
def collect_total_avail_cpg_values(btFeatures):
	cpg,cg = collect_total_avail_cpg_in_column(btFeatures,'combineString')
	if not btFeatures['directionality'].isnull().all():
		negFeatures = (btFeatures.loc[btFeatures['directionality']=='-'])
		nonFeatures = btFeatures.loc[btFeatures['directionality']!='-']
		negFeatures['reverseComplement']=negFeatures.apply(lambda row: reverse_complement_dictionary(row['combineString']),axis=1)
		negFeatures.drop(labels='combineString',axis=1,inplace=True)
		newnegFeatures = negFeatures.rename(columns={'reverseComplement':'combineString'})
		revFeatures = pd.concat([newnegFeatures,nonFeatures])
		revcpg,revcg = collect_total_avail_cpg_in_column(revFeatures,'combineString')
	else:
		revcpg,revcg=cpg,cg
		print 'There are no features with "-" directionality'
	return cpg,revcpg,cg,revcg

# make the original input data frame
def collect_input_data_frame(file):
	rangeFeatures = collect_element_coordinates(file)
	directionFeatures = assign_directionality_from_arg_or_boundary(rangeFeatures,file)
	return directionFeatures

# run the whole series of steps to make the data frame for each set of elements to plot
def run_whole_analysis_for_boundaries(pdfeatures,label):
	cpg,revcpg,cg,revcg = collect_total_avail_cpg_values(pdfeatures)
	allstream = collect_methylation_data_by_element(pdfeatures)
	allstream['group'] = label
	allstream['organization'] = 'unsorted'
	allstream['cpgsum'] = cpg
	allstream['cgsum'] = cg
	allreverse = modify_negative_directionality_elements(allstream)
	allreverse['group'] = label
	allreverse['organization'] = 'rcsorted'
	allreverse['cpgsum'] = revcpg
	allreverse['cgsum'] = revcg
	return allstream,allreverse

# separate the elements and the random regions
def seperate_elements_and_random(pdfeatures,column,value):
	element = pdfeatures[pdfeatures[column]==value]
	random  = pdfeatures[pdfeatures[column]!=value]
	return element,random

# get standard deviation, from ruth's random region script
def getPopSD(arObservedOverlaps):
	floatLen = float(len(arObservedOverlaps))
	floatMean = float(sum(arObservedOverlaps))/len(arObservedOverlaps)
	dSumOfSquares = sum([((float(number) - floatMean) ** 2) for number in arObservedOverlaps])
	dVariance = float(dSumOfSquares) / floatLen
	return math.sqrt(dVariance)

# ks test from ruth's random region script
def KSTest(aOverlapBP):
	mean = float(sum(aOverlapBP)) / len(aOverlapBP)
	sd = getPopSD(aOverlapBP)
	rvNormMatched = stats.norm.rvs(loc=mean, scale=sd, size=len(aOverlapBP))
	npArOverlapBP = np.array(aOverlapBP)
	ksStat, KsPval = stats.ks_2samp(npArOverlapBP, rvNormMatched)
	if KsPval <= 0.05:
		strKSresult = "No"
	else:
		strKSresult = "Yes"
	return ksStat, KsPval, strKSresult

# run ks test for normal distribution and choose appropriate stats test
def run_appropriate_test(element,random):
	ksStat,KsPval,strKSresult = KSTest(element)
	if strKSresult == 'Yes':
		statcoef,statpval = stats.ttest_ind(element,random)
		stattest = 'TTest'
		formatpval = '{:.01e}'.format(statpval)
	else:
		statcoef,statpval = stats.mannwhitneyu(element,random)
		stattest = 'MannWhitneyUTest'
		formatpval = '{:.01e}'.format(statpval)
	return formatpval,statcoef,stattest

# save panda
def save_panda(pdData, strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True)

# count by feature for graphing
def group_and_count_data_frame_by_column(pdfeatures,incol,outcol):
	methsort = pdfeatures.sort_values(by=['Tissue','group',incol],ascending=True)
	methsort['methgroupby'] = methsort.groupby(['Tissue',incol,'group'])[incol].transform('count')
	removedups = methsort.drop_duplicates(['group',incol,'Tissue','methgroupby'],keep='last') # only count each groupby object once!!
	removedups[outcol] = removedups.groupby(['group','Tissue','methgroupby',incol])['methgroupby'].transform('sum')
	return removedups

# format features for graphing
def format_data_frame_by_column(pdfeatures,countcol):
	new_index = range(0,num)
	pivotfeatures = pd.pivot_table(pdfeatures,index='methlocation',columns='Tissue',values=countcol)
	pivotfeatures.columns.name = None
	reindexfeatures = pivotfeatures.reindex(new_index,fill_value=0)
	reindexfeatures.fillna('0',inplace=True)
	reindexfeatures.index.name = None
	transposefeatures = reindexfeatures.T
	outfeatures = transposefeatures[transposefeatures.columns].astype(float)
	return outfeatures

# Make graphs for fangs
def graph_methylation(pdfeatures,filelabel):
	set_ploting_parameters()
	methfiles = [(file.split("-")[0]) for file in mFiles]
	methcount = Counter(methfiles)
	info = str(filelabel)
	sns.set_style('ticks')
	if reverseComplement:
		sorted = pdfeatures.loc[pdfeatures['organization']=='rcsorted']
		pp = PdfPages('Methylation_Heatmap_ReverseComplementSorted_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.pdf'.format(eFiles,stringName,filelabel,elementsize,binDir,periphery,methPerThreshlower,methPerThreshupper,methCovThresh))
		pstat = 'Statistic_Heatmap__ReverseComplementSorted_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.pdf'.format(eFiles,stringName,filelabel,elementsize,binDir,periphery,methPerThreshlower,methPerThreshupper,methCovThresh)
	else:
		sorted = pdfeatures.loc[pdfeatures['organization']=='unsorted']
		pp = PdfPages('Methylation_Heatmap__{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.pdf'.format(eFiles,stringName,filelabel,elementsize,binDir,periphery,methPerThreshlower,methPerThreshupper,methCovThresh))
		pstat = 'Statistic_Heatmap__{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.pdf'.format(eFiles,stringName,filelabel,elementsize,binDir,periphery,methPerThreshlower,methPerThreshupper,methCovThresh)
	plt.figure(figsize=(14,7))
	plt.suptitle(info,fontsize=16)
	removedups = group_and_count_data_frame_by_column(pdfeatures,'methlocation','countlocation')
	element,random = seperate_elements_and_random(removedups,'group','element')
	formatelement = format_data_frame_by_column(element,'countlocation')
	collectgroup = []
	for group in random['group'].unique():
		grouprandom = random[random['group']==group]
		formatrandom = format_data_frame_by_column(grouprandom,'countlocation')
		collectgroup.append(formatrandom)
	averagegroup = pd.concat([each.stack() for each in collectgroup],axis=1).apply(lambda x:x.mean(),axis=1).unstack()
	collectstats = []
	formatpval,statcoef,stattest = run_appropriate_test(formatelement.values.flatten(),averagegroup.values.flatten())
	statstable = pd.DataFrame(['loccount','whole set',formatpval,statcoef,stattest],index=['count set','comparison group','p value','coefficient','stats test'])
	collectstats.append(statstable)
	formatpvalelement,statcoefelement,stattestelement = run_appropriate_test(formatelement.loc[:,:plotLineLocationOne].values.flatten(),formatelement.loc[:,plotLineLocationTwo:].values.flatten())
	statstableelement = pd.DataFrame(['loccount','element',formatpvalelement,statcoefelement,stattestelement],index=['count set','comparison group','p value','coefficient','stats test'])
	collectstats.append(statstableelement)
	formatpvalrandom,statcoefrandom,stattestrandom = run_appropriate_test(averagegroup.loc[:,:plotLineLocationOne].values.flatten(),averagegroup.loc[:,plotLineLocationTwo:].values.flatten())
	statstablerandom = pd.DataFrame(['loccount','random',formatpvalrandom,statcoefrandom,stattestrandom],index=['count set','comparison group','p value','coefficient','stats test'])
	collectstats.append(statstablerandom)
	for tissue in element['Tissue'].unique():
		tissueelement = formatelement.loc[tissue]
		tissuerandom = formatrandom.loc[tissue]
		formatpval,statcoef,stattest = run_appropriate_test(tissueelement,tissuerandom)
		statstable = pd.DataFrame(['loccount',tissue,formatpval,statcoef,stattest],index=['count set','comparison group','p value','coefficient','stats test'])
		collectstats.append(statstable)
		formatpvalelement,statcoefelement,stattestelement = run_appropriate_test(tissueelement.loc[:,:plotLineLocationOne].values.flatten(),tissueelement.loc[:,plotLineLocationTwo:].values.flatten())
		statstableelement = pd.DataFrame(['loccount','{0}_element'.format(tissue),formatpvalelement,statcoefelement,stattestelement],index=['count set','comparison group','p value','coefficient','stats test'])
		collectstats.append(statstableelement)
		formatpvalrandom,statcoefrandom,stattestrandom = run_appropriate_test(tissuerandom.loc[:,:plotLineLocationOne].values.flatten(),tissuerandom.loc[:,plotLineLocationTwo:].values.flatten())
		statstablerandom = pd.DataFrame(['loccount','{0}_random'.format(tissue),formatpvalrandom,statcoefrandom,stattestrandom],index=['count set','comparison group','p value','coefficient','stats test'])
		collectstats.append(statstablerandom)
	catstat = pd.concat(collectstats,axis=1)
	catstat.reset_index(drop=True,inplace=True)
	save_panda(catstat.T,'{0}.txt'.format(pstat))
	gs = gridspec.GridSpec(1,1,height_ratios=[1],width_ratios=[1])
	gs.update(hspace=.8)
	ax0 = plt.subplot(gs[0,:])
	heatmap0 = sns.heatmap(formatelement,cmap='PuBu',ax=ax0,xticklabels=100)
	ax0.set_ylabel('Tissue',size=8)
	ax0.set_xlabel('Location',size=6)
	ax0.tick_params(labelsize=8)
	ylabels0 = formatelement.index
	ax0.set_yticklabels(ylabels0,minor=False,rotation=0)
	ax0.set_yticks(np.arange(formatelement.shape[0]) + 0.5, minor=False)
	ax0.set_title('Methylation Counts per Location - Elements',size=8)
	sns.despine()
	pp.savefig()
	gs = gridspec.GridSpec(1,1,height_ratios=[1],width_ratios=[1])
	gs.update(hspace=.8)
	ax1 = plt.subplot(gs[0,:])
	heatmap1 = sns.heatmap(averagegroup,cmap='RdPu',ax=ax1,xticklabels=100)
	ax1.set_ylabel('Tissue',size=8)
	ax1.set_xlabel('Location',size=6)
	ax1.tick_params(labelsize=8)
	ylabels1 = formatrandom.index
	ax1.set_yticklabels(ylabels1,minor=False,rotation=0)
	ax1.set_yticks(np.arange(formatrandom.shape[0]) + 0.5, minor=False)
	ax1.set_title('Methylation Counts per Location - Random Regions',size=8)
	sns.despine()
	plt.savefig(pp,format='pdf')
	pp.close()

# run the whole series of steps through to graphing
def run_whole_script_for_group(pdfeatures,rFiles,label):
	collect = []
	allstream,allreverse = run_whole_analysis_for_boundaries(pdfeatures,'element')
	collect.append(allstream)
	collect.append(allreverse)
	for randomFile in rFiles:
		print 'collecting data frame for {0}'.format(randomFile)
		randomFeatures = collect_input_data_frame(randomFile)
		if labelcolumn:
			typefeatures = (randomFeatures[randomFeatures['type'] == label])
		else:
			typefeatures = randomFeatures
		ranstream,ranreverse = run_whole_analysis_for_boundaries(typefeatures,'random{0}'.format(randomFile))
		collect.append(ranstream)
		collect.append(ranreverse)
		print 'collected data frame for {0}'.format(randomFile)
	concat = pd.concat(collect)
	graph_methylation(concat,'{0}'.format(label))

def main():
	args = get_args()
	set_global_variables(args)
	directionFeatures = collect_input_data_frame(eFiles)
	print 'collected data frame for {0}'.format(eFiles)
	if labelcolumn:
		typeList = directionFeatures['type'].unique()
		for type in typeList:
			typefeatures = (directionFeatures[directionFeatures['type'] == type])
			run_whole_script_for_group(typefeatures,rFiles,type)
	else:
		run_whole_script_for_group(directionFeatures,rFiles,'All')

if __name__ == "__main__":
	main()