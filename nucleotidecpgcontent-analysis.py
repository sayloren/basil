"""
Script to run mean CpG conent analysis

Wren Saylor
October 14 2017

To Do:
unittest
change directionality arg column
streamline and remove excess columns
.index len to break up large datasets

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
import math
from scipy.stats import chisquare
from scipy.interpolate import splrep, splev
import scipy.stats as ss
from scipy.stats import mstats
import time

# set command line arguments
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=str,help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("-r","--randomfile",type=argparse.FileType('rU'),help="A file containing a list of paths to the random regions equable with your elements to plot in contrast") # option to not use random at all, or generate randoms?
	parser.add_argument("-l", "--labelcolumn",type=int,help='column in the element file where the label (exonic, intergenic, intronic) is')
	parser.add_argument("-d", "--directionalitycolumn",type=int,help='column in the element file where directionality is, if not supplied will infer by AT content')
	parser.add_argument("-g","--genome",type=str, default="hg19.genome")
	parser.add_argument("-f","--fasta",type=str,default="hg19.fa")
	parser.add_argument("-t","--total",type=int,default="2200",help='total size of region to look at (region + flanks), should be an even number, suggested to be at least triple your element')
	parser.add_argument("-e","--element",type=int,help='size of your element (region without flanks), should be an even number, if not provided will use the smallest size of your input data')
	parser.add_argument("-i","--inset",type=int,default="50",help='size into your element from the boundaries, should be an even number')
	parser.add_argument("-w","--window",type=int,default="11",help='size of sliding window, should be an odd number, previous studies have used 11')
	parser.add_argument("-b","--bin",type=int,default="100",help='size of bins used to compare element ends and determine directionality')
	parser.add_argument('-s',"--stringname",type=str,help='string to add to the outfile name')
	parser.add_argument('-p',"--plotlinesize",type=int,default=2,help="size of the line to plot")
	parser.add_argument('-c',"--reversecomplement",action='store_true',help='if you want the reverse complement to be plotted')
	return parser.parse_args()

# set all args that will be used throughout the script
def set_global_variables(args):
	# Integer parameters
	global num
	global elementsize
	global inuce
	global window
	global binDir
	global halfwindow
	global fillX
	global eFiles
	global rFiles
	global labelcolumn
	global directionalitycolumn
	global sizeGenome
	global faGenome
	global nucList
	global stringName
	global reverseComplement
	global plotlinesize
	num = args.total
	elementsize = args.element
	inuce = args.inset
	window = args.window
	binDir = args.bin
	halfwindow = ((window/2)+1)
	fillX = range(0,(num-window))
	eFiles = args.efile
	rFiles = [line.strip() for line in args.randomfile]
	labelcolumn = args.labelcolumn
	directionalitycolumn = args.directionalitycolumn
	sizeGenome = args.genome
	faGenome = args.fasta
	nucList = ['CG','G','C']
	stringName = args.stringname
	reverseComplement = args.reversecomplement
	plotlinesize = args.plotlinesize

def set_ploting_parameters():
	global plotLineLocationOne # upstream element boundary
	global plotLineLocationTwo # downstream element boundary
	global plotLineLocationThree # upstream element inset
	global plotLineLocationFour # downstream element inset
	plotLineLocationOne = (((num-uce)/2)+(inuce-halfwindow))
	plotLineLocationTwo = (((num-uce)/2)+(uce-inuce-halfwindow))
	plotLineLocationThree = (((num-uce)/2)-halfwindow)
	plotLineLocationFour = (((num-uce)/2)+uce-halfwindow)

# get bt features
def get_bedtools_features(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

# get the correct range for fang evaluation
def collect_coordinates_for_element_positions(btFeatures):
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
	midFeatures['sCenter'] = midFeatures['middle'].astype(int) - (inregion/2)
	midFeatures['eCenter'] = midFeatures['middle'].astype(int) + (inregion/2)
	midFeatures['sEdge'] = midFeatures['start'].astype(int) + inuce
	midFeatures['eEdge'] = midFeatures['end'].astype(int) - inuce
	midFeatures['sBoundary'] = midFeatures['start'].astype(int) - flankSize
	midFeatures['eBoundary'] = midFeatures['end'].astype(int) + flankSize
	midFeatures=midFeatures.drop(['middle'],axis=1)
	if directionalitycolumn:
		midFeatures['directionality'] = midFeatures.loc[:,directionalitycolumn]
	if labelcolumn:
		midFeatures['type'] = midFeatures.loc[:,labelcolumn]
	midFeatures.insert(len(midFeatures.columns),'id',range(0,0+len(midFeatures)))
	print 'collected the coordinates to create the sequence strings'
	return midFeatures

# get the strings for sliding window regions
def get_fasta_for_element_coordinates(rangeFeatures):
	rangeFeatures['sBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','sBoundary','start']].values.astype(str).tolist()))
	rangeFeatures['sEdgeSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','start','sEdge']].values.astype(str).tolist()))
	rangeFeatures['MiddleSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','sCenter','eCenter']].values.astype(str).tolist()))
	rangeFeatures['eEdgeSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','eEdge','end']].values.astype(str).tolist()))
	rangeFeatures['eBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','end','eBoundary']].values.astype(str).tolist()))
	rangeFeatures['combineString'] = rangeFeatures['sBoundarySeq'].astype(str) + rangeFeatures['sEdgeSeq'].astype(str) + rangeFeatures['MiddleSeq'].astype(str) + rangeFeatures['eEdgeSeq'].astype(str) + rangeFeatures['eBoundarySeq'].astype(str)
	rangeFeatures['combineString'] = rangeFeatures['combineString'].str.upper()
	rangeFeatures=rangeFeatures.drop(['size','sBoundarySeq','sEdgeSeq','MiddleSeq','eEdgeSeq','eBoundarySeq','sBoundary','sEdge','sCenter','eCenter','eEdge','eBoundary'],axis=1)
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
	subsetFeatures = collect_coordinates_for_element_positions(btFeatures)
	rangeFeatures = get_fasta_for_element_coordinates(subsetFeatures)
	return rangeFeatures

# Do the comparison between boundaries to get + - or =
def calculate_nucleotides_at(element,size):
	start = element[:size]
	end = element[-size:]
	perSize = []
	perSize.append(eval('100*float(start.count("A") + start.count("a") + start.count("T") + start.count("t"))/len(start)'))
	perSize.append(eval('100*float(end.count("A") + end.count("a") + end.count("T") + end.count("t"))/len(end)'))
	return perSize

# Directionality, as inferred by comparing first and last n base pairs from input parameters
def compare_boundaries_size_n(element,size):
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
		print 'sorting the element boundaries by bin size {0}'.format(binDir)
	rangeFeatures=rangeFeatures.drop(['chr','start','end'],axis=1)
	return rangeFeatures

# run the sliding window for each nucleotide string
def run_sliding_window_for_each_nucleotide_string(features,label):
	outCollect = []
	for element,id in zip(features,label):
		outElement = {id: []}
		outList = {key:[] for key in nucList}
		n = num
		s = 1 # size to jump for sliding window
		start, end = 0, window
		while end < n:
			current = element[start:end]
			for key in nucList:
				percentage = float(100*current.count(key)/len(current))
				outList[key].append(percentage)
			start, end = start + s, end + s
		outElement[id].append(outList)
		outCollect.append(outElement)
	print 'ran the sliding window'
	return outCollect

# convert to a data frame with a list for each element x nucleotide string
def flatten_data_from_sliding_window(outCollect):
	outFlatten = pd.DataFrame()
	outIndex = []
	for d in outCollect:
		for key,values in d.items():
			outFlatten = outFlatten.append(values,ignore_index=True)
			outIndex.append(key)
	outFlatten.index = outIndex
	print 'flattend sliding window data to df'
	return outFlatten

# turn each list of element x nucleotide string into a separate df, within a larger df
def convert_sliding_window_to_dataframe(outFlatten):
	names = outFlatten.columns.tolist()
	collectNucDF = []
	for nuc in names:
		nuc = outFlatten[[nuc]]
		nuc.columns = ['temp']
		split = pd.DataFrame(nuc['temp'].values.tolist(),index=nuc.index)
		collectNucDF.append(split)
	print 'converted sliding window to df'
	return collectNucDF,names

# Put all the random regions for each string search into one df, for stat collection
def sliding_window_df_to_collect_all_random(collectRandom,allNames):
	dictNames = {key:[] for key in allNames}
	for random in collectRandom:
		for df,key in zip(random,dictNames):
			dictNames[key].append(df)
	catCollect = []
	for key,values in dictNames.items():
		cat = pd.concat(values)
		catCollect.append(cat)
	print 'combined random sliding windows'
	return catCollect

# Wrapper for running sliding window and converting it into easy to use format
def sliding_window_wrapper(features,label):
	outCollect = run_sliding_window_for_each_nucleotide_string(features,label)
	outFlatten=flatten_data_from_sliding_window(outCollect)
	outDataFrame,names = convert_sliding_window_to_dataframe(outFlatten)
	print 'retrieved sliding window data for nucleotides strings {0}'.format(names)
	return outDataFrame,names

# Get group, mean and standard deviation for AT
def collect_linear_two_nucleotides(dfWindow,names,nuc1):
	CpGNames = [names.index(i) for i in names if i == nuc1]
	CpGDataFrames = [dfWindow[i] for i in CpGNames]
	CpGconcat = pd.concat(CpGDataFrames,axis=1)
	CpGgroup = CpGconcat.groupby(CpGconcat.columns,axis=1).sum()
	CpGmean = CpGgroup.mean()
	CpGstd = CpGgroup.std()
	return CpGgroup, CpGmean, CpGstd

# save panda
def save_panda(pdData, strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True)

# Make graphs for fangs
def print_element_line_means(dfWindow,names,fileName,denseRandom):
	set_ploting_parameters()
	CpGgroup,CpGmean,CpGstd = collect_linear_two_nucleotides(dfWindow,names,'CG')
	Cgroup,Cmean,Cstd = collect_linear_two_nucleotides(dfWindow,names,'C')
	Ggroup,Gmean,Gstd = collect_linear_two_nucleotides(dfWindow,names,'G')
	totalCGmean = Cmean+Gmean
	normCpGmean = CpGmean/totalCGmean
	totalnumberelements = str(len(CpGgroup.index))
	ranCpGgroup,ranCpGmean,ranCpGstd = collect_linear_two_nucleotides(denseRandom,names,'CG')
	ranCgroup,ranCmean,ranCstd = collect_linear_two_nucleotides(denseRandom,names,'C')
	ranGgroup,ranGmean,ranGstd = collect_linear_two_nucleotides(denseRandom,names,'G')
	rantotalCGmean = ranCmean+ranGmean
	rannormCpGmean = ranCpGmean/rantotalCGmean
	datatable = pd.DataFrame([normCpGmean,rannormCpGmean],index=['Element','Random'])
	save_panda(datatable.T,'Data_CpGcontent_{0}.txt'.format(fileName))
	wilcoxonsignedrank = ss.wilcoxon(CpGmean,ranCpGmean)
	statstable = pd.DataFrame([wilcoxonsignedrank],
		columns=['statistic','pvalue'],
		index=['wsr'])
	save_panda(statstable,'Stats_CpGcontent_{0}.txt'.format(fileName))

# For type groups, separate the groups and run the analyses
def separate_dataframe_by_group(listgroup,directionFeatures,typecolumn,fileName):
	bool = (directionFeatures[directionFeatures[typecolumn] == listgroup])
	dfWindow,Names = sliding_window_wrapper(bool['combineString'],bool['id'])
	print 'separated df by type'
	return bool,dfWindow,Names

# Sliding window RCsorting
def sort_sliding_window_by_directionality(negStr,posStr):
	negStr['reverseComplement']=negStr.apply(lambda row: reverse_complement_dictionary(row['combineString']),axis=1)
	negDF,negNames=sliding_window_wrapper(negStr['reverseComplement'],negStr['id'])
	posDF,posNames=sliding_window_wrapper(posStr['combineString'],posStr['id'])
	compWindow=[]
	for x, y in zip(negDF,posDF):
		tempCat=pd.concat([x,y],axis=1)
		tempGroup=tempCat.groupby(tempCat.columns,axis=1).sum()
		compWindow.append(tempGroup)
	print 'ran sliding window for directionality sorted df'
	return compWindow,negNames

# Separate and sort by plus and minus orientation
def sort_elements_by_directionality(directionFeatures,columnCompare):
	negStr = (directionFeatures[(directionFeatures[columnCompare] == '-')])
	posStr = (directionFeatures[(directionFeatures[columnCompare] == '+')])
	compWindow, compNames = sort_sliding_window_by_directionality(negStr,posStr)
	return compWindow,compNames

# Get the empirical probability for each direction classification
def make_probabilites_for_direction(directionFeatures,probabilitycolumn):
	lenAll = float(len(directionFeatures.index))
	numPlus = (directionFeatures[probabilitycolumn] == '+').sum()/lenAll
	numMinus = (directionFeatures[probabilitycolumn] == '-').sum()/lenAll
	numEqual = (directionFeatures[probabilitycolumn] == '=').sum()/lenAll
	probOptions = [numMinus,numPlus,numEqual]
	print 'made probabilities for df: {0} for +, {1} for - and {2} for ='.format(numPlus,numMinus,numEqual)
	return probOptions

def main():
	args = get_args()
	set_global_variables(args)
	paramlabels = '{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(elementsize,inuce,num,binDir,window,eFiles,stringName)
	rangeFeatures = collect_element_coordinates(eFiles)
	directionFeatures = assign_directionality_from_arg_or_boundary(rangeFeatures,eFiles)
	dirOptions = ['-','+','=']
	probOptions = make_probabilites_for_direction(directionFeatures,'directionality')
	if labelcolumn:
		typeList = directionFeatures['type'].unique()
		for type in typeList:
			print 'Now running {0} elements'.format(type)
			bool = (rangeFeatures[rangeFeatures['type'] == type])
			if reverseComplement:
				boolWindow,boolNames = sort_elements_by_directionality(bool,'directionality')
			else:
				boolWindow,boolNames = sliding_window_wrapper(bool['combineString'],bool['id'])
			probOptionstype = make_probabilites_for_direction(bool,'directionality')
			spreadRandomtype,denseRandomtype=[],[]
			for r in rFiles:
				randomtypeFeatures = collect_element_coordinates(r)
				directiontypeRandom = assign_directionality_from_arg_or_boundary(randomtypeFeatures,r)
				if reverseComplement:
					typedirWindow,typedirNames = sort_elements_by_directionality(directiontypeRandom,'directionality')
				else:
					typedirWindow,typedirNames = sliding_window_wrapper(directiontypeRandom['combineString'],directiontypeRandom['id'])
# 				bool['randomDirectiontype'] = np.random.choice(dirOptions,len(bool.index),p=probOptionstype)
# 				typedirWindow,typedirNames = sort_elements_by_directionality(bool,'randomDirectiontype')
				spreadRandomtype.append(typedirWindow)
			denseRandomtype = sliding_window_df_to_collect_all_random(spreadRandomtype,boolNames)
			print_element_line_means(boolWindow,boolNames,'{0}_{1}'.format(type,paramlabels),denseRandomtype)
	else:
		if reverseComplement:
			allWindow,allNames = sort_elements_by_directionality(directionFeatures,'directionality')
		else:
			allWindow,allNames = sliding_window_wrapper(directionFeatures['combineString'],directionFeatures['id'])
		spreadRandom,denseRandom=[],[]
		for r in rFiles:
			randomFeatures = collect_element_coordinates(r)
			directionRandom = assign_directionality_from_arg_or_boundary(randomFeatures,r)
			if reverseComplement:
				randirWindow,randirNames = sort_elements_by_directionality(directionRandom,'directionality')
			else:
				randirWindow,randirNames = sliding_window_wrapper(directionRandom['combineString'],directionRandom['id'])
# 			directionFeatures['randomDirection'] = np.random.choice(dirOptions,len(directionFeatures.index),p=probOptions)
# 			randirWindow,randirNames = sort_elements_by_directionality(directionFeatures,'randomDirection')
			spreadRandom.append(randirWindow)
		denseRandom = sliding_window_df_to_collect_all_random(spreadRandom,allNames)
		print_element_line_means(allWindow,allNames,'all_{0}'.format(paramlabels),denseRandom)

if __name__ == "__main__":
	main()