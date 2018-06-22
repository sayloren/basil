"""
Script to run mean AT conent analysis

Wren Saylor
April 2018

To Do:
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
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=str,help='A file containing a list of paths to the element files with unique names separated by newlines')
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
	parser.add_argument('-p',"--plotlinesize",type=int,default=3,help='size of the line to plot')
	parser.add_argument("-n", "--numberrandomassignments",default=100,type=int,help='the number of times to to randomly assign direction, will only be used when "--reversecomplement" is "random"')
	parser.add_argument('-c',"--reversecomplement",action='store_true',help='if you want the reverse complement to be plotted')
	return parser.parse_args()

# set all args that will be used throughout the script
def set_global_variables(args):
	global num
	global elementsize
	global inuce
	global window
	global binDir
	global halfwindow
	global fillX
	global eFiles
	global labelcolumn
	global directionalitycolumn
	global sizeGenome
	global faGenome
	global nucList
	global stringName
	global randomassignments
	global plotlinesize
	global reverseComplement
	num = args.total
	elementsize = args.element
	inuce = args.inset
	window = args.window
	binDir = args.bin
	halfwindow = ((window/2)+1)
	fillX = range(0,(num-window))
	eFiles = args.efile
	labelcolumn = args.labelcolumn
	directionalitycolumn = args.directionalitycolumn
	sizeGenome = args.genome
	faGenome = args.fasta
	nucList = ['A','T']
	stringName = args.stringname
	randomassignments = args.numberrandomassignments
	plotlinesize = args.plotlinesize
	reverseComplement = args.reversecomplement
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
def collect_coordinates_for_element_positions(btFeatures):
	midFeatures = pd.read_table(btFeatures.fn, header=None)
	midFeatures['middle'] = midFeatures.loc[:,1:2].mean(axis=1).astype(int)
	midFeatures['size'] = midFeatures.loc[:,2].astype(int)-midFeatures.loc[:,1].astype(int)
	if directionalitycolumn:
		midFeatures['directionality'] = midFeatures.loc[:,directionalitycolumn]
	if labelcolumn:
		midFeatures['type'] = midFeatures.loc[:,labelcolumn]
	midFeatures.insert(len(midFeatures.columns),'id',range(0,0+len(midFeatures)))
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
	rangeFeatures['elementString'] = rangeFeatures['sEdgeSeq'].astype(str) + rangeFeatures['MiddleSeq'].astype(str) + rangeFeatures['eEdgeSeq'].astype(str)
	rangeFeatures['elementString'] = rangeFeatures['elementString'].str.upper()
	rangeFeatures['flankString'] = rangeFeatures['sBoundarySeq'].astype(str) + rangeFeatures['eBoundarySeq'].astype(str)
	rangeFeatures['flankString'] = rangeFeatures['flankString'].str.upper()
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
	perSize.append(eval('100*float(start.count("A") + start.count("T"))/len(start)'))
	perSize.append(eval('100*float(end.count("A") + end.count("T"))/len(end)'))
	return perSize

# Directionality, as inferred by comparing first and last n base pairs from input parameters
def compare_boundaries_size_n(element,size):
# 	if motifdirectionality:
# 		perSize = locate_boundary_with_motif(element,size,motifdirectionality)
# 	else:
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

# wrapper for running sliding window and converting it into easy to use format
def sliding_window_wrapper(features,label):
	outCollect = run_sliding_window_for_each_nucleotide_string(features,label)
	outFlatten=flatten_data_from_sliding_window(outCollect)
	outDataFrame,names = convert_sliding_window_to_dataframe(outFlatten)
	print 'retrieved sliding window data for nucleotides strings {0}'.format(names)
	return outDataFrame,names

# get group, mean and standard deviation for AT
def collect_sum_two_nucleotides(dfWindow,names,nuc1,nuc2):
	ATNames = [names.index(i) for i in names if nuc1 in i or nuc2 in i]
	ATDataFrames = [dfWindow[i] for i in ATNames]
	ATconcat = pd.concat(ATDataFrames,axis=1)
	ATgroup = ATconcat.groupby(ATconcat.columns,axis=1).sum()
	ATmean = ATgroup.mean()
	ATstd = ATgroup.std()
	print 'summed a and t for graph'
	return ATgroup, ATmean, ATstd

# save panda
def save_panda(pdData, strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True)

# graph
def graph_element_line_means_random_below(dfWindow,names,fileName,denseRandom): # Extra
	set_ploting_parameters()
	ATgroup,ATmean,ATstd = collect_sum_two_nucleotides(dfWindow,names,'A','T')
	totalnumberelements = str(len(ATgroup.index))
	info = str(fileName) + ', '+ totalnumberelements + ' - ' "UCES"
	sns.set_style('ticks')
	gs = gridspec.GridSpec(1,1,height_ratios=[1])#,1
	gs.update(hspace=.8)
	pp = PdfPages('Fangs_{0}.pdf'.format(fileName))
	plt.figure(figsize=(10,10))
	plt.suptitle(info,fontsize=10)
	ax0 = plt.subplot(gs[0,:])
# 	ax1 = plt.subplot(gs[1,:],sharex=ax0)
	ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(denseRandom,names,'A','T')
	upstreamelement = ATgroup.loc[:,plotLineLocationThree:centerelementpoint].mean()
	downstreamelement = ATgroup.loc[:,centerelementpoint:plotLineLocationFour].mean()
	upstreamrandom = ranATgroup.loc[:,plotLineLocationThree:centerelementpoint].mean()
	downstreamrandom = ranATgroup.loc[:,centerelementpoint:plotLineLocationFour].mean()
	wilcoxonsignedrank = ss.wilcoxon(upstreamelement,downstreamelement)
	wilcoxonsignedrankrandom = ss.wilcoxon(upstreamrandom,downstreamrandom)
	statstable = pd.DataFrame([wilcoxonsignedrank,wilcoxonsignedrankrandom],
		columns=['statistic','pvalue'],
		index=['wsr-element','wsr-random'])
	save_panda(statstable,'Stats_ATcontent_{0}.txt'.format(fileName))

	datatable = pd.DataFrame([ATmean,ranATmean],index=['Element','Random'])
	save_panda(datatable.T,'Data_ATcontent_{0}.txt'.format(fileName))

	ax0.set_ylabel('% AT Content',size=16)
	ax0.set_xlabel('Position (bp)',size=16)
	plt.xlim(0,num)
	ax0.plot(fillX,ranATmean,linewidth=plotlinesize,label='Random',color='#d0abb0')
	ax0.plot(fillX,ATmean,linewidth=plotlinesize,label='Element',color='#4d5c82')
# 	ax1.set_xlabel('Position (bp)',size=16)
# 	ax1.set_ylabel('% AT Content',size=16)
# 	plt.setp(ax1.get_xticklabels(),visible=True)
	subplots = [ax0]# ax1
	for plot in subplots:
# 		plot.axvline(x=plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#c5969d')
# 		plot.axvline(x=plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#c5969d')
# 		plot.axvline(x=plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#c5969d')
# 		plot.axvline(x=plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#c5969d')
# 		plot.set_yticks(plot.get_yticks()[::2])
		plot.tick_params(axis='both',which='major',labelsize=20)
# 		plot.hlines(y=62,xmin=20,xmax=31,linewidth=.5,color='#081d58',zorder=0)
# 		plot.text(32,62,'{0}bp sliding window'.format(window),size=12)
	sns.despine()
	pp.savefig()
	pp.close()

# sliding window RCsorting
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

# separate and sort by plus and minus orientation
def sort_elements_by_directionality(directionFeatures,columnCompare):
	negStr = (directionFeatures[(directionFeatures[columnCompare] == '-')])
	posStr = (directionFeatures[(directionFeatures[columnCompare] == '+')])
	compWindow, compNames = sort_sliding_window_by_directionality(negStr,posStr)
	return compWindow,compNames

# get the empirical probability for each direction classification
def make_probabilites_for_direction(directionFeatures,probabilitycolumn):
	lenAll = float(len(directionFeatures.index))
	numPlus = (directionFeatures[probabilitycolumn] == '+').sum()/lenAll
	numMinus = (directionFeatures[probabilitycolumn] == '-').sum()/lenAll
	numEqual = (directionFeatures[probabilitycolumn] == '=').sum()/lenAll
	probOptions = [numMinus,numPlus,numEqual]
	print 'made probabilities for df: {0} for +, {1} for - and {2} for ='.format(numPlus,numMinus,numEqual)
	return probOptions

def main():
	starttime = time.time()
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
			for i in range(randomassignments):
				bool['randomDirectiontype'] = np.random.choice(dirOptions,len(bool.index),p=probOptionstype)
				typedirWindow,typedirNames = sort_elements_by_directionality(bool,'randomDirectiontype')
				spreadRandomtype.append(typedirWindow)
			denseRandomtype = sliding_window_df_to_collect_all_random(spreadRandomtype,boolNames)
			graph_element_line_means_random_below(boolWindow,boolNames,'{0}_{1}'.format(type,paramlabels),denseRandomtype)
	else:
		if reverseComplement:
			allWindow,allNames = sort_elements_by_directionality(directionFeatures,'directionality')
		else:
			allWindow,allNames = sliding_window_wrapper(directionFeatures['combineString'],directionFeatures['id'])
		spreadRandom,denseRandom=[],[]
		for i in range(randomassignments):
			directionFeatures['randomDirection'] = np.random.choice(dirOptions,len(directionFeatures.index),p=probOptions)
			randirWindow,randirNames = sort_elements_by_directionality(directionFeatures,'randomDirection')
			spreadRandom.append(randirWindow)
		denseRandom = sliding_window_df_to_collect_all_random(spreadRandom,allNames)
		graph_element_line_means_random_below(allWindow,allNames,'all_{0}'.format(paramlabels),denseRandom)
	endtime = time.time()
	print 'total time elapsed is {0}'.format(endtime-starttime)

if __name__ == "__main__":
	main()