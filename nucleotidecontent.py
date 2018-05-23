"""
Script to run mean nucleotide conent analysis

Wren Saylor
April 2018

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
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=str,help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("-l", "--labelcolumn",type=int,help='column in the element file where the label (exonic, intergenic, intronic) is')
	parser.add_argument("-g","--genome",type=str, default="hg19.genome")
	parser.add_argument("-f","--fasta",type=str,default="hg19.fa")
	parser.add_argument("-t","--total",type=int,default="2200",help='total size of region to look at (region + flanks), should be an even number, suggested to be at least triple your element')
	parser.add_argument("-e","--element",type=int,help='size of your element (region without flanks), should be an even number, if not provided will use the smallest size of your input data')
	parser.add_argument("-i","--inset",type=int,default="50",help='size into your element from the boundaries, should be an even number')
	parser.add_argument("-w","--window",type=int,default="11",help='size of sliding window, should be an odd number, previous studies have used 11')
	parser.add_argument('-s',"--stringname",type=str,help='string to add to the outfile name')
	parser.add_argument('-p',"--plotlinesize",type=int,default=2,help='size of the line to plot')
	return parser.parse_args()

# set all args that will be used throughout the script
def set_global_variables(args):
	global num
	global elementsize
	global inuce
	global window
	global halfwindow
	global fillX
	global eFiles
	global labelcolumn
	global sizeGenome
	global faGenome
	global nucList
	global stringName
	global plotlinesize
	num = args.total
	elementsize = args.element
	inuce = args.inset
	window = args.window
	halfwindow = ((window/2)+1)
	fillX = range(0,(num-window))
	eFiles = args.efile
	labelcolumn = args.labelcolumn
	sizeGenome = args.genome
	faGenome = args.fasta
	nucList = ['A','C','G','T']
	stringName = args.stringname
	plotlinesize = args.plotlinesize
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

# wrapper for running sliding window and converting it into easy to use format
def sliding_window_wrapper(features,label):
	outCollect = run_sliding_window_for_each_nucleotide_string(features,label)
	outFlatten=flatten_data_from_sliding_window(outCollect)
	outDataFrame,names = convert_sliding_window_to_dataframe(outFlatten)
	print 'retrieved sliding window data for nucleotides strings {0}'.format(names)
	return outDataFrame,names

# get nucleotide mean
def get_nuc_group(dfwindow,names,nuc):
	dflocation = [names.index(i) for i in names if nuc in i]
	dfselect = [dfwindow[i] for i in dflocation]
	dfconcat = pd.concat(dfselect,axis=1)
	dfgroup = dfconcat.groupby(dfconcat.columns,axis=1).sum()
	totalnumberelements = str(len(dfgroup.index))
	dfmean = dfgroup.mean()
	return dfmean,totalnumberelements

# get means for all individual nucleotides
def collect_nucleotide_mean(dfwindow,names):
	Amean,Atotal = get_nuc_group(dfwindow,names,'A')
	Cmean,Ctotal = get_nuc_group(dfwindow,names,'C')
	Gmean,Gtotal = get_nuc_group(dfwindow,names,'G')
	Tmean,Ttotal = get_nuc_group(dfwindow,names,'T')
	return Amean,Cmean,Gmean,Tmean,Ttotal

# save panda
def save_panda(pdData, strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True)

# graph
def graph_element_line_means_random_below(dfWindow,names,fileName):
	set_ploting_parameters()
	Amean,Cmean,Gmean,Tmean,totalnumberelements = collect_nucleotide_mean(dfWindow,names)
	info = str(fileName) + ', '+ totalnumberelements + ' - ' "UCES"
	sns.set_style('ticks')
	gs = gridspec.GridSpec(1,1,height_ratios=[1])#,1
	gs.update(hspace=.8)
	pp = PdfPages('Nuc_{0}.pdf'.format(fileName))
	plt.figure(figsize=(14,7))
	plt.suptitle(info,fontsize=10)
	ax0 = plt.subplot(gs[0,:])
# 	upstreamelement = ATgroup.loc[:,plotLineLocationThree:centerelementpoint].mean()
# 	downstreamelement = ATgroup.loc[:,centerelementpoint:plotLineLocationFour].mean()
# 	wilcoxonsignedrank = ss.wilcoxon(upstreamelement,downstreamelement)
# 	wilcoxonsignedrankrandom = ss.wilcoxon(upstreamrandom,downstreamrandom)
# 	statstable = pd.DataFrame([wilcoxonsignedrank,wilcoxonsignedrankrandom]
# 		columns=['statistic','pvalue'],
# 		index=['wsr-element','wsr-random'])
# 	save_panda(statstable,'Stats_{0}.txt'.format(fileName))
	ax0.plot(fillX,Amean,linewidth=plotlinesize,label='A',color='#f4c88b')#ff0045
	ax0.plot(fillX,Cmean,linewidth=plotlinesize,label='C',color='#a47288')#db00cc
	ax0.plot(fillX,Gmean,linewidth=plotlinesize,label='G',color='#ae9ea1')#6b4dff
	ax0.plot(fillX,Tmean,linewidth=plotlinesize,label='T',color='#ff94b8')#ab00ff
	ax0.set_ylabel('% Nucleotide Content',size=16)
	ax0.set_xlabel('Position (bp)',size=16)
	ax0.legend()
	plt.xlim(0,num)
	subplots = [ax0]
	for plot in subplots:
		plot.axvline(x=plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#c5969d')
		plot.axvline(x=plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#c5969d')
		plot.axvline(x=plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#c5969d')
		plot.axvline(x=plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#c5969d')
		plot.set_yticks(plot.get_yticks()[::2])
		plot.tick_params(axis='both',which='major',labelsize=16)
		plot.hlines(y=34,xmin=20,xmax=31,linewidth=.5,color='#081d58',zorder=0)
		plot.text(32,34,'{0}bp sliding window'.format(window),size=12)
	sns.despine()
	pp.savefig()
	pp.close()

def main():
	args = get_args()
	set_global_variables(args)
	paramlabels = '{0}_{1}_{2}_{3}_{4}_{5}'.format(elementsize,inuce,num,window,eFiles,stringName)
	rangeFeatures = collect_element_coordinates(eFiles)
	if labelcolumn:
		typeList = rangeFeatures['type'].unique()
		for type in typeList:
			print 'Now running {0} elements'.format(type)
			bool = (rangeFeatures[rangeFeatures['type'] == type])
			boolWindow,boolNames = sliding_window_wrapper(bool['combineString'],bool['id'])
			graph_element_line_means_random_below(boolWindow,boolNames,'{0}_{1}'.format(type,paramlabels))
	else:
		allWindow,allNames = sliding_window_wrapper(rangeFeatures['combineString'],rangeFeatures['id'])
		graph_element_line_means_random_below(allWindow,allNames,'all_{0}'.format(paramlabels))

if __name__ == "__main__":
	main()
