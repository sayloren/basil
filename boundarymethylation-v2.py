"""
Script to run methylation analysis

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

To Do:
Get the number of C and G avaial to form CpG
Questions and appropriate stats
Other methylation questions
Which samples to use

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

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile",type=str,help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("-r","--randomfile",required=True,type=argparse.FileType('rU'),help="A file containing a list of paths to the random regions equable with your elements to plot in contrast")
	parser.add_argument("-m","--methylationfile",type=argparse.FileType('rU'),help="A file containing a list of paths to the methlylation bedfiles")

	parser.add_argument("-l", "--labelcolumn",type=int,help='column in the element file where the label (exonic, intergenic, intronic) is') # way to only read in certain entries, like only read in if 'intergenic'
	parser.add_argument("-d", "--directionalitycolumn",type=int,help='column in the element file where directionality is, if not supplied will infer by AT content')
	parser.add_argument("-i", "--idcolumn",type=int,help='column in the element file where the id is')

	parser.add_argument("-g","--genome",type=str, default="hg19.genome")
	parser.add_argument("-f","--fasta",type=str,default="hg19.fa")

	parser.add_argument("-p","--periphery",type=int,default="10",help='number bp from your boundary you want to include in the analysis')
	parser.add_argument("-b","--bin",type=int,default="100",help='size of bins used to compare element ends and determine directionality')
	parser.add_argument("-e","--element",type=int,help='size of your element (region without flanks), should be an even number, if not provided will use the smallest size of your input data')

	parser.add_argument("-mp","--methylationthresholdpercentage",type=int, default="10",help='size to threshold percentage methylation data')
	parser.add_argument("-mc","--methylationthresholdcoverage",type=int,default="10",help='size to threshold uncapped coverage of methylation data to send to percentage methylation, field often uses 10')

	parser.add_argument('-s',"--stringname",type=str,help='string to add to the outfile name')
	parser.add_argument('-c',"--reversecomplement",action='store_true',help='if you want the reverse complement to be plotted')
	return parser.parse_args()

# set all args that will be used throughout the script
def set_global_variables(args):
	global periphery
	global binDir
	global elementsize
	global methPerThresh
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
	periphery = args.periphery
	binDir = args.bin
	elementsize = args.element
	methPerThresh = args.methylationthresholdpercentage
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

# get bt features
def get_bedtools_features(strFileName):
	return pbt.BedTool(strFileName)

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
	rangeFeatures['upstreamsequence'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','startup','startdown']].values.astype(str).tolist()))
	rangeFeatures['upstreamsequence'] = rangeFeatures['upstreamsequence'].str.upper()
	rangeFeatures['downstreamsequence'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','endup','enddown']].values.astype(str).tolist()))
	rangeFeatures['downstreamsequence'] = rangeFeatures['downstreamsequence'].str.upper()
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
	pdmethThresh = (pdmethFeatures[(pdmethFeatures.loc[:,3] >= methCovThresh) & (pdmethFeatures.loc[:,4] >= methPerThresh)])
	btmethThresh = convert_panda_to_bed_format(pdmethThresh)
	print 'methylation coverage is being thresholded at {0} and percentage at {1}'.format(methCovThresh,methPerThresh)
	initiallength = len(pdmethFeatures.index)
	checklength = initiallength - len(pdmethThresh.index)
	print 'there were {0} out of {1} total removed by thresholding {2}'.format(checklength,initiallength,methName)
	return btmethThresh

# Intersect regions from the methylation data with element regions
def intersect_methylation_data(btFeatures,meFeature,meName,eStart,eEnd,sBound,eBound,column,boundary):
	initiallength = len(meFeature)
	intersectboundary = meFeature.intersect(btFeatures[['chr',eStart,eEnd,'id']].values.astype(str).tolist(),wb=True,wa=True)
	checklength = len(intersectboundary)
	print 'there were {0} intersections out of {1} in {2}'.format(checklength,initiallength,meName)
	if len(intersectboundary) != 0:
		pdmeth=convert_bedtool_to_panda(intersectboundary)
		pdmeth['int']=0
		pdmeth.columns = ['mchr','mstart','mstop','coverage','percentage','chr',sBound,eBound,'id','int']
		pdmeth['strand'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmeth[['mchr','mstart','mstop']].values.astype(str).tolist()))
		pdmeth['methlocation'] = pdmeth['int'].astype(int)+(pdmeth['mstart'].astype(int)-pdmeth[sBound].astype(int))
		outmeth = pdmeth[['chr','mstart','mstop',sBound,eBound,'id','percentage','methlocation','strand']]
		outmeth.columns = ['chr','mStart','mStop','eStart','eStop','id','percentage','methlocation','strand']
		outmeth['Tissue'] = meName.replace('.bed','')
		subfeatures = btFeatures[['id','directionality',column]]
		merge = pd.merge(outmeth,subfeatures,how='left',on='id')
		merge['methcount'] = len(merge.index)
		merge['boundary'] = boundary
		sortoutmeth = merge[['chr','id','directionality','methlocation','percentage','strand','methcount','boundary','Tissue']]
		sortoutmeth.columns = ['chr','id','directionality','methlocation','percentage','strand','methcount','boundary','Tissue']
	else:
		sortoutmeth = pd.DataFrame(np.nan,index=[0],columns=['chr','id','directionality','methlocation','percentage','strand','methcount'])
		sortoutmeth['boundary'] = boundary
		sortoutmeth['Tissue'] = meName.replace('.bed','')
	return sortoutmeth

# Run the analysis to extract methylation properties
def collect_methylation_data_by_element(pdfeatures):
	capture = []
	for methname in mFiles:
		mfeatures=get_bedtools_features(methname)
		pdthresh=threshold_methylation_data(mfeatures,methname)
		uintersect=intersect_methylation_data(pdfeatures,pdthresh,methname,'startup','startdown','sBoundary','sEdge','upstreamsequence','up stream')
		dintersect=intersect_methylation_data(pdfeatures,pdthresh,methname,'endup','enddown','eEdge','eBoundary','downstreamsequence','down stream')
		capture.append(uintersect)
		capture.append(dintersect)
	totalconcat = pd.concat(capture)
	totalconcat.reset_index(drop=True,inplace=True)
	return totalconcat

# Reverse those column that need to be reversed in case of - directionality
def negative_directionality_columns_to_modify(negFeatures):
	seqDict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
	originalRange = range(0,periphery*2)
	reverseRange = originalRange[::-1]
	rangeDict = dict(zip(originalRange,reverseRange))
	boundaryDict = {'up stream':'down stream','down stream':'up stream'}
	negFeatures['newmethlocation'] = negFeatures.loc[:,'methlocation'].map(rangeDict)
	negFeatures['newboundary'] = negFeatures.loc[:,'boundary'].map(boundaryDict)
	negFeatures['newstrand'] = negFeatures.loc[:,'strand'].map(seqDict)
	newnegFeatures = negFeatures[['chr','id','directionality','newmethlocation','percentage','newstrand','newboundary','Tissue','group']]
	newnegFeatures.columns = ['chr','id','directionality','methlocation','percentage','strand','boundary','Tissue','group']
	return newnegFeatures

# Separate on plus and minus orientation, rcsort and return methylation
def modify_negative_directionality_elements(btFeatures):
	if not btFeatures['directionality'].isnull().all():
		negFeatures = (btFeatures.loc[btFeatures['directionality']=='-'])
		nonFeatures = btFeatures.loc[btFeatures['directionality']!='-']
		newnegFeatures = negative_directionality_columns_to_modify(negFeatures)
		totalFeatures = pd.concat([newnegFeatures,nonFeatures])
		totalFeatures['methgroupby'] = totalFeatures.groupby(['Tissue','id','methlocation','boundary'])['Tissue'].transform('count')
		totalFeatures['methcount'] = totalFeatures.groupby(['methgroupby','Tissue','boundary'])['methgroupby'].transform('sum')
		outFeatures = totalFeatures[['chr','id','directionality','methlocation','percentage','strand','methcount','boundary','Tissue']]
	else:
		outFeatures = btFeatures
		print 'There are no features with "-" directionality, will use the original unsorted set for the analysis'
	return outFeatures

# get the total cpg in a column of strings
def collect_total_avail_cpg_in_column(btFeatures,column):
	btFeatures['cpgcount'] = btFeatures.apply(lambda row: float(row[column].count("CG")),axis=1)
	btFeatures['cgcount'] = btFeatures.apply(lambda row: (row[column].count("G")+row[column].count("C")),axis=1)
	return btFeatures['cpgcount'].sum(),btFeatures['cgcount'].sum()

# collect the values for up/down/rev up/rev down stream total cpg counts
def collect_total_avail_cpg_values(btFeatures):
	upcpg,upcg = collect_total_avail_cpg_in_column(btFeatures,'upstreamsequence')
	downcpg,downcg = collect_total_avail_cpg_in_column(btFeatures,'downstreamsequence')
	if not btFeatures['directionality'].isnull().all():
		negFeatures = (btFeatures.loc[btFeatures['directionality']=='-'])
		nonFeatures = btFeatures.loc[btFeatures['directionality']!='-']
		newnegFeatures = negFeatures.rename(columns={'upstreamsequence':'downstreamsequence','downstreamsequence':'upstreamsequence'})
		revFeatures = pd.concat([newnegFeatures,nonFeatures])
		revupcpg,revupcg = collect_total_avail_cpg_in_column(revFeatures,'upstreamsequence')
		revdowncpg,revdowncg = collect_total_avail_cpg_in_column(revFeatures,'downstreamsequence')
	else:
		revupcpg,revdowncpg,revupcg,revdowncg=upcpg,downcpg,upcg,downcg
		print 'There are no features with "-" directionality'
	return upcpg,downcpg,revupcpg,revdowncpg,upcg,downcg,revupcg,revdowncg

# make the original input data frame
def collect_input_data_frame(file):
	rangeFeatures = collect_element_coordinates(file)
	directionFeatures = assign_directionality_from_arg_or_boundary(rangeFeatures,file)
	return directionFeatures

# run the whole series of steps to make the data frame for each set of elements to plot
def run_whole_analysis_for_boundaries(pdfeatures,label):
	upcpg,downcpg,revupcpg,revdowncpg,upcg,downcg,revupcg,revdowncg = collect_total_avail_cpg_values(pdfeatures)
	allstream = collect_methylation_data_by_element(pdfeatures)
	allstream['group'] = label
	allstream['organization'] = 'unsorted'
	allstream['cpgsum'] = np.where(allstream['boundary'] == 'up stream',upcpg,downcpg)
	allstream['cgsum'] = np.where(allstream['boundary'] == 'up stream',upcg,downcg)
	allreverse = modify_negative_directionality_elements(allstream)
	allreverse['group'] = label
	allreverse['organization'] = 'rcsorted'
	allreverse['cpgsum'] = np.where(allreverse['boundary'] == 'up stream',revupcpg,revdowncpg)
	allstream['cgsum'] = np.where(allstream['boundary'] == 'up stream',revupcg,revdowncg)
	return allstream,allreverse

# remove duplicate counts by groupby params
def group_data_frame_by_column(pdfeatures,groupbycol,countcol):
	methsort = pdfeatures.sort_values(by=['group','boundary','Tissue'],ascending=True)
	methsort[countcol] = pdfeatures.groupby(['group','boundary','Tissue',groupbycol])[groupbycol].transform('count')
	methdup = methsort.drop_duplicates(['group','boundary','Tissue',groupbycol,countcol],keep='last')
	return methdup

# separate the elements and the random regions
def seperate_elements_and_random(pdfeatures,column,value):
	element = pdfeatures[pdfeatures[column]==value]
	random  = pdfeatures[pdfeatures[column]!=value]
	return element,random

# perform barplot for each set
def boxplot_params(removedups,yval,hval,pp,label):
	element,random = seperate_elements_and_random(removedups,'group','element')
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1],width_ratios=[1])
	gs.update(hspace=.8)
	ax0 = plt.subplot(gs[0,:])
	ax1 = plt.subplot(gs[1,:])
	sns.barplot(data=element,x='Tissue',y=yval,hue=hval,ax=ax0)
	sns.barplot(data=random,x='Tissue',y=yval,hue=hval,ax=ax1)
	ax0.set_title("UCEs")
	ax1.set_title("Random Regions")
	subplots = [ax0,ax1]
	for plot in subplots:
		plot.tick_params(axis='both',which='major',labelsize=16)
		plot.set_ylabel(label,size=16)
		plot.set_xlabel("Tissue",fontsize=16)
	sns.despine()
	plt.savefig(pp,format='pdf')

# Make graphs for fangs
def graph_boundary_methylation(pdfeatures,filelabel):
	info = str(filelabel)
	sns.set_style('ticks')
	if reverseComplement:
		sorted = pdfeatures.loc[pdfeatures['organization']=='rcsorted']
		pp = PdfPages('Methylation_ReverseComplementSorted_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.pdf'.format(eFiles,stringName,filelabel,elementsize,binDir,periphery,methPerThresh,methCovThresh))
	else:
		sorted = pdfeatures.loc[pdfeatures['organization']=='unsorted']
		pp = PdfPages('Methylation_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.pdf'.format(eFiles,stringName,filelabel,elementsize,binDir,periphery,methPerThresh,methCovThresh))
	plt.figure(figsize=(14,7))
	plt.suptitle(info,fontsize=16)
	surroundingfang = periphery*2
	sns.set_palette("Blues")
	# groupby and count different params
	removedupcpgper = group_data_frame_by_column(sorted,'methcount','methcounttest')
	removedupstrand = group_data_frame_by_column(sorted,'strand','strandcount')
	removedupdir = group_data_frame_by_column(sorted,'directionality','dircount')
	removeduploc = group_data_frame_by_column(sorted,'methlocation','loccount')
	removedupid = group_data_frame_by_column(sorted,'id','idcount')
	removedupavailcg = group_data_frame_by_column(sorted,'createdcpg','createdcpgtest')
	# graph cpgs methylated
# 	boxplot_params(removedupcpgper,'percpgmeth','boundary',pp,'% CpGs Methylated')
	boxplot_params(removedupcpgper,'methcount','boundary',pp,'Count CpGs Methylated')
	boxplot_params(removedupcpgper,'methcount','strand',pp,'Count CpGs Methylated')
	boxplot_params(removedupcpgper,'methcount','directionality',pp,'Count CpGs Methylated')
	# graph strand cpgs methylated
# 	boxplot_params(removedupstrand,'strandcount','boundary',pp,'Count Methylation Strand')
	boxplot_params(removedupstrand,'strandcount','strand',pp,'Count Methylation Strand')
# 	boxplot_params(removedupstrand,'strandcount','directionality',pp,'Count Methylation Strand')
	# graph percentage
	boxplot_params(sorted,'percentage','boundary',pp,'Count % Methylation')
	boxplot_params(sorted,'percentage','strand',pp,'Count % Methylation')
	boxplot_params(sorted,'percentage','directionality',pp,'Count % Methylation')
	# graph directionality
# 	boxplot_params(removedupdir,'dircount','boundary',pp,'Count Directionality')
# 	boxplot_params(removedupdir,'dircount','strand',pp,'Count Directionality')
	boxplot_params(removedupdir,'dircount','directionality',pp,'Count Directionality')
	# graph id
	boxplot_params(removedupid,'idcount','id',pp,'Count UCE ID Methylated') # needs work
	# graph location
	boxplot_params(removeduploc,'methlocation','loccount',pp,'Count Methylation Location') # needs work
	# graph avail c and g to cpg
	boxplot_params(removedupavailcg,'createdcpg','boundary',pp,'% C and G creating CpG')
	pp.close()

# run the whole series of steps through to graphing
def run_whole_script_for_group(pdfeatures,rFiles,label):
	collect = []
	allstream,allreverse = run_whole_analysis_for_boundaries(pdfeatures,'element')
	collect.append(allstream)
	collect.append(allreverse)
	for randomFile in rFiles:
		randomFeatures = collect_input_data_frame(randomFile)
		ranstream,ranreverse = run_whole_analysis_for_boundaries(randomFeatures,'random{0}'.format(randomFile))
		collect.append(ranstream)
		collect.append(ranreverse)
	concat = pd.concat(collect)
	concat['percpgmeth'] = (concat['methcount']/concat['cpgsum'])*100.0
	concat['createdcpg'] = (concat['cpgsum']/concat['cgsum'])*100.0
	concat.reset_index(drop=True,inplace=True)
	graph_boundary_methylation(concat,'{0}'.format(label))

def main():
	args = get_args()
	set_global_variables(args)
	directionFeatures = collect_input_data_frame(eFiles)
	run_whole_script_for_group(directionFeatures,rFiles,'All')
	if labelcolumn:
		typeList = directionFeatures['type'].unique()
		for type in typeList:
			typefeatures = (directionFeatures[directionFeatures['type'] == type])
			run_whole_script_for_group(typefeatures,rFiles,type)

if __name__ == "__main__":
	main()