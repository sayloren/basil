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
	parser.add_argument("-lc", "--labelcolumn",type=int,help='column in the element file where the label (exonic, intergenic, intronic) is') # way to only read in certain entries, like only read in if 'intergenic'
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
	parser.add_argument('-t',"--twographs",action='store_true',help='if you want to see the upstream and downstream boundaries separately')
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
	global reverseComplement
	stringName = args.stringname
	splitgraphs = args.twographs
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

# get the reverse complement
def reverse_complement_dictionary(sequence):
	seqDict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
	return "".join([seqDict[base] for base in reversed(sequence)])

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
# 		rangeFeatures=rangeFeatures.drop(['feature'],axis=1)
		numval = rangeFeatures['directionality'].value_counts()
	print '{0} in {1}'.format(numval,fileName)
	return rangeFeatures

# save file from bedtool
def save_bedtool_as_bedfile(btObject,strFilename):
	btObject.saveas(strFilename)

# convert bedtool to panda
def convert_bedtool_to_panda(btobject):
	save_bedtool_as_bedfile(btobject,'temp.bed')
	pdObject = pd.read_table(btobject.fn, header=None)
	return pdObject

# convert panda to bedtool
def convert_panda_to_bed_format(panda):
	arArFeatures = panda.values.tolist()
	btoutFeatures = get_bedtools_features(arArFeatures)
	return btoutFeatures

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
def intersect_methylation_data_by_element(btFeatures,meFeature,meName,eStart,eEnd,sBound,eBound):
	initiallength = len(meFeature)
	intersectboundary = meFeature.intersect(btFeatures[['chr',eStart,eEnd,'id']].values.astype(str).tolist(),wb=True,wa=True)
	checklength = len(intersectboundary)
	print 'there were {0} intersections out of {1} in {2}'.format(checklength,initiallength,meName)
	if len(intersectboundary) != 0:
		pdmeth=convert_bedtool_to_panda(intersectboundary)
		pdmeth['int']=0
		pdmeth.columns = ['mchr','mstart','mstop','coverage','percentage','chr',sBound,eBound,'id','int']
		pdmeth['cytosine'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmeth[['mchr','mstart','mstop']].values.astype(str).tolist()))
		pdmeth['methlocation'] = pdmeth['int'].astype(int)+(pdmeth['mstart'].astype(int)-pdmeth[sBound].astype(int))
		outmethylation = pdmeth[['chr','mstart','mstop',sBound,eBound,'id','percentage','methlocation','cytosine']]
		outmethylation.columns = ['chr','mStart','mStop','eStart','eStop','id','percentage','methlocation','cytosine']
		sortoutmethylation = outmethylation.sort_values(['methlocation'],ascending=True)
	else:
		sortoutmethylation = None
	return sortoutmethylation

# Combine directionality and id information to the methylation data, format
def make_methylation_data_frame_and_verify_cytosine(methdf,name,features,column,boundary):
	if methdf is not None:
		methdf['Tissue'] = name.replace('.bed','')
		subfeatures = features[['id','directionality',column]]
		merge = pd.merge(methdf,subfeatures,how='left',on='id')
		merge['methcount'] = len(merge.index)
		merge['boundary'] = boundary
		subMeth = merge[['chr','id','directionality','methlocation','percentage','cytosine','methcount','boundary','Tissue']]
		subMeth.columns = ['chr','id','directionality','methlocation','percentage','cytosine','methcount','boundary','Tissue']
	else:
		subMeth = pd.DataFrame(np.nan,index=[0],columns=['chr','id','directionality','methlocation','percentage','cytosine','methcount','boundary'])
		subMeth['Tissue'] = name.replace('.bed','')
	return subMeth

# Run the analysis to extract methylation properties
def collect_methylation_data_by_element(pdfeatures):
	capture = []
	for methname in mFiles:
		mfeatures=get_bedtools_features(methname)
		pdthresh=threshold_methylation_data(mfeatures,methname)
		uintersect=intersect_methylation_data_by_element(pdfeatures,pdthresh,methname,'startup','startdown','sBoundary','sEdge')
		dintersect=intersect_methylation_data_by_element(pdfeatures,pdthresh,methname,'endup','enddown','eEdge','eBoundary')
		umerge=make_methylation_data_frame_and_verify_cytosine(uintersect,methname,pdfeatures,'upstreamsequence','up stream')
		dmerge=make_methylation_data_frame_and_verify_cytosine(dintersect,methname,pdfeatures,'downstreamsequence','down stream')
		capture.append(umerge)
		capture.append(dmerge)
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
	negFeatures['newmethlocation'] = negFeatures.methlocation.map(rangeDict)
	negFeatures['newcytosine'] = negFeatures.cytosine.map(seqDict)
	negFeatures['newboundary'] = negFeatures.boundary.map(boundaryDict)
	newnegFeatures = negFeatures[['chr','id','directionality','newmethlocation','percentage','newcytosine','newboundary','Tissue','group']]
	newnegFeatures.columns = ['chr','id','directionality','methlocation','percentage','cytosine','boundary','Tissue','group']
	return newnegFeatures

# Separate on plus and minus orientation, rcsort and return methylation
def modify_negative_directionality_elements(btFeatures):
	negFeatures = btFeatures.loc[btFeatures['directionality']=='-']
	if not negFeatures.empty:
		nonFeatures = btFeatures.loc[btFeatures['directionality']!='-']
		newnegFeatures = negative_directionality_columns_to_modify(negFeatures)
		totalFeatures = pd.concat([newnegFeatures,nonFeatures])
		totalFeatures['methgroupby'] = totalFeatures.groupby(['Tissue','id','methlocation','boundary'])['Tissue'].transform('count')
		totalFeatures['methcount'] = totalFeatures.groupby(['methgroupby','Tissue','boundary'])['methgroupby'].transform('sum')
		outFeatures = totalFeatures[['chr','id','directionality','methlocation','percentage','cytosine','methcount','boundary','Tissue']]
	else:
		outFeatures = btFeatures
		print 'There are no features with "-" directionality'
	return outFeatures

# get the total cpg in a column of strings
def collect_total_avail_cpg_in_column(btFeatures,column):
	btFeatures['cpgcount'] = btFeatures.apply(lambda row: float(row[column].count("CG")),axis=1)
	cpgsum = btFeatures['cpgcount'].sum()
	return cpgsum

# collect the values for up/down/rev up/rev down stream total cpg counts
def collect_total_avail_cpg_values(btFeatures):
	upcpgcount = collect_total_avail_cpg_in_column(btFeatures,'upstreamsequence')
	downcpgcount = collect_total_avail_cpg_in_column(btFeatures,'downstreamsequence')
	negFeatures = btFeatures.loc[btFeatures['directionality']=='-']
	if not negFeatures.empty:
		nonFeatures = btFeatures.loc[btFeatures['directionality']!='-']
		newnegFeatures = negFeatures.rename(columns={'upstreamsequence':'downstreamsequence','downstreamsequence':'upstreamsequence'})
		revFeatures = pd.concat([newnegFeatures,nonFeatures])
		revupcpgcount = collect_total_avail_cpg_in_column(revFeatures,'upstreamsequence')
		revdowncpgcount = collect_total_avail_cpg_in_column(revFeatures,'downstreamsequence')
# 		stream['cpgsequencecountsumcombineboundary'] = upcpgcountsum + downcpgcountsum
	else:
		revupcpgcount,revdowncpgcount=upcpgcount,downcpgcount
		print 'There are no features with "-" directionality'
	return upcpgcount,downcpgcount,revupcpgcount,revdowncpgcount

# Make graphs for fangs
def graph_boundary_methylation(pdfeatures,fileName):
	info = str(fileName)
	sns.set_style('ticks')
	pp = PdfPages('Methylation_{0}.pdf'.format(fileName))
	plt.figure(figsize=(14,7))
	plt.suptitle(info,fontsize=10)
	surroundingfang = periphery*2

	if not splitgraphs:
		gs = gridspec.GridSpec(1,1,height_ratios=[1])
		gs.update(hspace=.8)
		
		ax0 = plt.subplot(gs[0,:])
		sns.barplot(data=pdfeatures,x='Tissue',y='percpgmeth',hue='organization',ax=ax0)#,palette={"Element":'#8ba6e9',"Random":'#bed0f4'}
		ax0.set_ylabel('% CpGs Methylated',size=16)
		ax0.tick_params(axis='both',which='major',labelsize=16)

	else:
		gs = gridspec.GridSpec(2,1,height_ratios=[1,1],width_ratios=[1])
		gs.update(hspace=.8)
		ax0 = plt.subplot(gs[0,:])
		ax1 = plt.subplot(gs[1,:])
		
# 		Elements = sortstreams[sortstreams['barcolor']=='Element']
# 		Randoms = sortstreams[sortstreams['barcolor']=='Random']
# 		
# 		sns.barplot(data=Elements,x='Tissue',y='percpgmeth',hue='boundary',ax=ax0)#,palette={"Upstream":'#8ba6e9',"Downstream":'#d7b7bc'}
# 		sns.barplot(data=Randoms,x='Tissue',y='percpgmeth',hue='boundary',ax=ax1)#,palette={"Upstream":'#bed0f4',"Downstream":'#f7e4e0'}
# 		ax0.tick_params(axis='both',which='major',labelsize=16)
# 		ax1.tick_params(axis='both',which='major',labelsize=16)
# 		ax0.set_ylabel('% CpGs Methylated',size=16)
# 		ax1.set_ylabel('% CpGs Methylated',size=16)

	sns.despine()
	pp.savefig()
	pp.close()

def main():
	# Get and set parameters
	args = get_args()
	set_global_variables(args)
	paramlabels = '{0}_{1}_{2}_{3}_{4}'.format(elementsize,periphery,binDir,eFiles,stringName)
	
	rangeFeatures = collect_element_coordinates(eFiles)
	directionFeatures = assign_directionality_from_arg_or_boundary(rangeFeatures,eFiles)
	
	upcpg,downcpg,revupcpg,revdowncpg = collect_total_avail_cpg_values(directionFeatures)
	
	collect = [] # initiate collections
	
	allstream = collect_methylation_data_by_element(directionFeatures)
	allstream['group'] = 'element'
	allstream['organization'] = 'unsorted'
	allstream['cpgsum'] = np.where(allstream['boundary'] == 'up stream',upcpg,downcpg)
	collect.append(allstream)
	
	if rFiles:
		for randomFile in rFiles:
			randomFeatures = collect_element_coordinates(randomFile)
			randirFeatures= assign_directionality_from_arg_or_boundary(randomFeatures,randomFile)
			ranupcpg,randowncpg,ranrevupcpg,ranrevdowncpg = collect_total_avail_cpg_values(randirFeatures)
			allrandom = collect_methylation_data_by_element(randirFeatures)
			allrandom['group'] = 'random{0}'.format(randomFile)
			allrandom['organization'] = 'unsorted'
			allrandom['cpgsum'] = np.where(allrandom['boundary'] == 'up stream',ranupcpg,randowncpg)
			collect.append(allrandom)
			
			if reverseComplement:
				randomrev = modify_negative_directionality_elements(randirFeatures)
				randomrev['group'] = 'random{0}'.format(randomFile)
				randomrev['organization'] = 'rcsorted'
				randomrev['cpgsum'] = np.where(randomrev['boundary'] == 'up stream',ranrevupcpg,ranrevdowncpg)
				collect.append(randomrev)
	
	if reverseComplement:
		allreverse = modify_negative_directionality_elements(allstream)
		allreverse['group'] = 'element'
		allreverse['organization'] = 'rcsorted'
		allreverse['cpgsum'] = np.where(allreverse['boundary'] == 'up stream',revupcpg,revdowncpg)
		collect.append(allreverse)
		
	concat = pd.concat(collect)
	concat['percpgmeth'] = concat['methcount']/concat['cpgsum']*100.0
	print concat
	graph_boundary_methylation(concat,'all_{0}'.format(paramlabels))

# 	if labelcolumn:
# 		typeList = directionFeatures['type'].unique()
# 		for type in typeList:
# 		typeBool,typeupstreammethylation,typedownstreammethylation = separate_dataframe_by_group(type,rangeFeatures,'type',eFiles)

if __name__ == "__main__":
	main()