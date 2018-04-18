"""
Script to run methylation analysis for flanks with heatmap

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
	parser.add_argument("-p","--periphery",type=int,default="20",help='number bp from your boundary you want to include in the analysis')
	parser.add_argument("-b","--bin",type=int,default="100",help='size of bins used to compare element ends and determine directionality')
	parser.add_argument("-e","--element",type=int,help='size of your element (region without flanks), should be an even number, if not provided will use the smallest size of your input data')
	parser.add_argument("-mu","--methylationthresholdpercentageupper",type=int, default="100",help='size to threshold percentage methylation data upper bound')
	parser.add_argument("-ml","--methylationthresholdpercentagelower",type=int, default="1",help='size to threshold percentage methylation data lower bound')
	parser.add_argument("-mc","--methylationthresholdcoverage",type=int,default="10",help='size to threshold uncapped coverage of methylation data to send to percentage methylation, field often uses 10')
	parser.add_argument("-ms","--combinemethylationsamples",action='store_true',help='whether to combine those samples with the same tissue/cell type or leave as seperate')
	parser.add_argument('-s',"--stringname",type=str,help='string to add to the outfile name')
	parser.add_argument('-c',"--reversecomplement",action='store_true',help='if you want the reverse complement to be plotted')
	parser.add_argument("-hb","--histogrambins",type=int,default="20",help='')
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
	global histogrambins
	periphery = args.periphery
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
	histogrambins = args.histogrambins

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
	pdmethThresh = (pdmethFeatures[(pdmethFeatures.loc[:,3] >= methCovThresh) & (pdmethFeatures.loc[:,4] <= methPerThreshupper) & (pdmethFeatures.loc[:,4] >= methPerThreshlower)])
	btmethThresh = convert_panda_to_bed_format(pdmethThresh)
	print 'methylation coverage is being thresholded at {0} and percentage between {1} and {2}'.format(methCovThresh,methPerThreshlower,methPerThreshupper)
	initiallength = len(pdmethFeatures.index)
	checklength = initiallength - len(pdmethThresh.index)
	print 'there were {0} out of {1} total removed by thresholding {2}'.format(checklength,initiallength,methName)
	return btmethThresh

# Intersect regions from the methylation data with element regions
def intersect_methylation_data(btFeatures,meFeature,meName,eStart,eEnd,column,boundary,subtractcolumn):
	initiallength = len(meFeature)
	intersectboundary = meFeature.intersect(btFeatures[['chr',eStart,eEnd,'id',subtractcolumn]].values.astype(str).tolist(),wb=True,wa=True)
	checklength = len(intersectboundary)
	print 'there were {0} intersections out of {1} in {2}'.format(checklength,initiallength,meName)
	if len(intersectboundary) != 0:
		pdmeth=convert_bedtool_to_panda(intersectboundary)
		pdmeth.columns = ['mchr','mstart','mstop','coverage','percentage','chr','estart','estop','id','subtract']
		pdmeth['strand'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmeth[['mchr','mstart','mstop']].values.astype(str).tolist()))
		pdmeth['strand'] = pdmeth['strand'].str.upper()
		pdmeth['methlocation'] = pdmeth['subtract'].astype(int)-pdmeth['mstart'].astype(int)
		outmeth = pdmeth[['chr','mstart','mstop','estart','estop','id','percentage','methlocation','strand']]
		outmeth['Tissue'] = meName.split(".")[0]
		subfeatures = btFeatures[['id','directionality',column]]
		merge = pd.merge(outmeth,subfeatures,how='left',on='id')
		merge['boundary'] = boundary
		sortoutmeth = merge[['chr','id','directionality','methlocation','percentage','strand','boundary','Tissue']]
	else:
		sortoutmeth = pd.DataFrame(np.nan,index=[0],columns=['chr','id','directionality','methlocation','percentage','strand'])
		sortoutmeth['boundary'] = boundary
		sortoutmeth['Tissue'] = meName.split(".")[0]
	return sortoutmeth

# Run the analysis to extract methylation properties
def collect_methylation_data_by_element(pdfeatures):
	capture = []
	for methname in mFiles:
		mfeatures=get_bedtools_features(methname)
		pdthresh=threshold_methylation_data(mfeatures,methname)
		uintersect=intersect_methylation_data(pdfeatures,pdthresh,methname,'startup','startdown','upstreamsequence','up stream','start')
		dintersect=intersect_methylation_data(pdfeatures,pdthresh,methname,'endup','enddown','downstreamsequence','down stream','end')
		uinvert=invert_location(uintersect)
		capture.append(uinvert)
		capture.append(dintersect)
	totalconcat = pd.concat(capture)
	totalconcat.reset_index(drop=True,inplace=True)
	return totalconcat

# invert methlocation column
def invert_location(pdfeatures):
	originalRange = range(-periphery,periphery+1)#range(0,periphery*2)
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
	boundaryDict = {'up stream':'down stream','down stream':'up stream'}
	invFeatures['newboundary'] = invFeatures.loc[:,'boundary'].map(boundaryDict)
	invFeatures['newstrand'] = invFeatures.loc[:,'strand'].map(seqDict)
	newnegFeatures = invFeatures[['chr','id','directionality','methlocation','percentage','newstrand','newboundary','Tissue','group']]
	newnegFeatures.columns = ['chr','id','directionality','methlocation','percentage','strand','boundary','Tissue','group']
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
	allreverse['cgsum'] = np.where(allreverse['boundary'] == 'up stream',revupcg,revdowncg)
	return allstream,allreverse

# count by feature for graphing
def group_and_count_data_frame_by_column(pdfeatures,incol,outcol):
	methsort = pdfeatures.sort_values(by=['Tissue','group',incol],ascending=True)
	methsort['methgroupby'] = methsort.groupby(['Tissue',incol,'group'])[incol].transform('count')
	removedups = methsort.drop_duplicates(['group',incol,'Tissue','methgroupby'],keep='last') # only count each groupby object once!!
	removedups[outcol] = removedups.groupby(['group','Tissue','methgroupby',incol])['methgroupby'].transform('sum')
	return removedups

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
	if (len(element.index) == 0) or (len(random.index) ==0):
		formatpval,statcoef,stattest='nan','nan','nan'
	else:
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

# plot params
def set_plot_params(removedups,xval,yval,hval,pp,setxlabel,whichplot,elementpalette,randompalette):
	element,random = seperate_elements_and_random(removedups,'group','element')
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1],width_ratios=[1])
	gs.update(hspace=.8)
	ax0 = plt.subplot(gs[0,:])
	ax1 = plt.subplot(gs[1,:])
	collectstats = []
	if whichplot == 'boxplot':
		sns.barplot(data=element,x=xval,y=yval,hue=hval,ax=ax0,palette=elementpalette)
		sns.barplot(data=random,x=xval,y=yval,hue=hval,ax=ax1,palette=randompalette)
		for label in ax0.get_xticklabels():
			label.set_rotation(15)
		for label in ax1.get_xticklabels():
			label.set_rotation(15)
		for bartype in element[hval].unique():
			typeelement = element[element[hval]==bartype]
			typerandom = random[random[hval]==bartype]
			typefillnaelement = typeelement[yval].fillna(0)
			typefillnarandom = typerandom[yval].fillna(0)
			formatpvaltype,statcoeftype,stattesttype = run_appropriate_test(typefillnaelement,typefillnarandom)
			statstable = pd.DataFrame([yval,bartype,formatpvaltype,statcoeftype,stattesttype],index=['count set','comparison group','p value','coefficient','stats test'])
			collectstats.append(statstable)
	else:
		for tissue in element[hval].unique():
			tissueelement = element[element[hval]==tissue]
			tissuerandom = random[random[hval]==tissue]
			tissueelement.dropna(axis=0,inplace=True)
			tissuerandom.dropna(axis=0,inplace=True)
			sns.distplot(tissueelement[xval],ax=ax0,label=tissue,bins=histogrambins)
			sns.distplot(tissuerandom[xval],ax=ax1,label=tissue,bins=histogrambins)#,norm_hist=False
		ax0.legend()
		ax1.legend()
	fillnaelement = element[yval].fillna(0)
	fillnarandom = random[yval].fillna(0)
	formatpval,statcoef,stattest = run_appropriate_test(fillnaelement,fillnarandom)
	type = 'whole set'
	statstable = pd.DataFrame([yval,type,formatpval,statcoef,stattest],index=['count set','comparison group','p value','coefficient','stats test'])
	collectstats.append(statstable)
	for tissue in element['Tissue'].unique():
		tissueelement = element[element['Tissue']==tissue]
		tissuerandom = random[random['Tissue']==tissue]
		tissuefillnaelement = tissueelement[yval].fillna(0)
		tissuefillnarandom = tissuerandom[yval].fillna(0)
		tissueformatpval,tissuestatcoef,tissuestattest = run_appropriate_test(tissuefillnaelement,tissuefillnarandom)
		statstable = pd.DataFrame([yval,tissue,tissueformatpval,tissuestatcoef,tissuestattest],index=['count set','comparison group','p value','coefficient','stats test'])
		collectstats.append(statstable)
	ax0.set_title("Ultraconserved Elements")
	ax1.set_title("Random Regions")
	subplots = [ax0,ax1]
	for plot in subplots:
		plot.tick_params(axis='both',which='major',labelsize=16)
		plot.set_ylabel('Count Methylation',size=16)
		plot.set_xlabel(setxlabel,fontsize=16)
	sns.despine()
	plt.savefig(pp,format='pdf')
	return collectstats

# Make graphs for fangs
def graph_boundary_methylation(pdfeatures,filelabel):
	methfiles = [(file.split("-")[0]) for file in mFiles]
	methcount = Counter(methfiles)
	info = str(filelabel)
	sns.set_style('ticks')
	sns.set_palette("husl",n_colors=len(pdfeatures['Tissue'].unique()))
	if reverseComplement:
		sorted = pdfeatures.loc[pdfeatures['organization']=='rcsorted']
		pp = PdfPages('Methylation_Boundary_ReverseComplementSorted_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.pdf'.format(eFiles,stringName,filelabel,elementsize,binDir,periphery,methPerThreshlower,methPerThreshupper,methCovThresh))
		pstat = 'Statistic_Boundary_ReverseComplementSorted_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.pdf'.format(eFiles,stringName,filelabel,elementsize,binDir,periphery,methPerThreshlower,methPerThreshupper,methCovThresh)
	else:
		sorted = pdfeatures.loc[pdfeatures['organization']=='unsorted']
		pp = PdfPages('Methylation_Boundary_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.pdf'.format(eFiles,stringName,filelabel,elementsize,binDir,periphery,methPerThreshlower,methPerThreshupper,methCovThresh))
		pstat = 'Statistic_Boundary_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.pdf'.format(eFiles,stringName,filelabel,elementsize,binDir,periphery,methPerThreshlower,methPerThreshupper,methCovThresh)
	plt.figure(figsize=(14,7))
	plt.suptitle(info,fontsize=16)
	removedupcpgper = group_and_count_data_frame_by_column(sorted,'boundary','boundarycount')
# 	removedupcpgper['percpgmeth'] = (removedupcpgper[outcol]/removedupcpgper['cpgsum'])*100.0
	removedupstrand = group_and_count_data_frame_by_column(sorted,'strand','strandcount')
	strandDict = {'C':'+','G':'-'}
	removedupstrand['strandedness'] = removedupstrand.loc[:,'strand'].map(strandDict)
	removedupdir = group_and_count_data_frame_by_column(sorted,'directionality','dircount')
	dirDict = {'+':'AT rich','-':'AT poor','=':'AT balanced'}
	removedupdir['ATcontent'] = removedupdir.loc[:,'directionality'].map(dirDict)
	removeduploc = group_and_count_data_frame_by_column(sorted,'methlocation','loccount')
	removeduppercentage = group_and_count_data_frame_by_column(sorted,'percentage','percount')
	collectstat = []
	boundarypaletteelement = {'up stream':'#8ba6e9','down stream':'#adc0ef'}
	boundarypaletterandom = {'up stream':'#c5969d','down stream':'#d6b5ba'}
	statboundary = set_plot_params(removedupcpgper,'Tissue','boundarycount','boundary',pp,'Tissue','boxplot',boundarypaletteelement,boundarypaletterandom)
	strandpaletteelement = {'+':'#8ba6e9','-':'#adc0ef'}
	strandpaletterandom = {'+':'#c5969d','-':'#d6b5ba'}
	statstrand = set_plot_params(removedupstrand,'Tissue','strandcount','strandedness',pp,'Tissue','boxplot',strandpaletteelement,strandpaletterandom)
	directionpaletteelement = {'AT rich':'#8ba6e9','AT poor':'#adc0ef','AT balanced':'#6174a3'}
	directionpaletterandom = {'AT rich':'#c5969d','AT poor':'#d6b5ba','AT balanced':'#89696d'}
	statdirection = set_plot_params(removedupdir,'Tissue','dircount','ATcontent',pp,'Tissue','boxplot',directionpaletteelement,directionpaletterandom)
	statlocation = set_plot_params(removeduploc,'methlocation','loccount','Tissue',pp,'Distance from Boundary','distplot','husl','husl') 
	statpercentage = set_plot_params(removeduppercentage,'percentage','percount','Tissue',pp,'Percentage','distplot','husl','husl')
	collectstat.extend(statboundary)
	collectstat.extend(statstrand)
	collectstat.extend(statdirection)
	collectstat.extend(statlocation)
	collectstat.extend(statpercentage)
	catstat = pd.concat(collectstat,axis=1)
	catstat.reset_index(drop=True,inplace=True)
	save_panda(catstat.T,'{0}.txt'.format(pstat))

# 	methsort['createdcpg'] = (methsort['cpgsum']/methsort['cgsum'])*100.0
# 	if combinesamples:
# 		methsort['newTissue'] = methsort['Tissue'].str.extract('(^.*)[-].*?$',expand=True)
# 		methsort.drop(labels='Tissue',axis=1,inplace=True)
# 		methsort.rename(columns={'newTissue':'Tissue'},inplace=True)
# 	removedupid = group_data_frame_by_column(methsort,'id','idcount')
# 	boxplot_params(removedupid,'idcount','id',pp,'Count UCE ID Methylated') # needs work
# 	removedupavailcg = group_data_frame_by_column(methsort,'createdcpg','createdcpgtest')
# 	boxplot_params(removedupavailcg,'createdcpg','boundary',pp,'% C and G creating CpG')
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
	graph_boundary_methylation(concat,'{0}'.format(label))

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