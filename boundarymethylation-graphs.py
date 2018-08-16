"""
Script to graph methylation analysis for flanks with heatmap

Wren Saylor
August 2018

Copyright 2018 Harvard University, Wu Lab

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
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import re

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile",type=str,help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("-hb","--histogrambins",type=int,default="20",help='')
	return parser.parse_args()

# set all args that will be used throughout the script
def set_global_variables(args):
	global eFiles
	global histogrambins
	eFiles = args.efile
	histogrambins = args.histogrambins
	outfilename = parse_file_name(eFiles)
	return outfilename

# splitting the filename into a dataframe
def parse_file_name(file):
	outfilename = re.search('Data_MethBoundary_(.*?).txt',file).group(1)
	return outfilename

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

# plot params
def set_plot_params(removedups,xval,yval,hval,pp,setxlabel,whichplot,elementpalette,randompalette):
	element,random = seperate_elements_and_random(removedups,'group','element')
	numberrandom = len(random.index)
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1],width_ratios=[1])
	gs.update(hspace=.8)
	if numberrandom != 0:
		plt.figure(figsize=(10,10))
		ax0 = plt.subplot(gs[0,:])
		ax1 = plt.subplot(gs[1,:])
		if whichplot == 'boxplot':
			sns.barplot(data=element,x=xval,y=yval,hue=hval,ax=ax0,palette=elementpalette,errwidth=1)
			sns.barplot(data=random,x=xval,y=yval,hue=hval,ax=ax1,palette=randompalette,errwidth=1)
			for label in ax0.get_xticklabels():
				label.set_rotation(15)
			for label in ax1.get_xticklabels():
				label.set_rotation(15)
			for bartype in element[hval].unique():
				typeelement = element[element[hval]==bartype]
				typerandom = random[random[hval]==bartype]
				typefillnaelement = typeelement[yval].fillna(0)
				typefillnarandom = typerandom[yval].fillna(0)
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
		for tissue in element['Tissue'].unique():
			tissueelement = element[element['Tissue']==tissue]
			tissuerandom = random[random['Tissue']==tissue]
			tissuefillnaelement = tissueelement[yval].fillna(0)
			tissuefillnarandom = tissuerandom[yval].fillna(0)
		ax0.set_title("Ultraconserved Elements")
		ax1.set_title("Random Regions")
		subplots = [ax0,ax1]
		for plot in subplots:
			plot.tick_params(axis='both',which='major',labelsize=16)
			plot.set_ylabel('Count Methylation',size=16)
			plot.set_xlabel(setxlabel,fontsize=16)
		sns.despine()
	else:
		plt.figure(figsize=(10,5))
		ax0 = plt.subplot(gs[:,:])
		if whichplot == 'boxplot':
			sns.barplot(data=element,x=xval,y=yval,hue=hval,ax=ax0,palette=elementpalette,errwidth=1)
			for label in ax0.get_xticklabels():
				label.set_rotation(15)
			for bartype in element[hval].unique():
				typeelement = element[element[hval]==bartype]
				typefillnaelement = typeelement[yval].fillna(0)
		else:
			for tissue in element[hval].unique():
				tissueelement = element[element[hval]==tissue]
				tissueelement.dropna(axis=0,inplace=True)
				sns.distplot(tissueelement[xval],ax=ax0,label=tissue,bins=histogrambins)
			ax0.legend()
		fillnaelement = element[yval].fillna(0)
		for tissue in element['Tissue'].unique():
			tissueelement = element[element['Tissue']==tissue]
			tissuefillnaelement = tissueelement[yval].fillna(0)
		ax0.set_title("Ultraconserved Elements")
		subplots = [ax0]
		for plot in subplots:
			plot.tick_params(axis='both',which='major',labelsize=16)
			plot.set_ylabel('Count Methylation',size=16)
			plot.set_xlabel(setxlabel,fontsize=16)
		sns.despine()

	plt.savefig(pp,format='pdf')

# Make graphs for fangs
def graph_boundary_methylation(sorted,filelabel):
	pp = PdfPages('{0}.pdf'.format(filelabel))
	sns.set_style('ticks')
	sns.set_palette("husl",n_colors=len(sorted['Tissue'].unique()))
	removedupcpgper = group_and_count_data_frame_by_column(sorted,'boundary','boundarycount')
	removedupstrand = group_and_count_data_frame_by_column(sorted,'strand','strandcount')
	strandDict = {'C':'+','G':'-'}
	removedupstrand['strandedness'] = removedupstrand.loc[:,'strand'].map(strandDict)
	removedupdir = group_and_count_data_frame_by_column(sorted,'directionality','dircount')
	dirDict = {'+':'AT rich','-':'AT poor','=':'AT balanced'}
	removedupdir['ATcontent'] = removedupdir.loc[:,'directionality'].map(dirDict)
	removeduploc = group_and_count_data_frame_by_column(sorted,'methlocation','loccount')
	removeduppercentage = group_and_count_data_frame_by_column(sorted,'percentage','percount')
	boundarypaletteelement = {'up stream':'#afbfe3','down stream':'#6f84ba'}
	boundarypaletterandom = {'up stream':'#d6b5ba','down stream':'#9d787d'}
	set_plot_params(removedupcpgper,'Tissue','boundarycount','boundary',pp,'Tissue','boxplot',boundarypaletteelement,boundarypaletterandom)
	strandpaletteelement = {'+':'#afbfe3','-':'#6f84ba'}
	strandpaletterandom = {'+':'#d6b5ba','-':'#9d787d'}
	set_plot_params(removedupstrand,'Tissue','strandcount','strandedness',pp,'Tissue','boxplot',strandpaletteelement,strandpaletterandom)
	directionpaletteelement = {'AT rich':'#afbfe3','AT poor':'#6f84ba','AT balanced':'#465475'}
	directionpaletterandom = {'AT rich':'#d6b5ba','AT poor':'#9d787d','AT balanced':'#624b4f'}
	set_plot_params(removedupdir,'Tissue','dircount','ATcontent',pp,'Tissue','boxplot',directionpaletteelement,directionpaletterandom)
	set_plot_params(removeduploc,'methlocation','loccount','Tissue',pp,'Distance from Boundary','distplot','husl','husl') 
	set_plot_params(removeduppercentage,'percentage','percount','Tissue',pp,'Percentage','distplot','husl','husl')
	pp.close()	

def main():
	args = get_args()
	outfilename = set_global_variables(args)
	df = pd.read_csv(eFiles,sep='\t',header=0,index_col=0)
	graph_boundary_methylation(df,outfilename)

if __name__ == "__main__":
	main()