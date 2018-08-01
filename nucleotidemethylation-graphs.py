"""
Script to run methylation analysis as heatmap graph

Wren Saylor
July 2018

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
	return parser.parse_args()

def set_global_variables(args):
	global eFiles
	eFiles = args.efile
	outfilename = parse_file_name(eFiles)
	return outfilename

# splitting the filename into a dataframe
def parse_file_name(file):
	outfilename = re.search('Data_MethFlanks_(.*?).txt',file).group(1)
	return outfilename

# Make graphs for fangs
def graph_methylation(outfilename,df):
	sns.set_style('ticks')
	plt.figure(figsize=(10,10))
	pp = PdfPages('Methylation_Heatmap_{}.pdf'.format(outfilename))
	
	readelement = df.filter(like='_Element')
	readelement.columns = [col.replace('_Element','')  for col in readelement.columns]
	readrandom = df.filter(like='_Random')
	readrandom.columns = [col.replace('_Random','')  for col in readrandom.columns]
	
	formatelement = readelement.T
	formatrandom = readrandom.T
	
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
	heatmap1 = sns.heatmap(formatrandom,cmap='RdPu',ax=ax1,xticklabels=100)
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

def main():
	args = get_args()
	outfilename = set_global_variables(args)
	df = pd.read_csv(eFiles,sep='\t',header=0,index_col=0)
	graph_methylation(outfilename,df)

if __name__ == "__main__":
	main()