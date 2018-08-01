"""
Script to run mean AT conent analysis graph

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

# set command line arguments
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=str,help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument('-p',"--plotlinesize",type=int,default=3,help='size of the line to plot')
	return parser.parse_args()

# set all args that will be used throughout the script
def set_global_variables(args):
	global eFiles
	global plotlinesize
	eFiles = args.efile
	plotlinesize = args.plotlinesize
	outfilename = parse_file_name(eFiles)
	return outfilename

# splitting the filename into a dataframe
def parse_file_name(file):
	outfilename = re.search('Data_ATcontent_(.*?).txt',file).group(1)
	return outfilename

# graph
def graph_element_line_means_random_below(outfilename,df):
	sns.set_style('ticks')
	gs = gridspec.GridSpec(2,1)
	gs.update(hspace=.8)
	pp = PdfPages('Fangs_{0}.pdf'.format(outfilename))
	plt.figure(figsize=(10,10))
	ATmean = df['Element']
	ranATmean = df['Random']
	ax0 = plt.subplot(gs[0,:])
	ax1 = plt.subplot(gs[1,:])
	ax0.set_ylabel('% AT Content',size=16)
	ax0.set_xlabel('Position (bp)',size=16)
	allnums = [int(s) for s in outfilename.split("_") if s.isdigit()]
	num = allnums[1]
	window = allnums[3]
	fillX = range(0,(num-window))
	ax0.plot(fillX,ATmean,linewidth=plotlinesize,label='Element',color='#4d5c82')
	ax1.plot(fillX,ranATmean,linewidth=plotlinesize,label='Random',color='#d0abb0')
	subplots = [ax0,ax1]
	for plot in subplots:
		plot.tick_params(axis='both',which='major',labelsize=20)
		plot.set_xlim(0,num)
	sns.despine()
	pp.savefig()
	pp.close()

def main():
# 	starttime = time.time()
	args = get_args()
	outfilename = set_global_variables(args)
	df = pd.read_csv(eFiles,sep='\t',header=0,index_col=0)
	graph_element_line_means_random_below(outfilename,df)

if __name__ == "__main__":
	main()