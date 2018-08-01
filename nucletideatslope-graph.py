"""
Script to run mean AT conent slope analysis

Wren Saylor
October 2017

Other: Analog and discrete signals, Fourier series, Spectral analysis, Fourier transform, Discrete Fourier Transform, Nyquist-Shannon sampling theorem, aliasing

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
from scipy.interpolate import splrep, splev
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

# Get smoothed mean, first and second derivatives
def collect_smoothed_lines(ATmean,fillX,window):
	halfwindow = ((window/2)+1)
	f=splrep(fillX,ATmean,k=3,s=50)
	smoothMean = splev(fillX,f)
	firstDer = splev(fillX,f,der=1)
	firstDer[0:halfwindow] = 0 # small edge effect
	firstDer[-halfwindow:] = 0 # small edge effect
	secondDer = splev(fillX,f,der=2)
	secondDer[0:window] = 0 # small edge effect
	secondDer[-window:] = 0 # small edge effect
	return smoothMean,firstDer,secondDer

# graph
def graph_element_line_means_random_below(outfilename,df):
	allnums = [int(s) for s in outfilename.split("_") if s.isdigit()]
	num = allnums[1]
	window = allnums[3]
	fillX = range(0,(num-window))
	sns.set_style('ticks')
	gs = gridspec.GridSpec(2,1)
	gs.update(hspace=.8)
	pp = PdfPages('Slope_{0}.pdf'.format(outfilename))
	plt.figure(figsize=(10,10))
	ATmean = df['Element']
	smoothMean,firstDer,secondDer = collect_smoothed_lines(ATmean,fillX,window)
	ranATmean = df['Random']
	ransmoothMean,ranfirstDer,ransecondDer = collect_smoothed_lines(ranATmean,fillX,window)
	ax0 = plt.subplot(gs[0,:])
	ax1 = plt.subplot(gs[1,:])
	ax0.set_ylabel('% AT Content',size=16)
	ax0.set_xlabel('Position (bp)',size=16)
	plt.xlim(0,num)
	ax1.plot(fillX,ransecondDer,linewidth=plotlinesize,label='Random',color='#d0abb0')
	ax0.plot(fillX,secondDer,linewidth=plotlinesize,label='Element',color='#4d5c82')
	subplots = [ax0,ax1]
	for plot in subplots:
		plot.tick_params(axis='both',which='major',labelsize=20)
	sns.despine()
	pp.savefig()
	pp.close()

def main():
	args = get_args()
	outfilename = set_global_variables(args)
	df = pd.read_csv(eFiles,sep='\t',header=0,index_col=0)
	graph_element_line_means_random_below(outfilename,df)

if __name__ == "__main__":
	main()