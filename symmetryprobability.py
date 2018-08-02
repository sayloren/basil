"""
Script to determine the liklihood of getting AT rich, poor and balanced elements

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
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy
import seaborn as sns
import math
import pybedtools as pbt
import pandas as pd
import numpy as np

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=str,help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("-r","--randomfile",type=argparse.FileType('rU'),help="A file containing a list of paths to the random regions equable with your elements to plot in contrast") # option to not use random at all, or generate randoms?
	parser.add_argument("-l", "--labelcolumn",type=int,help='column in the element file where the label (exonic, intergenic, intronic) is')
	parser.add_argument("-d", "--directionalitycolumn",type=int,help='column in the element file where directionality is, if not supplied will infer by AT content')
	parser.add_argument("-g","--genome",type=str, default="hg19.genome")
	parser.add_argument("-f","--fasta",type=str,default="hg19.fa")
	parser.add_argument("-b","--bin",type=int,default="100",help='size of bins used to compare element ends and determine directionality')
	return parser.parse_args()

# set all args that will be used throughout the script
def set_global_variables(args):
	global binDir
	global labelcolumn
	global directionalitycolumn
	global sizeGenome
	global faGenome
	global randomassignments
	global rFiles
	global eFiles
	binDir = args.bin
	labelcolumn = args.labelcolumn
	directionalitycolumn = args.directionalitycolumn
	sizeGenome = args.genome
	faGenome = args.fasta
	eFiles = args.efile
	rFiles = [line.strip() for line in args.randomfile]

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
	midFeatures=midFeatures.drop(['middle'],axis=1)
	return midFeatures

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
	return subsetFeatures

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

# get the empirical probability for each direction classification
def make_probabilites_for_direction(directionFeatures,probabilitycolumn,filename):
	lenAll = float(len(directionFeatures.index))
	numPlus = (directionFeatures[probabilitycolumn] == '+').sum()/lenAll
	numMinus = (directionFeatures[probabilitycolumn] == '-').sum()/lenAll
	numEqual = (directionFeatures[probabilitycolumn] == '=').sum()/lenAll
	probOptions = [filename,numPlus,numMinus,numEqual]
	return probOptions




# save panda
def save_panda(pdData, strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True)

def main():
	args = get_args()
	set_global_variables(args)
	collectpd = []
	rangeFeatures = collect_element_coordinates(eFiles)
	directionFeatures = assign_directionality_from_arg_or_boundary(rangeFeatures,eFiles)
	probFeatures = make_probabilites_for_direction(directionFeatures,'directionality',eFiles)
	collectpd.append(probFeatures)
	for r in rFiles:
		randomFeatures = collect_element_coordinates(r)
		directionRandom = assign_directionality_from_arg_or_boundary(randomFeatures,r)
		probRandom = make_probabilites_for_direction(directionRandom,'directionality',r)
		collectpd.append(probRandom)
	df = pd.DataFrame(collectpd,columns=['File','ATrich','ATpoor','ATbalanced'])
	paramlabels = '{0}_{1}'.format(binDir,eFiles)
	save_panda(df,'Orientation_{0}.txt'.format(paramlabels))



if __name__ == "__main__":
	main()
