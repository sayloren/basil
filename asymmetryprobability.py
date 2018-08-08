"""
Script to determine how the probability of getting boundaries with the same
value changes as the bin size of the element boundaries changes

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
	parser.add_argument("-e","--elementfile",type=argparse.FileType('rU'),help='A file containing a list of paths to the element files with unique names separated by newlines')
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
	global eFiles
	binDir = args.bin
	labelcolumn = args.labelcolumn
	directionalitycolumn = args.directionalitycolumn
	sizeGenome = args.genome
	faGenome = args.fasta
	if args.randomfile:
		global rFiles
		rFiles = [line.strip() for line in args.randomfile]
	eFiles = [line.strip() for line in args.elementfile]

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
	rangeFeatures['feature'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','start','end']].values.astype(str).tolist()))
	rangeFeatures['feature'] = rangeFeatures['feature'].str.upper()
	if not directionalitycolumn:
		rangeFeatures['directionality'] = rangeFeatures.apply(lambda row: (compare_boundaries_size_n(row['feature'],binDir)),axis=1)
	return rangeFeatures

# Run calculate_nucleotides_at for increasing bin sizes
def emperical_boundaries_for_increasing_bins_percentage(element,size):
	totalSteps=[]
	for i in np.arange(1,size*2):
		pairStep = calculate_nucleotides_at(element,i)
		totalSteps.append(pairStep)
	return totalSteps

# Run compare_boundaries_size_n over increasing bin sizes
def emperical_boundaries_for_increasing_bins_symbol(element,size):
	totalSteps=[]
	for i in np.arange(1,size*2):
		perSize = calculate_nucleotides_at(element,i)
		# give + - = depending on which side has larger AT content
		if perSize[0] > perSize[1]: outList = '+'
		if perSize[1] > perSize[0]: outList = '-'
		if perSize[1] == perSize[0]: outList = '='
		totalSteps.append(outList)
	return totalSteps

# Get the actual spread of = in the data over increasing bin size
def collect_emperical_boundary_comparisons(rangeFeatures,name):
	rangeAT = rangeFeatures.apply(lambda row: (emperical_boundaries_for_increasing_bins_percentage(row['feature'],binDir)),axis=1)
	pdrangeAT = pd.DataFrame(rangeAT.values.tolist())
	equalAT = rangeFeatures.apply(lambda row: (emperical_boundaries_for_increasing_bins_symbol(row['feature'],binDir)),axis=1)
	pdequalAT = pd.DataFrame(equalAT.values.tolist())
	countAT = pdequalAT.apply(pd.value_counts)
	equalCounts = countAT.loc['=']/len(pdequalAT.index)
	splitAT = pd.concat(dict([(row[0],row[1].apply(lambda y: pd.Series(y))) for row in pdrangeAT.iterrows()]),axis=1)
	outcollect = []
	for column in splitAT:
		outcollect.append(splitAT[column[0]])
	outcat = pd.concat(outcollect,axis=1)
	outcat /= 100 # convert to decimal
	pdBins = pd.concat([equalCounts],axis=1)
	pdBins.columns=[name]
	return pdBins

# Compute the theoretical probability
def run_boundary_probability_calculation(yrange,p):
	# https://math.stackexchange.com/questions/151810/probability-of-3-heads-in-10-coin-flips
	# need to include probability different than 0.5
	equal = []
	for y in yrange:
		totalPerm = pow(2,y)
		if totalPerm > 0:
			permsum = []
			for k in range (0, (y + 1)):
				permuationK = math.factorial(y)/(math.factorial(k)*(math.factorial(y -k)))#*pow(p,k)*pow((1-p),y-k)
				floatPForKHeads = float(permuationK)/float(totalPerm)
				floatPForKHeadsBothSides = floatPForKHeads * floatPForKHeads
				permsum.append(floatPForKHeadsBothSides)
			probabilitiesfloatPEqualSides = sum(permsum)
			equal.append(probabilitiesfloatPEqualSides)
		else:
			equal.append(probabilitiesfloatPEqualSides)
	return equal

# save panda
def save_panda(pdData, strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True)

# Graph the probability
def graph_equal_boundary_probability(features,paramlabels):
	yrange = numpy.arange(1,binDir*2)
	equal = run_boundary_probability_calculation(yrange,0.5)
	features['Theoretical'] = equal
	sns.set_style('ticks')
	sns.set_palette("husl",n_colors=8)
	plt.figure(figsize=(3.5,3.5))
	features.plot(yrange,linewidth=2,alpha=0.9)
	plt.xlabel('Bin Size',size=12)
	plt.ylabel('Probability',size=12)
	plt.title('Equal Boundary for {0} Bins'.format(binDir),size=16)
	plt.legend(loc=0,fontsize=6,labelspacing=0.05)
	plt.tight_layout()
	sns.despine()
	plt.savefig('Probability_{0}.pdf'.format(paramlabels))

def main():
	args = get_args()
	set_global_variables(args)
	collectFeatures = []
	collectName = []
	for e in eFiles:
		collectName.append(e)
		rangeFeatures = collect_element_coordinates(e)
		directionFeatures = assign_directionality_from_arg_or_boundary(rangeFeatures,e)
		binFeatures = collect_emperical_boundary_comparisons(directionFeatures,e)
		collectFeatures.append(binFeatures)
		if labelcolumn:
			typeList = directionFeatures['type'].unique()
			for type in typeList:
				typeFeatures = (directionFeatures[directionFeatures['type'] == type])
				binType = collect_emperical_boundary_comparisons(typeFeatures,'{0}_{1}'.format(type,e))
				collectFeatures.append(binType)
		concatFeatures = pd.concat(collectFeatures,axis=1)
	if args.randomfile:
		collectRandom = []
		for r in rFiles:
			randomFeatures = collect_element_coordinates(r)
			directionRandom = assign_directionality_from_arg_or_boundary(randomFeatures,r)
			binRandom = collect_emperical_boundary_comparisons(directionRandom,r)
			collectRandom.append(binRandom)	
		concatRandom = pd.concat(collectRandom,axis=1)
		aveRandom = concatRandom.mean(axis=1)
		concatFeatures['Random'] = aveRandom	
	paramlabels = '{0}_{1}'.format(binDir,eFiles)
	graph_equal_boundary_probability(concatFeatures,paramlabels)

if __name__ == "__main__":
	main()
