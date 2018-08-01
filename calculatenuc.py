"""
Script to get nucleotide content within element

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
import pybedtools as pbt

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file", type=str,help='A file with the bed coordinates for your elements')
	parser.add_argument("-l", "--labelcolumn",type=int,help='column in the element file where the label (exonic, intergenic, intronic) is') # way to only read in certain entries, like only read in if 'intergenic'
	parser.add_argument("-g","--genome",type=str, default="hg19.genome")
	parser.add_argument("-f","--fasta",type=str,default="hg19.fa")
	parser.add_argument("-p","--periphery",type=int,default="10",help='number bp from your boundary you want to include in the analysis')
	parser.add_argument('-s',"--stringname",type=str,help='string to add to the outfile name')
	return parser.parse_args()

def set_global_variables(args):
	global file
	global label
	global genome
	global fasta
	global boundary
	global stringname
	file = args.file
	label = args.labelcolumn
	genome = args.genome
	fasta = args.fasta
	boundary = args.periphery
	stringname = args.stringname

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
	midFeatures['startup'] = midFeatures.loc[:,1] - boundary
	midFeatures['startdown'] = midFeatures.loc[:,1] + boundary
	midFeatures['endup'] = midFeatures.loc[:,2] - boundary
	midFeatures['enddown'] = midFeatures.loc[:,2] + boundary
	if label:
		midFeatures['type'] = midFeatures.loc[:,label]
	return midFeatures

# used in get_fasta_for_element_coordinates to extract just the fasta strings
def get_just_fasta_sequence_for_feature(inFeature):
	seqFeature = inFeature.sequence(fi=fasta)
	outFeature = pd.read_table(seqFeature.seqfn)
	outSequence = outFeature[::2]
	return outSequence.reset_index(drop=True)

# get the strings for sliding window regions
def get_fasta_for_element_coordinates(rangeFeatures):
	rangeFeatures['feature'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','start','end']].values.astype(str).tolist()))
	rangeFeatures['feature'] = rangeFeatures['feature'].str.upper()
	rangeFeatures['upstreamsequence'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','startup','startdown']].values.astype(str).tolist()))
	rangeFeatures['upstreamsequence'] = rangeFeatures['upstreamsequence'].str.upper()
	rangeFeatures['downstreamsequence'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','endup','enddown']].values.astype(str).tolist()))
	rangeFeatures['downstreamsequence'] = rangeFeatures['downstreamsequence'].str.upper()
	return rangeFeatures

# get coordinates with flanking regions
def collect_element_coordinates(fileName):
	btFeatures = get_bedtools_features(fileName)
	subsetFeatures = collect_coordinates_for_element_positions(btFeatures)
	return get_fasta_for_element_coordinates(subsetFeatures)

# Get the percentage AT in the element
def percentage_at_for_element(features,label):
	collectAT = []
	collectCpG = []
	for f in features:
		collectAT.append(eval('100*float(f.count("A") + f.count("T"))/len(f)'))
		collectCpG.append(eval('100*float(f.count("CG"))/len(f)'))
	pdAT = pd.DataFrame(collectAT)
	pdCpG = pd.DataFrame(collectCpG)
	table = pd.DataFrame([pdAT.mean(),pdCpG.mean()],index=['%ATcontent','%CpGcontent'])
	table.columns = [label]
	return table.T

# save panda
def save_panda(pdData,strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True)

def main():
	args = get_args()
	set_global_variables(args)
	rangeFeatures = collect_element_coordinates(file)
	collect = []

	allwhole = percentage_at_for_element(rangeFeatures['feature'],'allwhole')
	collect.append(allwhole)
	allupboundary = percentage_at_for_element(rangeFeatures['upstreamsequence'],'allup{0}bp'.format(boundary))
	collect.append(allupboundary)
	alldownboundary = percentage_at_for_element(rangeFeatures['downstreamsequence'],'alldown{0}bp'.format(boundary))
	collect.append(alldownboundary)
	
	if label:
		typelist = rangeFeatures['type'].unique()
		for t in typelist:
			type = (rangeFeatures[rangeFeatures['type'] == t])
			typewhole = percentage_at_for_element(type['feature'],'{0}whole'.format(t))
			collect.append(typewhole)
			typeupboundary = percentage_at_for_element(type['upstreamsequence'],'{0}up{1}bp'.format(t,boundary))
			collect.append(typeupboundary)
			typedownboundary = percentage_at_for_element(type['downstreamsequence'],'{0}down{1}bp'.format(t,boundary))
			collect.append(typedownboundary)

	cat = pd.concat(collect)
	save_panda(cat,'nucleotide_content_{0}_{1}bpboundary.txt'.format(file,boundary))

if __name__ == "__main__":
	main()
