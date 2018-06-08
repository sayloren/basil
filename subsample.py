"""
Script to subsample from a bedfile

Wren Saylor
June 2018

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
import numpy as np

# set command line arguments
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file",type=argparse.FileType('rU'),help="A file containing a list of paths to the files to process")
	parser.add_argument("-s","--subsample",type=int,default="1000",help='the number of elements to subsample')# ratio or number?
	return parser.parse_args()

# get bt features
def get_bedtools_features(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

# convert bedtool to panda
def bedtool_to_panda(btobject):
	pdobject = pd.read_table(btobject.fn,header=None)
	return pdobject

# subsample random rows
def subsample_rows(pdobject,subsample):
	pdsubsample = pdobject.loc[np.random.choice(pdobject.index,subsample,replace=False)]
	return pdsubsample

# save panda
def save_panda(pdData, strFilename):
	pdData.to_csv(strFilename,sep='\t',index=False,header=False,columns=None)

def main():
	args = get_args()
	subsample = args.subsample
	files = [line.strip() for line in args.file]
	for file in files:
		btobject = get_bedtools_features(file)
		pdobject = bedtool_to_panda(btobject)
		pdsubsample = subsample_rows(pdobject,subsample)
		save_panda(pdsubsample,'{0}_subsample_{1}.bed'.format(file,subsample))

if __name__ == "__main__":
	main()