"""
Script to accept a bunch of bed files, and report GpC summary stats for each one. Requires a genome in fasta format
to be specified

Adapted from Ruth's GC_summary_stats4 with minimal changes

Inputs: List of bed file names
Genome in fasta format

Outputs: bedtools nuc output for each bed file
Summary stats for each bed file, written to a table

Ruth McCole + Wren Saylor
January 2018

"""

import argparse
import pandas as pd
import pybedtools as pbt

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file", type=argparse.FileType('rU'),help='A file containing a list of paths to the feature files you want to process, ''separated by newlines')
	parser.add_argument("-g", "--genome", type=str, default='hg19.fa',help='The full path of the fasta file where you want to draw the nucleotide sequences from')
	return parser.parse_args()

def getFileNames(args):
	arFilenames = [line.strip() for line in args.file]
	return arFilenames

def getFeatures(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

def saveBedTool(btObject, strFilename):
	btObject.saveas(strFilename)

def getNucForFilename(strFilename, strFastaPath):
	#Get a bedtool
	btFeatures = getFeatures(strFilename)
	#Run nuc
	btNucResult = btFeatures.nucleotide_content(fi=strFastaPath,pattern='CG',C=True)
	#save this result
	strNucFilename = 'Nuc_full_result_for_{0}.txt'.format(strFilename)
	saveBedTool(btNucResult, strNucFilename)
	return btNucResult

def getCpGColumnName(pandasNuc):
	arPandasNucColumns = list(pandasNuc)
	#print arPandasNucColumns
	for strPandascolumn in arPandasNucColumns:
		if 'patt_count' in strPandascolumn:
			#print strPandascolumn
			return strPandascolumn

def getTotalColumnName(pandasNuc):
	arPandasNucColumns = list(pandasNuc)
	#print arPandasNucColumns
	for strPandascolumn in arPandasNucColumns:
		if 'seq_len' in strPandascolumn:
			#print strPandascolumn
			return strPandascolumn

def getpdGC(btNucResult):
	#Make bedtool into pandas object
	pandasNuc = pd.read_table(btNucResult.fn)
	#Get the name of the column that holds the information for CpG counts
	strCpGColumn = getCpGColumnName(pandasNuc)
	# Get the name of the column that holds the information for total number of nucleotides in the seq
	strTotalColumn = getTotalColumnName(pandasNuc)
	#Get a pandas dataframe containing just CpG counts and total number of nuc in seq
	pdCpGPercentage = pandasNuc[strCpGColumn]/pandasNuc[strTotalColumn]
	return pdCpGPercentage

def getFeaturesAndGC(btNucResult):
	pandasNuc = pd.read_table(btNucResult.fn)
	strCpGColumn = getCpGColumnName(pandasNuc)
	strTotalColumn = getTotalColumnName(pandasNuc)
	pdFeatureCpG = pandasNuc[['#1_usercol','2_usercol','3_usercol',strCpGColumn,strTotalColumn]]
	pdFeatureCpG['14_pct_cpg'] = pdFeatureCpG[strCpGColumn]/pdFeatureCpG[strTotalColumn]
	pdFeatureCpG.drop(columns=[strCpGColumn,strTotalColumn],inplace=True)
	return pdFeatureCpG

def getSummaryStats(pdGC, btNucResult):
	intN = pdGC.count()
	intCoverage = btNucResult.total_coverage()
	floatMedian = pdGC.median()
	floatMean = pdGC.mean()
	float25thPercentile = pdGC.quantile(q=0.25)
	float75thPercentile = pdGC.quantile(q=0.75)
	floatMin = pdGC.min()
	floatMax = pdGC.max()
	arStats = [intN, intCoverage, floatMedian, floatMean, float25thPercentile, float75thPercentile, floatMin, floatMax]
	return arStats

def main():
	print 'The expected build is hg19, if your bed files are not hg19, you need to specify an appropriate genome in fasta format'
	args = get_args()
	strFastaPath = args.genome
	arFilenames = getFileNames(args)

	arArStats = []
	for strFilename in arFilenames:
		print 'Getting data for {0}'.format(strFilename)
		#Run bedtools Nuc
		btNucResult = getNucForFilename(strFilename, strFastaPath)

		#Get panda containing results
		pdGC = getpdGC(btNucResult)

		#Save GC contents
		strGCFilename = 'CpG_fractions_for_plots_{0}.txt'.format(strFilename)
		pdGC.to_csv(strGCFilename, sep='\t', index=False)

		#Save bed file with GC fraction for each feature
		pdFeatureGC = getFeaturesAndGC(btNucResult)
		strFeatureGCFilename = 'Feature_and_CpG_fraction_for_matrix_{0}'.format(strFilename)
		pdFeatureGC.to_csv(strFeatureGCFilename, sep='\t', index=False, header=False)

		#Get GC stats
		arStats = getSummaryStats(pdGC, btNucResult)
		arArStats.append(arStats)

	arStatColumns = ['NumObs', 'Coverage_bp', 'floatMedian', 'floatMean', 'float25thPercentile', 'float75thPercentile', 'floatMin', 'floatMax']
	#Produce a pandas dataframe of stats
	pdAllStats = pd.DataFrame(data=arArStats, index=arFilenames, columns=arStatColumns)
	#print pdAllStats
	pdAllStats.to_csv('CpGStats.txt', sep='\t')

if __name__ == "__main__":
	main()