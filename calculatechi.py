"""
Script to return tables 

Wren Sayor
July 13 2017

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
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import chisquare

# Make Methylation graphs
def create_chi_square_table(TableData,Title,ranTableData,ranTitle,fileName):
	sns.set_style('ticks')
	info = str(fileName)
# 	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Table_{0}.pdf'.format(fileName))
	plt.figure(figsize=(16,6))
	
	allData = pd.concat([TableData,ranTableData])#,keys=['elements','random']

	outchisquare = []
	for x,y in zip(TableData,ranTableData):
		outchisquare.append(chisquare(TableData[x],ranTableData[y]))
	pdchi = pd.DataFrame(outchisquare)
	pdchi.index = TableData.T.index
	
	allData.columns.name = None
	printTable = pd.concat([allData.T,pdchi],axis=1)

	gs = gridspec.GridSpec(1,1,height_ratios=[1])
	ax0 = plt.subplot(gs[0])
	ax0.set_frame_on(False)
	ax0.set_yticks([])
	ax0.set_xticks([])
	Table0 = ax0.table(cellText=[['']*2],colLabels=['Observed', 'Expected','Chi Square'],loc='center',bbox=[0, 0.6, 0.8, 0.1])
	Table1 = ax0.table(cellText=[['']],colLabels=['Chi Square'],loc='center', bbox=[0.8, 0.6, 0.2, 0.1])
	Table2 = ax0.table(cellText=printTable.values.round(3),rowLabels=printTable.index,colLabels=printTable.columns,cellLoc='center',rowLoc='center',loc='center', bbox=[0, 0.35, 1.0, 0.3])
	ax0.set_title('Chi Square Test for {0} and {1}'.format(Title, ranTitle),size=24)
	Table2.set_fontsize(22)
	plt.tight_layout()

	sns.despine()
	plt.savefig(pp, format='pdf')
	pp.close()

def main(TableData,Title,ranTableData,ranTitle,fileName):
	create_chi_square_table(TableData,Title,ranTableData,ranTitle,fileName)

if __name__ == "__main__":
	main()
