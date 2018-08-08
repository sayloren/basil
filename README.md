Brief description, dependent files can be found by using -h

asymmetryprobability.py

Takes in a file with a list of filenames, to return a plot the the probability of getting two strings of nucleotides that have the same number of AT content
If plotted with -e, each file will plot independently, if -r the average will be plotted
boundarymethylation-graphs.pyboundarymethylation-v1.pyboundarymethylation-v2-o2.pyboundarymethylation-v2.py

'v1' is a legacy file, depreciated
'v2' is the complete analysis of boundary methylation, 'o2' runs with out the graphs, the output of which is put into 'graphs'
calculatechi.py

Is not yet working
calculatenuc.py

Calculates the nucleotide content
cpgcontentcalculator.py

Calculates the cpg content
nucleotideatcontent-randomassingment-graphs.pynucleotideatcontent-randomassingment-o2.pynucleotideatcontent-randomassingment-two-graphs.pynucleotideatcontent-randomassingment.py'nucleotideatcontent-randomassingment.py' is the whole analysis for finding the AT content across a group of elements, and then performing the random assignment for -i iterations'two-graphs' plots the AT content for your elements and random files on two separate plots
'o2' runs with out the graphs, the output of which is put into 'graphs'nucleotidecontent-graphs.pynucleotidecontent-o2.pynucleotidecontent.py

Plots the nucleotide content. 'o2' runs with out the graphs, the output of which is put into 'graphs'nucleotidecpgcontent-graphs.pynucleotidecpgcontent-v1.pynucleotidecpgcontent-v2-o2.pynucleotidecpgcontent-v2.py

'v1' is a legacy file, depreciated
'v2' is the complete analysis for cog content, 'o2' runs with out the graphs, the output of which is put into 'graphs'
nucleotidemethylation-graphs.pynucleotidemethylation-o2.pynucleotidemethylation.py

'nucleotidemethylation.py' plots the counts of methylation across a group of elements. 'o2' runs with out the graphs, the output of which is put into 'graphs'nucletideatslope-graph.pynucletideatslope.py

'nucletideatslope.py' plots the slope of AT content for a group of elements. 'graph' can take in the 'Data*' file from the AT content random assignment script
subsample.py

Takes in a group of elements, and returns a random subsampling
symmetryprobability.py

Returns a matrix of the proabilities for a group of elements for being AT rich, poor, or balanced