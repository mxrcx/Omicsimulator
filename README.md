# Omicsimulator
Omicsimulator for Simulation of Gene Expression Data

## How to use it:
First install the package omicsimulator, then start the simulation with :
 
`Omicsimulator()`
 
You can add additional parameters, which are all optional, to optimize results:

* input.file (OPTIONAL) The location of a special formatted input file, containing gene expression relations, e.g. "inputDirectory/input.txt". Usally, the program creates this file automatically, using a KEGG pathway, so just specify the disease parameter.
* disease (OPTIONAL) The name of a specific cancer disease which specifies a KEGG pathway, e.g. "Breast cancer". It is used to extract the gene expression relations for the input file.
* output.directory (OPTIONAL) The directory of the output files, e.g. "specificPath/outputDirectory".
* sample.number (OPTIONAL) The number of samples used for the simulation, e.g. "10".
 
Example call:

`Omicsimulator(disease="Breast cancer")`
