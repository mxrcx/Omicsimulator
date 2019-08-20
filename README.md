# Omicsimulator
Omicsimulator for Simulation of Gene Expression Data

## How to use it:
First install the package omicsimulator, then start the simulation with :
 
`Omicsimulator()`
 
The standard case, without any parameters, just has the parameters as followed: disease as "Breast cancer", sample_number as 10, top_DEG_number as 100 and threshold_eQTls as 0.7.

These parameters are necessary to get specific results:
* `disease` The name of a specific cancer disease which specifies a KEGG pathway, e.g. "Breast cancer". It is used to extract the gene expression relations for the input file.
* `sample_number` The number of samples used for the simulation, e.g. "10".
* `top_DEG_number` The number of top expressed genes used to compare the simulation results, e.g. "100".
* `threshold_eQTls` Threshold of the number of eQTLs used for MAF file generation.

You can add optional parameters as well:
* `output_directory (OPTIONAL)` The directory of the output files, e.g. "specificPath/outputDirectory".
* `file_prefix (OPTIONAL)` The file name prefix of the output files.

The return value is a list, containing the DEA results of normal-tumor and normal-simulated DEA. Both matrices contain the gene, the fold change and the p-value.
 
## Example call:

`Omicsimulator(disease = disease, sample_number = 10, top_DEG_number = 1000, output_directory = "../output", file_prefix = "omicsimulatorResults_Breast cancer_10_1000", threshold_eQTls = 0.7)`
