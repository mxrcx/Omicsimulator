# R-Script - How to use Omicimulator (Sample script)

## Sample Call 1 - detailled Query
## result <- Omicsimulator(disease = disease, sample_number = 10, top_DEG_number = 1000, output_directory = "../output", file_prefix = "omicsimulatorResults_Breast cancer_10_1000")
##
## result contains a list which contains two matrices from DEA (normal-tumor and normal-simulated) with gene name, gene expression value, fold change and p-value

## Sample Call 2 - basic Query
## result <- Omicsimulator(disease = disease, sample_number = 10, top_DEG_number = 1000)
##
## Output folder and file prefix of the output files are standard values as seen in the docs.
## result contains a list which contains two matrices from DEA (normal-tumor and normal-simulated) with gene name, gene expression value, fold change and p-value

## Sample Call 3 - test Query
## result <- Omicsimulator()
##
## All input parameters are standard values as seen in the docs.
## result contains a list which contains two matrices from DEA (normal-tumor and normal-simulated) with gene name, gene expression value, fold change and p-value
