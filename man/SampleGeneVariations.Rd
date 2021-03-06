% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SampleGeneVariations.R
\name{SampleGeneVariations}
\alias{SampleGeneVariations}
\title{Generates a vcf file based on a given list of gene symbols}
\usage{
SampleGeneVariations(disease, output_directory, genes_variation,
  file_prefix)
}
\arguments{
\item{disease}{String; cancer name}

\item{output_directory}{String; name of the output directory}

\item{genes_variation}{List; genes with variations}

\item{file_prefix}{String; name of the results output file}
}
\value{
None
}
\description{
Gets gene range of every gene with variations, Searches for possible variations in this range, Writes vcf entry for first variation. Prints complete vcf file.
}
\examples{
SampleGeneVariations(genes.variation)
}
