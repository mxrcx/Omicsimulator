% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoadTCGAMatrix.R
\name{LoadTCGAMatrix}
\alias{LoadTCGAMatrix}
\title{Load the TCGA Matrix of specified disease}
\usage{
LoadTCGAMatrix(disease, sample_type, sample_number)
}
\arguments{
\item{disease}{String; cancer name}

\item{sample_type}{String; type of samples ('Primary solid Tumor' or 'Solid Tissue Normal')}

\item{sample_number}{Integer; the number of samples which are simulated}
}
\value{
tcga_matrix
}
\description{
Load the TCGA Matrix of specified disease
}
