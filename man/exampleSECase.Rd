\name{exampleSECase}
\alias{exampleSECase}
\docType{data}
\title{
Sample input data
}
\description{
A SummarizedExperiment object containing the input data for the case group.
It is part of the sample dataset exampleHummingbird.

}
\usage{exampleSECase}
\details{
It contains the two assays: abnormM, abnormUM and the CpG position 
information: pos.
}
\source{
Chen et al. (2017) Global misregulation of genes largely uncoupled to DNA
methylome epimutations characterizes a congenital overgrowth syndrome.
Scientific Reports 7, 12667

The raw FASTQ files of the WGBS experiment from this study are publicly
available at Gene Expression Omnibus (GEO) database with accession no. GSE93775.
}
\examples{
library(SummarizedExperiment)
data(exampleHummingbird)
exampleSECase
}
\keyword{datasets}
