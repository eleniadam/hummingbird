\name{normM}
\alias{normM}
\docType{data}
\title{
Sample matrix
}
\description{
A matrix containing the methylated read count data of the control group.
It is part of the sample dataset exampleHummingbird.
}
\usage{normM}
\details{
Each column of the matrix represents a replicate and each row represents 
a CpG position.
}
\source{
Chen et al. (2017) Global misregulation of genes largely uncoupled to DNA
methylome epimutations characterizes a congenital overgrowth syndrome.
Scientific Reports 7, 12667

The raw FASTQ files of the WGBS experiment from this study are publicly
available at Gene Expression Omnibus (GEO) database with accession no. GSE93775.
}
\examples{
data(exampleHummingbird)
normM
}
\keyword{datasets}
