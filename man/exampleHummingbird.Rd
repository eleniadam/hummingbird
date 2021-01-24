\name{exampleHummingbird}
\alias{exampleHummingbird}
\docType{data}
\title{
Sample dataset
}
\description{
Example of input data for the hummingbird package.

The sample dataset is partial data of chromosome 29 in the large offspring
syndrome (LOS) study described in Chen Z. et al (2017).
}
\usage{data("exampleHummingbird")}
\format{
    \describe{
    \item{experimentSEControl}{
A SummarizedExperiment object containing the input data for the control group:
The two assays: normM, normUM and the CpG position information: pos.
}
    \item{experimentSECase}{
A SummarizedExperiment object containing the input data for the case group:
The two assays: abnormM, abnormUM and the CpG position information: pos.
}
    \item{normM}{
A matrix containing the methylated read count data of the control group. Each
column of the matrix represents a replicate and each row represents a CpG
position.
}
    \item{normUM}{
A matrix containing the unmethylated read count data of the control group. Each
column of the matrix represents a replicate and each row represents a CpG
position.
}
    \item{abnormM}{
A matrix containing the methylated read count data of the case group. Each
column of the matrix represents a replicate and each row represents a CpG
position.
}
    \item{abnormUM}{
A matrix containing the unmethylated read count data of case group. Each
column of the matrix represents a replicate and each row represents a CpG
position.
}
    \item{pos}{
The CpG positions.
}
}
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
}
\keyword{datasets}
