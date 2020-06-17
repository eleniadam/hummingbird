\name{exampleHummingbird}
\alias{exampleHummingbird}
\docType{data}
\title{
Sample dataset
}
\description{
Example of input data of the hummingbird package. Data used as input to the hummingbird package: A SummarizedExperiment object containing four assays under the names normM, normUM, abnormM, abnormUM, the CpG positions under the name pos and the corresponding chromosome name.
}
\usage{data("exampleHummingbird")}
\format{
    \describe{
    \item{\code{exampleSE}}{
A SummarizedExperiment object containing the four assays and the position information.
}
    \item{normM}{
A matrix containing the methylated read count data of the normal group. Each column of the matrix represents a replicate and each row represents a CpG position.
}
    \item{normUM}{
A matrix containing the unmethylated read count data of the normal group. Each column of the matrix represents a replicate and each row represents a CpG position.
}
    \item{abnormM}{
A matrix containing the methylated read count data of the abnormal group. Each column of the matrix represents a replicate and each row represents a CpG position.
}
    \item{abnormUM}{
A matrix containing the unmethylated read count data of abnormal group. Each column of the matrix represents a replicate and each row represents a CpG position.
}
    \item{pos}{
The CpG positions.
}
}
}
\source{
The raw FASTQ files are publicly available at Gene Expression Omnibus, GEO accession no. GSE93775.
Chen, Z., Hagen, D. E., Ji, T., Elsik, C. G., and Rivera, R. M. (2017). Global misregulation of genes largely uncoupled to DNA methylome epimutations characterizes a congenital overgrowth syndrome. Scientific Reports 7, 12667.
}
\examples{
library(SummarizedExperiment)
data(exampleHummingbird)
assays(exampleSE)[["normM"]]
assays(exampleSE)[["normUM"]]
assays(exampleSE)[["abnormM"]]
assays(exampleSE)[["abnormUM"]]
rowRanges(exampleSE)
}
\keyword{datasets}
