\name{hummingbirdEM}
\alias{hummingbirdEM}
\title{
EM Algorithm for Fitting the Hidden Markov Model
}
\description{
This function reads in methylated and unmethylated read count data, transforms it into logarithm bin-wise data, sets up initial values and implements the EM algorithm to estimate HMM parameters and find the best sequence of hidden states based on model fitting.
}
\usage{
hummingbirdEM(normM, normUM, abnormM, abnormUM, pos, binSize)
}
\arguments{
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
The CpG positions. It must be in the structure of a matrix.
}
    \item{binSize}{
The size of a bin. Default value is: 40.
}

}

\value{
\item{obs }{ A matrix containing the transformed observations for each bin, distance between the current bin the bin ahead of it and the start and end positions of each bin.}
\item{normAbnorm }{A matrix containing the average methylation rate of abnormal and normal groups.}
}
\examples{
data(exampleHummingbird)
hmmbird1 <- hummingbirdEM(normM=m_normM, normUM=m_normUM, abnormM=m_abnormM, abnormUM=m_abnormUM, pos=m_pos, binSize=40)
}

