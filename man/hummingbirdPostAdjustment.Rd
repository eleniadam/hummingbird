\name{hummingbirdPostAdjustment}
\alias{hummingbirdPostAdjustment}
\title{
Post Adjustment algorithm for the output of the EM
}
\description{
This function adjusts HMM output such that each detected DMR has a minimum length and a minimum number of CpGs in each DMR.
}
\usage{
hummingbirdPostAdjustment(em, pos, minCpGs, minLength, maxGap)
}
\arguments{
    \item{em}{
The obs output of the hummingbirdEM function.
}
    \item{pos}{
The CpG positions. It must be in the structure of a matrix.
}
    \item{minCpGs}{
The minimum number of CpGs contained in a DMR. Default value: 10.
}
    \item{minLength}{
The minimum length of a DMR. Default value: 500.
}
    \item{maxGap}{
The maximum gap between any two CpGs. Default value: 300.
}
}
\value{
\item{obsPostAdj }{A matrix containing the transformed observations for each bin and the adjusted predictions for each bin (0-normal, 1-hyper, 2-hypo).}
\item{DMRs }{A matrix containing the detected regions based on the user-defined arguments. The columns of the DMRs matrix are: Region number, Start genomic position, End genomic position, Length, Direction, Number of CpGs.}
}

\examples{
hmmbird2 <- hummingbirdPostAdjustment(em=hmmbird1$obs, pos=m_pos, minCpGs=10, minLength=100, maxGap=300)
hmmbird2$DMRs
}

