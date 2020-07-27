\name{hummingbird-package}
\alias{hummingbird-package}
\alias{hummingbird}
\docType{package}
\title{
    A Bayesian Hidden Markov Model for the detection of differentially
    methylated regions
}
\description{
    A package for identifying differentially methylated regions (DMRs) between
    case and control groups using whole genome bisulfite sequencing (WGBS) or
    reduced representative bisulfite sequencing (RRBS) experiment data.

    The hummingbird package uses a Bayesian hidden Markov model (HMM) for
    detecting DMRs. It fits a Bayesian HMM for one chromosome at a time. The
    final output of hummingbird are the detected DMRs with start and end
    positions in a given chromosome, directions of the DMRs (hyper- or hypo-),
    and the numbers of CpGs in these DMRs.

    The hummingbird package implements the algorithm described in the
    publication below.
}
\details{
    The main functions of the package are:
    hummingbirdEM, hummingbirdPostAdjustment and hummingbirdGraph.
}
\author{
Eleni Adam, Tieming Ji, Desh Ranjan

Maintainer: Eleni Adam <eadam002@odu.edu>
}
\references{
    Ji (2019) A Bayesian hidden Markov model for detecting differentially
    methylated regions. Biometrics 75(2):663-673.
}
