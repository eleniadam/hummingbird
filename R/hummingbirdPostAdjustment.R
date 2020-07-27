hummingbirdPostAdjustment <- function(experimentInfoControl,
                                    experimentInfoCase, emInfo, minCpGs = 10,
                                    minLength = 500, maxGap = 300){

    ## Input
    obs <- matrix(c(emInfo$direction, emInfo$distance, start(ranges(emInfo)),
                    end(ranges(emInfo))), nrow=length(emInfo), ncol=4)

    ## Computation
    hmmbirdPA <- hummingbirdPostAdjustmentInternal(obs, matrix(rowRanges(
        experimentInfoControl)@ranges@pos), minCpGs, minLength, maxGap)

    ## Output
    GobsPostAdj <- GRanges(matrix(rowRanges(
        experimentInfoControl)@seqnames@values), IRanges(hmmbirdPA$obs[,3],
                                                        hmmbirdPA$obs[,4]),
        distance = hmmbirdPA$obs[,2], direction = hmmbirdPA$obs[,1])

    if(dim(hmmbirdPA$DMRs)[1] == 0){
        ## No DMRs
        gDMRs <- GRanges(matrix(rowRanges(
            experimentInfoControl)@seqnames@values), IRanges(0,0))
    }
    else{
        gDMRs <- GRanges(matrix(rowRanges(
            experimentInfoControl)@seqnames@values),
            IRanges(hmmbirdPA$DMRs[,2], hmmbirdPA$DMRs[,3]),
            length = hmmbirdPA$DMRs[,4], direction = hmmbirdPA$DMRs[,5],
            CpGs = hmmbirdPA$DMRs[,6])
    }
    outPA <- list(obsPostAdj=GobsPostAdj, DMRs=gDMRs)

    return(outPA)
}


