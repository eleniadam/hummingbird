hummingbirdGraph <- function(pos, normM, normUM, abnormM, abnormUM, dmrs, coord1, coord2){

    ## Data within the selected coordinates
    genomicPositions <- which(pos >= coord1 & pos <= coord2)
    cntPos <- NULL
    cpgPos <- NULL
    labelPos <- NULL
    k <- 1
    found <- 0
    for(cnt in genomicPositions){
        i <- 1
        while(i <= dim(dmrs)[1]){
            ## Within a DMR
            if( (pos[cnt] >= dmrs[i, 2]) & (pos[cnt] <= dmrs[i, 3]) ){
                cntPos[k] <- cnt
                cpgPos[k] <- pos[cnt]
                labelPos[k] <- dmrs[i, 5]
                found <- 1
                break;
            }
            i <- i + 1
        }
        ## Not in a DMR
        if(found == 0){
            cntPos[k] <- cnt
            cpgPos[k] <- pos[cnt]
            labelPos[k] <- 0
        }
        found <- 0
        k <- k + 1
    }

    ## Case and control
    norm_p <- NULL
    abnorm_p <- NULL
    k <- 1
    for(cnt in cntPos){
        norm_p[k] <- sum(normM[cnt,])/sum(normM[cnt,], normUM[cnt,])
        abnorm_p[k] <- sum(abnormM[cnt,])/sum(abnormM[cnt,], abnormUM[cnt,])
        k <- k + 1
    }

    ## Data presentation for the figures
    valuesDifference <- norm_p - abnorm_p
    valMin <- min(valuesDifference)
    valMax <- max(valuesDifference)
    fig2Limit <- (max(abs(valMax), abs(valMin))) + 0.1

    ## Figure 1 - Methylation Level
    par(mar=c(5, 5, 1, 1))
    plot(cpgPos, norm_p, pch=1, main="", xlab="", ylab="Methylation Level", ylim=c(0, 1))
    lines(cpgPos, abnorm_p, pch=19, type="p")
    legend("bottomleft", c("Norm", "Abnorm"), pch=c(1, 19), cex=.75)

    ## Figure 2 - Predictions
    par(mar=c(5, 5, 1, 1))
    plot(cpgPos, valuesDifference, pch=NA, main="", xlab="", ylab="Predictions", ylim=c(-fig2Limit, fig2Limit))
    abline(h=0, lty=4)
    text(cpgPos, valuesDifference, labels=labelPos)
    legend("bottomleft", c("0: No change", "1: Hyper Methylation", "2: Hypo Methylation"), cex=.75)

}


