#' getAltCircRNA
#' 
#' Given the coordinates of the backsplicing junction, builds possible alternative internal
#'  structures of a circRNA by permuting the introns and exons within boundaries.
#'
#' @param g The path to a gff file, or a GRanges object of (at least) exons with the 
#' attributes `transcript_id` and `type`.
#' @param circRNA The coordinates of the circRNA backsplice junction; either a GRanges 
#' containing a single range, or a string of the form `seqname:start-end`
#' @param circName The name of the circRNA, serves as a prefix for the names of the 
#' alternative forms
#' @param intronIncMaxSize The maximum size for an intron to be included.
#'
#' @return A BED-like data.frame
#' @export
getAltCircRNA <- function(g, circRNA="chr1:7110274-7110696", circName="circRNA", intronIncMaxSize=10000){
    library(GenomicRanges)
    if(is.character(g)) g <- rtracklayer::import(g)
    if(is(circRNA,"GRanges")){
        cr <- circRNA
    }else{
        x <- strsplit(gsub("-",":",circRNA),":",fixed=T)[[1]]
        cr <- GRanges(seqnames=x[[1]],IRanges(as.numeric(x[[2]]),as.numeric(x[[3]])))
    }
    
    # we identify the transcripts that have exons starting and ending at exactly these coordinates (+/-1 to account for notation difference)
    S <- function(x){ as.character(seqnames(x)) }
    t1 <- unique(as.character(g$transcript_id[which( g$type=="exon" & S(g)==S(cr) & start(g) %in% (-1:1 + start(cr)))]))
    t2 <- unique(as.character(g$transcript_id[which( g$type=="exon" & S(g)==S(cr) & end(g) %in% (-1:1 + end(cr)))]))
    tx <- intersect(t1,t2)
    if(length(tx)==0) return(NULL) # no matching tx
    print(cr)
    print(tx)
    
    g <- g[which(g$type=="exon" & g$transcript_id %in% tx),]
    library(gtools)
    
    # for each tx, we find the exons between the boundaries and run the exon combination function
    res <- unlist(lapply(split(g[,c()],g$transcript_id), cr, intronIncMaxSize=intronIncMaxSize, FUN=function(x, cr, intronIncMaxSize){
        x <- sort(x[findOverlaps(x,cr,type="within")@from])
        nbExons <- length(x)
        if(nbExons>8){
            # we'll do something else
            return(NULL)
        }else{
            x <- getExonCombinations(x)
            if(nbExons > 1 & intronIncMaxSize>0){
                x <- unlist(sapply(x,intronIncMaxSize=intronIncMaxSize,FUN=getIntronCombinations))
            }
            return(x)
        }
    }))

    # we convert to bed format
    res <- t(sapply(strsplit(res,",",fixed=T), FUN=function(blocks){
        x <- strsplit(blocks,"-",fixed=T)
        starts <- as.numeric(sapply(x,function(x) x[[1]]))
        ends <- as.numeric(sapply(x,function(x) x[[2]]))
        c(starts[1],ends[length(ends)],length(starts),paste(ends-starts,collapse=","), paste(starts-starts[1],collapse=","))
    }))
    n <- nrow(res)
    res <- data.frame( chr=rep(as.character(seqnames(cr)),n), 
                       start=res[,1],
                       end=res[,2],
                       name=paste(circName,1:n,sep="."), 
                       score=rep(0,n), 
                       strand=rep(as.character(strand(g)[1]),n), 
                       V7=rep(".",n), V8=rep(".",n), V9=rep(".",n), 
                       blockCounts=res[,3],
                       blockSizes=res[,4],
                       blockStarts=res[,5] 
                     )
    row.names(res) <- NULL
    return(res)
}

getExonCombinations <- function(g){
    if(length(g)<=2) return(paste(paste(start(g),end(g),sep="-"),collapse=","))
    g <- cbind(start(g),end(g))
    ex <- apply(g,1,collapse="-",FUN=paste)
    apply(permutations(2,nrow(g)-2,c(T,F),repeats.allowed=T),ex=ex,1,FUN=function(w,ex){
        paste(ex[c(1,1+which(w),length(ex))],collapse=",")
    })
}

getIntronCombinations <- function(tx, intronIncMaxSize=10000){
    exons <- strsplit(tx,",",fixed=T)[[1]]
    x <- strsplit(exons,"-",fixed=T)
    introns <- cbind(as.numeric(sapply(x[-length(x)],FUN=function(x)x[[2]])),
                     as.numeric(sapply(x[-1],FUN=function(x)x[[1]])))
    w <- which( (introns[,2]-introns[,1]+1)<intronIncMaxSize )
    if(length(w)==0) return(tx)
    x <- strsplit(tx,"-",fixed=T)[[1]]
    apply(permutations(2,length(w),c(T,F),repeats.allowed=T),w=w,tx=x, 1,FUN=function(x,w,tx){
        if(all(!x)) return(paste(tx,collapse="-"))
        return(paste(tx[-(w[which(x)]+1)],collapse="-"))
    })
}

# not used
# x is a matrix of the start/end coord, i in the current position
.oneStep <- function(i,e,curpath=""){
    # we add the current exon to the path
    w <- which(e[,2]==i)
    curpath <- paste0(curpath,ifelse(curpath=="","",","),paste(e[w,1],e[w,2],sep="-"))
    # we find the exons that start after the current position
    w <- which(e[,1]>=i)
    # if there is none, we're done and return the path
    if(length(w)==0)    return(curpath)
    # if there are some exons left downstream, we call the current function on them
    unlist(lapply(e[w,2],e=e[w,,drop=F],curpath,FUN=.oneStep))
}
