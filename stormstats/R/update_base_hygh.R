update_base_hygh <- function(x,y,verbose=TRUE){
  ## x: a list of objects from base_hydrograph
  ## y: review results showing the
  ##    base event to be combined with adjacent storm events
  helper <- function(hygh,line){
    ## hygh: a list containing peak, base and value data frame
    ##       output from hydrograph package
    ## line: a instruction from manual review of some baseflow events
    ## return: the same list with instruction incorporated
    #browser()
    base_ <- subset(hygh$base,index==line$base)
    peak_ <- subset(hygh$peak,index==line$storm1)
    storm_ <- subset(hygh$value,index==line$storm1)
    run_ <- peak_$run[1]
    lookup <- function(hygh_,run_,idx_){
      ## helper function to look up event ID
      ## based on the time stamp location
      event_ <- subset(hygh_$event,run==run_)
      loc_ <- with(event_,which(start<=idx_&(end+1)>=idx_))
      ## when there are multiple runs, choose the latter one
      loc_ <- loc_[length(loc_)]
      hygh_$event$event[loc_]
    }
    combine <- function(i1,a1){
      ## helper function to combine an interval
      ## with a couple of adjacent intervals
      ## with a simple error check to ensure these
      ## intervals are indeed adjacent
      ## return the combined interval
      a_ <- na.omit(rbind(i1,a1))
      l_ <- apply(a_,1,function(x) diff(range(x))+1)
      n_ <- diff(range(c(a_)))+1
      stopifnot(sum(l_)==n_)
      range(c(a_))
    }
    ## the base flow to append
    astart_ <- c(base_$iStart,NA,NA)
    aend_ <- c(base_$iStop,NA,NA)
    if(!is.na(line$base_2)){
      ## additional base flow to append
      base2_ <- subset(hygh$base,index==line$base_2)
      astart_[2] <- base2_$iStart
      aend_[2] <- base2_$iStop
    }
    if(!is.na(line$storm2)){
      ## additional storm flow to convert
      peak2_ <- subset(hygh$peak,index==line$storm2)
      storm2_ <- subset(hygh$value,index==line$storm2)
      astart_[3] <- storm2_$iStart
      aend_[3] <- storm2_$iEnd
    }
    
    if(line$rise==0){
      ## the interval of the last falling limb
      istart_ <- peak_[nrow(peak_),"iFallStart"]
      iend_ <- peak_[nrow(peak_),"iFallEnd"]
      ## combine with baseflow etc.
      out_ <- combine(c(istart_,iend_),cbind(astart_,aend_))
      ## update the corresponding limb
      peak_$iFallStart[nrow(peak_)] <- out_[1]
      peak_$iFallEnd[nrow(peak_)] <- out_[2]
      peak_$jFallStart[nrow(peak_)] <- lookup(hygh,run_,out_[1])
      peak_$jFallEnd[nrow(peak_)] <- lookup(hygh,run_,out_[2])
    }else{
      ## the interval of the first rising limb
      istart_ <- peak_[1,"iRiseStart"]
      iend_ <- peak_[1,"iRiseEnd"]
      ## combine with baseflow etc.
      out_ <- combine(c(istart_,iend_),cbind(astart_,aend_))
      ## update the corresponding limb
      peak_$iRiseStart[1] <- out_[1]
      peak_$iRiseEnd[1] <- out_[2]
      peak_$jRiseStart[1] <- lookup(hygh,run_,out_[1])
      peak_$jRiseEnd[1] <- lookup(hygh,run_,out_[2])
    }
    storm_$iStart <- min(peak_$iRiseStart)
    storm_$iEnd <- max(peak_$iFallEnd)
    
    ## Remove base interval
    hygh$base <- subset(hygh$base,index != line$base) 
    if(!is.na(line$base_2)){
      ## Remove adjacent base (if any)
      hygh$base <- subset(hygh$base, index != line$base_2)
    }
    
    if(!is.na(line$storm2)){
      ## Remove adjacent storm-flow (if any)
      hygh$peak <- subset(hygh$peak, index != line$storm2)
      hygh$value <- subset(hygh$value, index != line$storm2)
    }
    
    ## Update storm data
    within(hygh,{
      peak[peak$index==line$storm1,]=peak_
      value[value$index==line$storm1,]=storm_
    })
  }
  if(is.null(y)){
    y <- data.frame(station=NA,	base=NA,
                    storm1=NA,	storm2=NA,	base_2=NA,	
                    rise=NA)
  }
  #browser()
  for(k in 1:length(x)){
    #if(k==6) browser()
    cat("station=",k,"\n")
    x_ <- x[[k]]
    y_ <- subset(y,station==k)
    if(nrow(y_)==0) next
    for(i in 1:nrow(y_)){
      text_ <- paste0("append base ",y_$base[i])
      if(!is.na(y_$base_2[i]))
        text_ <- paste0(text_," and base ",y_$base_2[i])
      if(!is.na(y_$storm2[i])) 
        text_ <- paste0(text_," , and convert storm ",y_$storm2[i])
      if(y_$rise[i]==0){
        if(verbose) cat(text_, "to falling limb of storm",
                        y_$storm1[i],"\n")
      }
      if(y_$rise[i]==1){
        if(verbose) cat(text_,"to rising limb of storm",
                        y_$storm1[i],"\n")
      }
      x_ <- helper(hygh = x_,line = y_[i,])
    }
    x[[k]] <- x_
  }
  return(x)
}

dev_ <- function(){
  rm(list=ls());gc(T)
  load(file="~/Rscratch/update_base_hygh_debug.rData")
  source("update_base_hygh.R")
  names(x[[1]])
  str(x[[1]]$event)
  str(x[[1]]$peak)
  str(x[[1]]$value)
  str(x[[1]]$base)
  x2 <- update_base_hygh(x,y)
  subset(x[[1]]$base,index %in% c(31,210,338),select=c(1,3))
  subset(x2[[1]]$base,index %in% c(31,210,338),select=c(1,3))
  i <- 32
  subset(x[[1]]$value,index ==i,select=c(3,4))
  subset(x2[[1]]$value,index ==i,select=c(3,4))
  subset(x[[1]]$peak,index ==i,select=c(5,6,9,10))
  subset(x2[[1]]$peak,index ==i,select=c(5,6,9,10))
}