plotBase <- function(x,i,buffer=50){
  ##x: output from hydrograph function
  ##i: the base-flow event number (index)
  base_ <- x$base
  data <- x$data
  stopifnot(i<=nrow(base_))
  lb_ <- subset(base_,index==i)
  if(nrow(lb_)==0) return(1)
  #browser()
  break_ <- c(lb_$iStartDate[1],lb_$iStopDate[1])
  seq_ <- seq(lb_$iStart[1]-buffer,lb_$iStop[1]+buffer)
  tmp_ <- data[seq_,]
  
  ## find all base events in this window
  overlap_ <- rep(F,nrow(base_))
  for(i in 1:nrow(base_)){
    duration_ <- seq(base_$iStart[i],base_$iStop[i])
    overlap_[i] <- any(duration_ %in% seq_)
  }
  if(any(overlap_)){
    base_overlap_ <- base_[overlap_,]
    break2_ <- with(base_overlap_,c(iStartDate,iStopDate))
    break2_ <- break2_[!(break2_ %in% break_)]
    cat("overlap with baseflow events ",which(overlap_),"\n")
  }
  
  ## find all storm events in this window
  overlap2_ <- rep(F,nrow(x$value))
  for(i in 1:nrow(x$value)){
    duration_ <- seq(x$value$iStart[i],x$value$iEnd[i])
    overlap2_[i] <- any(duration_ %in% seq_)
  }
  if(any(overlap2_)){
    cat("overlap with storm events ",which(overlap2_),"\n")
    storm_overlap_ <- x$value[overlap2_,]
    break3_ <- with(storm_overlap_,data$dateTimeRound[c(iStart,iEnd)])
    break3_ <- break3_[!(break3_ %in% break_)]
  }
  
  y.xts <- xts(x=tmp_$y,order.by = tmp_$dateTimeRound)
  i.xts <- xts(x=tmp_$index,order.by=tmp_$dateTimeRound)
  define_event <- function(type_){
    r <- rep(0,nrow(tmp_))
    loc_ <- which(tmp_$type==type_)
    loc_ <- loc_-1
    loc_ <- loc_[loc_>0]
    #y_ <- c(tmp_$y[1],(tmp_$y[-1]+tmp_$y[-nrow(tmp_)])/2)
    r[loc_] <- tmp_$y[loc_]
    xts(x=r,order.by=tmp_$dateTimeRound)
  }
  rise.xts <- define_event(1)
  fall.xts <- define_event(2)
  null.xts <- define_event(0)
  dk.xts <- define_event(8)
  ts_ <- cbind(y.xts,rise.xts,fall.xts,null.xts,dk.xts,i.xts)
  plot_ <- dygraph(ts_) %>%
    dySeries("null.xts",color="green",label="null",
             stepPlot = TRUE, fillGraph = TRUE) %>%
    dySeries("rise.xts",color="blue",label="rise",
             stepPlot = TRUE,fillGraph = TRUE) %>%
    dySeries("fall.xts",color="red",label="fall",
             stepPlot = TRUE,fillGraph = TRUE) %>%
    dySeries("dk.xts",color="black",label="dk",
             stepPlot = TRUE,fillGraph = TRUE)%>%
    dySeries("i.xts",color="white",axis = "y2") %>%
    dyAxis("y2",valueRange = c(0,1))%>%dyRangeSelector()%>% 
    dyEvent(x=break_,strokePattern = "dashed",color = "red")
  
  if(any(overlap_)){
    plot_ <- plot_ %>% dyEvent(x=break2_,strokePattern = "dashed",
                               color="blue")
  }
  if(any(overlap2_)){
    plot_ <- plot_ %>% dyEvent(x=break3_,strokePattern = "solid",
                               color="blue")
  }
  plot_
}

example_ <- function(){
  rm(list=ls())
  dummy_ <- sapply(dir("."),source)
  load("../data/qdat.RData")
  ## delineate storm for the first 5000 discharge values
  ## use 5000 to save time.
  r_ <- hydrograph(qdat[1:5000,],verbose = FALSE)
  storm_review <- data.frame(station=1,storm=1,peak1=2,peak2=3)
  r1_ <- update_hydrograph(r_,storm_review)
  r2_ <- base_hydrograph(r1_,verbose = F)
  ## show the baseflow event with largest yield.
  plotBase(r2_,which.max(r2_$base$y))
  
}
