plotStorm <- function(x,i){
  ##x: output from hydrograph function
  ##i: the storm number (i.e. index)
  multi_ <- x$misc$multi
  data <- x$data
  x <- x$value
  if(!("index" %in% names(x))){
    x$index <- x$storm
  }
  lx <- subset(x,index==i)
  if(nrow(lx)==0) return(1)
  
  seq_ <- seq(lx$iStart[1],lx$iEnd[1]+1)
  tmp_ <- data[seq_,]
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
    dyAxis("y2",valueRange = c(0,1))%>%dyRangeSelector()
  
  ## check multiple peaks
  #browser()
  storm_ <- x$storm[i]
  mp_ <- subset(multi_,storm==storm_&is.na(eos)&!eor)
  if(nrow(mp_)>0){
    mp1_ <- subset(multi_,storm==storm_&is.na(eos)&!eor&eos2)
    break_ <- data$dateTimeRound[mp1_$iFallEnd]
    mp2_ <- subset(multi_,storm==storm_&is.na(eos)&!eor&!eos2)
    break2_ <- data$dateTimeRound[mp2_$iFallEnd]
    plot_ <- plot_ %>% dyEvent(x=break_,strokePattern = "solid") %>% 
      dyEvent(x=break2_)
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
  plotStorm(r_,1)
  
}