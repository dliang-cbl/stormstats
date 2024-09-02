base_hydrograph <- function(x,verbose=TRUE){
  ## x: output from update hydrograph
  ## value: a data frame of base-flow periods
  data_pad_ <- x$data
  
  ## identify runs of valid values
  run_ <- rle(x=is.na(data_pad_$y))
  run_start <- cumsum(c(1,run_$lengths[-length(run_$lengths)]))
  run_NA_start <- run_start[run_$values]
  run_NA_len <- run_$lengths[run_$values]
  run_NA_stop <- run_NA_start+run_NA_len-1
  run_NA_start_date <- data_pad_$dateTimeRound[run_NA_start]
  run_NA_stop_date <- data_pad_$dateTimeRound[run_NA_stop]
  if(length(run_NA_len)>0){
    run_NA_check_ <- rep(NA,length(run_NA_len))
    for(i in 1:length(run_NA_len)){
      run_NA_check_[i] <- all(is.na(
        data_pad_$y[run_NA_start[i]:run_NA_stop[i]]))
    }
    stopifnot(all(run_NA_check_))
    
    #browser()
    
    res_NA <- data.frame(
      iStart=run_NA_start,iStartDate=run_NA_start_date,
      iStop=run_NA_stop,iStopDate=run_NA_stop_date,
      len=run_NA_len,y=NA)## q: average flow/yield during this period
    
  }else{
    ## no missing value
    res_NA <- NULL
  }
  
  run_y_start <- run_start[!run_$values]
  run_y_len <- run_$lengths[!run_$values]
  run_y_stop <- run_y_start+run_y_len-1
  
  if(length(run_y_len)==0) stop("no data.")
  
  ## identify all time periods covered in "storms"
  idx_in_storms <- with(x$value,sequence(nvec=len,from=iStart))
  
  storage <- vector("list",length(run_y_len))
  for(i in 1:length(run_y_len)){
    ## go through each run
    ## replace the storm events there with NA
    ## because the storm does not involve gaps
    ## this should break the run of yiled into segments
    ## each segment will reprent the baseflow period
    stopifnot(all(!is.na(data_pad_$y[run_y_start[i]:run_y_stop[i]])))
    ## identify time periods of storms in this run
    lrun_idx <- run_y_start[i]:run_y_stop[i]
    lstorm_flag <- lrun_idx %in% idx_in_storms
    ## extract the yield values
    ly_ <- data_pad_$y[lrun_idx]
    ## extract the dates
    ldate_ <- data_pad_$dateTimeRound[lrun_idx]
    ## replace yield in the storms with NA
    ly_[lstorm_flag] <- NA
    ## the following routine just create the valid
    ## yields within this run
    lrun_ <- rle(x=!is.na(ly_))
    lrun_start <- cumsum(c(1,lrun_$lengths[-length(lrun_$lengths)]))
    lrun_start <- lrun_idx[lrun_start]
    lrun_base_start <- lrun_start[lrun_$values]
    lrun_base_len <- lrun_$lengths[lrun_$values]
    lrun_base_stop <- lrun_base_start+lrun_base_len-1
    lrun_base_start_date <- data_pad_$dateTimeRound[lrun_base_start]
    lrun_base_stop_date <- data_pad_$dateTimeRound[lrun_base_stop]
    if(verbose){
      cat("run ",i,"from ",format(ldate_[1],"%y-%m-%d %H:%M"),
          "to",format(ldate_[length(ldate_)],"%y-%m-%d %H:%M"),
          ":",sum(lrun_$values),"baseflow events.\n")
    }
    lq_ <- rep(NA,length(lrun_base_start))
    for(j in 1:length(lrun_base_start)){
      lq_[j] <- mean(data_pad_$y[lrun_base_start[j]:lrun_base_stop[j]])
    }
    storage[[i]] <- data.frame(
      iStart=lrun_base_start,iStartDate=lrun_base_start_date,
      iStop=lrun_base_stop,iStopDate=lrun_base_stop_date,
      len=lrun_base_len,y=lq_)## q: average flow/yield during this period
    
  }
  value_ <- do.call(rbind,storage)
  stopifnot(all(!is.na(value_$q)))
  value_ <- rbind(res_NA,value_)
  
  Checking_code <- function(){
    check_base_idx_ <- sequence(nvec=value_$len,from=value_$iStart)
    any(duplicated(check_base_idx_))
    check_storm_idx_ <- with(x$value,sequence(nvec=len,from=iRise))
    check_storm_idx_ <- unique(check_storm_idx_) ## remove overlapping start and end of multi-peak events.
    check_ <- sort(c(check_base_idx_,check_storm_idx_))
    table(diff(check_))
    dim(data_pad_)
  }
  x$base <- value_[order(value_$iStart),]
  x$base$index <- seq(1,nrow(value_))
  return(x)
}

example_ <- function(){
  rm(list=ls())
  library(xts);library(dygraphs)
  dummy_ <- sapply(dir("."),source)
  load("../data/qdat.RData")
  ## delineate storm for a subset of 5000 discharge values
  ## to save time.
  r_ <- hydrograph(qdat[10000+1:5000,],verbose = FALSE)
  storm_review <- data.frame(station=1,storm=5,peak1=1)
  r1_ <- update_hydrograph(r_,storm_review)
  ## prepare the base-flow events
  r2_ <- base_hydrograph(r1_,verbose = F)
  summary(r2_$base$y)
  ## base-flow has low water yield, this is good.
}
