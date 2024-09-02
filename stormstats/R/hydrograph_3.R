hydrograph <- function(
  data, abstol = c(0.02,0.10),reltol = c(0.01,0.05),k=2,
  q_base_th = 0.05, gap_rise=4*3, gap_fall = c(1,4)*4,
  rel_valley_th =c(0.5,0.7),
  gap = 3,   min_len = 12, alpha = 0.975, 
  peak_q_th =1,verbose=TRUE
){
  asp <- 8/5 ## determine angles of the changes internally
  ## documentation see hydrograph.docx
  ## remove duplicated time periods
  data <- subset(data,!duplicated(dateTimeRound))
  data$index <- seq(1,nrow(data))
  
  ## identify runs of missing values
  run_ <- rle(x=is.na(data$y))
  run_start <- cumsum(c(1,run_$lengths[-length(run_$lengths)]))
  run_NA_start <- run_start[run_$values]
  run_NA_len <- run_$lengths[run_$values]
  
  ## identify short runs of missing values
  short_run_NA_start <- run_NA_start[run_NA_len<=gap]
  short_run_NA_len <- run_NA_len[run_NA_len<=gap]
  
  if(length(short_run_NA_start)>0){
    ## if there are some short runs, linearly interpolate
    ## the missing yields.
    impute_NA_loc <- sequence(short_run_NA_len,from=short_run_NA_start)
    approx_yield <- approx(x=data$index,y=data$y,xout = impute_NA_loc)
    data$y[impute_NA_loc] <- approx_yield$y
  }
  
  ## separate out baseflow and quickflow (in terms of yield)
  bfi_out <- BFI(data$y,alpha = alpha,ReturnQbase = TRUE)
  data$base <- bfi_out$Qbase
  data$quick <- data$y - data$base
  
  #browser()
  ## estimate abstol, q_base_th and peak_q_th if not provided
  if(is.null(abstol)){
    base_y_ <- subset(data,base>0.9*y)$y
    base_delta_y <- 1*max(base_y_)
    abstol <- c(1,5)*base_delta_y
    if(verbose) cat("abstol estimated to be ",abstol,"\n")
  }
  
  if(is.null(q_base_th)){
    q_base_th <- 2*max(base_y_)
    if(verbose) cat("q_base_th estimated to be", q_base_th,"\n")
  }
  
  if(is.null(peak_q_th)){
    storm_y_ <- subset(data,base<0.1*y&y>q_base_th)$y
    peak_q_th <- 1.5*min(storm_y_)
    if(verbose) cat("peak_q_th estimated to be",peak_q_th,"\n")
  }
  
  # with(data,table(is.na(y),is.na(base)))
  ## 408 missing Baseflow, thus the change analysis
  ## has to base on the total flow data.
  
  ## compute the change in total flow per time unit
  ## (the baseflow change is not ideal due to NAs)
  data$delta <- c(NA,diff(data$y))
  
  ## identify runs of valid values
  run_ <- rle(x=is.na(data$y))
  run_start <- cumsum(c(1,run_$lengths[-length(run_$lengths)]))
  run_y_start <- run_start[!run_$values]
  run_y_len <- run_$lengths[!run_$values]
  
  ## consider segments of flow with minimum number of valid
  ## observations
  run_y_start <- run_y_start[run_y_len>=min_len]
  run_y_len <- run_y_len[run_y_len>=min_len]
  run_y_n <- length(run_y_len)
  
  ## type of hydrograph movement
  ## 0: no change
  ## 1: rising limb
  ## 2: falling limb
  ## 8: unknown change type
  data$type <- NA 
  
  decision_tree <- function(x,y){
    ## x: flow time series
    ## y: change of flow time series
    ## parameters: abstol, reltol,k defined in the functions
    ## value: 1= increase, 2=decrease, 8=no practical change
    ## use k to differentiate the thresholds for null 
    ## between the rises and the falls.
    
    delta_ <- y[-1] ## all valid changes
    base_ <- x[-length(x)] ## base of each change
    p_ <- delta_/base_ ## proportion of change
    r <- rep(NA,length(delta_)) ## results
    
    for(i in 1:length(r)){
      if(delta_[i]==0){
        r[i] <- 0
      }
      else if(delta_[i]>0){
        ## the change is potentially in the rising limb
        if(p_[i]<reltol[1]){
          if (delta_[i]<abstol[2]) {
            ## if percent change is small but actual change is median
            ## estimate as no change
            r[i] <- 0
          }else{
            ## if percent change is small, but actual change is large
            ## estimate the unknown change
            r[i] <- 8
          }
        }else if(p_[i]<reltol[2]){
          if(delta_[i]<abstol[1]){
            ## if relative is median, and actual change is small,
            ## estimate as no change ~0
            r[i] <- 0 
          }else if (delta_[i]<abstol[2]) {
            ## if percent changes and actual changes are median
            ## estimate as unknown change
            r[i] <- 8
          }else{
            ## if percent change is median, but actual change is large
            ## estimate the change
            r[i] <- 1 
          }
        }
        else{
          if(delta_[i]<abstol[1]){
            ## if relative change is large, but the actual change
            ## is small, estimate as no change ~0
            r[i] <- 0
          }else{
            ## if relative change is large, and the actual change
            ## is median or large, estimate as change
            r[i] <- 1 #2-as.integer(delta_[i]>0)
          }
        }
      }else{
        ## the change is potentially in the falling limb
        if(abs(p_[i])<reltol[1]/k){
          if (abs(delta_[i])<abstol[2]/k) {
            ## if percent change is small but actual change is median
            ## estimate as no change
            r[i] <- 0
          }else{
            ## if percent change is small, but actual change is large
            ## estimate the unknown change
            r[i] <- 8
          }
        }else if(abs(p_[i])<reltol[2]/k){
          if(abs(delta_[i])<abstol[1]/k){
            ## if relative is median, and actual change is small,
            ## estimate as no change ~0
            r[i] <- 0 
          }else if (abs(delta_[i])<abstol[2]/k) {
            ## if percent changes and actual changes are median
            ## estimate as unknown change
            r[i] <- 8
          }else{
            ## if percent change is median, but actual change is large
            ## estimate the change
            r[i] <- 2 #2-as.integer(delta_[i]>0)
          }
        }
        else{
          if(abs(delta_[i])<abstol[1]/k){
            ## if relative change is large, but the actual change
            ## is small, estimate as no change ~0
            r[i] <- 0
          }else{
            ## if relative change is large, and the actual change
            ## is median or large, estimate as change
            r[i] <- 2 #2-as.integer(delta_[i]>0)
          }
        }
      }
    }
    r
  }
  
  ## within each run of valid yield values
  if(verbose) cat("Identifying segments of limbs from valid yields.\n")
  storage <- vector("list",run_y_n)
  for(i in 1:run_y_n){
    #browser()
    #cat("i=",i,"\n")
    ## locate the run of valid yields
    run_y_loc_ <- seq(run_y_start[i],len=run_y_len[i])
    
    #browser()
    ## within the run identify the type of changes
    change_type_y_ <- decision_tree(
      x=data$y[run_y_loc_],
      y=data$delta[run_y_loc_]
    )
    ## store change
    data$type[run_y_loc_] <- c(NA,change_type_y_)
    
    ## estimate runs of changes
    run_change_y_ <- rle(x=change_type_y_)
    ## number of "event"s ~ unchange, rising or falling limbs
    n_events_ <- length(run_change_y_$lengths)
    ## locate the start of each "event"
    lstart_ <- cumsum(c(1,run_change_y_$lengths[-n_events_]))
    run_y_start_ <- run_y_loc_[lstart_]
    ## locate the end of each "event
    lend_ <- cumsum(run_change_y_$lengths)
    run_y_end_ <- run_y_loc_[lend_]
    ## estimate the duration of each "event"
    run_y_duration_ <- run_change_y_$lengths
    
    ## calculate the characteristics of each event
    ## cumulative total, quick and base flows
    ## average change in total flow and standard deviation
    yield_mat_ <- matrix(NA,length(run_y_start_),6)
    avg_helper <- function(x,...){
      if(length(x)>1){
        r <- mean((x[-length(x)]+x[-1])/2,...)
      }else{
        r <- mean(x,...)
      }
      return(r)
    }
    for(j in 1:length(run_y_start_)){
      #cat("j=",j,"\n")
      ## for each run of flow
      tmp_ <- seq(run_y_start_[j],run_y_end_[j])
      ## average of the total, quick and base flow
      ## measure of the size of the event
      yield_mat_[j,3] <- avg_helper(data$y[tmp_])
      yield_mat_[j,2] <- avg_helper(data$quick[tmp_],na.rm=T)
      yield_mat_[j,1] <- avg_helper(data$base[tmp_],na.rm=T)
      ## average change of the total flow
      ## measure of the change of the flow.
      local_yield_ <- na.omit(data$y[tmp_])
      if(length(local_yield_)>1){
        yield_mat_[j,4] <- mean(diff(local_yield_))
      }
      if(length(local_yield_)>2){
        yield_mat_[j,5] <- sd(diff(local_yield_))
      }
      yield_mat_[j,6] <- max(data$y[tmp_])
    }
    ## store results
    storage[[i]] <- data.frame(
      run=i,event=1:n_events_,
      start=run_y_start_,
      len=run_y_duration_,
      end=run_y_end_,
      type=run_change_y_$values,
      total=yield_mat_[,3],
      quick=yield_mat_[,2],
      base=yield_mat_[,1],
      mean=yield_mat_[,4],
      sd=yield_mat_[,5],
      peak=yield_mat_[,6]
    )
  }
  
  ## all segments   
  event <- do.call(rbind,storage)
  
  ## retrieve starting time of each event
  event <- merge(event,data[,c("index","dateTimeRound")],
                 by.x="start",by.y="index",all.x=T)
  #any(is.na(event$dateTimeRound))
  
  ## processing each runs of valid values  
  l <- split(event,event$run)
  
  ## store intermediate results
  storage <- vector("list",length(l))
  for(i in 1:length(l)){
    if(verbose) cat("processing run",i,"length=",sum(l[[1]]$len),".\n")
    #if(i==10) browser()
    ## all potential rising limb locations
    loc_ <- which(l[[i]]$type==1)
    ## remove locations at the end
    ## where the next event is unknown
    loc_ <- loc_[(loc_+1)<=nrow(l[[i]])]
    
    if(length(loc_)==0){
      if(verbose)  cat("No potential rising limbs in this run.\n")
      next
    }
    
    if(verbose) cat("Found",length(loc_),"potential rising limbs.\n")
    
    ## start with the first potential rising limb in the run
    curr_ <- loc_[1]
    
    ## temporary storage
    storms_lst <- vector("list",length(loc_))
    j <- 1  ## peak number
    
    ##browser()
    update_base_flow <- TRUE ## flag
    while(any(loc_>=curr_)){
      ## work through the potential rising limbs
      
      ## record the entering flow
      if(update_base_flow){
        flow0_ <- data$y[l[[i]]$start[curr_]-1]
      }
      
      start_ <- curr_ ## start of the rising limb
      
      ## read next event until reaching the end
      ## of run or a falling limb
      if(verbose) cat("rising limb",j,"starts at event",curr_,"\n")
      next_type_ <- l[[i]]$type[curr_+1]
      if(!is.na(next_type_)&(next_type_==2)){
        if(verbose) cat("  event at ", curr_+1," type =",next_type_,
                        ".\n")
      }
      while(next_type_ %in% c(0,1,8)){
        if(verbose) cat("  event at ", curr_+1," type =",next_type_,
                        "continue.\n")
        curr_ <- curr_+1
        if(next_type_==0){
          if(l[[i]]$len[curr_]>gap_rise){
            if(verbose) cat("  null event length=",l[[i]]$len[curr_],
                            "skip to next rising limb.\n")
            #browser()
            next_type_ <- 99
            next
          }
        }
        next_type_ <- l[[i]]$type[curr_+1]
        if(!is.na(next_type_)&(next_type_==2)){
          if(verbose) cat("  event at ", curr_+1," type =",next_type_,
                          ".\n")
        }
      }
      
      if(is.na(next_type_)){
        if(verbose) cat("end of run",i,".\n")
        next
      }
      
      if(next_type_==99){
        ## a null between rising limbs that is long,
        ## likely a false null.
        if(any(loc_>curr_)){
          curr_ <- min(loc_[loc_>curr_])
        }else{
          curr_ <- max(loc_)+1
        }
        
        if(j>1){
          ## not the first event
          ## if there are pending unknown EOS status?
          ## flush this as a EOS before skipping.
          if(is.na(end_of_storm)){
            if(verbose) cat("Reset EOS at event ",j-1,"from NA to TRUE.\n")
            storms_lst[[j-1]]$eos <- TRUE 
          }
          
          ## if the previous event is not EOS,
          ## flush it as a EOS before skipping
          if(!is.na(end_of_storm)&(!end_of_storm)){
            if(verbose) cat("Reset EOS at event ",j-1,"from FALSE to TRUE.\n")
            storms_lst[[j-1]]$eos <- TRUE 
          }
        }
        
        next
      }
      
      end_ <- curr_ ## end of the rising limb
      
      
      if(verbose) cat("rising limb",j,"from event",start_,"to ",end_,
                      ",base Q=",flow0_,"\n")
      
      ## move to the falling limb  
      curr_ <- curr_+1
      
      ## record the start of the falling limb
      start2_ <- curr_
      if(verbose) cat("falling limb",j,"starts at event",curr_,"\n")
      
      ## read the next event until reaching the end of run 
      ## or another rising limb
      next_type_ <- l[[i]]$type[curr_+1]
      if(!is.na(next_type_)&(next_type_==1)){
        if(verbose) cat("  event at ", curr_+1," type =",next_type_,
                        ".\n")
      }
      while(next_type_ %in% c(0,2,8)){
        if(verbose) cat("  event at ", curr_+1," type =",next_type_,
                        "continue.\n")
        curr_ <- curr_+1
        next_type_ <- l[[i]]$type[curr_+1]
        if(!is.na(next_type_)&(next_type_==1)){
          if(verbose) cat("  event at ", curr_+1," type =",next_type_,
                          ".\n")
        }
      }
      
      ## record the ending of the falling limb
      end2_ <- curr_
      if(verbose) cat("falling limb",j,"from event",start2_,"to",end2_,".\n")
      
      ## extract flow of the falling limb
      falling_limb_event_ <- l[[i]][seq(start2_,end2_),]
      len_ <- sum(falling_limb_event_$len)
      idx_ <- with(l[[i]],seq(start[start2_],len=len_+1))
      y_falling_limb_ <- data$y[idx_]
      min_ <- min(y_falling_limb_)
      
      ## determine falling limb
      th__ <- max(q_base_th,flow0_,na.rm = T)
      if(min_<th__){
        ## if falling limb back to the precedent flow condition
        end_of_storm <- TRUE
        idx_ <- idx_[1:min(which(y_falling_limb_<=th__))]
        if(length(idx_)>gap_fall[2]){
          idx_ <- idx_[1:gap_fall[2]]
          if(verbose) cat("falling limb",j,"truncated to baseflow and length at",
                          length(idx_),"EOS.\n")
        }else{
          if(verbose) cat("falling limb",j,"truncated to baseflow at",
                          length(idx_),"EOS.\n")
        }
        update_base_flow <- TRUE
      }else{
        if(len_<gap_fall[1]){
          ## if falling limb is short (w.r.t. gap_fall[1])
          end_of_storm <- FALSE
          if(verbose) cat("falling limb",j," not an EOS.\n")
        }
        else if(len_<gap_fall[2]){
          ## if falling limb is median (w.r.t. gap_fall)
          end_of_storm <- NA
          ## unknown whether it's end of storm
          if(verbose) cat("falling limb",j," may be an EOS.\n")
        }
        else{
          ## if falling limb is long (w.r.t. gap_fall)
          end_of_storm <- TRUE
          idx_ <- idx_[1:gap_fall[2]]
          if(verbose) cat("falling limb",j,"truncated by length",
                          length(idx_),"EOS.\n")
        }
        update_base_flow <- FALSE
      }
      
      ## characteristics of the falling limb
      if(!is.na(end_of_storm) &(end_of_storm)){
        ## truncation of the rising limb occurred
        ## re-read the falling limb.
        y_falling_limb_ <- data$y[idx_]
        min_ <- min(y_falling_limb_)
        len_ <- length(idx_)
      }
      
      ## whether this is the end of run
      eor_flag <- is.na(next_type_)
      
      if(eor_flag){
        if(!is.na(end_of_storm)&(!end_of_storm)){
          if(verbose) cat("skipping partial falling limb",j,".\n")
          next
        }
      }
      
      ## falling limb characteristics
      q_ <- quantile(y_falling_limb_,probs = 0.10)
      total_ <- avg_helper(y_falling_limb_)
      mean_ <- sd_ <- NA
      if(length(y_falling_limb_)>1){
        mean_ <- mean(diff(y_falling_limb_))
      }
      if(length(y_falling_limb_)>2){
        sd_ <- sd(diff(y_falling_limb_))
      }
      
      
      ## record this event
      if(verbose) cat("falling limb",j,"len=",len_,"min=",min_,
                      "change=",mean_,"\n")
      
      ## extract location of the rising limb
      rising_limb_event_ <- l[[i]][seq(start_,end_),]
      idx2_ <- with(rising_limb_event_,
                    seq(from=start[1],len=sum(len)+1))
      
      ## characteristics of the rising limb
      y_rising_limb_ <- data$y[idx2_]
      total2_ <- avg_helper(y_rising_limb_)
      mean2_ <- sd2_ <- NA
      if(length(y_rising_limb_)>1){
        mean2_ <- mean(diff(y_rising_limb_))
      }
      if(length(y_rising_limb_)>2){
        sd2_ <- sd(diff(y_rising_limb_))
      }
      
      ## baseflow
      base_ <- data$base[c(idx_,idx2_)]
      mean3_ <- avg_helper(base_,na.rm=T)
      mean4_ <- avg_helper(c(y_rising_limb_,y_falling_limb_))
      
      
      
      
      if(!eor_flag){
        curr_ <- curr_+1
      }
      
      if(eor_flag){
        if(verbose) cat("end of run ",i,"\n")
      }
      
      ## if this is the last event and is a rising limb
      if(curr_==nrow(l[[i]])){
        eor_flag <- TRUE
        if(verbose) cat("end of run ",i,"\n")
      }
      
      storms_lst[[j]] <- data.frame(
        run=i,peak=j,
        iRiseStart=idx2_[1],iRiseEnd=idx2_[length(idx2_)],
        jRiseStart=start_,jRiseEnd=end_,
        iFallStart=idx_[1],iFallEnd=idx_[length(idx_)],
        jFallStart=start2_,jFallEnd=end2_,
        eos=end_of_storm,
        base0=flow0_, nFall=len_,minFall=min_,q10Fall=q_,
        qFall=total_,meanDeltaFall=mean_,sdDeltaFall=sd_,
        qRise=total2_,meanDeltaRise=mean2_,sdDeltaRise=sd2_,
        qBase=mean3_,q=mean4_,eor=eor_flag,
        qinit=y_rising_limb_[1],qpeak=y_falling_limb_[1],
        qend=y_falling_limb_[length(y_falling_limb_)]
      )
      ## work on the next event
      j <- j+1
    }
    
    storage[[i]] <- do.call(rbind,storms_lst)
  }
  
  peak_ <- do.call(rbind,storage)
  
  #with(peak_,table(eor,eos))
  
  ## calculate storm number by runs
  currS_ <- 1
  lst2 <- split(peak_,peak_$run)
  for(i in 1:length(lst2)){
    storm_id_ <- rep(NA,nrow(lst2[[i]]))
    for(j in 1:nrow(lst2[[i]])){
      storm_id_[j] <- currS_
      if(!(is.na(lst2[[i]]$eos[j]))&(lst2[[i]]$eos[j])){
        ## if current falling limb is an end-of-storm
        currS_ <- currS_+1
      }
      
      ## last falling limb is not sure a storm
      ## treat as one
      unknown_last_ <- with(lst2[[i]],is.na(eos[j])&(eor[j]))
      if(unknown_last_){
        currS_ <- currS_+1
      }
    }
    lst2[[i]]$storm <- storm_id_
  }
  peak_ <- do.call(rbind,lst2)
  
  ## To identify singleton rise/fall events
  
  ## all rising and falling limbs
  all_rise_fall_ <- subset(event,type %in% c(1,2,8))
  
  ## all rising and falling limbs in the storms
  lst3 <- split(all_rise_fall_,all_rise_fall_$run)
  for(i in 1:length(lst3)){
    work2_ <- subset(peak_,run==lst3[[i]]$run[1])
    lst3[[i]]$single <- TRUE
    if(nrow(work2_)>0){
      jrise_ <- with(work2_,sequence(nvec=jRiseEnd-jRiseStart+1,
                                     from=jRiseStart))
      jfall_ <- with(work2_,sequence(nvec=jFallEnd-jFallStart+1,
                                     from=jFallStart))
      lst3[[i]]$single <- !(lst3[[i]]$event %in% c(jrise_,jfall_))
    }
  }
  all_rise_fall_ <- do.call(rbind,lst3)
  
  ## singleton events not covered by the identified storms
  singleton_ <- subset(all_rise_fall_,single)
  #summary(singleton_$peak)
  
  singleton_ <- subset(singleton_,peak > peak_q_th) 
  if(nrow(singleton_)){
    warnings(nrow(singleton_)," singleton events identified")
  }
  
  ## summarize each storm
  summarize_storm <- function(elmt){
    storm_ <- elmt$storm[1]     ## storm id
    sub_ <- 1  ## sub-storm id
    start3_ <- min(elmt$iRiseStart) ## starting point
    end3_ <- max(elmt$iFallEnd) ## ending point
    run_ <- elmt$run[1]  ## run of valid yield values
    start4_ <- min(elmt$jRiseStart) ## starting run
    end4_ <- max(elmt$jFallEnd) ## ending run
    base_ <- elmt$base0[1] ## base-flow estimate
    dur_ <- end3_-start3_+1 ## duration (in unit of time frequency)
    w_ <- with(elmt,iFallEnd-iRiseStart+1) 
    q_ <- sum(elmt$q*w_)/sum(w_) ## average field
    qbase_ <- sum(elmt$qBase*w_)/sum(w_) ## average base yield
    ## maximum distance between peaks
    dist_ <- 0
    if(nrow(elmt)>1){
      dist_ <- max(diff(elmt$iRiseEnd))
    }
    ## ignore the unknown peak at the end of run
    tmp_ <- subset(elmt,!eor)
    nevent_ <- sum(is.na(tmp_$eos))
    
    data.frame(
      storm=storm_,sub=sub_,iStart=start3_,iEnd=end3_,
      npeak=nrow(elmt),len=dur_,dist=dist_,nevent=nevent_,
      run=run_,jStart=start4_,jEnd=end4_,
      q=q_,base=qbase_,init=base_
    )
    
  }
  #browser()
  ## intermediate storage for storms
  lst4 <- split(peak_,peak_$storm)
  storage2 <- vector("list",length(lst4))
  storage3 <- vector("list",length(lst4))
  for(i in 1:length(lst4)){
    #cat("i=",i,"\n")
    storage2[[i]] <- summarize_storm(lst4[[i]])
    ## ignore the unknown peak at the end of run
    tmp_ <- subset(lst4[[i]],!eor)
    nevent_ <- sum(is.na(tmp_$eos))
    
    if(nevent_>0) {
      #cat("possible multi-peak storm=",storm_,"\n")
      #browser()
      
      ## storage
      storage3[[i]] <- lst4[[i]]
      storage3[[i]]$storm <- lst4[[i]]$storm[1]
      
      ## X-axis width in time units
      nrise_ <- with(lst4[[i]],iRiseEnd-iRiseStart+1)
      nfall_ <- with(lst4[[i]],nFall)
      ntotal_ <- nrise_+nfall_
      X_ <- sum(ntotal_)
      
      ## Y-axis height in yield units
      idx_ <- with(lst4[[i]],sequence(ntotal_,from=iRiseStart))
      local_mp_yield_ <- data$y[idx_]
      Y_ <- max(local_mp_yield_)-min(local_mp_yield_)
      
      ## constant to change yield unit to graph unit
      k_ <- (X_/Y_)*(1/asp)
      
      angle_helper <- function(idx,k){
        ## idx: location of the limb
        ## return: the magnitude of the angle of the limb and x-axis
        y_ <- data$y[idx]
        dy_ <- y_[length(idx)]-y_[1]
        dx_ <- length(idx)
        if(dy_<0){
          return(-atan(k*(dy_/dx_)))
        }else{
          return(atan(k*(dy_/dx_)))
        }
      }
      
      ## Rising limb angles 
      angle_rise_ <- rep(NA,nrow(lst4[[i]]))
      for(j in 1:length(angle_rise_)){
        rloc_ <- with(lst4[[i]],seq(iRiseStart[j],iRiseEnd[j]+1))
        angle_rise_[j] <- angle_helper(rloc_,k_)
      }
      
      ## falling limb angles
      angle_fall_ <- rep(NA,nrow(lst4[[i]]))
      for(j in 1:length(angle_fall_)){
        floc_ <- with(lst4[[i]],seq(iFallStart[j]-1,iFallEnd[j]))
        angle_fall_[j] <- angle_helper(floc_,k_)
      }
      
      storage3[[i]]$angle_rise <- angle_rise_
      storage3[[i]]$angle_fall <- angle_fall_
      
      ## relative change of angles from the previous rising limbs
      ## to the next rising limb
      rel_change_angle_ <- rep(NA,nrow(lst4[[i]]))
      for(j in 1:(length(rel_change_angle_)-1)){
        rel_change_angle_[j] <- log(angle_rise_[j+1]/angle_rise_[j])
      }
      storage3[[i]]$delta_angle_rise <- rel_change_angle_
      
      ## relative magnitudes of each valley versus the
      ## previous previous peaks
      rel_prev_peak_ <- rep(NA,nrow(lst4[[i]]))
      for(j in 1:length(rel_prev_peak_)){
        qpeak_ <- with(lst4[[i]],max(qpeak[1:j]))
        rel_prev_peak_[j] <- lst4[[i]]$qend[j]/qpeak_
      }
      ## relative magnitudes of each valley versus the
      ## following peaks
      rel_next_peak_ <- rep(NA,nrow(lst4[[i]]))
      for(j in 1:(length(rel_next_peak_)-1)){
        qpeak2_ <- with(lst4[[i]],max(qpeak[(j+1):length(rel_next_peak_)]))
        rel_next_peak_[j] <- lst4[[i]]$qend[j]/qpeak2_
      }
      
      decision_tree_multi <- function(x,y,z){
        ## x: valley divided by previous peak
        ## y: valley divided by next peak
        ## z: change of angles between current and next rising limb
        ## value: preliminary reclassification of eos
        ##      TRUE= eos, FALSE= not eos
        if(y>rel_valley_th[2]){
          ## next peak is shallow, this is not an EOS regardless
          ## of the previous peak
          return(FALSE)
        }else if(y<rel_valley_th[1]){
          if(x<rel_valley_th[1]){
            ## next peak is tall, and previous peak is also tall
            ## an EOS
            return(TRUE)
          }else if(x<rel_valley_th[2]){
            ## next peak is tall, and previous peak is neither
            ## tall or shallow
            return(abs(z)>0.3)
          }else{
            ## next peak is tall, and previous peak is shallow
            return(TRUE)
          }
        }
        else{
          if(x<rel_valley_th[1]){
            ## next peak is un-determined
            ## and previous peak is tall, not an EOS
            return(FALSE)
          }else if(x<rel_valley_th[2]){
            ## next peak is un-determined
            ## previous peak is also un-determined,
            return(abs(z)>0.3)
          }
          else{
            ## next peak is un-determined,
            ## previous peak is shallow
            return(TRUE)
          }
        }
      }
      #browser()
      eos_ <- lst4[[i]]$eos
      for(j in 1:nrow(lst4[[i]])){
        if(lst4[[i]]$eor[j]){
          next
        }
        if(is.na(eos_[j])){
          eos_[j] <- decision_tree_multi(
            rel_prev_peak_[j],rel_next_peak_[j],
            rel_change_angle_[j])
        }
      }
      storage3[[i]]$rel_prev_peak <- rel_prev_peak_
      storage3[[i]]$rel_next_peak <- rel_next_peak_
      storage3[[i]]$eos2 <- eos_
    }
    

  }
  multi_peak <- do.call(rbind,storage3)
  
  misc <- list(singleton=singleton_,multi=multi_peak)
  storm_ <- do.call(rbind,storage2)
  storm_$index <- storm_$storm
  peak_$sub <- 1
  peak_$index <- peak_$storm
  list(value=storm_,event=event,peak=peak_,data=data,misc=misc)
}
example_ <- function(){
  rm(list=ls())
  dummy_ <- sapply(dir("."),source)
  load("../data/qdat.RData")
  r_ <- hydrograph(qdat[1:5000,],verbose = FALSE)
}

