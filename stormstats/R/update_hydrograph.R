update_hydrograph <- function(x,y){
  ##x: output from hydrograph or a list of output from hydrograph 
  ##   if a list, it must be un-named.
  ##y: manual review to update automatically assigned storm boundaries
  ##   a data frame with station, storm id and peaks to flip
  ## value: updated output incorporating the results
  ##        either a list of hydrograph objects or a single object.
  
  is_single <- !is.null(names(x))
  if(is_single){
    if((names(x)[1]=="value")&(names(x)[2]=="event")&
       (names(x)[3]=="peak")&(names(x)[4]=="data")&
       (names(x)[5]=="misc")){
      is_single <- TRUE
      ## this is a single station
      x <- list(x)
    }else{
      warning("input argument is wrong")
      return(x)
    }
  }
  for(k in 1:length(x)){
    cat("station ",k,"\n")
    y_ <- subset(y,station==k)
    x_ <- x[[k]]
    ## correct selected automated determined boundaries of multi-peak
    ## storms
    for(i in 1:nrow(y_)){
      #browser()
      line_ <- with(x_$misc$multi,which(storm==y_$storm[i]&is.na(eos)))
      col_id <- grep("peak",names(y_))
      peak_ <- na.omit(unlist(y_[i,col_id]))
      cat("revision ",i,"storm",y_$storm[i],"n=",length(line_),
          "peak=",peak_,"\n")
      if(any(peak_>length(line_))){
        warning("wrong peaks found, not updating.")
      }else{
        est_ <- x_$misc$multi[line_[peak_],"eos2"]
        x_$misc$multi[line_[peak_],"eos2"] <- !est_
      }
    }
    
    #browser()
    ## create unique sub storm ids
    multi_ <- x_$misc$multi
    lst3 <- split(multi_,multi_$storm)
    storage2 <- vector("list",length(lst3))
    for(i in 1:length(lst3)){
      currS_ <- 1
      storm_id_ <- rep(NA,nrow(lst3[[i]]))
      for(j in 1:nrow(lst3[[i]])){
        storm_id_[j] <- currS_
        if(lst3[[i]]$eos2[j]){
          currS_ <- currS_+1
        }
      }
      lst3[[i]]$sub <- storm_id_
    }
    multi_ <- do.call(rbind,lst3)
    
    ## update Peak data frame to possibly divide 
    ## multi-peak events
    peak_ <- x_$peak
    peak_$sub <- NULL
    peak_$index <- NULL
    peak_update_ <- multi_[,c("run","peak","storm","sub")]
    peak2_ <- merge(peak_,peak_update_,by=c("run","peak"),
                    suffixes = c("","_b"),all.x = T)
    #with(peak2_,summary(storm-storm_b))
    peak2_$storm_b <- NULL
    peak2_$sub[is.na(peak2_$sub)] <- 1
    
    ## summarize each storm
    summarize_storm <- function(elmt){
      storm_ <- elmt$storm[1]     ## storm id
      sub_ <- elmt$sub[1]  ## sub-storm id
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
    lst4 <- split(multi_,list(storm=multi_$storm,sub=multi_$sub),drop=T)
    storage <- lapply(lst4,summarize_storm)
    # for(i in 1:length(lst4)){
    #   #cat("i=",i,"\n")
    #   storm_ <- lst4[[i]]$storm[1]     ## storm id
    #   sub_ <- lst4[[i]]$sub[1]  ## sub-storm id
    #   start3_ <- min(lst4[[i]]$iRiseStart) ## starting point
    #   end3_ <- max(lst4[[i]]$iFallEnd) ## ending point
    #   run_ <- lst4[[i]]$run[1]  ## run of valid yield values
    #   start4_ <- min(lst4[[i]]$jRiseStart) ## starting run
    #   end4_ <- max(lst4[[i]]$jFallEnd) ## ending run
    #   base_ <- lst4[[i]]$base0[1] ## base-flow estimate
    #   dur_ <- end3_-start3_+1 ## duration (in unit of time frequency)
    #   w_ <- with(lst4[[i]],iFallEnd-iRiseStart+1) 
    #   q_ <- sum(lst4[[i]]$q*w_)/sum(w_) ## average field
    #   qbase_ <- sum(lst4[[i]]$qBase*w_)/sum(w_) ## average base yield
    #   ## maximum distance between peaks
    #   dist_ <- 0
    #   if(nrow(lst4[[i]])>1){
    #     dist_ <- max(diff(lst4[[i]]$iRiseEnd))
    #   }
    #   ## ignore the unknown peak at the end of run
    #   tmp_ <- subset(lst4[[i]],!eor)
    #   nevent_ <- sum(is.na(tmp_$eos))
    #   
    #   storage[[i]] <- data.frame(
    #     storm=storm_,sub=sub_,iRise=start3_,iFall=end3_,
    #     npeak=nrow(lst4[[i]]),len=dur_,dist=dist_,nevent=nevent_,
    #     run=run_,jRise=start4_,jFall=end4_,
    #     q=q_,base=qbase_,init=base_
    #   )
    # }
    
    ## remove multi-events from the original results
    value_ <- subset(x_$value,!(storm %in% unique(multi_$storm)))
    value_$index <- NULL
    #value_$sub <- 1 
    
    ## prepare new results
    value_2 <- rbind(value_,do.call(rbind,storage))
    value_2 <- value_2[order(value_2$storm,value_2$sub),]
    value_2$index <- seq(1,nrow(value_2))
    
    ## update storm id in peak data frame
    peak_2 <- merge(peak2_,value_2[,c("storm","sub","index")])
    #any(is.na(peak_2$index))
    
    ## store results
    x_$value <- value_2
    x_$peak <- peak_2
    x[[k]] <- x_
  }
  if(is_single){
    return(x[[1]])
  }
  x
}
example_ <- function(){
  rm(list=ls())
  library(xts);library(dygraphs)
  dummy_ <- sapply(dir("."),source)
  load("../data/qdat.RData")
  ## delineate storm for the first 5000 discharge values
  ## use 5000 to save time.
  r_ <- hydrograph(qdat[1:5000,],verbose = FALSE)
  plotStorm(r_,1)
  ## suppose a reviewer decided to 
  ## change the second and third delineations
  storm_review <- data.frame(station=1,storm=1,peak1=2,peak2=3)
  r1_ <- update_hydrograph(r_,storm_review)
  plotStorm(r1_,2)
  ## note the second storm was cut out
  
}
