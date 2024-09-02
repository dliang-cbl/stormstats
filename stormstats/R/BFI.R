BFI <-
function(Q, alpha=0.925, passes=3, ReturnQbase=TRUE, n.reflect=30) {
  
  ############ Define Functions #########################  
  
  ############# first pass ################################
  
  # This function takes a vector of surface flows and returns a list containing
  # baseflows Qbase and quick flows Qquick
  #
  # Note: Qquick can be negative (if you are planning to use to the Qquick output other than
  # for input to the next stage of filtering, these negative values will need to be set to zero)
  #
  
  
  FirstPass <- function(Q,a) {  
    Qf1 <-  vector('numeric', length=length(Q))
    
    Qf1[1] <- Q[1] # set the initial value equal to the surface flow  
    
    for(i in 2:length(Q)) {
      Qf1[i] <- a*Qf1[i-1]+0.5*(1+a)*(Q[i]-Q[i-1])
    }
    
    Qb1 <- ifelse(Qf1 > 0, Q-Qf1,Q )
    return(data.frame(Qquick = Qf1, Qbase=Qb1))  
    
  }
  
  
  
  
  
  ##################### Backward pass #########################
  
  # This function takes a vector of quick flows and a vector of baseflows and retuns a list containing
  # baseflows and quickflows after a backward pass of the Lyne and Hollick Filter
  #
  # Note: Qquick can be negative (if you are planning to use to the Qquick out put other than
  # for input to the next stage of filtering, these negative values will need to be set to zer)
  #
  
  BackwardPass <- function(Q,a) {
    
    # expecting a dataframe, Q that contains baseflow and quick flow
    
    Qq <- Q$Qquick
    Qb <- Q$Qbase
    num.rows <- nrow(Q)
    
    Qf2 <- vector('numeric', length=num.rows)
    
    # set the last value (which is the first one the algorithm encounters equal to
    # the base flow from the previous pass through
    
    Qf2[num.rows] <- Qb[num.rows]
    
    for(i in (num.rows-1):1) {
      Qf2[i] <- a*Qf2[i+1]+0.5*(1+a)*(Qb[i]-Qb[i+1])
      
    }
    
    Qb2 <- ifelse(Qf2 >0, Qb-Qf2,Qb )
    
    return(data.frame(Qquick = Qf2, Qbase=Qb2))
  }
  
  
  ##################### forward pass #########################
  
  # This function takes a vector of quick flows and a vector of baseflows and retuns a list containing
  # baseflows and quickflows after a forward pass of the Lyne and Hollick Filter
  #
  # Note: Qquick can be negative (if you are planning to use to the Qquick out put other than
  # for input to the next stage of filtering, these negative values will need to be set to zero)
  #
  
  ForwardPass <- function(Q,a) {  
    
    Qq <- Q$Qquick
    Qb <- Q$Qbase
    num.rows <- length(Qq)
    
    Qf2 <- vector('numeric', length=num.rows)
    
    # set the first value equal to the first baseflow from the previous time step
    
    Qf2[1] <- Qb[1]
    for(i in 2:num.rows) {
      Qf2[i] <- a*Qf2[i-1]+0.5*(1+a)*(Qb[i]-Qb[i-1])
      
    }
    Qb2 <- ifelse(Qf2 >0, Qb-Qf2,Qb )
    
    return(data.frame(Qquick = Qf2, Qbase=Qb2))
  }
  
  
  
  ################## BFI.calc #######################################
  
  BFI.calc <- function(Q, alpha, passes, n.reflect) {
    
    # we reflect the first n.reflect values and the last n.reflect values.  this is to 
    # get rid of 'warm up' problems
    
    # Save the input data  
    Qin <- Q  
    
    # create a vector that includes the flow and reflection at either end
    
    Q.reflect <- vector(mode='numeric', length=length(Q) + 2*n.reflect)
    Q.reflect[1:n.reflect] <- Q[(n.reflect+1):2]
    Q.reflect[(n.reflect+1):(n.reflect + length(Q))] <- Q
    Q.reflect[(n.reflect+length(Q)+1):(length(Q)+2*n.reflect)] <- Q[(length(Q)-1):(length(Q)-n.reflect)]
    
    
    # run the first pass followed by a number of backwards/forwards passes
    Q1 <- FirstPass(Q.reflect,alpha) # 1
    
    # how many backwards/forward passes to we need
    n.pass <- round((passes-1)/2)
    
    BackwardPass(Q1,alpha)
    
    for(i in 1:n.pass){
      Q1 <- ForwardPass(BackwardPass(Q1,alpha),alpha)
    }
    
    
    ################# end of passes  ##############################
    
    
    # chop the first and last n.reflect values off so we are
    # back with the original number of data points
    # Return the Baseflow vector
    
    Qbase <- Q1$Qbase[(n.reflect+1):(length(Q1$Qbase)-n.reflect)]
    Qbase[Qbase < 0] <- 0  # set any values less than zero to zero
    
    #  Qbase <- Q1$Qbase # Use this if you want the reflected flow to be included
    
    BFI.value <- sum(Qbase)/sum(Qin)
    
    return(list(BFI = BFI.value, Qbase = Qbase))
  }
  
  
  
  ############## Flow.index ###########################
  
  # this function is used to deal with cases where there are missing values
  
  Flow.index <- function(Q, n.reflect){
    # Description 
    
    # Compute indexes to the start and end of
    # sequences of flow where there are at least n.reflect consecutive non-missing values
    
    # Usage
    # Flow.index(Q, len=31)
    
    # Arguements
    #
    # Q     a vector of flows
    # len   length of the sequence of non-missing values
    
    # Details
    # Based on p32 Spector(2008) Data manipulation with R. Springer
    
    # Value
    #
    # List with indicies to Q specifying start and end of sequences of non-missing flows >= to len
    
    Flow.rle <- rle(is.na(Q))
    
    # Check if there is any non-missing flow sequences of length n.reflect or greater (stop if not)
    
    if(!any(!Flow.rle$values & Flow.rle$lengths > n.reflect)) stop('Must have at least 31 consecutive non-missing flow values \n')
    
    #ends
    index <- which(!Flow.rle$values & Flow.rle$lengths > n.reflect )
    ends <- cumsum(Flow.rle$lengths)[index]
    
    #starts
    newindex <- ifelse(index > 1, index-1,0)
    starts <- cumsum(Flow.rle$lengths)[newindex] +1
    if(0 %in% newindex) starts =c(1,starts)
    
    #   list(starts=starts,ends=ends)
    
    return(list(starts=starts,ends=ends))
  }
  
  
  #
  ################ Main Routine ##########################################
  
  
  
  
  # A small amount of error checking
  # passes must be an odd number greater than or equal to 3
  
  if(passes %% 2 == 0 | passes < 3) stop('passes must be odd and greater than 2')
  
  # alpha must be between zero and 0
  
  if(alpha < 0 | alpha >=1 ) stop('alpha must be between zero and one \n') 
  
  # We need more than n.reflect values
  if(length(Q) <= n.reflect ) stop(paste('n.reflect must be <= length(Q)', ' \n'))
  
  
  # if there are not any missing  values, just cal BFI.Calc and return the BFI
  if(!any(is.na(Q))) {
    BFI.data <- BFI.calc(Q,  alpha, passes, n.reflect) 
    BFI.out=BFI.data$BFI
    FractionUsed = 1
    Qbase=BFI.data$Qbase
  }
  else { # there are missing values
    # get the indexes that point to the start and end of segments of flow longer
    # than n.reflect non-missing valuess
    seg.index <- Flow.index(Q, n.reflect)
    
    # number of segments
    seg.num <- length(seg.index$starts)
    
    
    
    
    # create vectors to store results
    w <- numeric(seg.num) # vector to store weights
    Q.seg.BFI <- numeric(seg.num) # vector to store BFI values
    Q.seg.Qbase <- rep(NA, times=length(Q))
    
    # Loop through all segments calculating the BFI for each segment
    # and getting the weight for each segment (i.e. the number of flow values)
    for(i in 1:seg.num){
      Q.seg <- Q[seg.index$starts[i]:seg.index$ends[i]]
      w[i] <- length(Q.seg) # weights are the number of values in each segment
      BFI.data <- BFI.calc(Q.seg, alpha, passes, n.reflect)
      Q.seg.BFI[i] <- BFI.data$BFI # store the BFI associated with each segment
      Q.seg.Qbase[seg.index$starts[i]:seg.index$ends[i]] <- BFI.data$Qbase
    }
    
    BFI.out <- weighted.mean(Q.seg.BFI,w) # Weighted average of the BFIs
    
    # Portion of data used (only count non missing values)
    FractionUsed <- sum(w)/sum(!is.na(Q))
    
    Qbase=Q.seg.Qbase
    
    
    
  }
  
  if(ReturnQbase){
    return(list(BFI = BFI.out, alpha=alpha, FractionUsed = FractionUsed, Qbase=Qbase))
  } else {
    return(list(BFI = BFI.out, alpha=alpha, FractionUsed = FractionUsed))
  }
  
}

example_ <- function(){
  rm(list=ls())
  dummy_ <- sapply(dir("."),source)
  load("../data/qdat.RData")
  r <- BFI(qdat$y[1:1000])
  str(r)
}
