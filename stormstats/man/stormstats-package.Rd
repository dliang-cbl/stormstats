\name{stormstats-package}
\alias{stormstats-package}
\alias{stormstats}
\docType{package}
\title{
\packageTitle{stormstats}
}
\description{
\packageDescription{stormstats}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{stormstats}
\packageIndices{stormstats}
A collection of tools to delineate storms from discharge time series. The main function is \code{hydrograph}.
}
\author{
\packageAuthor{stormstats}

Maintainer: \packageMaintainer{stormstats}
}
%\references{
%}
\keyword{ storm delineation }
%\seealso{
%}
\examples{
data(qdat)
## delineate storm for the first 5000 discharge values
## use 5000 to save time.
r_ <- hydrograph(qdat[1:5000,],verbose = FALSE)
plotStorm(r_,1)
## Cut of a seperate storm at break 2,
## while remove the cut at break 3.
storm_review <- data.frame(station=1,storm=1,peak1=2,peak2=3)
r1_ <- update_hydrograph(r_,storm_review)
plotStorm(r1_,2)
plotStorm(r1_,3)
## Note a second storm 2,
## and storm 3 included more falling limb.
## Now calculate baseflow
r2_ <- base_hydrograph(r1_,verbose = FALSE)
## show the baseflow event with largest yield.
plotBase(r2_,which.max(r2_$base$y))
## suppose we want to combine the center baseflow event with 
## the preceding storm flow
base_review <- data.frame(station=1,base=3,storm1=5,storm2=NA,base_2=NA,rise=0)
r3_ <- update_base_hygh(list(r2_),base_review)
## storm 5 before the update
plotStorm(r2_,5)
## storm 5 after the update, note the baseflow 3
## was combined with storm 5.
plotStorm(r3_[[1]],5)

}
