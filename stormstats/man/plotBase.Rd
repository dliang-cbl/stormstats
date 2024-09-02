\name{plotBase}
\alias{plotBase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a base flow event}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function generates time series of base-flow. Colors are used to delineate the rising and falling limbs, and the intermediate periods.
}
\usage{
plotBase(x, i, buffer = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A \code{list} structured according to the function \code{hydrograph}
}
  \item{i}{
%%     ~~Describe \code{i} here~~
An integer indicating the storm number from the list \code{x}

}
  \item{buffer}{
%%     ~~Describe \code{buffer} here~~
An integer denoting the buffer (in the time unit of input flow time series) to include around the baseflow event.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function identified the baseflow and stormflow within the \code{buffer} of the baseflow event. The console outputs all the identified event identification numbers. Baseflows are denoted by dashed vertical lines on the time series plot, stormflows are denoted by solid vertical lines. The center baseflow is red while the adjacent events are shown in blue.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
This function returns a \code{dygraph} object, which can be visualized in R. 
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
%%  ~~who you are~~
Dong Liang (dliang@umces.edu)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
data(qdat)
r_ <- hydrograph(qdat[1:5000,],verbose = FALSE)
storm_review <- data.frame(station=1,storm=1,peak1=2,peak2=3)
r1_ <- update_hydrograph(r_,storm_review)
r2_ <- base_hydrograph(r1_,verbose = FALSE)
## show the baseflow event with largest yield.
plotBase(r2_,which.max(r2_$base$y))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
