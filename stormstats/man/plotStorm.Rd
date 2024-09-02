\name{plotStorm}
\alias{plotStorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a storm}
\description{
This function generates a time series plot of a storm. Colors denote the rising and falling limbs of the storm, as well as intermediate periods between the limbs.
}
\usage{
plotStorm(x, i)
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Colors denotes rising and falling limbs. There are two other kinds of runs. The null run is a run with minor changes in discharge. The unknown run is a run with possible rising or falling limb.
}
\value{
This function returns a \code{dygraph} object, which can be visualized in R. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
## delineate storm for the first 5000 discharge values
## use 5000 to save time.
r_ <- hydrograph(qdat[1:5000,],verbose = FALSE)
plotStorm(r_,1)
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
