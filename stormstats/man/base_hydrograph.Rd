\name{base_hydrograph}
\alias{base_hydrograph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Baseflow delineation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function delineates baseflows between storms events. 
}
\usage{
base_hydrograph(x, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A \code{list} structured according to \code{hydrograph}
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
Logical indicating whether to show debug texts.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function divides the time series except those included in the storm delineations. The time series data must be padded with missing values to avoid gaps in between adjacent time stamps.

This function may be used to identify incorrectly classifed storms. If the discharge in those baseflow exceeds certain threshold, it should be manually reviewed to determine whether to combine this event with an adjacent storm flow, and which one to combine. The \code{plotBase} function can be used to aid the review process.
}
\value{
A \code{list} structured according to function \code{hydrograph}, with an entry named \code{base} with each row defines a potential baseflow event, containing the following variables.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
	\item{index}{an ID variable of each base flow}
	\item{y}{This is the mean of the raw values within the baseflow period. A missing value of y denotes that no measurement is recorded during the period, which could be either true baseflow or missed storm flows.}
	\item{Start,Stop}{The Start and Stop variables indicate the position of the baseflow event,which are expressed in terms of the indices in the raw time series or the time stamps}
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
## delineate storm for a subset of 5000 discharge values
## to save time.
r_ <- hydrograph(qdat[10000+1:5000,],verbose = FALSE)
storm_review <- data.frame(station=1,storm=5,peak1=1)
r1_ <- update_hydrograph(r_,storm_review)
## prepare the base-flow events
r2_ <- base_hydrograph(r1_,verbose = FALSE)
summary(r2_$base$y)
## base-flow has low water yield, this is good.
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
