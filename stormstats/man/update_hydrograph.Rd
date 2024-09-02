\name{update_hydrograph}
\alias{update_hydrograph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Update storm delineation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function incorporates manual review into the automated storm delineation. This function thus allows semi-automated storm delineation.
}
\usage{
update_hydrograph(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
output from \code{hydrograph}, or a list of output from \code{hydrograph}. If it is a list, it must be un-named.
}
  \item{y}{
%%     ~~Describe \code{verbose} here~~
a \code{data.frame} coding manual review results, and to update automatically assigned storm boundaries. The variables include station ID, storm ID, and the peak ID within each storm. Each automated decision to end a storm (or not) is flipped after manual review.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Line types are used to visualize the automated delineation results for complex storms with multiple peaks. A solid line denotes the end of a storm peak. A dashed line denotes part of a rising or falling limb. Manul review results can be stored in a spreadsheet before loading it into R.

}
\value{
A list of objects structured according to \code{hydrograph}, each item of the list is one station. The following variables in the data frames \code{value} and \code{peak} from each \code{hydrograph} object are updated according to the manual review.

\item{storm}{the multi-peak storm iD}
\item{sub}{the ID for individual event within the storm}
\item{index}{the ID for the individual event}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(qdat)
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
