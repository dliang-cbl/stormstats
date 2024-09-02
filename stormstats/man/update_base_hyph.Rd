\name{update_base_hygh}
\alias{update_base_hygh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Update storm delineation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function incorporates manual review of baseflow events into the automated storm delineation. This function thus allows semi-automated storm delineation. Specifically high flow events in the automatically delineated baseflow events are combined into the adjacent storm events.
}
\usage{
update_base_hygh(x,y,verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A \code{list} of lists, each list is structured according to the \code{hydrograph} function.

}
  \item{y}{
%%     ~~Describe \code{verbose} here~~
A \code{data.frame} of review results showing the baseflow events to be combined with adjacent storm events.
}
  \item{verbose}{
  logical, whether to print debug texts.
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
The \code{data.frame} object with each row indicates an action on one baseflow event. The action reads combine \code{base} and possibly \code{base2}, and possibly convert storm \code{storm2}, into the \code{rise} limb of storm \code{storm1}.
\itemize{
	\item{station: an integer of station id}
	\item{base: an integer of the baseflow event id}
	\item{storm1: the target storm ID into which the baseflow should be appended}
	\item{rise: The target limb of the target storm, 0 indicates falling limb, 1 indicates rising limb}
	\item{base_2: sometimes an adjacent baseflow event should also be combined, this variable denotes the additional baseflow event id}. 
	\item{storm2: sometimes an adjacent storm event should be converted into a limb of the target storm event. This variable denotes the additional stormflow event ID.}
}
}
\value{
A \code{list} of objects containing \code{value} and \code{base}, shows the reviewed storm and baseflow delineation results. Each item of the \code{list} denotes one station. The item is structured according to the \code{hydrograph}
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
r_ <- hydrograph(qdat[1:5000,],verbose = FALSE)
storm_review <- data.frame(station=1,storm=1,peak1=2,peak2=3)
r1_ <- update_hydrograph(r_,storm_review)
r2_ <- base_hydrograph(r1_,verbose = FALSE)
## show the baseflow event with largest yield.
plotBase(r2_,which.max(r2_$base$y))
## support we want to combine the center baseflow event with 
## the preceding storm flow
base_review <- data.frame(station=1,base=3,storm1=5,storm2=NA,base_2=NA,rise=0)
r3_ <- update_base_hygh(list(r2_),base_review)
## storm 5 before the update
plotStorm(r2_,5)
## storm 5 after the update, note the baseflow 3
## was combined with storm 5.
plotStorm(r3_[[1]],5)
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
