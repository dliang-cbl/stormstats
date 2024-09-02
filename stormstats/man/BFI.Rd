\name{BFI}
\alias{BFI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Baseflow seperation using the Lyne and Hollick filter
}
\description{
This function uses the Lyne and Hollick digitial filter method to calculate baseflow.  Input is a vector of streamflows.
}
\usage{
BFI(Q, alpha = 0.925, passes = 3, ReturnQbase = TRUE, n.reflect = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{
vector of stream flow e.g. daily flows
}
  \item{alpha}{
filter parameter
}
  \item{passes}{
number of times the filter is passed over the data
}
  \item{ReturnQbase}{
logical, should the baseflows be returned (TRUE) or just the Baseflow Index (FALSE)
}
  \item{n.reflect}{
number of values to reflect at the start and end of Q
}
}
\details{
The Baseflow index is the ratio of the volume of baseflow to the total flow volume. The Lyne and Hollick filter is used to calculate the baseflow at each time step. The baseflows are then summed and divided by the sum of the input flows

Missing data

Missing data are handled as follows: any continuous segments of flow data longer than n.reflect values with no missing values are dealt with separately. BFIs are calculated for each segment.The BFI for the whole record is the weighted average of the BFI for each segment. Weights are the lengths of each segment.

The FractionUsed is report which is the proportion of flow values used to calculate the BFI.  Where there are no missing values, the fraction used is 1, where there are missing values the FractionUsed may be less than the number of non-missing values because of the need for continuous flow segments to be longer than n.reflect.

}
\value{
If ReturnQbase=FALSE, the default, then a list with the following elements
\item{BFI}{Baseflow index} 
\item{alpha}{the input filter parameter}
\item{FractionUsed}{the fraction of flow values used to calculate the BFI}

If ReturnQbase=TRUE, the baseflow values are returned as well as a list element named Qbase.

}
\references{
Ladson, A. R., R. Brown, B. Neal and R. Nathan (2013) A standard approach to baseflow separation using the Lyne and Hollick filter. Australian Journal of Water Resources 17(1): 173-180.

Lyne, V. D. and M. Hollick (1979) Stochastic time-variable rainfall-runoff modelling. Hydrology and Water Resources Symposium, Perth, Institution of Engineers, Australia.

Nathan, R. J. and T. A. McMahon (1990). "Evaluation of automated techniques for base flow and recession analysis." Water Resources Research 26: 1465-1473.

Grayson, R. B., R. M. Argent, R. J. Nathan, T. A. McMahon and R. G. Mein (1996). Hydrological Recipes. Melbourne, Cooperative Research Centre for Catchment Hydrology.
}
\author{
  Tony Ladson(tony.ladson@gmail.com)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{EcoHydRology::BaseflowSeparation}
}
\examples{
data(qdat)
r <- BFI(qdat$y[1:1000])
str(r)
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
