\name{hydrograph}
\alias{hydrograph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Storm delineation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function delineates storms using minimum change criterion. Each storm needs to be manually reviewed. 
}
\usage{
hydrograph(data, abstol = c(0.02, 0.1), reltol = c(0.01, 0.05), k = 2,
  q_base_th = 0.05, gap_rise = 4 * 3, gap_fall = c(1, 4) * 4,
  rel_valley_th = c(0.5, 0.7), gap = 3, min_len = 12, alpha = 0.975,
  peak_q_th = 1, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
A \code{data.frame} containing the variables \code{dateTimeRound} (time stamp) and \code{y} (the water yield). Missing values are allowed in \code{y}.
}
  \item{abstol}{
%%     ~~Describe \code{abstol} here~~
the absolute yield change thresholds to detect change, first mainly for detecting changes during baseflow period for small change, and second mainly for detecting changes during multi-peak events.
}
  \item{reltol}{
%%     ~~Describe \code{reltol} here~~
the relative yield change thresholds to detect changes, first mainly for detecting changes during baseflow period for small change, and second mainly for detecting changes during multi-peak events.
}
  \item{k}{
%%     ~~Describe \code{k} here~~
a constant (larger than 1) for the ratio between the slopes of the rising limb over the falling limb. Default is 2, When classifying changes, use a more specific criterion (smaller threshold) to identify null period when the changes are negative, rather than positive.
}
  \item{q_base_th}{
%%     ~~Describe \code{q_base_th} here~~
the absolute yield threshold below which to consider as baseflow, mainly used to determine the end of the falling limb.
}
  \item{gap_rise}{
%%     ~~Describe \code{gap_rise} here~~
a null denotes a run of discharge that does not exceed the threshold to be a rising or falling limb. This is a constant for the null between rises, all runs of length below this number are considered true null. Default is 12, that is around 3 hours.
}
  \item{gap_fall}{
%%     ~~Describe \code{gap_fall} here~~
the thresholds, the first is to consider falling limbs as for sure in between peaks of a single event, the second is to cut the falling limb after peak discharge, mainly to consider the next storm a discrete event and to determine the end of the falling limb, the second default 4 hours chosen based on Miller et al. (2021) Hydrological Processes. 2021;35:e14287, the first is to identify potential multiple events lumped as a single one.
}
  \item{rel_valley_th}{
%%     ~~Describe \code{rel_valley_th} here~~
when the valley is less than the first threshold of the peak, it is considered deep, when the valley is greater than the second threshold of the peak, it is considered shallow, otherwise, it's un-determined, mainly to classify undetermined end of storm events.
}
  \item{gap}{
%%     ~~Describe \code{gap} here~~
the yield time series contain missing values, this parameter determines the short NA values to linearly interpolate from adjacent values, in the frequency of the original time series, for 15 minutes, this implies interpolation of a gap at most 45 minutes.
}
  \item{min_len}{
%%     ~~Describe \code{min_len} here~~
 the minimum number of non-missing flow per run to start looking for events. Due to missing values, the yield time series is broken into runs, if the run is shorter than this number, it will be skipped.
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
the parameter for the Lynne-Hollick filter to filter out base-flow.  
}
  \item{peak_q_th}{
%%     ~~Describe \code{peak_q_th} here~~
the minimum threshold of Q that indicates a storm (in yield unit), default was chosen based on Smith et al. (2015) Water Resources Research, 51(3),1552-1569.
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
logical indicating whether to print debugging texts.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Missing values of short duration are filled-in using linear interpolation. The extent of interpolation is controlled by the function argument \code{gap}.

This function may be used to identify incorrectly classifed complext storms. Those storms tend to show multiple peaks and valleys. The algorithm automatically determine whether a valley represents the end of an event, or just a period in-between the same limb. These events should be manually reviewed to determine whether to accept the end of event decisions. The \code{plotStorm} function can be used to aid the review process.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...event=event,peak=peak_,data=data,misc
This function returns a \code{list}
\item{event}{A \code{data.frame} with each row defines consecutive values of positive, negative, or null (close to zero) values in the raw discharge data.
  \itemize{
    \item run: which run of non-missing values this event belongs to.
    \item start, end: the starting and ending position of this event in the run.Start implies the change (whether rise or fall) occurred between the current and the next position. Likewise end implies that the current run ended between the current and the next position.
    \item type: indicates whether this run is positive (1), negative (2), too close to zero(0), or otherwise undefined (8).
	\item data frame also contains water flow/yield statistics within each event such as total, quick flow, base flow, mean flow, standard deviation, and peak.
  }
}
\item{peak}{A \code{data.frame} with each row defines a peak of discharge which should consist of a rising limb and a falling limb. Consecutive events defined earlier make up a peak.
	\itemize{
		\item run: which run of non-missing value this peak belongs to.
		\item peak: an ID variable for this peak.
		\item storm: a preliminary assignment of the storm, which will be manually reviewed and updated later.
		\item sub: a preliminary assignment of individual sub-storm with a storm.
		\item index: a preliminary index for a sub-storm.
		\item The RiseStart and RiseEnd define the position of the rising limb. The FallStart and Fallend define the position of the falling limb. The prefixes denote whether the position is expressed in terms of the raw time series (prefix i), or the event (prefix j).
		\item eos: whether this peak satisfies the end-of-storm criteria in the algorithm. Missing eos indicates undecided end-of-storm status that requires manual review. 
		\item eor: whether this peak is ended by a missing value.
		\item data frame also contain statistics of the flow during the rising and falling limbs. 
	}
}
\item{value}{A \code{data.frame} with each row defines a storm. Consecutive peaks defined earlier make up a storm. 
	\itemize{
		\item run: which run of non-missing value this peak belongs to.
		\item storm: a preliminary ID for a multi-peak storm, which was copied from the peak data frame.
		\item sub: a preliminary ID for a sub-storm within a multi-peak storm, this was also copied from the peak data frame.
		\item index: a preliminary unique ID for each sub-storm, this was also copied from the peak data frame.
		\item iStart and iEnd: define the starting and ending positions of the storm in terms of the index in the raw time series.
		\item len: the duration of the storm expressed in the frequency of the input time series.
		\item npeak: the number of peaks included in the storm.
		\item nevent: the number of unknown eos assignments in this sub-storm.
		\item this data frame also contains the statistics of the storm such as the q the mean flow, base the mean base flow, and init the initial flow.
	}
}
\item{data}{A \code{data.frame} of the input data, with added variable \code{type} showing the kind of change in discharge.}
\item{misc}{A \code{list} of entry \code{singleton}, which is a \code{data.frame} showing singleton events with large discharge and yet was not included in any storm, and \code{multi}, which is a \code{data.frame} showing storm events with multiple peaks and valleys.}
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
