\name{qdat}
\alias{qdat}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Example discharge data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset is downloaded from NWIS. This data contain the discharge at 15 minute frequency from the station NWIS#01589290. Missing values are maintained in the time series. There is no gap in the time stamps.
}
\usage{data("qdat")}
\format{
  A data frame with 595070 observations on the following 5 variables.
  \describe{
    \item{\code{dateTimeRound}}{a POSIXct, the rounded time stamp}
    \item{\code{n}}{a numeric vector, the number of raw discharge record within each rounded time period}
    \item{\code{q}}{a numeric vector, the median of raw discharge within each rounded time period, cubit feet per second}
    \item{\code{site_no}}{a character vector, the NWIS site number}
    \item{\code{y}}{a numeric vector, the discharge normalized by watershed drainage area, cubit meter per second per meter squared }
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The raw data were rounded to every 15 minutes, with missing values padded so that the timestamp does not contain any gap.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
https://waterdata.usgs.gov/usa/nwis/uv?site_no=01589290
}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
data(qdat)
}
\keyword{datasets}
