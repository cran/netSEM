% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.netSEM.R
\name{summary.netSEM}
\alias{summary.netSEM}
\title{Summary of net-SEM}
\usage{
\method{summary}{netSEM}(object, ...)
}
\arguments{
\item{object}{An object of class "netSEM", the returned list from netSEMm().}

\item{...}{A S3 generic/method consistency.}
}
\value{
NULL. A summary of data and fitting result is printed on screen.
}
\description{
summarize netSEM result
}
\details{
summary.netSEM gives a summary about the net-SEM analysis.
}
\examples{
data(acrylic)
ans <- netSEMm(acrylic)
summary(ans)
}
