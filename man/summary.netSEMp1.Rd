% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.netSEMp1.R
\name{summary.netSEMp1}
\alias{summary.netSEMp1}
\title{Summary of netSEMp1}
\usage{
\method{summary}{netSEMp1}(object, ...)
}
\arguments{
\item{object}{An object of class "netSEM", the returned list from netSEMp1().}

\item{...}{A S3 generic/method consistency.}
}
\value{
NULL. A summary of data and fitting result is printed on screen.
}
\description{
summarize netSEMp1 result
}
\details{
summary.netSEMp1 gives a summary about the netSEMp1 analysis.
}
\examples{
data(acrylic)
ans <- netSEMp1(acrylic)
summary(ans)
}
