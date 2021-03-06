% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccm_summary.R
\name{ccm_summary}
\alias{ccm_summary}
\title{ccm_summary}
\usage{
ccm_summary(ccm_output, predtype = "rho")
}
\arguments{
\item{ccm_output}{output of the convergent cross mapping.}

\item{predtype}{A character indicating which metric is used as predictive skill. Defaults to "rho". Can be "rho", "mae", or "rmse".}
}
\value{
list of mean, lower and upper confidence intervals of ccm iterations for the predictive skill used.
}
\description{
Computes the summary statistics of the convergent cross mapping output.
}
\keyword{ccm}
\keyword{rEDM,}
