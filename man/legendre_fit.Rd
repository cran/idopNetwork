% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{legendre_fit}
\alias{legendre_fit}
\title{generate curve based on legendre polynomials}
\usage{
legendre_fit(par, x)
}
\arguments{
\item{par}{vector of legendre polynomials coefficients}

\item{x}{vector equal to the x value for legendre polynomials(in this case times)}
}
\value{
the polynomials value
}
\description{
generate curve based on legendre polynomials
}
\examples{
legendre_fit(rep(1,5),1:14)
}
