% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{get_legendre_matrix}
\alias{get_legendre_matrix}
\title{generate legendre matrix}
\usage{
get_legendre_matrix(x, legendre_order)
}
\arguments{
\item{x}{vector equal to the x value for legendre polynomials(in this case times)}

\item{legendre_order}{the order of legendre polynomials}
}
\value{
the polynomials value of each order
}
\description{
generate legendre matrix
}
\examples{
get_legendre_matrix(1:14,4)
}
