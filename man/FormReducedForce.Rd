% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microTimoshenkoBeam.R
\name{FormReducedForce}
\alias{FormReducedForce}
\title{Reduced force vector.}
\usage{
FormReducedForce(forcevector)
}
\arguments{
\item{forcevector}{A vector of applied nodal loads (force,moments etc)
supplied in the form c(F1,M1, etc)}
}
\value{
A column matrix of the reduced force vector
}
\description{
Forms the reduced stiffness matrix after
         applying the boundary condition(s).
}
