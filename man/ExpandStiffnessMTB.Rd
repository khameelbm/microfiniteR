% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microTimoshenkoBeam.R
\name{ExpandStiffnessMTB}
\alias{ExpandStiffnessMTB}
\title{Expanded stiffness matrix of a microscale Timoshenko beam}
\usage{
ExpandStiffnessMTB(total_DOF, eMatrix, i, j)
}
\arguments{
\item{total_DOF}{Total degree of freedom in a
connected system of beams.}

\item{eMatrix}{Element's stiffness matrix.}

\item{i}{First node of an element.}

\item{j}{Second node of an element.}
}
\value{
The expanded stiffness matrix
                     a microscale Timoshenko beam element.
}
\description{
Accepts the element's stiffness matrix and returns the expanded form
of the matrix to conform with the global degree of freedom.
}
