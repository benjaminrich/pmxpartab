% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmxpartab.R
\name{parse_parameter_description}
\alias{parse_parameter_description}
\title{Parse a parameter description string}
\usage{
parse_parameter_description(string)
}
\arguments{
\item{string}{A \code{character}, starting with a name, followed by an optional
comma separated list of key-value pairs that can be parsed as R code. See
Examples.}
}
\value{
A named \code{list}.
}
\description{
Parse a parameter description string
}
\examples{
# Example 1: all elements present
x <- "CL, label='Clearance', units='L/h', trans=exp, type='Structural'"
parse_parameter_description(x)

# Example 2: Some elements missing (trans), will take default value (NULL)
x <- "CL, label='Clearance', units='L/h', type='Structural'"
parse_parameter_description(x)

# Example 3: Only the name is given
x <- "CL"
parse_parameter_description(x)

# Example 4: positional arguments
x <- "CL, 'Clearance', 'L/h', type='Structural'"
parse_parameter_description(x)
}
