\name{Gen_Data}
\alias{Gen_Data}
\title{
Function to generate type II proregressive data from Nakgami distribution
}
\description{
Thi function uses Aggarwala and Balakrishnan (1998) to genarte type II progressive sample from Nakagami distribution.
}
\usage{
Gen_Data(nu,xi, R)
}
\arguments{
  \item{nu,xi}{(shape,scale) parameters of Nakagami distribution}
  \item{R}{The vector of removals scheme}
}
\value{
returns a vector consist of failure times from first-failure type I progressive sample  Nakgami distribution}
\references{
Aggarwala, R. and Balakrishnan, N. (1998). Some properties of progressive censored order statistics from arbitrary and uniform distributions with applications to inference and simulation. Journal of Statistical Planning and Inference, 70:35–49.}
\author{
Riyadh Al-Mosawi
}
\examples{
nu=1
xi=1
R=c(0,0,0,0,10)
Gen_Data(nu,xi, R) 
}
