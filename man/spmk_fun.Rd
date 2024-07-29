\name{spmk_fun}
\alias{spmk_fun}
\title{Compute Process Capability Index (PCI) using Spmk function of Nakagami distribution}
\description{Compute Process Capability Index (PCI) using Spmk function of Nakagami distribution}
\usage{spmk_fun(para, l, u, t, m, s)}
\arguments{
  \item{para}{: c(shape,scae) parameters}
  \item{l}{: The lower specification limit}
  \item{u}{: The upper specification limit}
  \item{t}{: The target value}
  \item{m}{: mean of the process}
  \item{s}{: standrad devation of the process}
}
\value{a numeic value represents the Spmk function}
\author{Riyadh Al-Mosawi}
\examples{spmk_fun(c(1,1), 1, 10, 1, 1, 1)}
