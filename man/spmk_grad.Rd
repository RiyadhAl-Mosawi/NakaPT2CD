\name{spmk_grad}
\alias{spmk_grad}
\title{Compute the gradient of Process Capability Index (PCI) using Spmk function of Nakagami distribution}
\description{Compute the gradient of Process Capability Index (PCI) using Spmk function of Nakagami distribution}
\usage{spmk_grad(para, l, u, t, m, s)}
\arguments{
  \item{para}{: c(shape,scae) parameters}
  \item{l}{: The lower specification limit}
  \item{u}{: The upper specification limit}
  \item{t}{: The target value}
  \item{m}{: mean of the process}
  \item{s}{: standrad devation of the process}
}
\value{a vector represents the gradient of Spmk function of nakagami distribution}
\author{Riyadh Al-Mosawi}
\examples{spmk_grad(c(1,1), 1, 10, 1, 1, 1)}
