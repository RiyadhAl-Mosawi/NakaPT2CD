\name{Estim}
\alias{Estim}
\title{Compute the MLE (or MPSE) using Nelder-Mead algorithm of Nakagami distribution under type II progressive censoring data}
\description{
Compute the MLE (or MPSE) using Nelder-Mead algorithm of Nakagami distribution under type II progressive censoring data. Information matrix using observed information principle is also computed.
}
\usage{
Estim(True_Par,X,R,l,u, t,gm,para,type,lw,up, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{True_Par}{(shape,scale,PCI) True values parameters of Nakagami distribution}
  \item{X}{Failure times}
  \item{R}{The vector of removals scheme}
  \item{l}{Lower specification}
  \item{u}{upper specification}
  \item{T}{Target}
  \item{gm}{gamma parameter}
  \item{type}{"LK" for MLE and "PS" for MPSE}
}
\value{
 returns a list which includes
  \item{par}{estimated values}
  \item{value}{value of log-likelihood function at the estimated value}
  \item{inform}{infomation matrix using observed information principle}
}
\author{
Riyadh Al-Mosawi
}
\examples{
Estim(True_Par,X,R,l,u, t,gm,para,type,lw,up, method)
{
    .Call(`_NDLibraryFunction_Estim`,X,R,l,u, t,gm,para,type,lw,up, method)
  }
}
