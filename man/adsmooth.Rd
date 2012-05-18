\name{adsmooth}
\alias{adsmooth}
\title{Adaptive Kernel density estimation}
\description{
Generates kernel density estimates of data defined on polygonal regions. Automatically
selects the global bandwidth (if necessary) via cross validation, and has boundary bias correction performed via an
adaptive boundary kernel. Bias is of order h0^4 in the interior, and of order h0^2 for points on or near the boundary.
}
\usage{
adsmooth(points, boundary, h0=0, kernel="gaussian", type="adaptivecorrected", nx=50, ny=50, plot=FALSE)
}
\arguments{
\item{points}{a vector of x,y coordinates.}
\item{boundary}{a polygon boundary file (list of vertices of boundary)}
\item{h0}{the (global) smoothing bandwidth to use. If h0=0 it is chosen automatically via cross validation}
\item{kernel}{a character string specifying the type of kernel to use.  Choices are "gaussian" or "biweight"}
\item{type}{a character string specifying the type of estimate used. Choices are "fixed", "fixedcorrected", "fixedcorrected0", "adaptive", "adaptivecorrected", and "adaptivecorrected0"}
\item{nx}{the grid resolution in the x direction at which the estimate is evaluated}
\item{ny}{the grid resolution in the y direction at which the estimate is evaluated}
\item{plot}{whether or not adsmooth should plot the resulting estimate}
}
\value{
\item{p}{the kernel density estimate, in (x,y,z) coordinates}
}
\note{
Large data sets (eg > 1000 points) will take a reasonable length of time (more than 30 seconds) to generate the
"adaptivecorrected" estimate if cross validation is to be used to generate the global bandwidth.

The "corrected0" estimates are estimates formed by dividing each kernel by the amount of weight contained within the boundary. "corrected" estimates use a linear boundary kernel.
}
\examples{
## load in some data from the splancs library
library(splancs)
data(southlancs)
## compute adaptive kernel estimate and plot directly
p <- adsmooth(southlancs.pts[southlancs.cc == 0,],southlancs.bdy, plot=TRUE)
## compute the adaptive kernel without plotting, on a refined scale
p <- adsmooth(southlancs.pts[southlancs.cc == 0,],southlancs.bdy,nx=100,ny=100)
## plot the boundary
polymap(southlancs.bdy,border="grey")
## plot as a contour plot
contour(p, nlevels=10, add=TRUE,drawlabels=FALSE)
## overlay the case points
pointmap(southlancs.pts[southlancs.cc == 1,], add=TRUE, pch=19, col="green")
title("Density map using adaptive smoothing on controls.  Cases in green.")
}
\author{
Jonathan Marshall \email{J.C.Marshall@massey.ac.nz}
}
\keyword{smooth}