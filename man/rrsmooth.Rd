\name{rrsmooth}
\alias{rrsmooth}
\title{Relative Risk Adaptive Kernel density estimation}
\description{
Generates relative risk estimates of case-control data defined on polygonal regions.  Uses adaptive, boundary corrected kernel density estimates. Automatically
selects the global bandwidth (if necessary) via cross validation, and has boundary bias correction performed via an
adaptive boundary kernel. Bias is of order h0^4 in the interior, and of order h0^2 for points on or near the boundary.
}
\usage{
rrsmooth(case_points, control_points, boundary, h0=0, delta0=0, logarithmic=TRUE, kernel="gaussian", type="adaptivecorrected", nx=50, ny=50, plot=FALSE)
}
\arguments{
\item{case_points}{a vector of x,y coordinates.}
\item{control_points}{a vector of x,y coordinates.}
\item{boundary}{a polygon boundary file (list of vertices of boundary)}
\item{h0}{the (global) smoothing bandwidth to use. If h0=0 it is chosen automatically via cross validation}
\item{delta0}{the correction added to the kernel density estimates of the numerator and denominator to flatten the density in regions where no data resides. If delta0=0 it is chosen automatically as max(g)/20, where g is the density estimate of the controls.}
\item{logarithmic}{whether we should create an estimate for log(f/g) or f/g, where f and g are estimates of cases and controls respectively. Defaults to FALSE.}
\item{kernel}{a character string specifying the type of kernel to use.  Choices are "gaussian" or "biweight"}
\item{type}{a character string specifying the type of estimate used. Choices are "fixed", "fixedcorrected", "fixedcorrected0", "adaptive", "adaptivecorrected", and "adaptivecorrected0"}
\item{nx}{the grid resolution in the x direction at which the estimate is evaluated}
\item{ny}{the grid resolution in the y direction at which the estimate is evaluated}
\item{plot}{whether or not adsmooth should plot the resulting estimate}
}
\value{
\item{p}{the relative risk estimate, in (x,y,z) coordinates}
}
\note{
Large data sets (eg > 1000 points) will take a reasonable length of time (more than 30 seconds) to generate the
"adaptivecorrected" estimate if cross validation is to be used to generate the global bandwidth.

The "corrected0" estimates are estimates formed by dividing each kernel by the amount of weight contained within the boundary. "corrected" estimates use a linear boundary kernel.
}
\examples{
## load in some data from the splancs library
library(adsmooth)
library(splancs)
data(southlancs)
## compute adaptive relative risk estimate and plot directly
p = rrsmooth(southlancs[southlancs.cc==1,1:2], southlancs[southlancs.cc==0,1:2], southlancs.bdy, h0=3000)
## plot the boundary
polymap(southlancs.bdy,border="grey")
## plot as a contour plot
contour(p, nlevels=10, add=TRUE,drawlabels=FALSE)
## overlay the case and control points
pointmap(southlancs.pts[southlancs.cc == 0,], add=TRUE, pch=19, col="green")
pointmap(southlancs.pts[southlancs.cc == 1,], add=TRUE, pch=19, col="red")
title("Relative risk map using adaptive smoothing.  Cases in red, controls in green")
}
\author{
Jonathan Marshall \email{J.C.Marshall@massey.ac.nz}
}
\keyword{smooth}