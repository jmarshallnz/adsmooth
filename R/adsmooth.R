estimatemise <- function(h, points, triangulation, nx, ny, grid, numberkernel, numbertype)
{
  returned_data = .C('estimate_mise_R', num_data = as.integer(nrow(points)), data = as.double(t(points)), kernel = as.integer(numberkernel), bandwidth = as.double(h), num_gridX = as.integer(nx), num_gridY = as.integer(ny), grid = as.double(grid), boundarysize = as.integer(length(triangulation) / 6), boundary = as.double(triangulation), type = as.integer(numbertype), mise = double(1), PACKAGE="adsmooth")

  return(returned_data$mise)
}

estimate_rr_mise <- function(h, d, points_f, points_g, triangulation, nx, ny, grid, numberkernel, numbertype)
{
  returned_data = .C('estimate_rr_mise_R', num_data_f = as.integer(nrow(points_f)), data_f = as.double(t(points_f)), num_data_g = as.integer(nrow(points_g)), data_g = as.double(t(points_g)), kernel = as.integer(numberkernel), bandwidth = as.double(h), delta = as.double(d), num_gridX = as.integer(nx), num_gridY = as.integer(ny), grid = as.double(grid), boundarysize = as.integer(length(triangulation) / 6), boundary = as.double(triangulation), type = as.integer(numbertype), mise = double(1), PACKAGE="adsmooth")

  return(returned_data$mise)
}

gettype <- function(type)
{
  # type of estimate
  numbertype <- -1
  if (type == "adaptivecorrected")
    numbertype <- 6
  if (type == "adaptivecorrected0")
    numbertype <- 5
  if (type == "adaptive")
    numbertype <- 4
  if (type == "fixedcorrected")
    numbertype <- 2
  if (type == "fixedcorrected0")
    numbertype <- 1
  if (type == "fixed")
    numbertype <- 0
  if (numbertype < 0)
  {
    print('Invalid type, defaulting to fixed')
    numbertype = 0
  }
  return(numbertype)
}

getkernel <- function(type)
{
  # kernel to use
  numberkernel <- -1
  if (type == "gaussian")
    numberkernel <- 2
  if (type == "biweight")
    numberkernel <- 1
  if (numberkernel < 0)
  {
    print('Invalid Kernel, defaulting to Gaussian')
    numberkernel = 2
  }
  return(numberkernel)
}

minimizemise <- function(points, boundary, type="adaptivecorrected", kernel="gaussian", nx=50, ny=50)
{
  numbertype <- gettype(type)
  numberkernel <- getkernel(kernel)

  # We start by triangulating the boundary
  tri <- t(triangulate(as(boundary, "gpc.poly")))

  # compute an appropriate grid
  minmax <- bbox(boundary)
  gx <- seq(minmax[1,1],minmax[1,2],length=nx)
  gy <- seq(minmax[2,1],minmax[2,2],length=ny)
  tempgrid <- expand.grid(y=gy,x=gx)
  grid <- matrix(NA,2,nx*ny)
  grid[1,] = tempgrid$x
  grid[2,] = tempgrid$y

  # generate a suitable starting value for minimization of mise for h
  A <- min(sd(points), IQR(points[,1])/1.34, IQR(points[,2])/1.34)
  h_start <- 0.9 * A * (length(points)/2)^(-1/5)

  # optimize the MISE using cross-validation
  h <- optimize(estimatemise, interval = c(h_start * 0.2,h_start * 3.0), points = points, triangulation = tri, nx = nx, ny = ny, grid = grid, numbertype = numbertype, numberkernel = numberkernel)

  cat('bandwidth = ', h$minimum, '\n')

  return(h$minimum)
}

minimize_rr_mise <- function(points_f, points_g, delta, boundary, type="adaptivecorrected", kernel="gaussian", nx=50, ny=50)
{
  numbertype <- gettype(type)
  numberkernel <- getkernel(kernel)

  # We start by triangulating the boundary
  tri <- t(triangulate(as(boundary, "gpc.poly")))

  # compute an appropriate grid
  minmax <- bbox(boundary)
  gx <- seq(minmax[1,1],minmax[1,2],length=nx)
  gy <- seq(minmax[2,1],minmax[2,2],length=ny)
  tempgrid <- expand.grid(y=gy,x=gx)
  grid <- matrix(NA,2,nx*ny)
  grid[1,] = tempgrid$x
  grid[2,] = tempgrid$y

  # generate a suitable starting value for minimization of mise for h
  A <- min(sd(points_g), IQR(points_g[,1])/1.34, IQR(points_g[,2])/1.34)
  h_start <- 0.9 * A * (length(points_g)/2)^(-1/5)

  if (numberkernel == 1)
    h_start <- h_start * sqrt(8)

  # optimize the MISE using cross-validation
  h <- optimize(estimate_rr_mise, interval = c(h_start * 0.2,h_start * 3.0), points_f = points_f, points_g = points_g, d = delta, triangulation = tri, nx = nx, ny = ny, grid = grid, numbertype = numbertype, numberkernel = numberkernel)

  cat('bandwidth = ', h$minimum, '\n')

  return(h$minimum)
}

adsmooth <- function(points, boundary, h0=0, type="adaptivecorrected", kernel="gaussian", nx=50, ny=50, plot=FALSE)
{
  # type of estimate
  numbertype <- gettype(type)
  numberkernel <- getkernel(kernel)

  # We start by triangulating the boundary
  tri <- triangulate(as(boundary, "gpc.poly"))

  # compute an appropriate grid
  minmax <- bbox(boundary)
  gx <- seq(minmax[1,1],minmax[1,2],length=nx)
  gy <- seq(minmax[2,1],minmax[2,2],length=ny)
  tempgrid <- expand.grid(y=gy,x=gx)
  grid <- matrix(NA,2,nx*ny)
  grid[1,] = tempgrid$x
  grid[2,] = tempgrid$y

  if (h0 == 0)
  {
    # need to generate our h0 by cross-validation
    h0 <- minimizemise(points, boundary, type, kernel, nx, ny)
  }

  # compute the estimate
  returned_data = .C('construct_pdf_R', num_data = as.integer(length(points) / 2), data = as.double(t(points)), kernel = as.integer(numberkernel), bandwidth = as.double(h0), num_gridX = as.integer(nx), num_gridY = as.integer(ny), grid = as.double(grid), boundarysize = as.integer(length(tri) / 6), boundary = as.double(t(tri)), type = as.integer(numbertype), pdf = double(nx*ny), PACKAGE="adsmooth")

  out <- list(x=gx,y=gy,z=matrix(returned_data$pdf,nx,ny,byrow=TRUE))

  if (plot)
  {
    filled.contour(out, nlevels=20, col=topo.colors(20))
  }
  return(out)
}

rrsmooth <- function(points_f, points_g, boundary, h0=0, delta0=0, type="adaptivecorrected", kernel="gaussian", logarithmic=FALSE, nx=50, ny=50, plot=FALSE)
{
  # type of estimate
  numbertype <- gettype(type)
  if (logarithmic)
    numbertype = numbertype + 8

  numberkernel <- getkernel(kernel)

  # We start by triangulating the boundary
  tri <- triangulate(as(boundary, "gpc.poly"))

  # compute an appropriate grid
  minmax <- bbox(boundary)
  gx <- seq(minmax[1,1],minmax[1,2],length=nx)
  gy <- seq(minmax[2,1],minmax[2,2],length=ny)
  tempgrid <- expand.grid(y=gy,x=gx)
  grid <- matrix(NA,2,nx*ny)
  grid[1,] = tempgrid$x
  grid[2,] = tempgrid$y

  if (h0 == 0)
  {
    # need to generate our h0 by cross-validation
    h0 <- minimize_rr_mise(points_f, points_g, delta = delta0, boundary, type, kernel, nx, ny)
  }

  # compute the estimate - TODO: This is non-optimal as 
  data_f = .C('construct_pdf_R', num_data = as.integer(length(points_f) / 2), data = as.double(t(points_f)), kernel = as.integer(numberkernel), bandwidth = as.double(h0), num_gridX = as.integer(nx), num_gridY = as.integer(ny), grid = as.double(grid), boundarysize = as.integer(length(tri) / 6), boundary = as.double(t(tri)), type = as.integer(numbertype), pdf = double(nx*ny), PACKAGE="adsmooth")
  data_g = .C('construct_pdf_R', num_data = as.integer(length(points_g) / 2), data = as.double(t(points_g)), kernel = as.integer(numberkernel), bandwidth = as.double(h0), num_gridX = as.integer(nx), num_gridY = as.integer(ny), grid = as.double(grid), boundarysize = as.integer(length(tri) / 6), boundary = as.double(t(tri)), type = as.integer(numbertype), pdf = double(nx*ny), PACKAGE="adsmooth")


  if (delta0 == 0)
  {
    # TODO: need to optimize delta0 by some fancy method
    returned_data = .C('choose_delta_R', num_data_f = as.integer(nrow(points_f)), data_f = as.double(t(points_f)), num_data_g = as.integer(nrow(points_g)), data_g = as.double(t(points_g)), kernel = as.integer(numberkernel), bandwidth = as.double(h0), num_gridX = as.integer(nx), num_gridY = as.integer(ny), grid = as.double(grid), boundarysize = as.integer(length(tri) / 6), boundary = as.double(t(tri)), type = as.integer(numbertype), delta = double(1), PACKAGE="adsmooth")
    delta0 <- returned_data$delta
    cat('delta chosen = ',delta0,'\n')
  }

  rr <- (data_f$pdf + delta0) / (data_g$pdf + delta0)
  if (logarithmic)
    rr <- log(rr)

  out <- list(x=gx,y=gy,z=matrix(rr,nx,ny,byrow=TRUE))

  if (plot)
  {
    filled.contour(out, nlevels=20, col=topo.colors(20))
  }
  return(out)
}
