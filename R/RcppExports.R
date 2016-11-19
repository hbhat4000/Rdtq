.integrandmat <- function(xvec,yvec,h,driftfun,difffun)
{
  drifty = driftfun(yvec)
  diffy = difffun(yvec)
  mysd = diffy*sqrt(h)
  myvar = diffy*diffy*h

  c0 = 1/sqrt(2*pi)
  c0mod = c0/mysd
  dx = xvec[2] - xvec[1]
  bign = length(xvec)

  kvec = c(0)
  mydiags = list(NULL)

  # main diagonal
  maindiag = exp(-(h/2)*(drifty*drifty)/(diffy*diffy))*c0mod
  refsum = sum(abs(maindiag))*dx
  mydiags[[1]] = maindiag

  # super diagonals
  done = 0
  nk = 1
  mytol = 2.2e-16
  while (done==0)
  {
    curdiag = kvec[nk] + 1
    mymean = curdiag*dx + drifty*h
    thisdiag = exp(-mymean*mymean/(2*myvar))*c0mod
    thisdiag = thisdiag[-c(1:curdiag)]
    if ((nk==1) || (sum(abs(thisdiag)) > mytol*refsum))
    {
      mydiags[[nk+1]] = thisdiag
      kvec = c(kvec,curdiag)
      nk = nk + 1
    }
    else
    {
      done = 1
    }
  }
  kvec = c(kvec,-kvec[2:nk])
  # sub diagonals
  for (i in c((nk+1):length(kvec)))
  {
    curdiag = kvec[i]
    mymean = curdiag*dx + drifty*h
    thisdiag = exp(-mymean*mymean/(2*myvar))*c0mod
    thisdiag = thisdiag[-c((bign+curdiag+1):bign)]
    mydiags[[i]] = thisdiag
  }
  out = Matrix::bandSparse(bign,k=kvec,diagonals=mydiags)
  return(out)
}

.dtqsparse <- function(h, k=NULL, a=NULL, b=NULL, bigm, init, fT, drift, diffusion)
{
  numsteps = ceiling(T/h)

  if (is.null(k))
    xvec = seq(from=a,to=b,length.out=bigm)
  else
    xvec = k*c(-bigm:bigm)

  if (length(init)==1)
  {
    driftinit = drift(init)
    diffinit2 = diffusion(init)^2
    fs = exp(-(xvec-init-driftinit*h)^2/(2*diffinit2*h))/sqrt(2*pi*diffinit2*h)
    approxpdf = Matrix::Matrix(fs,ncol=1)
    startnum = 2
  }
  else
  {
    approxpdf = Matrix::Matrix(init,ncol=1)
    startnum = 1
  }
  A = k*(.integrandmat(xvec,xvec,h,driftfun=drift,difffun=diffusion))
  for (i in c(startnum:numsteps)) approxpdf = A %*% approxpdf
  return(list(xvec=xvec,pdf=as.numeric(approxpdf)))
}

#' Density Tracking by Quadrature (DTQ)
#'
#' \code{rdtq} implements DTQ algorithms to compute the probability density
#' function of a stochastic differential equation with
#' user-specified drift and diffusion functions.
#'
#' Consider the stochastic differential equation (SDE)
#'
#' \deqn{ dX(t) = f(X(t)) dt + g(X(t)) dW(t) }
#'
#' where \eqn{W(t)} is standard Brownian motion, \eqn{f} is the drift
#' function, and \eqn{g} is the diffusion function. Let \eqn{p(x,t)}
#' denote the time-dependent probability density function (PDF) of \eqn{X(t)};
#' then \code{rdtq} computes \eqn{p(x,T)} for a fixed time \eqn{T}.
#'
#' Note that the PDF is computed on a spatial grid that can be specified in
#' one of two ways:
#' \enumerate{
#'  \item Specify a real, positive value \eqn{k} and a positive integer
#' \eqn{M} = \code{bigm}. In this case, the PDF will be computed on the grid
#' \eqn{x_j = j k} where \eqn{j = -M, -M+1, ..., M-1, M}. In total,
#' there will be \eqn{2M+1} grid points.
#'  \item Specify a real, positive integer \eqn{M} and a computational domain
#' \eqn{[a,b]}. In this case, there will be exactly \eqn{M} equispaced grid
#' points. The grid spacing will be \eqn{k = (b-a)/(M-1)}.
#' }
#'
#' @param h Time step size, a positive numeric scalar.
#' @param k Spatial grid spacing, a positive numeric scalar.
#'  (Must be specified if \code{a} and \code{b} are not specified.)
#' @param bigm If k is specified, then bigm is a positive integer such that
#'  -\code{bigm*k} and \code{bigm*k} are, respectively, the minimum and maximum
#'  grid points.  If \code{a} and \code{b} are specified, then \code{bigm} is
#'  the total number of grid points.  Note that the fractional part of
#'  \code{bigm} is ignored, and that \code{floor(bigm)} must be at least 2.
#' @param a Left boundary, a numeric scalar.
#'  (Must be specified if \code{k} is not specified.)
#' @param b Right boundary, a numeric scalar.
#'  (Must be specified if \code{k} is not specified.)
#' @param init A numeric scalar indicating either a fixed initial condition of
#'  the form \eqn{X(0)}=\code{init}, or a numeric vector giving the PDF at time
#'  \eqn{t=0}. In the latter case, the vector must have the same size as the
#'  spatial grid.
#' @param fT The final time, a positive numeric scalar. The computation
#'  assumes an initial time of \eqn{t}=0 and then computes the PDF at time
#'  \eqn{t}=\code{fT}.
#' @param drift When the user chooses the \code{method="cpp"} algorithm,
#' this should be a pointer to a drift function that is implemented in C++
#' using \code{Rcpp}.  In our C++ code, we define the type \code{funcPtr}
#' using the following code:
#'
#' \code{typedef double (*funcPtr)(const double& x);}
#'
#' We expect the drift function to be a C++ function, implemented using Rcpp, of
#' type \code{XPtr<funcPtr>}.  See the first example below.
#'
#' When the user chooses the \code{method="sparse"} algorithm, this should be an
#' R function that takes as input a numeric vector of values.
#' The function should return as output a numeric vector containing the drift
#' function evaluated at each element of the input vector.
#' See the second example below.
#'
#' @param diffusion When the user chooses the \code{method="cpp"} algorithm,
#' this should be a pointer to a diffusion function that is implemented in C++
#' using \code{Rcpp}.  All of the details are analogous to that of the \code{drift}
#' function described above.
#'
#' When the user chooses the \code{method="sparse"} algorithm, this should be an
#' R function that takes as input a numeric vector of values.
#' The function should return as output a numeric vector containing the diffusion
#' function evaluated at each element of the input vector.
#' See the second example below.
#'
#' @param thresh This is an optional numeric scalar parameter that is only used
#' for the \code{method="cpp"} algorithm.  When the DTQ summand drops below
#' code{thresh}, the algorithm stops summing, even if it has not summed over
#' all grid points.  The default value of this parameter is zero, indicating that
#' the full DTQ sum is evaluated.  Setting this parameter to a small positive
#' value such as \eqn{2.2 \times 10^{-16}} can result in a substantial speed up
#' for computations on large spatial grids, especially when \eqn{h} is also small.
#' @param method A string that indicates which DTQ algorithm to use.
#' There are two choices:
#' \describe{
#'   \item{\code{"cpp"}}{This DTQ method is implemented in C++.  No matrices are formed;
#'   the method is highly conservative in its usage of memory.  For sufficiently small \code{h} and \code{k},
#'   it is necessary to use this method.  This method also allows for approximate
#'   evaluation of the DTQ algorithm by setting a positive threshold parameter.}
#'   \item{\code{"sparse"}}{This DTQ method is implemented in R using sparse matrices from
#'   the Matrix package.  The method uses more memory than the \code{"cpp"} method, but may be
#'   faster for larger values of \code{h} and \code{k}.  This is the default method.}
#' }
#' @return The output consists of a list with two elements:
#' \describe{
#'   \item{\code{xvec}}{a numeric vector that contains the spatial grid}
#'   \item{\code{pdf}}{a numeric vector that contains the PDF evaluated at the grid points}
#' }
#' @seealso H. S. Bhat and R. W. M. A. Madushani,
#' "Density Tracking by Quadrature for Stochastic Differential Equations,"
#' arXiv:1610.09572 [stat.CO], \url{http://bit.ly/2fbNsp5}
#'
#' @examples
#' # Example 1:
#' # Define the drift function f(x) = -x and diffusion function g(x) = 1
#' # using C++ code:
#' require(Rcpp)
#' sourceCpp(code = '#include <Rcpp.h>
#' using namespace Rcpp;
#' double drift(double& x)
#' {
#'   return(-x);
#' }
#' double diff(double& x)
#' {
#'   return(1.0);
#' }
#' typedef double (*funcPtr)(double& x);
#' // [[Rcpp::export]]
#' XPtr<funcPtr> driftXPtr()
#' {
#'   return(XPtr<funcPtr>(new funcPtr(&drift)));
#' }
#' // [[Rcpp::export]]
#' XPtr<funcPtr> diffXPtr()
#' {
#'   return(XPtr<funcPtr>(new funcPtr(&diff)));
#' }')
#' # Solve for the PDF (at final time fT=1) of the SDE with drift f,
#' # diffusion g, and deterministic initial condition X(0) = 0.
#' # First we solve using the grid specified by k and bigm.
#' # Then we solve using the grid specified by a, b, and bigm.
#' # We then check that we get the same PDF either way.
#' k = 0.01
#' M = 250
#' test1 = rdtq(h=0.1,k,bigm=M,init=0,fT=1,
#'              drift=driftXPtr(),diffusion=diffXPtr(),method="cpp")
#' test2 = rdtq(h=0.1,a=-2.5,b=2.5,bigm=501,init=0,fT=1,
#'              drift=driftXPtr(),diffusion=diffXPtr(),method="cpp")
#' print(k*sum(abs(test1$pdf-test2$pdf)))
#'
#' # Example 2:
#' # We again use the drift function f(x) = -x and diffusion function g(x) = 1.
#' # This time, we use the method="sparse" version of DTQ.
#' # This requires us to define the drift and diffusion functions in R:
#' mydrift = function(x) { -x }
#' mydiff = function(x) { rep(1,length(x)) }
#' test = rdtq(h=0.1,k=0.01,bigm=250,init=0,fT=1,
#'             drift=mydrift,diffusion=mydiff,method="sparse")
#' plot(test$xvec,test$pdf,type='l')
rdtq <- function(h, k=NULL, bigm, a=NULL, b=NULL, init, fT,
                 drift=NULL, diffusion=NULL, thresh=0,
                 method="sparse") {
  if ((method != "cpp") && (method != "sparse"))
    stop("Invalid method.")
  if (is.null(drift) || is.null(diffusion))
    stop("Both the cpp and sparse methods require non-NULL drift and diffusion.")
  if ((! is.numeric(h)) || (h <= 0))
    stop("Step size h must be numeric and positive.")
  if ((! is.numeric(fT)) || (fT <= 0))
    stop("Final time T must be numeric and positive.")
  if (! is.numeric(bigm))
    stop("Detected non-numeric value of bigm.")
  if (!is.numeric(init))
    stop("Deteched non-numeric value of init.")
  bigmout = floor(bigm)
  if (bigmout <= 1)
    stop("floor(bigm) must be at least 2.")

  if (is.null(a) || is.null(b))
  {
    print("Calling dtq with grid specified via k and bigm.")
    if ((! is.numeric(k)) || (k <= 0))
      stop("Grid spacing k must be numeric and positive.")
    linit = length(init)
    if ((linit != 1) && (linit != (bigmout+1)))
      stop("Initial condition init must either be a scalar or a vector with the size (bigm+1), the same size as the spatial grid.")
    if (method=="cpp")
    {
      print("Using cpp method.")
      return(.Call('Rdtq_rdtq', PACKAGE = 'Rdtq', h, k, bigmout, init, fT, drift, diffusion, thresh))
    }
    if (method=="sparse")
    {
      print("Using sparse method.")
      return(.dtqsparse(h=h, k=k, bigm=bigmout, init=init, fT=fT, drift=drift, diffusion=diffusion))
    }
  }
  else
  {
    print("Calling dtq with grid specified via a, b, and bigm.")
    if (! is.numeric(a))
      stop("Left boundary a must be numeric.")
    if (! is.numeric(b))
      stop("Right boundary b must be numeric.")
    linit = length(init)
    if ((linit != 1) && (linit != bigmout))
      stop("Initial condition init must either be a scalar or a vector with the size bigm, the same size as the spatial grid.")
    if (method=="cpp")
    {
      print("Using cpp method.")
      return(.Call('Rdtq_rdtqgrid', PACKAGE = 'Rdtq', h, a, b, bigmout, init, fT, drift, diffusion, thresh))
    }
    if (method=="sparse")
    {
      print("Using sparse method.")
      return(.dtqsparse(h=h, a=a, b=b, bigm=bigmout, init=init, fT=fT, drift=drift, diffusion=diffusion))
    }
  }
}

#' Study DTQ Convergence
#'
#' \code{studydtqconv} facilitates the generation of convergence plots, i.e., plots
#' where one studies the error (in various norms) as a function of the time step
#' \eqn{h}; the error is computed as the difference between the exact PDF and the
#' approximate PDF computed via the DTQ method.
#'
#' @seealso H. S. Bhat and R. W. M. A. Madushani,
#' "Density Tracking by Quadrature for Stochastic Differential Equations,"
#' arXiv:1610.09572 [stat.CO], \url{http://bit.ly/2fbNsp5}
#'
#' @param method This must be a string, either "cpp" or "sparse", that indicates
#'  which algorithm to use.  See the parameter of the same
#'  name in the \code{rdtq} function.
#' @param drift if \code{method="cpp"}, then this should be a pointer to a
#'  drift function implemented in C++ using \code{Rcpp}.
#'  If \code{method="sparse"}, then this should be the name of an R function.
#'  For further details, see the description of the \code{drift} parameter for
#'  the \code{rdtq} function.
#' @param diffusion if \code{method="cpp"}, then this should be a pointer to a
#'  diffusion function implemented in C++ using \code{Rcpp}.
#'  If \code{method="sparse"}, then this should be the name of an R function.
#'  For further details, see the description of the \code{diffusion} parameter for
#'  the \code{rdtq} function.
#' @param exact an R function that accepts two arguments, \code{xvec} and
#'  \code{T}, and returns the exact probability density function \eqn{p(x,T)}
#'  for each \eqn{x} in \code{xvec}.
#' @param hseq a numeric vector of values of \eqn{h}, the time step, to use
#'  for the computation of the DTQ solution.  Note that \code{hseq} and
#'  \code{kseq} must have the same lengths.
#' @param kseq a numeric vector of values of \eqn{k}, the grid spacing, to use
#'  for the computation of the DTQ solution.  Note that \code{hseq} and
#'  \code{kseq} must have the same lengths.
#' @param Mseq a numeric vector of integer values of \eqn{M}.  For each
#'   corresponding value of \eqn{k}, the spatial grid will cover the domain
#'   \eqn{[-y_M, y_M]} where \eqn{y_M = Mk}.  This corresponds to the parameter
#'   \code{bigm} in the \code{rdtq} function.  Note that \code{kseq} and
#'   \code{Mseq} must have the same lengths.
#' @param init a scalar initial condition.
#' @param fT a positive numeric scalar giving the final time at which to compare
#'  the exact and DTQ solutions.
#' @param thresh This optional positive scalar is only used when \code{method="cpp"}.
#'  See the parameter of the same name in the
#'  \code{rdtq} function.
#' @return The function returns the errors between the DTQ and exact solutions
#'  indexed by the corresponding value of hseq.  The errors are returned in
#'  the \eqn{L^1} norm, \eqn{L^\infty} norm, and the Kolmogorov-Smirnov norm.
#'  The errors are returned in the form of a data frame.
#' @examples
#' # In this example, we will study the convergence of the DTQ method
#' # for the SDE with drift f(x) = x/2 + (1 + x^2)^(1/2) and
#' # diffusion g(x) = (1 + x^2)^(1/2).
#'
#' library(Rdtq)
#' library(Rcpp)
#'
#' # implement the drift and diffusion functions using C++
#' sourceCpp(code = '#include <Rcpp.h>
#'           using namespace Rcpp;
#'           double drift(double& x) { return(0.5*x + sqrt(x*x + 1.0)); }
#'           double diff(double& x) { return(sqrt(x*x + 1.0)); }
#'           typedef double (*funcPtr)(double& x);
#'           // [[Rcpp::export]]
#'           XPtr<funcPtr> driftXPtr() { return(XPtr<funcPtr>(new funcPtr(&drift))); }
#'           // [[Rcpp::export]]
#'           XPtr<funcPtr> diffXPtr() { return(XPtr<funcPtr>(new funcPtr(&diff))); }')
#'
#' # implement the drift and diffusion functions using R
#' mydrift = function(y)
#' {
#'   return(0.5*y + sqrt(y^2 + 1))
#' }
#' mydiff = function(y)
#' {
#'   return(sqrt(y^2 + 1))
#' }
#'
#' # implement the exact solution at time t, i.e.,
#' # the analytical formula for the pdf p(x,t)
#' exactsol = function(xvec,t)
#' {
#'   transx = asinh(xvec) - t
#'   prefac = (1 + xvec^2)^(-1/2)
#'   z = prefac*dnorm(x=transx)
#'   return(z)
#' }
#'
#' # define the sequence of parameters that will be used to study convergence
#' hseq = c(0.5,0.2,0.1,0.05,0.02,0.01)
#' kseq = hseq^(0.55)
#' Mseq = ceiling(5*(-log(hseq))/kseq)
#'
#' # we will use the method="sparse" code for the three largest values in hseq,
#' # and then switch to the method="cpp" code for the three smallest values
#' firstpart = c(1:3)
#' errpart1 = studydtqconv(method="sparse",drift=mydrift,diffusion=mydiff,exact=exactsol,
#'                         hseq[firstpart],kseq[firstpart],Mseq[firstpart],
#'                         init=0,fT=1)
#' errpart2 = studydtqconv(method="cpp",drift=driftXPtr(),diffusion=diffXPtr(),exact=exactsol,
#'                         hseq[-firstpart],kseq[-firstpart],Mseq[-firstpart],
#'                         init=0,fT=1,thresh=1e-16)
#'
#' # now we will put everything together into one data frame
#' mydat = rbind(errpart1,errpart2)
#'
#' # we plot the convergence diagram, on a log-log scale, using ggplot2
#' library(ggplot2)
#' library(scales)
#' myplot = ggplot(data=mydat, aes(x=x,y=y,group=norm,color=norm))
#' myplot = myplot + theme_bw() + theme(plot.background = element_rect(fill='white'))
#' myxticks = sort(10^(round(log(hseq)/log(10)*10)/10))
#' rawyticks = round(log(mydat$y)/log(10)*10)/10
#' rawyticks = round(seq(from=min(rawyticks),to=max(rawyticks),length.out=length(myxticks))*1)/1
#' myyticks = unique(10^rawyticks)
#' myplot = myplot + scale_x_log10(breaks = hseq)
#' myplot = myplot + theme(axis.text.x = element_text(angle=90,hjust=1))
#' myplot = myplot + scale_y_log10(breaks = myyticks, labels = trans_format("log10", math_format(10^.x)))
#' myplot = myplot + labs(x="h (temporal step size)", y="error")
#' myplot = myplot + geom_line() + geom_point()
#'
#' # save the plot to a pdf (portable document format) file
#' ggsave(filename="example.pdf", plot=myplot, width=5, height=4)
studydtqconv <- function(method, drift, diffusion, exact, hseq, kseq, Mseq, init, fT, thresh=0)
{
  if (!is.numeric(hseq))
    stop("hseq must be a numeric vector.")
  if (!is.numeric(kseq))
    stop("kseq must be a numeric vector.")
  if (!is.numeric(Mseq))
    stop("Mseq must be a numeric vector.")
  if (!is.numeric(init))
    stop("init must be numeric.")
  if (length(init) > 1)
    stop("init must be a scalar.")
  if ((length(hseq) != length(kseq)) || (length(hseq) != length(Mseq)))
    stop("hseq, kseq, and Mseq must have the same lengths.")

  nt = length(hseq)
  l1errors = numeric(length=nt)
  linferrors = numeric(length=nt)
  kserrors = numeric(length=nt)
  for (i in c(1:nt))
  {
    if (method=="cpp")
      dtqsol = rdtq(h=hseq[i],k=kseq[i],bigm=Mseq[i],
                    init=init,fT=fT,drift=drift,diffusion=diffusion,
                    method="cpp",thresh=thresh)
     if (method=="sparse")
       dtqsol = rdtq(h=hseq[i],k=kseq[i],bigm=Mseq[i],
                     init=init,fT=fT,drift=drift,diffusion=diffusion,
                     method="sparse")
    exactsol = exact(dtqsol$xvec,fT)
    l1errors[i] = kseq[i]*sum(abs(dtqsol$pdf - exactsol))
    linferrors[i] = max(abs(dtqsol$pdf - exactsol))
    kserrors[i] = max(abs(kseq[i]*(cumsum(dtqsol$pdf) - cumsum(exactsol))))
  }

  mydat = data.frame(x=hseq,y=l1errors,norm="L1")
  mydat = rbind(mydat,data.frame(x=hseq,y=linferrors,norm="Linf"))
  mydat = rbind(mydat,data.frame(x=hseq,y=kserrors,norm="K-S"))

  return(mydat)
}

