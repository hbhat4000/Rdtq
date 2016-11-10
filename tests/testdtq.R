require(Rcpp)
sourceCpp(code = '#include <Rcpp.h>
  using namespace Rcpp;
  double drift(double& x)
  {
    return(-x);
  }
  double diff(double& x)
  {
    return(1.0);
  }
  typedef double (*funcPtr)(double& x);
  // [[Rcpp::export]]
  XPtr<funcPtr> driftXPtr()
  {
    return(XPtr<funcPtr>(new funcPtr(&drift)));
  }
  // [[Rcpp::export]]
  XPtr<funcPtr> diffXPtr()
  {
    return(XPtr<funcPtr>(new funcPtr(&diff)));
  }'
)

require(Rdtq)
k=0.01
M=250
test=rdtq(0.1,k,M,init=0,T=1,drift=driftXPtr(),diffusion=diffXPtr())
test2=rdtqgrid(0.1,-2.5,2.5,501,init=0,T=1,drift=driftXPtr(),diffusion=diffXPtr())

