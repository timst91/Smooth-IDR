# sanity check for local bandwidth smooth IDR (computational time ~10 min)

n=100


X=runif(n,0,10)
y=rgamma(n,shape=sqrt(X),scale=min(max(X,2),8))

x_test=6

intf=smooth_IDR_density_h_opt(y,X,nu=nu,x_test=x_test,y_test=10,h_init=2.5,
                              normalize = TRUE)
intf=intf$integral

int=integrate(function(x){
  smooth_IDR_density_h_opt(y,X,nu=nu,x_test=x_test,y_test=x,h_init=2.5,
                           normalize = TRUE)$density},
          lower=0,upper=20,subdivisions = 10000,rel.tol = 1e-4)

diff=smooth_IDR_CDF_h_opt(y,X,nu=nu,x_test=x_test,y_test=20,h_init=2.5)$cdf-
  smooth_IDR_CDF_h_opt(y,X,nu=nu,x_test=x_test,y_test=0,h_init =2.5)$cdf


abs(int-diff)



# sanity check for global bandwidth smooth IDR


n=100


X=runif(n,0,10)
y=rgamma(n,shape=sqrt(X),scale=min(max(X,2),8))

x_test=6


integrate(function(x){smooth_IDR_density(y,X,nu=nu,x_test=x_test,y_test=x,h=2.5)},
          lower=0,upper=20,subdivisions = 1000)

smooth_IDR_CDF(y,X,nu=nu,x_test=x_test,y_test=20,h=2.5)-
  smooth_IDR_CDF(y,X,nu=nu,x_test=x_test,y_test=0,h=2.5)


