function yhat = contrastadaptationfitone (beta, x)
%beta = [sigma, mu, maximumFiring]

p = normcdf (x, beta(2), beta(1));
yhat =  p*beta(3);