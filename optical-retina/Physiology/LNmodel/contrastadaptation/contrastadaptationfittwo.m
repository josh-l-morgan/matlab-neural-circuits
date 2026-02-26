function yhat = contrastadaptationfittwo (beta, passInput)
%beta = [sensitivity maintainedDrive]

p = normcdf (  beta(1)*(passInput.data + beta(2)), passInput.mu, passInput.sigma);
yhat = p * passInput.maximumFiring;