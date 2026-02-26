function yhat = inputoutputfitMP (beta, passInput)
%beta = [sensitivity maintainedDrive]
p = normcdf ( (beta(1) * (passInput.x + beta(2))), 0, 1);
yhat = p * passInput.maximumRate + passInput.minimumRate;