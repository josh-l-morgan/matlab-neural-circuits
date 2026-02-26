function [yhat] = timefieldfit2(beta,x)
% beta = [pOne tOne n pTwo tTwo]

pOne = beta(1); %amplitude of +lobe
tOne = beta(2); %peaktime of +lobe
n = beta(3);    

yhat = pOne * (x./tOne).^n .* exp(-n*(x./(tOne-1)));