function [yhat] = timefieldfit(beta,x)
% beta = [pOne tOne n pTwo tTwo]

pOne = beta(1); %amplitude of +lobe
tOne = beta(2); %peaktime of +lobe
n = beta(3);    
pTwo = beta(4); %amplitude of -lobe
tTwo = beta(5); %peaktime of -lobe

yhat = pOne * (x./tOne).^n .* exp(-n*(x./(tOne-1))) -...
    pTwo * (x./tTwo).^n .* exp(-n*(x./(tTwo-1)));