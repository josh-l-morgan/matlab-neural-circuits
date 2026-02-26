function[coefVarRmsError rmsError] = rmsClust(testDistribution)

%%Coefficient of variation of the root mean square error of a distribution
%%is used to measure clustering in a distribution.   In this case
%%clustering is considered to be distance from each measure looking like
%%the average. Equivilant to rms(error)/mean



Error = (testDistribution-mean(testDistribution));

rmsError = sqrt(mean((Error.^2)));

coefVarRmsError = rmsError/mean(testDistribution);



