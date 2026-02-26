function[errResults] = conError(con,predCon,useMat);

%%calculate  the error of a predicted connectivity graph


testCon.con = con;
testCon.pred = predCon;
testCon.use = useMat;

%}
clear errResults
%% Get matricies
con = testCon.con;

if ~exist('useMat','var')
    %%use all entries by default
    useMat = (con*0+1)>0;
else
    useMat = testCon.use>0;
end

if ~exist('predCon','var')
    %%Predict the mean for every entry if no prediction is given
    predCon = con*0;
    predCon(useMat) = mean(con(useMat));
end

useCon = con(useMat);
usePred = predCon(useMat);


%% Find error

difCon = useCon - usePred;
varCon = sum(difCon.^2)/(length(difCon)-1);

other.synErrRange = [sum(floor(abs(difCon))) ...
    round(sum(abs(difCon))) ...
    sum(ceil(abs(difCon)))];

%%Error rates

synErr = round(useCon) - round(usePred);
synErrorRate = sum(abs(synErr))/(sum(round(useCon)) + sum(round(usePred)))* 100;

misses = round(usePred)-round(useCon);
misses(misses<0) = 0;
missNum = sum(misses);
hitRate = (sum(round(usePred))-missNum)/sum(round(usePred))*100;

errResults.hitRate = hitRate; %percent of predicted synapses that landed in the right place.
errResults.synErrorPC = synErrorRate;


%% Binary accuracy
binCon = round(useCon)>0;
binPred = round(usePred)>0;

falsePositives = sum((binPred-binCon)>0);
falseNegatives = sum((binCon - binPred)>0);
errorRate = (falsePositives + falseNegatives ) / (sum(binCon) + sum(binPred)) * 100;

other.binary.falsePositives = falsePositives;
other.binary.falseNegatives = falseNegatives;
other.binary.errorRate = errorRate;
errResults.binaryErrorPC = errorRate;



%% Errors (unused)
errCon = useCon - mean(useCon);
errPred = usePred - mean(usePred);
sum(errCon .* errPred)/(length(errCon)-1);


%% Calculate covariance of syn and prediction relative to 0
var0Con = sum(useCon.^2)/(length(useCon)-1);
var0Pred = sum(usePred.^2)/(length(usePred)-1);
synCov0 = (sum(useCon .* usePred)/(length(useCon)-1)); %covariance relative to 0 (not mean)
synCoef0 = synCov0 / sqrt(var0Con * var0Pred);

rootCov0 = sqrt(synCov0);

errResults.synCoef0 = synCoef0;
other.synCov0 = synCov0;


%% Standard coefficient of variation
synCovarianceMatrix = cov([useCon usePred]);
[R,P,RLO,RUP]=corrcoef([useCon usePred]);

other.covariance = synCovarianceMatrix;

errResults.corrcoef = R(1,end);
errResults.other.R = R;
errResults.other.P = P;
errResults.other.RLO = RLO;
errResults.other.RUP = RUP;



%% Kruskal–Wallis one-way analysis of variance

errResults.other = other;



%% Record basic Stats


errResults.sampleN = length(useCon);
errResults.realSynNumber = sum(round(useCon));
errResults.predSynNumber = sum(round(usePred));
errResults.realBinNumber = sum(round(useCon)>0);
errResults.predBinNumber = sum(round(usePred)>0);
errResults.realCVRMSE = cvrmse(useCon);
errResults.predCVRMSE = cvrmse(usePred);

errResults.other.testCon = testCon;



