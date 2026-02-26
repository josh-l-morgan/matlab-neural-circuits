function [pABgivenS pAxpBgivenS] ...
    = independence(spikes, indA, indB, nValues)

nFrames = length(spikes(indA).logicTrain);
logicAB = spikes(indA).logicTrain .* spikes(indB).logicTrain;
pAxpB = spikes(indA).probabilityTrain .* spikes(indB).probabilityTrain;
pAxpBs = sort(pAxpB);

pAxpBLimit = zeros(nValues, 2);
for j=1:nValues
    pAxpBLimit(j,:) = ...
        [pAxpBs(round((j-1)*nFrames/nValues)+1),...
        pAxpBs(round(j*nFrames/nValues))];
end
% pAxpBgivenS = mean(pAxpBLimit, 2);
pABgivenS = zeros(nValues, 1);

for j=1:nValues
    positionVector = ...
        find( (pAxpB >= pAxpBLimit(j,1)) & (pAxpB < pAxpBLimit(j,2)) );
    if positionVector(end) == nFrames
        positionVector(end)=[];
    else
    end
    pAxpBgivenS(j) = median(pAxpB(positionVector));
    pABgivenS(j) = mean(logicAB(positionVector));
    clear positionVector
end