function logicalTrain = spiketorate (spikeTrain, frameTrain)
% transfers spike timestamps into rates
nFrames = length(frameTrain);
logicalTrain = zeros(nFrames,1);
for i=1:nFrames-1
    
        logicalTrain(i) =length(spikeTrain(spikeTrain > frameTrain(i)...
            & spikeTrain <= frameTrain(i+1)));

end

        