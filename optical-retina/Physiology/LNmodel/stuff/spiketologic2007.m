function logicalTrain = spiketologic2007 (spikeTrain, frameTrain)
nFrames = length(frameTrain);
logicalTrain = zeros(nFrames,1);
for i=1:nFrames-1
        logicalTrain(i) = any(spikeTrain(spikeTrain > frameTrain(i)...
            & spikeTrain <= frameTrain(i+1)));
end

        