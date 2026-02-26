function rateTrain = spiketorate (spikeTrain, frameTrain, binSize)
% transfers spike timestamps into rates

nFrames = length(frameTrain);
rateTrain = zeros(nFrames-1,1);
for i=1:nFrames-1
    
        rateTrain(i) =length(spikeTrain(spikeTrain > frameTrain(i)...
            & spikeTrain <= frameTrain(i+1))) / binSize;

end

        