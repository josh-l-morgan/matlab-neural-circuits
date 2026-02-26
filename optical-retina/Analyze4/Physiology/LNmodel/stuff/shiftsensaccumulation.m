[m nExperiments] = size(ONBSFAST);
for k=1:nExperiments
    if exist ('drive','var')
        nextEntry = length(drive);
        nCells = length(ONBSFAST(k).spikes.channel);
        for i=1:nCells
            sens(nextEntry+i) = ONBSFAST(k).spikes.highContrastFit(i,1); %#ok<AGROW>
            drive(nextEntry+i) = ONBSFAST(k).spikes.highContrastFit(i,2); %#ok<AGROW>
        end
    else
        nCells = length(ONBSFAST(k).spikes.channel);
        for i=1:nCells
            sens(i) = ONBSFAST(k).spikes.highContrastFit(i,1); %#ok<AGROW>
            drive(i) = ONBSFAST(k).spikes.highContrastFit(i,2); %#ok<AGROW>
        end
    end
end

        
