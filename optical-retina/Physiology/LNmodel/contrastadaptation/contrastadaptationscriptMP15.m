%% LOADING FILES
[fileName directoryName] = uigetfile('*.txt','Select .txt file');
cd(directoryName);
sp = nexread2006(fileName);
uiopen('load')
cd('C:\Program Files\MATLAB\R2006b\work')
clear D* fileName directoryName

analogStart = sp.data(find((diff(sp.data(:,2)) > 59)...
    & (diff(sp.data(:,2)) < 61),1, 'first'), 2);
analogStop = sp.data(find((diff(sp.data(:,1)) > 59)...
    & (diff(sp.data(:,1)) < 61),1, 'last')+1, 1);
sp.channels

%% ALIGN TIME AND DELETE CHANNELS
channelsToBeDeleted = [1 2 3 7 9 13 14 16];
stretchFactor = (Cadapt_S(end) - Cadapt_S(1)) / (analogStop - analogStart);
frameTrain = round((Cadapt_S - Cadapt_S(1))*1000);
timeDifference = diff(frameTrain,1,2);
repeatTime = mean(timeDifference(:));
spikes.channel = sp.channels;
spikes.channel(channelsToBeDeleted) = [];
[data] = deletecolumns (sp.data, channelsToBeDeleted);
spikes.channel

%% INTENSITY MATRIX, RESAMPLING TO 10MS (from 40ms hence the factor 4)
frameWindow = 100;
[m nChannels] = size(data);
rawIntensity = restoreintensity(Cadapt_int, .35, .64);        %this restores the input from output matrix (Cadapt_int) of intensities
intensity = rawIntensity - mean(rawIntensity(:));             %expressing stimulus entries as deviations from the mean
[M N] =size(intensity);
fourM = 4*M;
shortIndexVector = (1:M*N);
intensityVector = reshape(intensity, 1, M*N);
highIntensity = intensityVector;
highIntensity(mod(shortIndexVector, M) <= M/12 | mod(shortIndexVector, M) > M/3) =[];
lowIntensity = intensityVector;
lowIntensity(mod(shortIndexVector, M) <= 7*M/12 | mod(shortIndexVector, M) > 10*M/12) =[];
indexVector = (1:fourM*N);
intensityResampled = [intensityVector; intensityVector; intensityVector; intensityVector];
intensityMatrix = zeros(fourM*N, frameWindow);

for i=1:(fourM*N)
    if i >= frameWindow
        intensityMatrix(i,:) = intensityResampled(i - frameWindow + 1 : i);
    else
    end
end

%% COMPUTING STA, STE AND INPUT-OUTPUT RELATION
nValues = 50; %number of steps between upper and lower Limit of generator Signals in the input(X-axis) vector. The number of bins will be nValues-1.
g = waitbar(0, 'channels...');

for i=1:nChannels %this is what is done to each channel serially
    if (data(end,i) > 0)
        spikeTrainUnaligned = data(:,i);
    else
        lastPosition =  find (data(:,i) >= 0, 1, 'last');
        spikeTrainUnaligned = data(1:lastPosition,i);
    end

    spikesOfInterest = find ((analogStart < spikeTrainUnaligned)...
        & (spikeTrainUnaligned < analogStop));
    spikeTrain = round((spikeTrainUnaligned(spikesOfInterest(1)...
        : spikesOfInterest(end)) - analogStart) * stretchFactor * 1000);


    %the spikes during HIGH contrast stimulation
    highContrastSpikeTrain =  spikeTrain;
    highContrastSpikeTrain( mod(highContrastSpikeTrain, repeatTime) <= repeatTime/12 |...
        mod(highContrastSpikeTrain, repeatTime) > repeatTime/3 ) = [];
    highContrastIndex = zeros (1,length(highContrastSpikeTrain));

    for j=1:length(highContrastSpikeTrain)
        lastFrame = find((highContrastSpikeTrain(j) - frameTrain >= 0), 1, 'last');
        highContrastIndex(j) = ( (lastFrame-1)*4 ...
            + ceil((highContrastSpikeTrain(j) - frameTrain(lastFrame))/10) );
    end

    highContrastSTE = intensityMatrix(highContrastIndex,:);
    spikes.highContrastSTA(i,:) = mean(highContrastSTE);
    highContrastGeneratorPre = intensityMatrix * (spikes.highContrastSTA(i,:))';
    highContrastGeneratorPre(mod(indexVector, fourM) <= fourM/12 |...
        mod(indexVector, fourM) > fourM/3) =[];
    highContrastRatio = var(highIntensity)/var(highContrastGeneratorPre);
    spikes.highContrastSTAnormalized(i,:) = sqrt(highContrastRatio) * spikes.highContrastSTA(i,:);
    highContrastGeneratorSignal = sqrt(highContrastRatio) * highContrastGeneratorPre;
    highContrastGeneratorSorted = sort(highContrastGeneratorSignal);
    highContrastBoundaries = zeros(nValues, 2);
    for j=1:nValues
        highContrastBoundaries(j,:) = ...
            [highContrastGeneratorSorted(round((j-1)*length(highContrastGeneratorSignal)/nValues)+1), ...
            highContrastGeneratorSorted(round(j*length(highContrastGeneratorSignal)/nValues))];
    end
    %     highContrastInput = mean(highContrastBoundaries, 2);
    highContrastInput = zeros(nValues, 1);
    highContrastOutput = zeros(nValues, 1);
    highContrastFrameTrainVector = reshape(frameTrain, 1, M*N);
    highContrastFrameTime = mean(diff(highContrastFrameTrainVector));
    highContrastFrameTrainVector(mod(shortIndexVector, M) <= M/12 | mod(shortIndexVector, M) > M/3) =[];
    highContrastFrameTrain = [highContrastFrameTrainVector;...
        highContrastFrameTrainVector + highContrastFrameTime/4;...
        highContrastFrameTrainVector + highContrastFrameTime/2;...
        highContrastFrameTrainVector + 3*highContrastFrameTime/4];
    for j=1:nValues
        positionVector = find( (highContrastGeneratorSignal >= highContrastBoundaries(j,1))...
            & (highContrastGeneratorSignal < highContrastBoundaries(j,2)) );
        if positionVector(end) == numel(highContrastFrameTrain)
            positionVector(end)=[];
        else
        end
        highContrastInput(j) = median(highContrastGeneratorSignal(positionVector));
        momentaryOut = zeros(length(positionVector),1);
        for k=1:length(positionVector)
            momentaryOut(k) = length(find(highContrastSpikeTrain > highContrastFrameTrain(positionVector(k))...
                & highContrastSpikeTrain <= highContrastFrameTrain(positionVector(k)+1)));
        end
        highContrastOutput(j) = sum(momentaryOut)/length(positionVector);
        clear positionVector
    end
    muSeed = highContrastInput(find(highContrastOutput <=(max(highContrastOutput)...
        + min(highContrastOutput))/2, 1, 'last'));

    passInput.x = highContrastInput;
    passInput.maximumRate = max(highContrastOutput);
    passInput.minimumRate = min(highContrastOutput);

    highContrastFit0 = [10, -muSeed];
    %highContrastFit0 = [sensitivity, maintainedDrive]
    scalingOptions = statset ('MaxIter', 10000);
    [highContrastFit] = nlinfit(passInput, highContrastOutput,...
        @inputoutputfitMP, highContrastFit0, scalingOptions);
    spikes.highContrastInput(i,:) = highContrastInput;
    spikes.highContrastOutput(i,:) = highContrastOutput;
    spikes.highContrastFit(i,:) = [highContrastFit,...
        passInput.maximumRate, passInput.minimumRate];

    clear highContrast*

    %the spikes during LOW contrast stimulation
    lowContrastSpikeTrain =  spikeTrain;
    lowContrastSpikeTrain( mod(lowContrastSpikeTrain, repeatTime) <= 7*repeatTime/12 |...
        mod(lowContrastSpikeTrain, repeatTime) > 10*repeatTime/12 ) = [];
    lowContrastIndex = zeros (1,length(lowContrastSpikeTrain));

    for j=1:length(lowContrastSpikeTrain)
        lastFrame = find((lowContrastSpikeTrain(j) - frameTrain >= 0), 1, 'last');
        lowContrastIndex(j) = ( (lastFrame-1)*4 ...
            + ceil((lowContrastSpikeTrain(j) - frameTrain(lastFrame))/10) );
    end

    lowContrastSTE = intensityMatrix(lowContrastIndex,:);
    spikes.lowContrastSTA(i,:) = mean(lowContrastSTE);
    lowContrastGeneratorPre = intensityMatrix * (spikes.lowContrastSTA(i,:))';
    lowContrastGeneratorPre(mod(indexVector, fourM) <= 7*fourM/12 |...
        mod(indexVector, fourM) > 10*fourM/12) =[];
    lowContrastRatio = var(lowIntensity)/var(lowContrastGeneratorPre);
    spikes.lowContrastSTAnormalized(i,:) = sqrt(lowContrastRatio) * spikes.lowContrastSTA(i,:);
    lowContrastGeneratorSignal = sqrt(lowContrastRatio) * lowContrastGeneratorPre;
    lowContrastGeneratorSorted = sort(lowContrastGeneratorSignal);
    lowContrastBoundaries = zeros(nValues, 2);
    for j=1:nValues
        lowContrastBoundaries(j,:) = ...
            [lowContrastGeneratorSorted(round((j-1)*length(lowContrastGeneratorSignal)/nValues)+1), ...
            lowContrastGeneratorSorted(round(j*length(lowContrastGeneratorSignal)/nValues))];
    end
    lowContrastInput = zeros(nValues, 1);
    lowContrastOutput = zeros(nValues, 1);
    lowContrastFrameTrainVector = reshape(frameTrain, 1, M*N);
    lowContrastFrameTime = mean(diff(lowContrastFrameTrainVector));
    lowContrastFrameTrainVector(mod(shortIndexVector, M) <= 7*M/12 | mod(shortIndexVector, M) > 10*M/12) =[];
    lowContrastFrameTrain = [lowContrastFrameTrainVector;...
        lowContrastFrameTrainVector + lowContrastFrameTime/4;...
        lowContrastFrameTrainVector + lowContrastFrameTime/2;...
        lowContrastFrameTrainVector + 3*lowContrastFrameTime/4];
    for j=1:nValues
        positionVector = find( (lowContrastGeneratorSignal >= lowContrastBoundaries(j,1))...
            & (lowContrastGeneratorSignal < lowContrastBoundaries(j,2)) );
        if positionVector(end) == numel(lowContrastFrameTrain)
            positionVector(end)=[];
        else
        end
        lowContrastInput(j) = median(lowContrastGeneratorSignal(positionVector));
        momentaryOut = zeros(length(positionVector),1);
        for k=1:length(positionVector)
            momentaryOut(k) = length(find(lowContrastSpikeTrain > lowContrastFrameTrain(positionVector(k))...
                & lowContrastSpikeTrain <= lowContrastFrameTrain(positionVector(k)+1)));
        end
        lowContrastOutput(j) = sum(momentaryOut)/length(positionVector);
        clear positionVector
    end
    muSeed = lowContrastInput(find(lowContrastOutput <=(max(lowContrastOutput)...
        + min(lowContrastOutput))/2, 1, 'last'));
    passInput.x = lowContrastInput;
    lowContrastFit0 = [10, -muSeed];
    %lowContrastFit0 =[sensitivity, maintainedDrive]
    scalingOptions = statset ('MaxIter', 10000);
    [lowContrastFit] = nlinfit(passInput, lowContrastOutput,...
        @inputoutputfitMP, lowContrastFit0, scalingOptions);
    spikes.lowContrastInput(i,:) = lowContrastInput;
    spikes.lowContrastOutput(i,:) = lowContrastOutput;
    spikes.lowContrastFit(i,:) = [lowContrastFit,...
        passInput.maximumRate, passInput.minimumRate];

    clear lowContrast*

    waitbar(i/nChannels, g)
end
close(g)