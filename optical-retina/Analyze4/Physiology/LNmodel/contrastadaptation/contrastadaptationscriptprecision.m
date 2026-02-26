%% LOADING FILES
[fileName directoryName] = uigetfile('*.txt','Select .txt file');
cd(directoryName);
sp = nexread2006(fileName);
uiopen('load')
cd('F:\Program Files\MATLAB\R2006a\work')
clear D* fileName directoryName

analogStart = sp.data(find((diff(sp.data(:,2)) > 59) & (diff(sp.data(:,2)) < 61),1, 'first'), 2);
analogStop = sp.data(find((diff(sp.data(:,1)) > 59) & (diff(sp.data(:,1)) < 61),1, 'last')+1, 1);
sp.channels

%% ALIGN TIME AND DELETE CHANNELS
channelsToBeDeleted = [1 2 5 13]; %user-defined parameter
stretchFactor = (Cadapt_S(end) - Cadapt_S(1)) / (analogStop - analogStart);
frameTrain = round((Cadapt_S - Cadapt_S(1))*1000);
timeDifference = diff(frameTrain,1,2);
repeatTime = mean(timeDifference(:));
channel = sp.channels;
channel(channelsToBeDeleted) = [];
[data] = deletecolumns (sp.data, channelsToBeDeleted);
channel

%% INTENSITY MATRIX, RESAMPLING TO 10MS (from 40ms hence the factor 4)
frameWindow = 50;
[m nChannels] = size(data);
rawIntensity = restoreintensity(Cadapt_int, .35, .64);        %this restores the input from output matrix (Cadapt_int) of intensities
intensity = rawIntensity - mean(rawIntensity(:));             %expressing stimulus entries as deviations from the mean
[M N] =size(intensity);
fourM = 4*M;
shortIndexVector = (1:M*N);
indexVector = (1:fourM*N);
intensityVector = reshape(intensity, 1, M*N);
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

    spikesOfInterest = find ((analogStart < spikeTrainUnaligned) & (spikeTrainUnaligned < analogStop));
    spikeTrain = round((spikeTrainUnaligned(spikesOfInterest(1) : spikesOfInterest(end)) - analogStart) * stretchFactor * 1000);


    %the spikes during high contrast stimulation
    highContrastSpikeTrain =  spikeTrain;
    highContrastSpikeTrain( mod(highContrastSpikeTrain, repeatTime) <= repeatTime/6 |...
        mod(highContrastSpikeTrain, repeatTime) > repeatTime/2 ) = [];
    highContrastIndex = zeros (1,length(highContrastSpikeTrain));

    h = waitbar(0, 'high contrast spikes...');
    for j=1:length(highContrastSpikeTrain)
        lastFrame = find((highContrastSpikeTrain(j) - frameTrain >= 0), 1, 'last');
        highContrastIndex(j) = ( (lastFrame-1)*4 ...
            + ceil((highContrastSpikeTrain(j) - frameTrain(lastFrame))/10) );
        waitbar(j/length(highContrastSpikeTrain), h)
    end
    close(h)

    highContrastSTE = intensityMatrix(highContrastIndex,:);
    spikes(i).highContrastSTE = highContrastSTE;
    spikes(i).highContrastSTA = mean(highContrastSTE);
    highContrastSTAtranspose = (mean(highContrastSTE))';
    highContrastGeneratorSignal = intensityMatrix * highContrastSTAtranspose;
    highContrastGeneratorSignal(mod(indexVector, fourM) <= fourM/6 |...
        mod(indexVector, fourM) > fourM/2) =[];
    highContrastGeneratorSorted = sort(highContrastGeneratorSignal);
    highContrastBoundaries = zeros(nValues, 2);
    for j=1:nValues
        highContrastBoundaries(j,:) = ...
            [highContrastGeneratorSorted(round((j-1)*length(highContrastGeneratorSignal)/nValues)+1), ...
            highContrastGeneratorSorted(round(j*length(highContrastGeneratorSignal)/nValues))];
    end
    highContrastInput = mean(highContrastBoundaries, 2);
    highContrastOutput = zeros(length(highContrastInput), 1);
    highContrastFrameTrainVector = reshape(frameTrain, 1, M*N);
    highContrastFrameTime = mean(diff(highContrastFrameTrainVector));
    highContrastFrameTrainVector(mod(shortIndexVector, M) <= M/6 | mod(shortIndexVector, M) > M/2) =[];
    highContrastFrameTrain = [highContrastFrameTrainVector;...
        highContrastFrameTrainVector + highContrastFrameTime/4;...
        highContrastFrameTrainVector + highContrastFrameTime/2;...
        highContrastFrameTrainVector + 3*highContrastFrameTime/4];
    h = waitbar(0, 'high contrast input-output...');
    for j=1:nValues
        positionVector = find( (highContrastGeneratorSignal >= highContrastBoundaries(j,1))...
            & (highContrastGeneratorSignal < highContrastBoundaries(j,2)) );
        if positionVector(end) == numel(highContrastFrameTrain)
            positionVector(end)=[];
        else
        end
        momentaryOut = zeros(length(positionVector),1);
        for k=1:length(positionVector)
            momentaryOut(k) = length(find(highContrastSpikeTrain > highContrastFrameTrain(positionVector(k))...
                & highContrastSpikeTrain <= highContrastFrameTrain(positionVector(k)+1)));
        end
        highContrastOutput(j) = sum(momentaryOut)/length(positionVector);
        waitbar(j/(nValues), h)
        clear positionVector
    end
    spikes(i).highContrastInput = highContrastInput;
    spikes(i).highContrastOutput = highContrastOutput;
    close(h)
    clear high*

    %the spikes during low contrast stimulation
    lowContrastSpikeTrain =  spikeTrain;
    lowContrastSpikeTrain( mod(lowContrastSpikeTrain, repeatTime) <= 4*repeatTime/6) = [];
    lowContrastIndex = zeros (1,length(lowContrastSpikeTrain));

    h = waitbar(0, 'low contrast spikes...');
    for j=1:length(lowContrastSpikeTrain)
        lastFrame = find((lowContrastSpikeTrain(j) - frameTrain >= 0), 1, 'last');
        lowContrastIndex(j) = ( (lastFrame-1)*4 ...
            + ceil((lowContrastSpikeTrain(j) - frameTrain(lastFrame))/10) );
        waitbar(j/length(lowContrastSpikeTrain), h)
    end
    close(h)

    lowContrastSTE = intensityMatrix(lowContrastIndex,:);
    spikes(i).lowContrastSTE = lowContrastSTE;
    spikes(i).lowContrastSTA = mean(lowContrastSTE);
    lowContrastSTAtranspose = (mean(lowContrastSTE))';
    lowContrastGeneratorSignal = intensityMatrix * lowContrastSTAtranspose;
    lowContrastGeneratorSignal(mod(indexVector, fourM) < 4*fourM/6) =[];
    lowContrastGeneratorSorted = sort(lowContrastGeneratorSignal);
    lowContrastBoundaries = zeros(nValues, 2);
    for j=1:nValues
        lowContrastBoundaries(j,:) = ...
            [lowContrastGeneratorSorted(round((j-1)*length(lowContrastGeneratorSignal)/nValues)+1), ...
            lowContrastGeneratorSorted(round(j*length(lowContrastGeneratorSignal)/nValues))];
    end
    lowContrastInput = mean(lowContrastBoundaries, 2);
    lowContrastOutput = zeros(length(lowContrastInput), 1);
    lowContrastFrameTrainVector = reshape(frameTrain, 1, M*N);
    lowContrastFrameTime = mean(diff(lowContrastFrameTrainVector));
    lowContrastFrameTrainVector(mod(shortIndexVector, M) < 4*M/6) =[];
    lowContrastFrameTrain = [lowContrastFrameTrainVector;...
        lowContrastFrameTrainVector + lowContrastFrameTime/4;...
        lowContrastFrameTrainVector + lowContrastFrameTime/2;...
        lowContrastFrameTrainVector + 3*lowContrastFrameTime/4];

    h = waitbar(0, 'low contrast input-output...');
    for j=1:nValues
        positionVector = find( (lowContrastGeneratorSignal >= lowContrastBoundaries(j,1))...
            & (lowContrastGeneratorSignal < lowContrastBoundaries(j,2)) );
        if positionVector(end) == numel(lowContrastFrameTrain)
            positionVector(end)=[];
        else
        end
        momentaryOut = zeros(length(positionVector),1);
        for k=1:length(positionVector)
            momentaryOut(k) = length(find(lowContrastSpikeTrain > lowContrastFrameTrain(positionVector(k))...
                & lowContrastSpikeTrain <= lowContrastFrameTrain(positionVector(k)+1)));
        end
        lowContrastOutput(j) = sum(momentaryOut)/length(positionVector);
        waitbar(j/(nValues), h)
        clear positionVector
    end
    spikes(i).lowContrastInput = lowContrastInput;
    spikes(i).lowContrastOutput = lowContrastOutput;
    spikes(i).channel = channel(i);
    close(h)
    clear low*
    waitbar(i/nChannels, g)
end
close(g)