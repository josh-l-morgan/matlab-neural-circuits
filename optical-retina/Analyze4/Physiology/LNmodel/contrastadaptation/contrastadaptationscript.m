%% LOADING FILES
[fileName directoryName] = uigetfile('*.txt','Select .txt file');
cd(directoryName);
sp = nexread2006(fileName);
uiopen('load')
cd('F:\Program Files\MATLAB\R2006a\work')
clear D* fileName directoryName

% analogStart = sp.data(find((diff(sp.data(:,2)) > 59) & (diff(sp.data(:,2)) < 61),1, 'first'), 2);
% analogStop = sp.data(find((diff(sp.data(:,1)) > 59) & (diff(sp.data(:,1)) < 61),1, 'last')+1, 1);
sp.channels

%% ALIGN TIME AND DELETE CHANNELS
% channelsToBeDeleted = [1 2 8 11 18 20 22 23 24 28 29 30 32 35]; %user-defined parameter
channelsToBeDeleted = [1 3];
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
intensityVector = reshape(intensity, 1, M*N);
highIntensityEarly = intensityVector;
highIntensityEarly(mod(shortIndexVector, M) >= M/12) =[];
highIntensityLate = intensityVector;
highIntensityLate(mod(shortIndexVector, M) <= M/3 | mod(shortIndexVector, M) > M/2) =[];
lowIntensityEarly = intensityVector;
lowIntensityEarly(mod(shortIndexVector, M) <= M/2 | mod(shortIndexVector, M) > 7*M/12) =[];
lowIntensityLate = intensityVector;
lowIntensityLate(mod(shortIndexVector, M) <= 2*M/3) =[];
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

    spikesOfInterest = find ((analogStart < spikeTrainUnaligned) & (spikeTrainUnaligned < analogStop));
    spikeTrain = round((spikeTrainUnaligned(spikesOfInterest(1) : spikesOfInterest(end)) - analogStart) * stretchFactor * 1000);


    %the spikes EARLY during HIGH contrast stimulation
    highContrastEarlySpikeTrain =  spikeTrain;
    highContrastEarlySpikeTrain( mod(highContrastEarlySpikeTrain, repeatTime) >= repeatTime/12 ) = [];
    highContrastEarlyIndex = zeros (1,length(highContrastEarlySpikeTrain));

    h = waitbar(0, 'early high contrast spikes...');
    for j=1:length(highContrastEarlySpikeTrain)
        lastFrame = find((highContrastEarlySpikeTrain(j) - frameTrain >= 0), 1, 'last');
        highContrastEarlyIndex(j) = ( (lastFrame-1)*4 ...
            + ceil((highContrastEarlySpikeTrain(j) - frameTrain(lastFrame))/10) );
        waitbar(j/length(highContrastEarlySpikeTrain), h)
    end
    close(h)

    spikes(i).highContrastEarlySTE = intensityMatrix(highContrastEarlyIndex,:);
    spikes(i).highContrastEarlySTA = mean(spikes(i).highContrastEarlySTE);
    highContrastEarlyGeneratorPre = intensityMatrix * (spikes(i).highContrastEarlySTA)';
    highContrastEarlyGeneratorPre(mod(indexVector, fourM) >= fourM/12) =[];
    highContrastEarlyRatio = var(highIntensityEarly)/var(highContrastEarlyGeneratorPre);
    spikes(i).highContrastEarlySTAnormalized = sqrt(highContrastEarlyRatio) * spikes(i).highContrastEarlySTA;
    highContrastEarlyGeneratorSignal = sqrt(highContrastEarlyRatio) * highContrastEarlyGeneratorPre;
    highContrastEarlyGeneratorSorted = sort(highContrastEarlyGeneratorSignal);
    highContrastEarlyBoundaries = zeros(nValues, 2);
    for j=1:nValues
        highContrastEarlyBoundaries(j,:) = ...
            [highContrastEarlyGeneratorSorted(round((j-1)*length(highContrastEarlyGeneratorSignal)/nValues)+1), ...
            highContrastEarlyGeneratorSorted(round(j*length(highContrastEarlyGeneratorSignal)/nValues))];
    end
    highContrastEarlyInput = mean(highContrastEarlyBoundaries, 2);
    highContrastEarlyOutput = zeros(length(highContrastEarlyInput), 1);
    highContrastEarlyFrameTrainVector = reshape(frameTrain, 1, M*N);
    highContrastEarlyFrameTime = mean(diff(highContrastEarlyFrameTrainVector));
    highContrastEarlyFrameTrainVector(mod(shortIndexVector, M) >= M/12) =[];
    highContrastEarlyFrameTrain = [highContrastEarlyFrameTrainVector;...
        highContrastEarlyFrameTrainVector + highContrastEarlyFrameTime/4;...
        highContrastEarlyFrameTrainVector + highContrastEarlyFrameTime/2;...
        highContrastEarlyFrameTrainVector + 3*highContrastEarlyFrameTime/4];
    h = waitbar(0, 'early high contrast input-output...');
    for j=1:nValues
        positionVector = find( (highContrastEarlyGeneratorSignal >= highContrastEarlyBoundaries(j,1))...
            & (highContrastEarlyGeneratorSignal < highContrastEarlyBoundaries(j,2)) );
        if positionVector(end) == numel(highContrastEarlyFrameTrain)
            positionVector(end)=[];
        else
        end
        momentaryOut = zeros(length(positionVector),1);
        for k=1:length(positionVector)
            momentaryOut(k) = length(find(highContrastEarlySpikeTrain > highContrastEarlyFrameTrain(positionVector(k))...
                & highContrastEarlySpikeTrain <= highContrastEarlyFrameTrain(positionVector(k)+1)));
        end
        highContrastEarlyOutput(j) = sum(momentaryOut)/length(positionVector);
        waitbar(j/(nValues), h)
        clear positionVector
    end
    spikes(i).highContrastEarlyInput = highContrastEarlyInput;
    spikes(i).highContrastEarlyOutput = highContrastEarlyOutput;
    close(h)
    
    
    %the spikes LATE during HIGH contrast stimulation
    highContrastLateSpikeTrain =  spikeTrain;
    highContrastLateSpikeTrain( mod(highContrastLateSpikeTrain, repeatTime) <= repeatTime/3 |...
        mod(highContrastLateSpikeTrain, repeatTime) > repeatTime/2 ) = [];
    highContrastLateIndex = zeros (1,length(highContrastLateSpikeTrain));

    h = waitbar(0, 'late high contrast spikes...');
    for j=1:length(highContrastLateSpikeTrain)
        lastFrame = find((highContrastLateSpikeTrain(j) - frameTrain >= 0), 1, 'last');
        highContrastLateIndex(j) = ( (lastFrame-1)*4 ...
            + ceil((highContrastLateSpikeTrain(j) - frameTrain(lastFrame))/10) );
        waitbar(j/length(highContrastLateSpikeTrain), h)
    end
    close(h)

    spikes(i).highContrastLateSTE = intensityMatrix(highContrastLateIndex,:);
    spikes(i).highContrastLateSTA = mean(spikes(i).highContrastLateSTE);
    highContrastLateGeneratorPre = intensityMatrix * (spikes(i).highContrastLateSTA)';
    highContrastLateGeneratorPre(mod(indexVector, fourM) <= fourM/3 |...
        mod(indexVector, fourM) > fourM/2) =[];
    highContrastLateRatio = var(highIntensityLate)/var(highContrastLateGeneratorPre);
    spikes(i).highContrastLateSTAnormalized = sqrt(highContrastLateRatio) * spikes(i).highContrastLateSTA;
    highContrastLateGeneratorSignal = sqrt(highContrastLateRatio) * highContrastLateGeneratorPre;
    highContrastLateGeneratorSorted = sort(highContrastLateGeneratorSignal);
    highContrastLateBoundaries = zeros(nValues, 2);
    for j=1:nValues
        highContrastLateBoundaries(j,:) = ...
            [highContrastLateGeneratorSorted(round((j-1)*length(highContrastLateGeneratorSignal)/nValues)+1), ...
            highContrastLateGeneratorSorted(round(j*length(highContrastLateGeneratorSignal)/nValues))];
    end
    highContrastLateInput = mean(highContrastLateBoundaries, 2);
    highContrastLateOutput = zeros(length(highContrastLateInput), 1);
    highContrastLateFrameTrainVector = reshape(frameTrain, 1, M*N);
    highContrastLateFrameTime = mean(diff(highContrastLateFrameTrainVector));
    highContrastLateFrameTrainVector(mod(shortIndexVector, M) <= M/3 | mod(shortIndexVector, M) > M/2) =[];
    highContrastLateFrameTrain = [highContrastLateFrameTrainVector;...
        highContrastLateFrameTrainVector + highContrastLateFrameTime/4;...
        highContrastLateFrameTrainVector + highContrastLateFrameTime/2;...
        highContrastLateFrameTrainVector + 3*highContrastLateFrameTime/4];
    h = waitbar(0, 'late high contrast input-output...');
    for j=1:nValues
        positionVector = find( (highContrastLateGeneratorSignal >= highContrastLateBoundaries(j,1))...
            & (highContrastLateGeneratorSignal < highContrastLateBoundaries(j,2)) );
        if positionVector(end) == numel(highContrastLateFrameTrain)
            positionVector(end)=[];
        else
        end
        momentaryOut = zeros(length(positionVector),1);
        for k=1:length(positionVector)
            momentaryOut(k) = length(find(highContrastLateSpikeTrain > highContrastLateFrameTrain(positionVector(k))...
                & highContrastLateSpikeTrain <= highContrastLateFrameTrain(positionVector(k)+1)));
        end
        highContrastLateOutput(j) = sum(momentaryOut)/length(positionVector);
        waitbar(j/(nValues), h)
        clear positionVector
    end
    spikes(i).highContrastLateInput = highContrastLateInput;
    spikes(i).highContrastLateOutput = highContrastLateOutput;
    close(h)
    clear highContrast*
    
    %the spikes EARLY during LOW contrast stimulation
    lowContrastEarlySpikeTrain =  spikeTrain;
    lowContrastEarlySpikeTrain( mod(lowContrastEarlySpikeTrain, repeatTime) <= repeatTime/2 |...
        mod(lowContrastEarlySpikeTrain, repeatTime) > 7*repeatTime/12 ) = [];
    lowContrastEarlyIndex = zeros (1,length(lowContrastEarlySpikeTrain));

    h = waitbar(0, 'early low contrast spikes...');
    for j=1:length(lowContrastEarlySpikeTrain)
        lastFrame = find((lowContrastEarlySpikeTrain(j) - frameTrain >= 0), 1, 'last');
        lowContrastEarlyIndex(j) = ( (lastFrame-1)*4 ...
            + ceil((lowContrastEarlySpikeTrain(j) - frameTrain(lastFrame))/10) );
        waitbar(j/length(lowContrastEarlySpikeTrain), h)
    end
    close(h)

    spikes(i).lowContrastEarlySTE = intensityMatrix(lowContrastEarlyIndex,:);
    spikes(i).lowContrastEarlySTA = mean(spikes(i).lowContrastEarlySTE);
    lowContrastEarlyGeneratorPre = intensityMatrix * (spikes(i).lowContrastEarlySTA)';
    lowContrastEarlyGeneratorPre(mod(indexVector, fourM) <= fourM/2 |...
        mod(indexVector, fourM) > 7*fourM/12) =[];
    lowContrastEarlyRatio = var(lowIntensityEarly)/var(lowContrastEarlyGeneratorPre);
    spikes(i).lowContrastEarlySTAnormalized = sqrt(lowContrastEarlyRatio) * spikes(i).lowContrastEarlySTA;
    lowContrastEarlyGeneratorSignal = sqrt(lowContrastEarlyRatio) * lowContrastEarlyGeneratorPre;
    lowContrastEarlyGeneratorSorted = sort(lowContrastEarlyGeneratorSignal);
    lowContrastEarlyBoundaries = zeros(nValues, 2);
    for j=1:nValues
        lowContrastEarlyBoundaries(j,:) = ...
            [lowContrastEarlyGeneratorSorted(round((j-1)*length(lowContrastEarlyGeneratorSignal)/nValues)+1), ...
            lowContrastEarlyGeneratorSorted(round(j*length(lowContrastEarlyGeneratorSignal)/nValues))];
    end
    lowContrastEarlyInput = mean(lowContrastEarlyBoundaries, 2);
    lowContrastEarlyOutput = zeros(length(lowContrastEarlyInput), 1);
    lowContrastEarlyFrameTrainVector = reshape(frameTrain, 1, M*N);
    lowContrastEarlyFrameTime = mean(diff(lowContrastEarlyFrameTrainVector));
    lowContrastEarlyFrameTrainVector(mod(shortIndexVector, M) <= M/2 | mod(shortIndexVector, M) > 7*M/12) =[];
    lowContrastEarlyFrameTrain = [lowContrastEarlyFrameTrainVector;...
        lowContrastEarlyFrameTrainVector + lowContrastEarlyFrameTime/4;...
        lowContrastEarlyFrameTrainVector + lowContrastEarlyFrameTime/2;...
        lowContrastEarlyFrameTrainVector + 3*lowContrastEarlyFrameTime/4];
    h = waitbar(0, 'early low contrast input-output...');
    for j=1:nValues
        positionVector = find( (lowContrastEarlyGeneratorSignal >= lowContrastEarlyBoundaries(j,1))...
            & (lowContrastEarlyGeneratorSignal < lowContrastEarlyBoundaries(j,2)) );
        if positionVector(end) == numel(lowContrastEarlyFrameTrain)
            positionVector(end)=[];
        else
        end
        momentaryOut = zeros(length(positionVector),1);
        for k=1:length(positionVector)
            momentaryOut(k) = length(find(lowContrastEarlySpikeTrain > lowContrastEarlyFrameTrain(positionVector(k))...
                & lowContrastEarlySpikeTrain <= lowContrastEarlyFrameTrain(positionVector(k)+1)));
        end
        lowContrastEarlyOutput(j) = sum(momentaryOut)/length(positionVector);
        waitbar(j/(nValues), h)
        clear positionVector
    end
    spikes(i).lowContrastEarlyInput = lowContrastEarlyInput;
    spikes(i).lowContrastEarlyOutput = lowContrastEarlyOutput;
    close(h)


    %the LATE spikes during LOW contrast stimulation
    lowContrastLateSpikeTrain =  spikeTrain;
    lowContrastLateSpikeTrain( mod(lowContrastLateSpikeTrain, repeatTime) <= 2*repeatTime/3) = [];
    lowContrastLateIndex = zeros (1,length(lowContrastLateSpikeTrain));

    h = waitbar(0, 'late low contrast spikes...');
    for j=1:length(lowContrastLateSpikeTrain)
        lastFrame = find((lowContrastLateSpikeTrain(j) - frameTrain >= 0), 1, 'last');
        lowContrastLateIndex(j) = ( (lastFrame-1)*4 ...
            + ceil((lowContrastLateSpikeTrain(j) - frameTrain(lastFrame))/10) );
        waitbar(j/length(lowContrastLateSpikeTrain), h)
    end
    close(h)

    spikes(i).lowContrastLateSTE = intensityMatrix(lowContrastLateIndex,:);
    spikes(i).lowContrastLateSTA = mean(spikes(i).lowContrastLateSTE);
    lowContrastLateGeneratorPre = intensityMatrix * (spikes(i).lowContrastLateSTA)';
    lowContrastLateGeneratorPre(mod(indexVector, fourM) < 2*fourM/3) =[];
    lowContrastLateRatio = var(lowIntensityLate)/var(lowContrastLateGeneratorPre);
    spikes(i).lowContrastLateSTAnormalized = sqrt(lowContrastLateRatio) * spikes(i).lowContrastLateSTA;
    lowContrastLateGeneratorSignal = sqrt(lowContrastLateRatio) * lowContrastLateGeneratorPre;
    lowContrastLateGeneratorSorted = sort(lowContrastLateGeneratorSignal);
    lowContrastLateBoundaries = zeros(nValues, 2);
    for j=1:nValues
        lowContrastLateBoundaries(j,:) = ...
            [lowContrastLateGeneratorSorted(round((j-1)*length(lowContrastLateGeneratorSignal)/nValues)+1), ...
            lowContrastLateGeneratorSorted(round(j*length(lowContrastLateGeneratorSignal)/nValues))];
    end
    lowContrastLateInput = mean(lowContrastLateBoundaries, 2);
    lowContrastLateOutput = zeros(length(lowContrastLateInput), 1);
    lowContrastLateFrameTrainVector = reshape(frameTrain, 1, M*N);
    lowContrastLateFrameTime = mean(diff(lowContrastLateFrameTrainVector));
    lowContrastLateFrameTrainVector(mod(shortIndexVector, M) < 2*M/3) =[];
    lowContrastLateFrameTrain = [lowContrastLateFrameTrainVector;...
        lowContrastLateFrameTrainVector + lowContrastLateFrameTime/4;...
        lowContrastLateFrameTrainVector + lowContrastLateFrameTime/2;...
        lowContrastLateFrameTrainVector + 3*lowContrastLateFrameTime/4];

    h = waitbar(0, 'late low contrast input-output...');
    for j=1:nValues
        positionVector = find( (lowContrastLateGeneratorSignal >= lowContrastLateBoundaries(j,1))...
            & (lowContrastLateGeneratorSignal < lowContrastLateBoundaries(j,2)) );
        if positionVector(end) == numel(lowContrastLateFrameTrain)
            positionVector(end)=[];
        else
        end
        momentaryOut = zeros(length(positionVector),1);
        for k=1:length(positionVector)
            momentaryOut(k) = length(find(lowContrastLateSpikeTrain > lowContrastLateFrameTrain(positionVector(k))...
                & lowContrastLateSpikeTrain <= lowContrastLateFrameTrain(positionVector(k)+1)));
        end
        lowContrastLateOutput(j) = sum(momentaryOut)/length(positionVector);
        waitbar(j/(nValues), h)
        clear positionVector
    end
    spikes(i).lowContrastLateInput = lowContrastLateInput;
    spikes(i).lowContrastLateOutput = lowContrastLateOutput;
    spikes(i).channel = channel(i);
    close(h)
    clear lowContrast*
    
    % calculating some parameters of contrastadaptation
    spikes = contrastadaptationparameters(spikes, i);
    
    waitbar(i/nChannels, g)
end
close(g)