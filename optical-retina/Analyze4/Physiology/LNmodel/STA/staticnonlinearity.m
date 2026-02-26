%% LOADING FILES
[fileName directoryName] = uigetfile('*.txt','Select .txt file');
cd(directoryName);
sp = nexread2006(fileName);
uiopen('load')
cd('C:\Program Files\MATLAB\R2006a\work')
clear Cadapt* DIO SD Nrepeats Nframes Nepisodes ...
    DOnOff* fileName directoryName

for i=1:2
    lastPosition = find (sp.data(:,i) >= 0, 1, 'last');
    analog(:,i) = sp.data(1:lastPosition, i);
end
analog = fliplr(analog);
sp.channels

%% ALIGN TIME AND DELETE CHANNELS
channelsToBeDeleted = [1 2 3 4 5 7 8]; %user-defined parameter
mcStart = analog(27,1);
mcStop = analog(27,2);
stretchFactor = (DGaussian_S(end) - DGaussian_S(1)) / (mcStop - mcStart);
frameTrain = round((DGaussian_S - DGaussian_S(1))*1000);
channel = sp.channels;
channel(channelsToBeDeleted) = [];
[data] = deletecolumns (sp.data, channelsToBeDeleted);
channel
%% COMPUTING STA AND GENERATOR SIGNAL
nTimeBins = 13;
nValues = 50;
nRows = 60;
nColumns = 80;
nFrames = length(frameTrain);
SD = 0.175;
[m nChannels] = size(data);
tic
g = waitbar(0, 'channels...');
for i=1:nChannels %this is what is done to each channel serially

    if (data(end,i) > 0)
        spikeTrainUnaligned = data(:,i);
    else
        lastPosition =  find (data(:,i) >= 0, 1, 'last');
        spikeTrainUnaligned = data(1:lastPosition,i);
    end

    spikesOfInterest = find ((mcStart < spikeTrainUnaligned) ...
        & (spikeTrainUnaligned < mcStop));

    spikeTrain = round((spikeTrainUnaligned(spikesOfInterest(1) ...
        : spikesOfInterest(end)) - mcStart) * stretchFactor * 1000);

    correlationMatrix = zeros(nRows * nColumns, nTimeBins);
    momentaryMatrix = zeros(nRows * nColumns, nTimeBins);

    h = waitbar(0, 'spikes...');
    for j=1:length(spikeTrain)
        %this is what is done to each spike of a given channel serially
        lastFrame = find((spikeTrain(j) - frameTrain >= 0), 1, 'last');
        firstFrame = lastFrame - nTimeBins +1;
        framePosition = firstFrame : lastFrame;

        for k=1:nTimeBins
            randn('state', framePosition(k));
            checkerboard = SD * randn(nRows * nColumns,1);
            checkerboard(checkerboard > 0.5) = 0.5;
            checkerboard(checkerboard < -0.5) = -0.5;
            momentaryMatrix(:,k) = checkerboard;
        end
        correlationMatrix = correlationMatrix + momentaryMatrix;
        waitbar(j/length(spikeTrain), h)
    end

    STA = reshape(correlationMatrix, nRows * nColumns * nTimeBins,1)...
        /length(spikeTrain);
    close (h)

    generatorSignal = zeros(nFrames,1);
    recentStimulus = zeros(nRows * nColumns, nTimeBins);

    h = waitbar(0, 'generator signal...');
    for j=nTimeBins-1:nFrames
        %this is what is done to each frame serially
        for k=1:nTimeBins
            randn('state',j - nTimeBins + k)
            checkerboard = SD * randn(nRows * nColumns,1);
            checkerboard(checkerboard > 0.5) = 0.5;
            checkerboard(checkerboard < -0.5) = -0.5;
            recentStimulus(:,k) = checkerboard;
        end
        generatorSignal(j) = ...
            reshape(recentStimulus, 1, nRows * nColumns * nTimeBins) * STA;
        waitbar(j/(nFrames), h)
    end
    close(h)
    generatorSorted = sort(generatorSignal);
    generatorBoundaries = zeros(nValues, 2);
    for j=1:nValues
        generatorBoundaries(j,:) = ...
            [generatorSorted(round((j-1)*nFrames/nValues)+1), ...
            generatorSorted(round(j*nFrames/nValues))];
    end
    input = mean(generatorBoundaries, 2);
    output = zeros(nValues, 1);

    h = waitbar(0, 'input-output...');
    for j=1:nValues
        positionVector = ...
            find( (generatorSignal >= generatorBoundaries(j,1))...
            & (generatorSignal < generatorBoundaries(j,2)) );
        if positionVector(end) == nFrames
            positionVector(end)=[];
        else
        end
        momentaryOut = zeros(length(positionVector),1);
        for k=1:length(positionVector)
            momentaryOut(k) = ...
                length(find(spikeTrain > frameTrain(positionVector(k))...
                & spikeTrain <= frameTrain(positionVector(k)+1)));
        end
        output(j) = sum(momentaryOut)/length(positionVector);
        spikes(i).position(j).positionVector = positionVector;
        waitbar(j/(nValues), h)
        clear positionVector
    end
    spikes(i).input = input;
    spikes(i).output = output;
    close(h)
    
    spikes(i).logicalTrain = spiketologic (spikeTrain, frameTrain);
    
    waitbar(i/nChannels, g)
end
close(g)
toc





