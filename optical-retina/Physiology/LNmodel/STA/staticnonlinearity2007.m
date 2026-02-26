%% LOADING FILES
clear; clc;
[fileName pathName] = uigetfile('*.txt','Select TXT_ file');
sp = nexread2006([pathName fileName]);
clear fileName pathName

[fileName pathName] = uigetfile('*.mat', 'Select STIM_ file');
load([pathName fileName])
clear fileName pathName

[fileName pathName] = uiputfile('.*mat', 'Save Results As LN2007_');

[xlsFileName xlsPathName] = uiputfile('.*xls', 'Write to XLS_');

analog = sp.data(1:50,[1 2]);
% for i=1:2
%     lastPosition = find (sp.data(:,i) >= 0, 1, 'last');
%     analog(:,i) = sp.data(1:lastPosition, i);
% end
analog = fliplr(analog);
sp.channels

%% ALIGN TIME AND DELETE CHANNELS
channelsToBeDeleted = [1 2 5:15]; %user-defined parameter
mcStart = analog(27,1);
mcStop = analog(27,2);
gaussnoiseTime = DGaussian_S;
stretchFactor = (gaussnoiseTime(end) - gaussnoiseTime(1)) / (mcStop - mcStart);
frameTrain = round((gaussnoiseTime - gaussnoiseTime(1))*1000);
channel = sp.channels;
channel(channelsToBeDeleted) = [];
[data] = deletecolumns (sp.data, channelsToBeDeleted);
channel
clear sp

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
        & (spikeTrainUnaligned < (mcStop - ((mcStop -mcStart)/2))));

    spikeTrain = round((spikeTrainUnaligned(spikesOfInterest(1) ...
        : spikesOfInterest(end)) - mcStart) * stretchFactor * 1000);
    spikeTrain(spikeTrain < 1000)=[];
    %     logicTrain = spiketologic (spikeTrain, frameTrain);

    correlationMatrix = zeros(nRows * nColumns, nTimeBins);
    momentaryMatrix = zeros(nRows * nColumns, nTimeBins);

    for j=1:length(spikeTrain)
        %this is what is done to each spike of a given channel serially
        lastFrame = find((spikeTrain(j) - frameTrain >= 0), 1, 'last');
        firstFrame = lastFrame - nTimeBins +1;
        framePosition = firstFrame : lastFrame;

        for k=1:nTimeBins
            %             randn('state', stateM(:,framePosition(k)));
            randn('state', framePosition(k));
            checkerboard = SD * randn(nRows * nColumns,1);
            checkerboard(1)=0;
            checkerboard(checkerboard > 0.5) = 0.5;
            checkerboard(checkerboard < -0.5) = -0.5;
            momentaryMatrix(:,k) = checkerboard;
        end
        correlationMatrix = correlationMatrix + momentaryMatrix;
    end

    STA = reshape(correlationMatrix, nRows * nColumns * nTimeBins,1)...
        /length(spikeTrain);
    firstFrameSecondHalf = round(nFrames/2);
    nFramesSecondHalf = length(firstFrameSecondHalf:nFrames);
    generatorSignal = zeros(nFramesSecondHalf,1);
    recentStimulus = zeros(nRows * nColumns, nTimeBins);
    %     for j=nTimeBins-1:nFrames
    for j=firstFrameSecondHalf:nFrames
        %this is what is done to each frame serially
        for k=1:nTimeBins
            %             randn('state',stateM(:,j - nTimeBins + k))
            randn('state',j - nTimeBins + k)
            checkerboard = SD * randn(nRows * nColumns,1);
            checkerboard(1)=0;
            checkerboard(checkerboard > 0.5) = 0.5;
            checkerboard(checkerboard < -0.5) = -0.5;
            recentStimulus(:,k) = checkerboard;
        end
        generatorSignal(j - firstFrameSecondHalf + 1) = ...
            reshape(recentStimulus, 1, nRows * nColumns * nTimeBins) * STA;
    end
    ratio = SD^2 / var(generatorSignal);
    STA = STA * sqrt(ratio);
    generatorSignal = sqrt(ratio) * generatorSignal;
    generatorSorted = sort(generatorSignal);
    generatorBoundaries = zeros(nValues, 2);
    for j=1:nValues
        generatorBoundaries(j,:) = ...
            [generatorSorted(round((j-1)*nFramesSecondHalf/nValues)+1), ...
            generatorSorted(round(j*nFramesSecondHalf/nValues))];
    end
    input = zeros(nValues, 1);
    %     outputProbability = zeros(nValues, 1);
    outputRate = zeros(nValues, 1);

    for j=1:nValues
        positionVector = ...
            find( (generatorSignal >= generatorBoundaries(j,1))...
            & (generatorSignal < generatorBoundaries(j,2)) );
        if positionVector(end) == nFrames
            positionVector(end)=[];
        else
        end
        input(j) = median(generatorSignal(positionVector));
        %         outputProbability(j) = mean(logicTrain(positionVector));
        momentaryOut = zeros(length(positionVector),1);
        for k=1:length(positionVector)
            momentaryOut(k) = ...
                length(find(spikeTrain > frameTrain(positionVector(k))...
                & spikeTrain <= frameTrain(positionVector(k)+1)));
        end
        countToRate = 1000/mean(diff(frameTrain));
        outputRate(j) = (sum(momentaryOut)/length(positionVector)) * countToRate;
        clear positionVector
    end

    %     muSeed = input(find(outputProbability <=(max(outputProbability)...
    %         + min(outputProbability))/2, 1, 'last'));
    %
    %     probabilityFit0 = ...
    %         [10, -muSeed, max(outputProbability), min(outputProbability)];
    %     %probabilityFit0 =[sensitivity, maintainedDrive, maximum, minimum]
    %     scalingOptions = statset ('MaxIter', 10000);
    %     [probabilityFit] = nlinfit(input, outputProbability, @inputoutputfit,...
    %         probabilityFit0, scalingOptions);
    %     p = normcdf ( (probabilityFit(1) * generatorSignal + probabilityFit(2)),...
    %         0, 1);
    %     probabilityTrain =  p * probabilityFit(3) + probabilityFit(4);

    muSeed = input(find(outputRate <=(max(outputRate)...
        + min(outputRate))/2, 1, 'last'));
    scalingOptions = statset ('MaxIter', 10000);
    rateFit0 = [10, -muSeed];
    %rateFit0 =[sensitivity, maintainedDrive, maximum, minimum]
    passInput.x = input;
    passInput.maximumRate = max(outputRate);
    passInput.minimumRate = min(outputRate);
    [rateFit] = nlinfit(passInput, outputRate, @inputoutputfit,...
        rateFit0, scalingOptions);



    if i == 1
        spikes.channel = channel;
        spikes.STA = STA';
        spikes.input = input';
        spikes.generatorSignal = generatorSignal';
        spikes.outputRate = outputRate';
        spikes.rateFit = rateFit;
        save([pathName fileName], 'spikes');
        clear spike* momentaryMatrix correlationMatrix pass*
    elseif i == nChannels
        load([pathName fileName])
        spikes.STA(i,:) = STA';
        spikes.input(i,:) = input';
        spikes.generatorSignal(i,:) = generatorSignal';
        spikes.outputRate(i,:) = outputRate';
        spikes.rateFit(i,:) = rateFit;
        save([pathName fileName], 'spikes');
        export = [cell(channel'), num2cell(spikes.rateFit),...
            num2cell(spikes.outputRate(:,end))];
        xlswrite([xlsPathName xlsFileName], export)
        clear spike* momentaryMatrix correlationMatrix pass*
    else
        load([pathName fileName])
        spikes.STA(i,:) = STA';
        spikes.input(i,:) = input';
        spikes.generatorSignal(i,:) = generatorSignal';
        spikes.outputRate(i,:) = outputRate';
        spikes.rateFit(i,:) = rateFit;
        save([pathName fileName], 'spikes');
        clear spike* momentaryMatrix correlationMatrix pass*
    end
    waitbar(i/nChannels, g)
end
close(g)
toc





