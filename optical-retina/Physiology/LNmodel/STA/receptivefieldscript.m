%% LOADING FILES
[fileName directoryName] = uigetfile('*.txt','Select .txt file');
cd(directoryName);
sp = nexread2006(fileName);
uiopen('load')
cd('F:\Program Files\MATLAB\R2006a\work')
clear Cadapt* DIO SD Nrepeats Nframes Nepisodes DOnOff* fileName directoryName

for i=1:2
    lastPosition = find (sp.data(:,i) >= 0, 1, 'last');
    analog(:,i) = sp.data(1:lastPosition, i);
end
analog = fliplr(analog);
sp.channels

%% ALIGN TIME AND DELETE CHANNELS
channelsToBeDeleted = [1 2 3]; %user-defined parameter
mcStart = analog(32,1);
mcStop = analog(32,2);
stretchFactor = (DGaussian_S(end) - DGaussian_S(1)) / (mcStop - mcStart);
frameTrain = round((DGaussian_S - DGaussian_S(1))*1000);
channel = sp.channels;
channel(channelsToBeDeleted) = [];
[data] = deletecolumns (sp.data, channelsToBeDeleted);
channel
%% SPIKES OF INTEREST
timeWindow = 511; %temporal ROI around spike in ms
rowsOfCheckers = 60;
columnsOfCheckers = 80;
standardDeviation = 0.175;
meanGray = 0.5;
[m nChannels] = size(data);


g = waitbar(0, 'channels...');
for i=1:nChannels %this is what is done to each channel serially

    if (data(end,i) > 0)
        spikeTrainUnaligned = data(:,i);
    else
        lastPosition =  find (data(:,i) >= 0, 1, 'last');
        spikeTrainUnaligned = data(1:lastPosition,i);
    end

    spikesOfInterest = find ((mcStart < spikeTrainUnaligned) & (spikeTrainUnaligned < mcStop));
    spikeTrain = round((spikeTrainUnaligned(spikesOfInterest(1):spikesOfInterest(end)) - mcStart) * stretchFactor * 1000);
    correlationMatrix = zeros(rowsOfCheckers * columnsOfCheckers, timeWindow);
    momentaryMatrix = zeros(rowsOfCheckers * columnsOfCheckers, timeWindow);

    h = waitbar(0, 'spikes...');
    for j=1:length(spikeTrain) %this is what is done to each spike of a given channel serially
        left = find((spikeTrain(j) - 500 - frameTrain >= 0), 1, 'last');
        right = find((spikeTrain(j) + 10 - frameTrain >= 0), 1, 'last');
        frameTrainExcerpt = frameTrain (left : right);
        firstColumn = zeros (length(frameTrainExcerpt),1);

        for k=1:length(frameTrainExcerpt)
            framePosition = find (frameTrain == frameTrainExcerpt(k));
            randn('state', framePosition);
            checkerboard = standardDeviation * randn(rowsOfCheckers,columnsOfCheckers) + meanGray;
            checkerboard(checkerboard > 1) = 1;
            checkerboard(checkerboard < 0) = 0;
            intensity(:,k) = reshape(checkerboard',rowsOfCheckers*columnsOfCheckers,1);
            firstColumn(k) = (500 - (spikeTrain(j) - frameTrainExcerpt(k)));
        end
        
        for k=1:length(frameTrainExcerpt)
            if (k==1)
                for l=1:length(1:(firstColumn(k+1)-1))
                    momentaryMatrix(:,l) = intensity(:,k);
                end

            elseif ((k > 1) && (k < length(frameTrainExcerpt)))
                for l=1:length(firstColumn(k):firstColumn(k+1)-1)
                    momentaryMatrix(:,firstColumn(k)-1+l) = intensity(:,k);
                end

            else
                for l=1:length(firstColumn(k):timeWindow)
                    momentaryMatrix(:,firstColumn(k)-1+l) = intensity(:,k);
                end
            end
        end
        correlationMatrix = correlationMatrix + momentaryMatrix;
        waitbar(j/length(spikeTrain), h)
    end
    STA = correlationMatrix / length(spikeTrain);
    close (h)
    if i == 1
        spikes = struct('channel',{},'STA',{});
        spikes(i).channel = channel(i);
        spikes(i).STA = STA;
        save STA_P1moVnnIpp_0220a spikes; %replace trial by name of the file to be saved
        clear spike* momentaryMatrix correlationMatrix
    else
        load STA_P1moVnnIpp_0220a
        spikes(i).channel = channel(i);
        spikes(i).STA = STA;
        save STA_P1moVnnIpp_0220a spikes; %replace trial by name of the file to be saved
        clear spike* momentaryMatrix correlationMatrix
    end
    waitbar(i/nChannels, g)
end
close (g)
