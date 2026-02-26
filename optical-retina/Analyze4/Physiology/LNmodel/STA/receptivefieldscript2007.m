%% LOADING FILES
[fileName pathName] = uigetfile('*.txt','Select TXT_ file');
sp = nexread2006([pathName fileName]);
clear fileName pathName

[fileName pathName] = uigetfile('*mat', 'Select STIM_ file');
load([pathName fileName])
clear fileName pathName

[fileName pathName] = uiputfile('*mat', 'Save Results As');

for i=1:2
    lastPosition = find (sp.data(:,i) >= 0, 1, 'last');
    analog(:,i) = sp.data(1:lastPosition, i);
end
analog = fliplr(analog);
sp.channels

%% ALIGN TIME AND DELETE CHANNELS
channelsToBeDeleted = [1 2]; %user-defined parameter
mcStart = analog(27,1);
mcStop = analog(27,2);
stretchFactor = (gaussnoiseTime(end) - gaussnoiseTime(1)) / (mcStop - mcStart);
frameTrain = round((gaussnoiseTime - gaussnoiseTime(1))*1000);
channel = sp.channels;
channel(channelsToBeDeleted) = [];
[data] = deletecolumns (sp.data, channelsToBeDeleted);
channel

%% SPIKES OF INTEREST
timeBeforeSpike = 800;  %time (ms) before spike for which average is computed
timeAfterSpike = 10;    %time (ms) after spike for average is computed
timeWindow = timeBeforeSpike + 1 + timeAfterSpike; %t
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

    for j=1:length(spikeTrain) %this is what is done to each spike of a given channel serially
        left = find((spikeTrain(j) - timeBeforeSpike - frameTrain >= 0), 1, 'last');
        right = find((spikeTrain(j) + timeAfterSpike - frameTrain >= 0), 1, 'last');
        frameTrainExcerpt = frameTrain (left : right);
        firstColumn = zeros (length(frameTrainExcerpt),1);

        for k=1:length(frameTrainExcerpt)
            framePosition = find (frameTrain == frameTrainExcerpt(k));
            randn('state', framePosition);
            checkerboard = standardDeviation * randn(rowsOfCheckers,columnsOfCheckers) + meanGray;
            checkerboard(checkerboard > 1) = 1;
            checkerboard(checkerboard < 0) = 0;
            intensity(:,k) = reshape(checkerboard',rowsOfCheckers*columnsOfCheckers,1);
            firstColumn(k) = (800 - (spikeTrain(j) - frameTrainExcerpt(k)));
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
    end
    STA = correlationMatrix / length(spikeTrain);
    if i == 1
        spikes = struct('channel',{},'STA',{});
        spikes(i).channel = channel(i);
        spikes(i).STA = STA;
        save([pathName fileName], 'spikes');
        clear spike* momentaryMatrix correlationMatrix
    else
        load([pathName fileName])
        spikes(i).channel = channel(i);
        spikes(i).STA = STA;
        save([pathName fileName], 'spikes');
        clear spike* momentaryMatrix correlationMatrix
    end
    waitbar(i/nChannels, g)
end
close (g)
