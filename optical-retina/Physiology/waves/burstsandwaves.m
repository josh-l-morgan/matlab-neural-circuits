%% LOADING FILE AND SETTING UP FILES AND DIRECTORIES FOR SAVING RESULTS
clear; clc;
[fileName pathName] = uigetfile('*.txt','Select TXT_ file',...
    '\\128.208.64.36\wonglab\Daniel\waves_Data\TXT');
sp = nexread2006([pathName fileName]);
% clear fileName pathName
% [fileName pathName] = uiputfile('*.mat', 'Save WAVE_ file',...
%     '\\128.208.64.36\wonglab\Daniel\waves_Data\WAVE');
% [xlsFileName xlsPathName] = uiputfile('*.xls', 'Save STATS_ file',...
%     '\\128.208.64.36\wonglab\Daniel\waves_Data\XLS\STATS');
sp.channels

%% USER-DEFINED INPUTS
prompt = {'max. spike separation (s)', 'min. # of spikes/burst',...
    'max. burst separation (s)', 'surprise threshold', ...
    'channels in patch 1 []', 'channnels in patch 2[]', 'episode (s) []'};
title = 'Burst and Wave Identification'; nLines = 1;
defaultAnswers = {'0.4', '4', '5', '4', '3:', '', ''};
answer = inputdlg(prompt, title, nLines, defaultAnswers);
maxSpikeDist = str2double(answer{1});
minSpikeNum = str2double(answer{2});
maxBurstDist = str2double(answer{3});
surpriseThreshold = str2double(answer{4});
channelsPatch1 = str2num(answer{5}); %#ok<ST2NM>
nChannelsPatch1 = length(channelsPatch1);
channelsPatch2 = str2num(answer{6}); %#ok<ST2NM>
nChannelsPatch2 = length(channelsPatch2);
episode = str2num(answer{7}); %#ok<ST2NM>
channels = sp.channels([channelsPatch1 channelsPatch2]);
nChannels = length(channels);
data = sp.data(:,[channelsPatch1 channelsPatch2]);

%% BURST IDENTIFICATION AND FUSION OF BURST WIHTIN CHANNELS
averageRate = zeros(nChannels,1);
for i = 1:nChannels
    clear burstData waveData
    currentData = data(:,i);
    spikeTrain = currentData(currentData >= episode(1) &...
        currentData <= episode(2));
    averageRate(i) = length(spikeTrain) / (episode(2) - episode(1));
    diffTrain = diff(spikeTrain);
    possibleStart = diffTrain < maxSpikeDist;
    j=1; % this defines the counter for fusing spikes to bursts

    %%% first identify bursts wihtin channel(i) %%%
    while j < length(possibleStart)
        if possibleStart(j) == 1
            nStart = j;
            while possibleStart(j) == 1 && j < length(possibleStart)
                j = j+1;
                nStop = j;
                currentNum = nStop - nStart + 1;
            end
            if currentNum >= minSpikeNum
                averageNum = averageRate(i) * (spikeTrain(nStop) - ...
                    spikeTrain(nStart));
                P = exp(-averageNum) * averageNum^currentNum / ...
                    factorial(currentNum);
                surprise = -log10(P);
                currentRate = currentNum / (spikeTrain(nStop) - ...
                    spikeTrain(nStart));
                if surprise > surpriseThreshold
                    if exist('burstData', 'var')
                        [nRows nCols] = size(burstData);
                        burstData(nRows+1,:) = [spikeTrain(nStart),...
                            spikeTrain(nStop), currentNum, currentRate];
                    else
                        burstData = [spikeTrain(nStart),...
                            spikeTrain(nStop), currentNum, currentRate];
                    end
                else
                end
            else
            end
        else
            j = j+1;
        end
    end

    %%% calculate some burst parameters for channel(i) %%%
    burstWidth = mean(burstData(:,2) - burstData(:,1));
    burstFraction = sum(burstData(:,2) - burstData(:,1)) / ...
        (episode(2) - episode(1));
    spikesPerBurst = mean(burstData(:,3));
    fractionSpikesInBurst = sum(burstData(:,3))/length(spikeTrain);
    rateDuringBurst = mean(burstData(:,4));

    %%% now fuse bursts into waves wihtin channel(i) %%%
    realignBursts = [burstData(1:end-1,2), burstData(2:end,1)];
    burstDiff = diff(realignBursts,1,2);
    possibleWaveStart = burstDiff < maxBurstDist;
    k=1; % this initiates the counter for fusing bursts to waves
    while k < length(possibleWaveStart)
        if possibleWaveStart(k) == 1
            nWaveStart = k;
            while possibleWaveStart(k) == 1 && k < length(possibleWaveStart)
                k = k+1;
                nWaveStop = k;
                currentWaveNum = nWaveStop - nWaveStart +1;
                % # of bursts in the current wave
            end
            currentSpikeNum = sum(burstData(nWaveStart:nWaveStop,3));
            currentSpikeRate = currentSpikeNum /...
                (burstData(nWaveStop,2) - burstData(nWaveStart,1));
            currentWave = [burstData(nWaveStart,1),...
                burstData(nWaveStop,2), currentWaveNum,...
                currentSpikeNum, currentSpikeRate];
        else
            currentWaveNum = 1;
            currentWave = [burstData(k,1), burstData(k,2),...
                currentWaveNum, burstData(k,3), burstData(k,4)];
        end
        k = k+1;
        if exist('waveData', 'var')
            [nRows nCols] = size(waveData);
            waveData(nRows+1,:) = currentWave;
        else
            waveData = currentWave;
        end
    end
    waves(i).channel = channels(i); %#ok<AGROW>
    waves(i).spikeTrain = spikeTrain; %#ok<AGROW>
    waves(i).burstData = burstData; %#ok<AGROW>
    waves(i).waveData = waveData; %#ok<AGROW>
    burstStats(i,:)  = [burstWidth, spikesPerBurst, rateDuringBurst,...
        burstFraction, ];  %#ok<AGROW>
%     burstStats(i,:)  = [burstWidth, spikesPerBurst, rateDuringBurst,...
%         fractionSpikesInBurst, burstFraction, ];  %#ok<AGROW>
    %     burstStats(i,:)  = [fractionSpikesInBurst];  %#ok<AGROW>
end

%% VERIFYING WAVES BY COMPARING ACTIVITY OF ALL UNITS IN A PATCH

for i = 1:nChannels
    clear toBeRemoved
    if i <=nChannelsPatch1
        %         burstRate = mean(burstRatePatch1);
        channelsToCompare = 1:nChannelsPatch1;
        channelsToCompare(channelsToCompare == i) =[];
    else
        %         burstRate = mean(burstRatePatch2);
        channelsToCompare = nChannelsPatch1+1 :...
            nChannelsPatch1+nChannelsPatch2;
        channelsToCompare(channelsToCompare == i) =[];
    end
    [nWaves nCols] = size(waves(i).waveData);
    for k = 1:nWaves
        waveStart = waves(i).waveData(k,1);
        waveStop = waves(i).waveData(k,2);
        spikesOnOtherChannels = zeros(length(channelsToCompare),1);
        for l = 1:length(channelsToCompare)
            m = channelsToCompare(l);
            spikesOnOtherChannels(l) = numel(waves(m).spikeTrain...
                (waves(m).spikeTrain >= waveStart &...
                waves(m).spikeTrain <= waveStop));
        end
        actualSpikes = sum(spikesOnOtherChannels);
        expectedSpikes = (episode(2) - episode(1)) *...
            sum(averageRate(channelsToCompare));
        pSpikes = exp(-expectedSpikes) * expectedSpikes^actualSpikes / ...
            factorial(actualSpikes);
        waveSurprise = -log10(pSpikes);

        if waveSurprise <= surpriseThreshold
            if exist('toBeRemoved','var')
                nEntries = length(toBeRemoved);
                toBeRemoved(nEntries+1) = k;
            else
                toBeRemoved = k;
            end
        else
        end
    end
    if exist('toBeRemoved','var')
        waves(i).waveData(toBeRemoved,:) = []; %#ok<AGROW>
    else
    end
    %%% now calculate some wave parameters %%%
    waveWidth = mean(waves(i).waveData(:,2) - waves(i).waveData(:,1));
    waveFraction = sum(waves(i).waveData(:,2) - waves(i).waveData(:,1)) /...
        (episode(2) - episode(1));
    realignWaves = [waves(i).waveData(1:end-1,2),...
        waves(i).waveData(2:end,1)];
    waveInterval = mean(diff(realignWaves,1,2));
    burstsPerWave = mean(waves(i).waveData(:,3));
    spikesPerWave = mean(waves(i).waveData(:,4));
    spikeRateInWave = mean(waves(i).waveData(:,5));
    totalNumSpikes = length(waves(i).spikeTrain);
    nSpikesInWave = sum(waves(i).waveData(:,4));
    nSpikesOutWave = totalNumSpikes - nSpikesInWave;
    fractionSpikesInWave = nSpikesInWave / totalNumSpikes;
    fractionSpikesOutWave = nSpikesOutWave / totalNumSpikes;
    timeOutWave = (episode(2) - episode(1)) - ...
        sum((waves(i).waveData(:,2) - waves(i).waveData(:,1)));
    spikeRateOutWave = nSpikesOutWave / timeOutWave;
    waveStats(i,:) = [waveWidth, waveInterval, burstsPerWave,...
        spikesPerWave, spikeRateInWave, fractionSpikesInWave,...
        spikeRateOutWave, fractionSpikesOutWave, waveFraction]; %#ok<AGROW>
end
% save([pathName fileName],'waves', 'burstStats', 'waveStats', 'averageRate')
save('\\128.208.64.36\wonglab\Daniel\waves_Data\WAVE\WAVE_test.mat',...
    'waves', 'burstStats', 'waveStats', 'averageRate')

%% EXPORTING BURST AND WAVE STATS TO EXCEL
export = [cell(channels'), num2cell(burstStats), num2cell(waveStats),...
    num2cell(averageRate)];
% export = [cell(channels'), num2cell(burstStats)];
% xlswrite([xlsPathName xlsFileName], export)
xlswrite('\\128.208.64.36\wonglab\Daniel\waves_Data\XLS\STATS\test.xls',...
    export)

%% PLOT RESULTS OF BURST AND WAVE IDENTIFICATION FOR QUALITY CONTROL
burstplot(waves, episode, channels)

