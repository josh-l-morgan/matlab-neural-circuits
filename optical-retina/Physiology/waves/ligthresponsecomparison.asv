%% LOADING TXT FILE
clear; clc;
[fileName pathName] = uigetfile('*.txt','Select TXT_ file',...
    '\\128.208.64.36\wonglab\Daniel\waves_Data\TXT');
sp = nexread2006([pathName fileName]);
sp.channels

%% INPUTDLG
prompt = {'bin size (s)', 'channels []', '0 = OFF 1 = ON []'...
    'ctrl episode (s) []', 'drug episode (s) []'};
title = 'light response comparison';
defAns = {'0.05', '3:', '', '', ''};
answer = inputdlg(prompt,title,1,defAns);

binSize = str2double(answer{1});
channels = str2num(answer{2}); %#ok<ST2NM>
responseType = str2num(answer{3}); %#ok<ST2NM>
ctrlEpisode = str2num(answer{4}); %#ok<ST2NM>
drugEpisode = str2num(answer{5}); %#ok<ST2NM>
nChannels = length(channels);

%% DETERMINE SPIKE RATES
rates = zeros(nChannels,2);
for h=1:2
    onTime = sp.data(:,1);
    offTime = sp.data(:,2);
    if h==1
        episode = ctrlEpisode;
    else
        episode = drugEpisode;
    end
    onTime = onTime(onTime > episode(1) & onTime < episode(2));
    offTime = offTime(offTime > episode(1) & offTime < episode(2));
    nTrials = length(onTime);
    off = offTime(offTime > onTime(round(nTrials/2)) &...
        offTime < onTime(round(nTrials/2)+1));
    onDuration = off - onTime(round(nTrials/2));
    for i=1:nChannels
        spikes = sp.data(:,channels(i));
        spikes(spikes < episode(1) | spikes > episode(2)) = [];
        for j=1:nTrials
            if j < nTrials
                spikes(spikes >= onTime(j) & spikes < onTime(j+1)) =...
                    spikes(spikes >= onTime(j) & spikes < onTime(j+1))...
                    - onTime(j);
            else
                spikes(spikes >= onTime(j)) = spikes(spikes >= onTime(j))...
                    - onTime(j);
            end
        end
        trialLength = round(mean(diff(onTime)));
        peristimulusTime = 0:binSize:trialLength;
        firstOff = find(peristimulusTime > onDuration, 1, 'first');
        nBins= length(peristimulusTime)-1;
        peristimulusRate = zeros(nBins,1);
        for k=1:nBins
            peristimulusRate(k) = length(spikes(spikes > peristimulusTime(k)...
                & spikes <= peristimulusTime(k+1))) / (binSize * nTrials) ;
        end
        if responseType(i) == 1
            rates(i,h) = max(peristimulusRate(1:firstOff-1));
        else
            rates(i,h) = max(peristimulusRate(firstOff:end));
        end
    end
end

% plot(rates(responseType==1,1), rates(responseType==1,2),'ro')
% hold on
% plot(rates(responseType==0,1), rates(responseType==0,2),'ko')

onRates = rates(responseType==1,:);
offRates = rates(responseType==0,:);
