clear; clc; clf;
[fileName pathName] = uigetfile('*.txt','Select TXT_ file',...
    '\\128.208.64.36\wonglab\Daniel\waves_Data\TXT');
sp = nexread2006([pathName fileName]);
sp.channels
%% INPUTDLG
prompt = {'channels []'};
title = 'Select Channels';
answer = inputdlg(prompt, title);
channel = str2num(answer{1});
nChannels = length(channel);

%% DON'T KNOW
firstTrial = 23; %4 
nTrials = 15;
onTimer = sp.data(firstTrial:firstTrial+nTrials,2);
for i=1:nChannels
    spikesOfChannel = sp.data(:,channel(i));
    spikesOfChannel(spikesOfChannel<0)=[];
    h = subplot(nChannels,1,i);
    for j=1:nTrials
        spikesInTrial = spikesOfChannel;
        spikesInTrial(spikesInTrial < onTimer(j) |...
            spikesInTrial > onTimer(j+1))=[];
        nSpikes = length(spikesInTrial);
        for k=1:nSpikes
           plot(h,[spikesInTrial(k)-onTimer(j) spikesInTrial(k)-onTimer(j)], [j-1 j-0.2],'k');
           hold on
        end
        set(h,'xlim',[0 8],'ylim',[-0.2 nTrials])
    end
end

