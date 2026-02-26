clear; clc; 
% clf;
[fileName pathName] = uigetfile('*.txt','Select TXT_ file',...
    '\\128.208.64.36\wonglab\Daniel\waves_Data\TXT');
sp = nexread2006([pathName fileName]);
sp.channels
%% INPUTDLG
prompt = {'channels []'};
title = 'Select Channels';
answer = inputdlg(prompt, title);
channel = str2num(answer{1}); %#ok<ST2NM>
nChannels = length(channel);

%% DON'T KNOW
firstTrial = 23; % 4 22
nTrials = 15;
firstSpike = zeros(nTrials, nChannels);
offTimer = sp.data(firstTrial:firstTrial+nTrials,1);
for i=1:nChannels
    spikesOfChannel = sp.data(:,channel(i));
    spikesOfChannel(spikesOfChannel<0)=[];
    for j=1:nTrials
        spikesInTrial = spikesOfChannel;
        spikesInTrial(spikesInTrial < offTimer(j)) = [];
        firstSpike(j,i) = spikesInTrial(1)-offTimer(j); 
    end
end

