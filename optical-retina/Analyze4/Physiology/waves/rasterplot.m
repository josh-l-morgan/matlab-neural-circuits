%% LOADING FILE AND SETTING UP FILES AND DIRECTORIES FOR SAVING RESULTS
clear; clc; clf;
[fileName pathName] = uigetfile('*.txt','Select TXT_ file',...
    '\\128.208.64.36\wonglab\Daniel\waves_Data\TXT');
sp = nexread2006([pathName fileName]);
sp.channels

%% INPUTDLG TO SPECIFY TIME AND CHANNELS
prompt = {'episode', 'channels in order to be plotted'};
title = 'Rasterplot'; nLines =1;
answer = inputdlg(prompt, title, nLines);
episode = str2num(answer{1}); %#ok<ST2NM>
channelOrder = str2num(answer{2}); %#ok<ST2NM>
channels = sp.channels(channelOrder);

%% PLOT RASTERS
nChannels = length(channelOrder);
hAxes = axes;
for i = 1:nChannels
    spikeTrain = sp.data(:,channelOrder(i));
    spikeTrain(spikeTrain < episode(1) |...
        spikeTrain > episode(2)) = [];
    nSpikes = length(spikeTrain);
    for j = 1:nSpikes
        plot(hAxes, [spikeTrain(j) spikeTrain(j)], [i i+0.8],'k');
        hold on
    end
end

set(    hAxes,...
    'YTickMode', 'manual',...
    'YTickLabelMode', 'manual',...
    'YTick', (1:nChannels),...
    'YTickLabel',channels');

axis([episode(1) episode(2) 0.8 nChannels+1]);

