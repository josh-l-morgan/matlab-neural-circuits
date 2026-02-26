%spiketimesxls writes channel names and spike times to an excel file

%% LOADING FILE
clear; clc;
[fileName pathName] = uigetfile('*.txt','Select TXT_ file');
sp = nexread2006([pathName fileName]);
clear fileName pathName
[fileName pathName] = uiputfile('*.xls','Save .xls file as');
sp.channels

%% USER DEFINED INPUT
prompt = {'episode time (s) []', 'channels []'};
title = 'spike time export';
defAns = {'', '3:'};
answer = inputdlg(prompt,title,1,defAns);

episode = str2num(answer{1}); %#ok<ST2NM>
channels = str2num(answer{2}); %#ok<ST2NM>
nChannels = length(channels);
