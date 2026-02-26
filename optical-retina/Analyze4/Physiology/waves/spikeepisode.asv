%spikeepisode reads txt file and save a mat file containing a structure of

%% LOADING FILE
clear; clc;
[fileName pathName] = uigetfile('*.txt','Select TXT_ file');
sp = nexread2006([pathName fileName]);
clear fileName pathName
[fileName pathName] = uiputfile('*.mat','Save .mat file as');
sp.channels

%% USER DEFINED INPUT
prompt = {'episode time (s) []', 'channels []', '0 = OFF 1 = ON []'};
title = 'spike time export';
defAns = {'', '3:',''};
answer = inputdlg(prompt,title,1,defAns);

episode = str2num(answer{1}); %#ok<ST2NM>
channels = str2num(answer{2}); %#ok<ST2NM>
responseType = str2num(answer{3}); %#ok<ST2NM>
nChannels = length(channels);

%% DO THAT THING

for i = 1:nChannels
    clear currentData currentResponse
    currentData = sp.data(:,channels(i));
    currentData(currentData < episode(1) | currentData > episode(2)) = [];
    spikes(i).channel = sp.channels(channels(i)); %#ok<AGROW>
    if responseType(i) == 0
        spikes(i).responseType = 'OFF'; %#ok<AGROW>
    else
        spikes(i).responseType = 'ON'; %#ok<AGROW>
    end
    spikes(i).data = currentData; %#ok<AGROW>
end

save([pathName fileName], 'spikes');
