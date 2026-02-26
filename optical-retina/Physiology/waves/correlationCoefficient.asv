%% LOADING FILE
clear; clc;
[fileName pathName] = uigetfile('*.txt','Select TXT_ file');
sp = nexread2006([pathName fileName]);
clear fileName pathName
sp.channels

%% USER DEFINED INPUT
prompt = {'start of cross-correlation (s)', 'end of cross-correlation (s)',...
    'bin size (ms)', 'index of channels to compare([])'};
title = 'Cross-correlation';
nLines = 1;
answer = inputdlg(prompt,title,nLines);

corrStart = str2double(answer{1});
corrStop = str2double(answer{2});
binSize = str2double(answer{3})/1000; %changes units to seconds
channelsToCompare = str2num(answer{4});%#ok<ST2NM>
nChannels = length(channelsToCompare);
%% SPIKETRAIN TO LOGIC
if rem((corrStop-corrStart)/binSize,...
        round((corrStop-corrStart)/binSize)) == 0
    binTrain = corrStart:binSize:corrStop;
    nBins = length(binTrain);
else
    nBins = ceil( (corrStop-corrStart)/binSize )
    corrStop = corrStart + nBins * binSize;
    binTrain = corrStart:binSize:corrStop;
end

rateTrain = zeros(nBins-1, nChannels);

for i = 1:length(channelsToCompare)
    rateTrain(:,i) = spiketorate (sp.data(:,channelsToCompare(i)),...
        binTrain, binSize);
end

%% CORRELATION COEFFICIENTS
[R P] = corrcoef(rateTrain);
clc

%% OUTPUT

channels = sp.channels(channelsToCompare);
disp(sp.channels(channelsToCompare))
disp(R)
disp(mean(rateTrain))


