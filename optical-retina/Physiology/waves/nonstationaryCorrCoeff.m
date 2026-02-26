%% LOADING FILE
clear; clc;
[fileName pathName] = uigetfile('*.txt','Select TXT_ file');
sp = nexread2006([pathName fileName]);
clear fileName pathName
[fileName pathName] = uiputfile('*.xls','Save .xls file as'); 
sp.channels

%% USER DEFINED INPUT
prompt = {'start of cross-correlation (s)', 'end of cross-correlation (s)',...
    'bin size (ms)', 'nBins per average',...
    'index of channels to compare([])'};
title = 'Cross-correlation';
nLines = 1;
answer = inputdlg(prompt,title,nLines);

corrStart = str2double(answer{1});
corrStop = str2double(answer{2});
binSize = str2double(answer{3})/1000; %changes units to seconds
nBinsPerAverage = str2double(answer{4}); %# of bins per averaging window
channelsToCompare = str2num(answer{5});%#ok<ST2NM>
nChannels = length(channelsToCompare);
channels = sp.channels(channelsToCompare);
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
averageTrain = zeros(nBins-1, nChannels);
tic
for i = 1:length(channelsToCompare)
    rateTrain(:,i) = spiketorate (sp.data(:,channelsToCompare(i)),...
        binTrain, binSize);
    averageTrain(:,i) = ...
        slidingwindowaverage (nBinsPerAverage, rateTrain(:,i));
end
toc
clear sp*
%% CORRELATION COEFFICIENTS
[R P] = corrcoef(rateTrain);

A = rateTrain - averageTrain;
B = (A'*A) / (nBins-1);
C = zeros(nChannels);

for m = 1:nChannels
    for n = 1:nChannels
        if n <= m
            C(m,n) = B(m,n) / sqrt(B(m,m)*B(n,n));
        else
        end
    end
end

%% OUTPUT
for i = 1:nChannels
    for j = 1:nChannels - i
        if i==1
            cellA(j) = {char(channels(i))}; %#ok<AGROW>
            cellB(j) = {char(channels(i+j))}; %#ok<AGROW>
            Corr (j) = C(j+i,i); %#ok<AGROW>
        else
            if j==1
                lastEntry = length(cellA);
            else
            end
            cellA(lastEntry + j) = {char(channels(i))}; %#ok<AGROW>
            cellB(lastEntry + j) = {char(channels(i+j))}; %#ok<AGROW>
            Corr (lastEntry + j) = C(j+i,i); %#ok<AGROW>
        end
    end
end

export = [cellA' cellB' num2cell(Corr')];
xlswrite([pathName fileName], export)
