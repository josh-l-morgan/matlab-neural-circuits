%% LOADING FILE
clear; clc;
[fileName pathName] = uigetfile('*.txt','Select TXT_ file');
sp = nexread2006([pathName fileName]);
clear fileName pathName
[fileName pathName] = uiputfile('*.xls','Save .xls file as');
sp.channels

%% USER DEFINED INPUT
prompt = {'start of cross-correlation (s)', 'end of cross-correlation (s)',...
    'bin size (ms)', '# bins per average', 'correlation lag (ms)',...
    'index of channels to compare([])'};
title = 'Cross-correlation';
nLines = 1;
defAns = {'', '', '50', '100', '1000',''};
answer = inputdlg(prompt,title,nLines,defAns);

corrStart = str2double(answer{1});
corrStop = str2double(answer{2});
binSize = str2double(answer{3})/1000; %changes units to seconds
nBinsPerAverage = str2double(answer{4}); %# of bins per averaging window
maxLag = ceil(str2double(answer{5}))/(binSize*1000);
channelsToCompare = str2num(answer{6});%#ok<ST2NM>
nChannels = length(channelsToCompare);
channels = sp.channels(channelsToCompare);
%% SPIKETRAIN TO LOGIC
if rem((corrStop-corrStart)/binSize,...
        round((corrStop-corrStart)/binSize)) == 0
    binTrain = corrStart:binSize:corrStop;
    nBins = length(binTrain);
else
    nBins = ceil( (corrStop-corrStart)/binSize );
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


%% XCORR
[bC, blags] =  xcorr(rateTrain,maxLag,'none');
zeroC = bC(maxLag+1,:);
oneC = mean(bC([1 end],:));
BSI = (zeroC - oneC)./(zeroC + oneC);

normalizedTrain = rateTrain - averageTrain;
[C, lags] = xcorr(normalizedTrain, maxLag, 'unbiased');
[M, N] = size(C);
normCorr = zeros(M,N);
V = var(normalizedTrain);
for i = 1:N
    [p q] = ind2sub([nChannels nChannels],i);
    normCorr(:,i) = C(:,i)./sqrt(V(p)*V(q));
end
pZero = (normCorr(maxLag+1,:));
pOne = mean(normCorr([1 end],:));

%% CORRELATION COEFFICIENTS
% A = rateTrain - averageTrain;
% B = (A'*A) / (nBins-1);
% D = zeros(nChannels);
%
% for m = 1:nChannels
%     for n = 1:nChannels
%         if n <= m
%             D(m,n) = B(m,n) / sqrt(B(m,m)*B(n,n));
%         else
%         end
%     end
% end

%% OUTPUT
for i = 1:nChannels
    for j = 1:nChannels - i
        if i==1
            cellA(j) = {char(channels(i))}; %#ok<AGROW>
            cellB(j) = {char(channels(i+j))}; %#ok<AGROW>
            %             Corr (j) = D(j+i,i); %#ok<AGROW>
            Burst(j) = BSI(i+j); %#ok<AGROW>
            cZero(j) = pZero(i+j); %#ok<AGROW>
            cOne(j) = pOne(i+j); %#ok<AGROW>
        else
            if j==1
                lastEntry = length(cellA);
            else
            end
            cellA(lastEntry + j) = {char(channels(i))}; %#ok<AGROW>
            cellB(lastEntry + j) = {char(channels(i+j))}; %#ok<AGROW>
            %             Corr(lastEntry + j) = D(j+i,i); %#ok<AGROW>
            Burst(lastEntry + j) = BSI((i-1)*nChannels+i+j); %#ok<AGROW>
            cZero(lastEntry + j) = pZero((i-1)*nChannels+i+j);  %#ok<AGROW>
            cOne(lastEntry + j) = pOne((i-1)*nChannels+i+j);  %#ok<AGROW>
        end
    end
end

export = [cellA' cellB' num2cell(cZero') num2cell(cOne') num2cell(Burst')];
xlswrite([pathName fileName], export)
