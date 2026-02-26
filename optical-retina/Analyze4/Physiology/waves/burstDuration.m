%% IMPORTING FILE
clear; clc;
% [fileName pathName] = uigetfile('*.xls','Select ACORR_ file',...
%     '\\128.208.64.36\wonglab\Daniel\waves_Data\XLS\Autocorrelation');
% importAutocorrelation([pathName fileName])
importAutocorrelation('\\128.208.64.36\wonglab\Daniel\waves_Data\XLS\Autocorrelation\AUTO_P12.xls')
clear fileName pathName
%% BURST DURATION AS WHM OF AUTOCORRELOGRAMM
[nBins nChannels] = size(data);
maxima = max(data);
up = zeros(nChannels,1);
down = zeros(nChannels,1);
for i = 1:nChannels
    up(i) = find(data(:,i) >= maxima(i)/2, 1,'first');
    down(i) = find(data(:,i) >= maxima(i)/2, 1,'last');
end
burstWidth = (down - up) * 0.015;
%% EXPORT RESULTS TO EXCEL
export = [cell(textdata'), num2cell(burstWidth)];
% [fileName pathName] = uiputfile('*.xls','Save as',...
%     '\\128.208.64.36\wonglab\Daniel\waves_Data\XLS\Autocorrelation');
% xlswrite([pathName fileName], export)
xlswrite('\\128.208.64.36\wonglab\Daniel\waves_Data\XLS\Autocorrelation\STAT_P12.xls', export)


