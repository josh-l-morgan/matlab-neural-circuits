%% LOADING FILES AND GETTING NUMBER OF CHANNELS
[fileName,directoryName] = uigetfile;
load([directoryName fileName])
cd('F:\Program Files\MATLAB\R2006a\work')
s = whos('spikes');
nChannels = s.size(2);


%% ACCUMULATING SIGMA
if exist ('lowContrastEarly', 'var')
    [m n] = size(lowContrastEarly);
    lowContrastEarly = [lowContrastEarly; zeros(nChannels,3)];
    for i=1:nChannels
        lowContrastEarly(i + m, 1:2) = spikes(i).lowContrastEarlyParameters;
        if abs(max(spikes(i).highContrastLateSTA)) > abs(min(spikes(i).highContrastLateSTA)) %ON cells
            lowContrastEarly(i + m, 3) = 1;
        else
            lowContrastEarly(i + m, 3) = 0;
        end
    end
    for i=1:nChannels
        lowContrastLate(i + m, :) = spikes(i).lowContrastLateParameters;
    end
    for i=1:nChannels
        highContrastEarly(i + m, :) = spikes(i).highContrastEarlyParameters;
    end
else
    lowContrastEarly = zeros(nChannels, 3);
    for i=1:nChannels
        lowContrastEarly(i, 1:2) = spikes(i).lowContrastEarlyParameters;
        if abs(max(spikes(i).highContrastLateSTA)) > abs(min(spikes(i).highContrastLateSTA)) %ON cells
            lowContrastEarly(i, 3) = 1;
        else
            lowContrastEarly(i, 3) = 0;
        end
    end
    for i=1:nChannels
        lowContrastLate(i, :) = spikes(i).lowContrastLateParameters;
    end
    for i=1:nChannels
        highContrastEarly(i, :) = spikes(i).highContrastEarlyParameters;
    end
end
onOffDecider = lowContrastEarly(:,3);
onCells = find(onOffDecider == 1);
offCells = find(onOffDecider == 0);
clear spikes

% use on off information via e.g.
% lowContrastEarlyON = lowContrastEarly(onCells, :);