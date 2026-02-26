%% LOADING FILES AND GETTING NUMBER OF CHANNELS
[fileName,directoryName] = uigetfile;
load([directoryName fileName])
cd('F:\Program Files\MATLAB\R2006a\work')

if exist ('kernelMatrix', 'var')
    [nChannels m] = size(kernelMatrix);
else
    s = whos('spikes');
    nChannels = s.size(2);
end

timeWindow = (-500:10);
figure(2)
hold on
colorBank = {'r', 'g', 'b', 'c', 'm', 'y', 'k'};
colorSelector = randperm(nChannels);

if exist ('kernelMatrix', 'var')
    for i=1:nChannels
        plot ( timeWindow, kernelMatrix(i,:), colorBank{(mod(colorSelector(i),length(colorBank))+1)} )
    end
else
    for i=1:nChannels
        kernel = spikes(i).time.kernel;
        if abs(max(kernel)) > abs(min(kernel))
            kernel = kernel/max(kernel);
        else
            kernel = kernel/abs(min(kernel));
        end
        plot ( timeWindow, kernel, colorBank{(mod(colorSelector(i),length(colorBank))+1)} )
    end
end

axis tight
xlabel('Time to Spike (ms)')
ylabel('Relative STA Intensity')
title(fileName)
