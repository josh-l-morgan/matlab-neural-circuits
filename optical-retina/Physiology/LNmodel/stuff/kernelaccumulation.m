%% LOADING FILES AND GETTING NUMBER OF CHANNELS
[fileName,directoryName] = uigetfile;
load([directoryName fileName])
cd('F:\Program Files\MATLAB\R2006a\work')
s = whos('spikes');
nChannels = s.size(2);
timeWindow = (-500:10);


%% ACCUMULATING KERNELS
if exist ('kernelMatrix', 'var')
    [m n] = size(kernelMatrix);
    for i=1:nChannels
        kernel = spikes(i).time.kernel;
        if abs(max(kernel)) > abs(min(kernel))
            kernel = kernel/max(kernel);
        else
            kernel = kernel/abs(min(kernel));
        end
        kernelMatrix(i + m,:) = kernel;
    end
else
    kernelMatrix = zeros(nChannels, length(timeWindow));
    for i=1:nChannels
        kernel = spikes(i).time.kernel;
        if abs(max(kernel)) > abs(min(kernel))
            kernel = kernel/max(kernel);
        else
            kernel = kernel/abs(min(kernel));
        end
        kernelMatrix(i,:) = kernel;
    end
end
kernelMatrix = exciseRows(kernelMatrix);
clear spikes

