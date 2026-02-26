function [averageOut] = slidingwindowaverage (windowSize, rawIn)
%[averageOut] = slidingwindowaverage (windowSize, rawIn)
%This function computes a sliding window average of signals which are
%contained in the columns of the input Matrix rawIn. The output is a matrix
%of the same dimensions (averageOut).
% See also: CALCIUMSIGNAL

[nObservations nVariables] = size(rawIn);
averageOut = zeros(nObservations, nVariables);
halfWindow = floor(windowSize/2);

for i=1:nObservations
    if ((i - halfWindow) < 1) && ((i + halfWindow) <= nObservations)
        averageOut(i,:) = mean(rawIn(1 : windowSize, :));

    elseif ((i - halfWindow) >= 1) && ((i + halfWindow) <= nObservations)
        averageOut(i,:) = mean(rawIn(i - halfWindow : i + halfWindow, :));

    elseif ((i - halfWindow) >= 1) && ((i + halfWindow) > nObservations)
        averageOut(i,:) = mean(rawIn(nObservations - windowSize : nObservations,:));

    else
        averageOut(i,:) = mean(rawIn);
    end
end