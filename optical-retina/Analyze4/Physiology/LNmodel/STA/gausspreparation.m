function [z beta xy] = gausspreparation (STA)
%[z beta xy] = gausspreparation (STA)
% this computes the input to be used in 2D Gaussian fitting of spatial STA
% maps. It takes the STA matrix calculated byt the receptive field script
% as an input. The output is z a vector of intesity values to be fitted a
% matrix xy which contains the x and y coordinates of these z values and
% beta a row vector of initial parameters for the 2D Gaussian fit.
% See also: MASTERGAUSS, GAUSSIAN2D

%% FIRST-PIXEL ANOMALY CORRECTION
[pixel timeWindow] = size(STA);
meanGray = 0.5;
STA(1,:) = meanGray*ones(1,timeWindow);

%% PEAK FRAME NORMALIZATION
peakTime = find (std(STA) == max(std(STA)));
peakFrame = mean(STA(:,peakTime-1:peakTime+1),2);
m = min(peakFrame);
M = max(peakFrame);

if ( mean(STA(:))- min(STA(:))  <  max(STA(:))- mean(STA(:)) ) %for ON cells
    peakColumn = (peakFrame - m) / (M - m);
    disp ('ON')

elseif ( mean(STA(:))- min(STA(:))  >  max(STA(:))- mean(STA(:)) )%for OFF cells
    peakColumn = 1-((peakFrame - m) / (M - m));
    disp ('OFF')
else
    disp('ERROR')
end

peakMatrix = reshape(peakColumn, 80, 60)';

%% OUTPUT COMPUTATION AND ASSIGNMENT
% beta computation
s = size(peakMatrix);
IND = (1:pixel);
[yCoordinate xCoordinate] = ind2sub(s,IND);
[yOffset xOffset] = find(peakMatrix == max(peakMatrix(:)));
background = mean(peakMatrix(:));

% beta assignment
amplitude =1;
theta = pi;
sigmaX = 4;
sigmaY = 4;

% output
z = peakMatrix(:);
beta = [amplitude xOffset yOffset theta sigmaX sigmaY background];
xy = [xCoordinate' yCoordinate'];
