

%SPN = GetMyDir
SPN = 'F:\gxI_13-15\Aligned Stack\';
dSPN = dir([SPN '*.tif']);
inams = cat(1,dSPN.name);
TPN = [SPN(1:end-1) '_crop'];
for i = 1
    I = imread
