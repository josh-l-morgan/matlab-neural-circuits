SPN = 'C:\Users\joshm\Documents\myWork\myPublications\LGNs1\figures\draft5_revision\pngBig\';
%SPN = 'C:\Users\joshm\Documents\myWork\myPublications\LGNs1\figures\draft4\ai\figResources\';

TPN = [SPN(1:end-1) 'Small\'];
if ~exist(TPN),mkdir(TPN),end;


maxPix = 100000000;

sourceDir = dir([SPN '*.png']);
iNams = {sourceDir.name};

for i = 1:length(iNams)
    
    I = imread([SPN iNams{i}]);
    Imin = min(I,[],3);
    [y x] = find(Imin<255);
    Icrop = I(min(y):max(y),min(x):max(x),:);
    [ys xs zs] = size(Icrop);
    pixNum = ys*xs;
    scalePix = maxPix/ pixNum;
    if scalePix < 1
       Ire = imresize(Icrop,sqrt(scalePix)) ;
    else
        Ire = Icrop;
    end
    
    imwrite(Ire,[TPN iNams{i}]);
    
    
end