


sourceDir = 'C:\Users\joshm\Documents\myWork\myPublications\LGNs1\figures\movies\rotation\combineSeedAssociation\conCat2\';
targDir = 'C:\Users\joshm\Documents\myWork\myPublications\LGNs1\figures\movies\rotation\combineSeedAssociation\cropI1\';
if ~exist(targDir,'dir'),mkdir(targDir),end

xStart = 363;
xStop = 2058;
yStart = 312;
yStop = 1782;

countI = 0;
iNams = {};
clear readNames writeNames

ims = dir([sourceDir '*.png']);
iNams = cat(1,{ims.name});
    for i = 1:length(iNams);
        countI = countI + 1;
        readNames{countI} = [sourceDir iNams{i}];
        writeNames{countI} = [targDir sprintf('i_%05.0f.png',countI)];
        
end

parfor i = 1:length(readNames)
    I = imread(readNames{i});
     If = I(yStart:yStop,xStart:xStop,:);   
        imwrite(If,writeNames{i});

end



