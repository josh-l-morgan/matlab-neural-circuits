


clear sourceDir
sourceDir{1} = 'C:\Users\joshm\Documents\myWork\myPublications\LGNs1\figures\movies\rotation\combineSeedAssociation\seed4fullSeed\';
sourceDir{2} = 'C:\Users\joshm\Documents\myWork\myPublications\LGNs1\figures\movies\rotation\combineSeedAssociation\seed4fullRGC_trans\';
sourceDir{3} = 'C:\Users\joshm\Documents\myWork\myPublications\LGNs1\figures\movies\rotation\combineSeedAssociation\seed4fullTCR_saturate\';
targDir = 'C:\Users\joshm\Documents\myWork\myPublications\LGNs1\figures\movies\rotation\combineSeedAssociation\conCat2\';
if ~exist(targDir,'dir'),mkdir(targDir),end
fullSize = [2574 2030];


countI = 0;
iNams = {};
clear readNames writeNames
for d = 1:length(sourceDir)
    ims = dir([sourceDir{d} '*.png']);
    iNams = cat(1,{ims.name});
    for i = 1:length(iNams);
        countI = countI + 1;
        readNames{countI} = [sourceDir{d} iNams{i}];
        writeNames{countI} = [targDir sprintf('i_%05.0f.png',countI)];
        
    end
end

parfor i = 1:length(readNames)
    I = imread(readNames{i});
        
        [ys xs cs] = size(I);
        startY = fix(fullSize(2)-ys)/2;
        startX = fix(fullSize(1)-xs)/2;
        If = zeros(fullSize(2), fullSize(1),3,'uint8');
        If(startY+1:startY+ys,startX+1:startX+xs,:) = I;
        imwrite(If,writeNames{i});

end



