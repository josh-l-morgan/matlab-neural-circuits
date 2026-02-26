


%%Make stack of all object Sums.

colormap gray(256)

clear all
OPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export\'
TPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\matout\sumStack\'
%OPN = GetMyDir
%TPN = GetMyDir
SPN = 'D:\LGNs1\segmentation\VAST\Joshm\rayAlign1Synapse\eportSyn\'
%IPN = 'D:\LGNs1\HP_processing\rayAlign01\images\'


OPN = 'D:\LGNs1\segmentation\VAST\Joshm\1-8_cutoutAligned\export\'
TPN =  'D:\LGNs1\segmentation\VAST\Joshm\1-8_cutoutAligned\export_sum\'

if ~exist(TPN,'dir'),mkdir(TPN),end

aspectRat = 30/4;

dOPN = dir([OPN '*.txt']);
fileName = [OPN dOPN(1).name];


%fileName = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\ColorNames.txt'
ids = readVastColors(fileName);



allO = readVast(OPN);
obNum = max(allO(:));
if obNum<256
    allO = uint8(allO);
elseif obNum<(2^16)
    allO = uint16(allO);
else
    allO = single(allO);
end



%allS = uint8(readVast(SPN));



%% Draw sums
omax = max(allO(:));
idNums = cell2mat(ids(:,1));

[y x z] = size(allO);
newZ = z * aspectRat;
'Starting image extraction'
numBadName = 0;
parfor o = 1:omax
    
    targ = find(idNums==o);
    newName = sprintf('%sid%04.0f_%s.tif',TPN,o,ids{targ,2});
    if ~exist(newName,'file')
        
        
        sum1 = squeeze(sum(allO==o,1))';
        sum1 = imresize(sum1,[newZ x],'nearest');
        sum2 = squeeze(sum(allO==o,2));
        sum2 = imresize(sum2,[y newZ],'nearest');
        sum3 = sum(allO==o,3);
        
        
        %%tweak histo
        gammaCor = .7;
        minVal = 40;
        sum3 = sum3.^gammaCor;
        sum3 = sum3 * (256-minVal)/max(sum3(:));
        sum3(sum3>0) = sum3(sum3>0)+minVal;
        
        sum2 = sum2.^gammaCor;
        sum2 = sum2 * (256-minVal)/max(sum2(:));
        sum2(sum2>0) = sum2(sum2>0)+minVal;
        
        sum1 = sum1.^gammaCor;
        sum1 = sum1 * (256-minVal)/max(sum1(:));
        sum1(sum1>0) = sum1(sum1>0)+minVal;
        
        allSum = sum3;
        allSum(1:size(sum2,1),size(sum3,2)+1:size(sum3,2)+size(sum2,2))= sum2;
        allSum(size(sum3,1)+1:size(sum3,1)+size(sum1,1),1:size(sum1,2)) = sum1;
        allSum = uint8(allSum);
        
        %%write file
        
        numDone = length(dir([TPN '*.tif']))
        
        try
            imwrite(allSum,newName)
        catch err
            numBadName = numDone+fix(rand*100);
            badName = sprintf('%s%id04.0f_badName.tif',TPN,o);
            imwrite(allSum,badName)
            
        end
        
        disp(sprintf('Finished object %d of %d',numDone,omax))
        image(allSum)
        
        pause(.01)
    end
    
end


















