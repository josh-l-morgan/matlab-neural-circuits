

clear all
OPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export\'
%TPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\matOut\'
OPN = GetMyDir
TPN = GetMyDir
%SPN = 'D:\LGNs1\segmentation\VAST\Joshm\rayAlign1Synapse\eportSyn\'
%IPN = 'D:\LGNs1\HP_processing\rayAlign01\images\'


dOPN = dir(OPN);
for i = 1:length(dOPN)
    if sum(regexp(dOPN(i).name,'.txt'))
        fileName = [OPN dOPN(i).name];
    end
end
        
%fileName = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\ColorNames.txt'
ids = readVastColors(fileName);


%%

displayName = 'Segment'
foundCell = 0;
displayID = [];
for i = 1:length(ids)
    nam = ids{i,2};
    ID = ids{i,1};
    if ID>0;
    if ~sum(regexp(nam,displayName))
        foundCell = foundCell+1;
        displayID(foundCell) = ID;
    end
    end
    
end

%% find images
SPN = OPN;
dSPN = dir(SPN);
iNam = {};
for i =1:length(dSPN)
    if sum(regexp(dSPN(i).name,'.tif')) |...
        sum(regexp(dSPN(i).name,'.png'));
       iNam{length(iNam)+1,1} = dSPN(i).name;
    end
end

%% colormap
cmap = hsv(length(displayID))
cmap2 = zeros(length(displayID)*2,3);
cmap2(displayID+1,:) = cmap;


planes = length(iNam);
cmap2 = cmap2*300/planes;

imageInfo = imfinfo([SPN iNam{1}]);
ys = imageInfo.Height;
xs = imageInfo.Width;
Ic = zeros(ys,xs,3);
cTemp = zeros(ys,xs);
IcSum = Ic;

for i = 1:planes
    disp(sprintf('reading %d of %d.',i,planes))
        I = imread([SPN iNam{i}]);
        cTemp(:) = cmap2(I+1,1);
        Ic(:,:,1) = cTemp;
         cTemp(:) = cmap2(I+1,2);
        Ic(:,:,2) = cTemp;
         cTemp(:) = cmap2(I+1,3);
        Ic(:,:,3) = cTemp;
        IcSum = IcSum + Ic;
end


image(uint8(IcSum))
imwrite(uint8(IcSum),[TPN 'IcSum.tif'])