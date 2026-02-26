

clear all
OPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export\'
%TPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\matOut\'
% OPN = GetMyDir
% TPN = GetMyDir
%SPN = 'D:\LGNs1\segmentation\VAST\Joshm\rayAlign1Synapse\eportSyn\'
%IPN = 'D:\LGNs1\HP_processing\rayAlign01\images\'
%readWindow = {[3679	6616]	[5298	9585]}
%readWindow = {[1 1532] [1 2796]}

OPN = 'D:\LGNs1\segmentation\VAST\Joshm\1-8_cutoutAligned\export\';

dOPN = dir(OPN);
for i = 1:length(dOPN)
    if sum(regexp(dOPN(i).name,'.txt'))
        fileName = [OPN dOPN(i).name];
    end
end
        
%fileName = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\ColorNames.txt'
ids = readVastColors(fileName);


%%


typeLab = [];
for i = 1:length(ids)
    nam = ids{i,2}
    ID = ids{i,1}
    if ID>0;
    if sum(regexp(lower(nam),'background'))
        typeLab(ID) = 0;
    elseif sum(regexp(lower(nam),'rgc'))
        typeLab(ID) = 2;
    elseif sum(regexp(lower(nam),'lin'))
        typeLab(ID) = 2;
    elseif sum(regexp(lower(nam),'tcp'))
        typeLab(ID) = 1;
    elseif sum(regexp(lower(nam),'gli'))
        typeLab(ID) = 1;
    elseif sum(regexp(lower(nam),'cb'))
        typeLab(ID) = 1;
    elseif sum(strcmp(nam,'Segment 1'))
        typeLab(ID) = 1;
    else
        typeLab(ID) = 2;
    end
    end
    
end

%% Use list 
%useStrings = {'cb 1' 'rgc 12'	'rgc 13'	'rgc10'	'rgc 11'}
useStrings = ids(:,2);
%useStrings = {'cb 3'}
useList = typeLab * 0;
for i = 1:length(ids)
    nam = ids{i,2}
    ID = ids{i,1}
   
    if ID>0;
        foundStrings = regexp(lower(nam),lower(useStrings))
    if sum([foundStrings{:}])
       useList(ID) = 1;
    else
        useList(ID) = 0;
    end
    end
    
    
end

%%
% displayName = 'Segment'
% foundCell = 0;
% displayID = [];
% for i = 1:length(ids)
%     nam = ids{i,2};
%     ID = ids{i,1};
%     if ID>0;
%     if ~sum(regexp(nam,displayName))
%         foundCell = foundCell+1;
%         displayID(foundCell) = ID;
%     end
%     end
%     
% end

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
planes = length(iNam);

%% colormap

cmap = josh2col(typeLab);
filterMap = repmat(useList',[1 3]);
filterMap = cat(1,[0 0 0],filterMap);

cmap = cmap.*filterMap;

%%
imageInfo = imfinfo([SPN iNam{1}]);
readWindow{2}(2) = imageInfo.Width;
readWindow{1}(2) = imageInfo.Height;

Ic = zeros(readWindow{1}(2)-readWindow{1}(1)+1,readWindow{2}(2)-readWindow{2}(1)+1,3);


ys = size(Ic,1);
xs = size(Ic,2);
Ic = zeros(ys,xs,3);
cTemp = zeros(ys,xs);
IcSum = Ic;

for i = 1:planes
    disp(sprintf('reading %d of %d.',i,planes))
        I = imread([SPN iNam{i}]);
%         cTemp = cTemp+double(I(readWindow{1}(1):readWindow{1}(2),...
%              readWindow{2}(1):readWindow{2}(2)));
%         
        cTemp(:) = cmap(I(readWindow{1}(1):readWindow{1}(2),...
            readWindow{2}(1):readWindow{2}(2))+1,1);
        Ic(:,:,1) = cTemp;
                cTemp(:) = cmap(I(readWindow{1}(1):readWindow{1}(2),...
            readWindow{2}(1):readWindow{2}(2))+1,2);
        Ic(:,:,2) = cTemp;
                 cTemp(:) = cmap(I(readWindow{1}(1):readWindow{1}(2),...
            readWindow{2}(1):readWindow{2}(2))+1,3);
        Ic(:,:,3) = cTemp;
        IcSum = IcSum + Ic;
end


image(uint8(IcSum*2))
imwrite(uint8(IcSum),[TPN 'IcSumType7.tif'])


