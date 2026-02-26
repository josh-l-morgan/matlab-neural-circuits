

clear all
OPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export\'
TPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\matOut\'
SPN = 'D:\LGNs1\segmentation\VAST\Joshm\rayAlign1Synapse\eportSyn\'
IPN = 'D:\LGNs1\HP_processing\rayAlign01\images\'

fileName = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\ColorNames.txt'
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

allS = uint8(readVast(SPN));

I = uint8(readVast(IPN));



%%


typeLab = [];
for i = 1:length(ids)
    nam = ids{i,2}
    ID = ids{i,1}
    if ID>0;
    if sum(regexp(lower(nam),'background'))
        typeLab(ID) = 0;
    elseif sum(regexp(lower(nam),'rgc'))
        typeLab(ID) = 1;
    elseif sum(regexp(lower(nam),'lin'))
        typeLab(ID) = 2;
    elseif sum(regexp(lower(nam),'tcp'))
        typeLab(ID) = 3;
    elseif sum(regexp(lower(nam),'gli'))
        typeLab(ID) = 4;
    else
        typeLab(ID) = 5;
    end
    end
    
end



typeKey = [' RGC,       Inhibitory,       thalamocortical,       glia,       other'];
[sortType typePos] = sort(typeLab);

cells = double(1:max(allO(:)));


%% color Segs

cmap = hsv(length(cells));
cmap = cat(1,[0 0 0 ],cmap);
typeCol = cmap;
typeSeg = allO*0;
[ys xs zs ] = size(I);
typeI = zeros(ys, xs, 3 ,zs,'uint8');
planeI = typeI(:,:,1,1);
for i = 1:size(typeI,4)
    for c = 1:3
        planeI(:) = typeCol(allO(:,:,i)+1,c)*80;
        %image(planeI*100),pause
        typeI(:,:,c,i) = planeI;
    end
end
%%
imageDir = [TPN 'colorSeg\'];
mkdir(imageDir)
for i = 1:size(typeI,4)
    tempI =(repmat(I(:,:,i),[1 1 3])-80)*2 -  typeI(:,:,:,i);
    image(tempI),pause(.01)
    %imwrite(tempI,[imageDir sprintf('colorType_%04.0f.tif',i)])
end
