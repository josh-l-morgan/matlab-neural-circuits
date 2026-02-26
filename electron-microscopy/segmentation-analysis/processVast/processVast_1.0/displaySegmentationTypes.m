

clear all
OPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export\'
TPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\matOut\'
SPN = 'D:\LGNs1\segmentation\VAST\Joshm\rayAlign1Synapse\exportSyn\'
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


[typeLab typeKey] = getTypes(ids);
[sortType typePos] = sort(typeLab);

cells = double(1:max(allO(:)));


%% color types
typeCol = [0 0 0; 0 1 0 ; 1 0 0; 0 0 1; .75 .75 0; 0 .75 .75]
typeCol = [0 0 0; 0 1 0 ; 1 0 0; 0 0 1; .75 .75 0; 0 .75 .75]
typeCol = [0 0 0; 1 0 1; 0 1 1; 1 1 0; 0 0 5; 1 0 0]; 

typeCol = [0 0 0; 0 1 0 ; 1 0 0; 0 0 1; .75 .75 0; .25 0 .25; .25 0 .25;  .25 0 .25;.25 0 .25;]
typeCol = 1-typeCol;

typeSeg = allO*0;
typeSeg(allO>0) = typeLab(allO(allO>0));
[ys xs zs ] = size(I);
typeI = zeros(ys, xs, 3 ,zs,'uint8');
planeI = typeI(:,:,1,1);
synI = typeI;
SE = strel('disk',10);
for i = 1:size(typeI,4)
    sprintf('drawing %d of %d',i,size(typeI,4))
    for c = 1:3
        planeI(:) = typeCol(typeSeg(:,:,i)+1,c);
        planeS = allS(:,:,i);
        perimS = bwperim(planeS);
        perimS = imdilate(perimS,SE);
        if c == 1;
           
             planeS = planeS-uint8(perimS);
        else
             planeS = perimS*1000;
           
            
        end
        %image(planeI*100),pause
        typeI(:,:,c,i) = planeI;
        synI(:,:,c,i) = planeS;
    end
end


synI2 = synI;

synI2(:,:,2,:) = synI(:,:,1,:)*1000; 
synI2(:,:,3,:) = synI(:,:,1,:)*1000; 
synI2(:,:,1,:) = synI(:,:,2,:); 




%%
imageDir = [TPN 'colorTypeSyn2\'];
mkdir(imageDir)
countFrame = 0;
for i = 1:size(typeI,4)
    countFrame = countFrame + 1;
    tempI =((repmat(double(I(:,:,i)),[1 1 3])-60)*1.6 -  double(typeI(:,:,:,i))*40);
    %(repmat(double(I(:,:,i)),[1 1 3])-40) +  double(typeI(:,:,:,i))*40;
    image(uint8(tempI)),pause(.01)
    imwrite(uint8(tempI),[imageDir sprintf('colorType_%04.0f.tif',countFrame)])
end
for i = size(typeI,4)-1:-1:1
    countFrame = countFrame + 1;
    tempI =(repmat(double(I(:,:,i)),[1 1 3])-60)*1.6 -  double(typeI(:,:,:,i))*40 + double(synI2(:,:,:,i));
    image(uint8(tempI)),pause(.01)
    imwrite(uint8(tempI),[imageDir sprintf('colorType_%04.0f.tif',countFrame)])
end

