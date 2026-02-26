

colormap(gray(256))
SPN = GetMyDir;
TPN = [SPN(1:end-1) '_filtBlood3' '\'];
mkdir(TPN)
%% read names
dSPN = dir(SPN);
iNam = {};
for i  =   1:length(dSPN)
    nam = dSPN(i).name;
    if sum(regexp(nam,'.tif'))
        iNam{length(iNam)+1} = nam;
    end
end

%%

minArea = 100;
maxArea = 10000;

parfor i= 1:length(iNam)
    %%
    disp(sprintf('processing image %d of %d',i,length(iNam)))
    
    I = double(imread([SPN iNam{i}]));
    %I = double(imresize(I,0.25,'bicubic'));
    I= I - min(I(:));
    I = I * 256/max(I(:));
    
    colI = zeros(size(I,1),size(I,2),3,'uint8');
    colI(:,:,3) = I;
    
    useVals = double(I((I>0) & (I<200)));
    modeI = mode(useVals);
    I(I<modeI) = modeI;
    sdev = std(useVals);
    
    %% Find vessles
    vtI = I > (220);
    vProps = regionprops(vtI,'Area','PixelIdxList');
    vArea = [vProps.Area];
    goodV = find((vArea<10000)&(vArea>4));
    
    %find outside
    outV = find(vArea>10000);
    outI = I * 0;
    for v = 1:length(outV)
        outI(vProps(outV(v)).PixelIdxList) = 1;
    end
    
%     %fill in outsinde
%     oProps = regionprops(outI,'BoundingBox','FilledImage');
%     outI2 = outI;
%     for o = 1:length(oProps)
%         bb = round(oProps(o).BoundingBox);
%         outI(1:bb(4),1:bb(3)) = oProps(o).FilledImage;
%     end
%     
    outI = outI>0;
    vI = vtI*0;
    for g = 1:length(goodV)
        vI(vProps(goodV(g)).PixelIdxList)=1;
    end
    
    vI(outI) = 0;
    colI(:,:,1) = vI*10000;
    image(vI*1000)
    %% find cell bodies
    
    fI = mexHatV2(I,1,10);
    image(fI*10)
    fdev= std(fI(:))
    tI = fI>mode(fI(:))+fdev*1;
    
    
    
%     SE = strel('disk',2);
%     cI = imdilate(vI,SE);
%     tI(cI>0) = 0;
%     cI = imerode(cI,SE);


    image(tI*100)


    %%
    tProps = regionprops(tI,'Area','PixelIdxList');
    Areas = [tProps.Area];
    goodArea = find((Areas>=50) & (Areas<=300));
    
    goodI = tI*0;
    for g = 1:length(goodArea)
        goodI(tProps(goodArea(g)).PixelIdxList)=1;
    end
    image(goodI * 1000)
    
    colI(:,:,2) = goodI*1000;
    %image(colI),pause(.1)
    %%
    %
    %   % Refine cell bodies
    %
    %  cbProps = regionprops(goodI>0,'Area','PixelIdxList',...
    %      'BoundingBox','ConvexImage');
    %
    %  cbI= I * 0;
    %  for c = 1:length(cbProps)
    %     bb = round(cbProps(c).BoundingBox)
    %      cbI(bb(2):bb(4),bb(1):bb(3)) = cbProps(c).ConvexImage;
    %
    %
    %  end
    
    
    
    %%
    
    image(colI),pause(.01)
    imwrite(colI,[TPN iNam{i}]);
end

