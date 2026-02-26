

colormap(gray(256))
SPN = 'D:\LGNs1\testAlign\s8Sample_filt4\';
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

vessleRange = [200 100000];
minArea = 100;
maxArea = 10000;
cbFilt = [2 20];
cbThresh = .2; %number of standard deviations away from mode
cbArea = [500 10000];

%%
for i= 1:length(iNam)
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
    vtI = I > (240);
    vProps = regionprops(vtI,'Area','PixelIdxList');
    vArea = [vProps.Area];
    goodV = find((vArea<vessleRange(2))&(vArea>vessleRange(1)));
    
    %find outside
    outV = find(vArea>vessleRange(2));
    outI = I * 0;
    for v = 1:length(outV)
        outI(vProps(outV(v)).PixelIdxList) = 1;
    end
    %%
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
    
    fI = mexHatV2(I,cbFilt(1),cbFilt(2));
    image(fI*10)
    fdev= std(fI(:))
    tI = fI>mode(fI(:))+fdev*cbThresh;
    
%     SE = strel('disk',2);
%     cI = imdilate(vI,SE);
%     tI(cI>0) = 0;
%     cI = imerode(cI,SE);

    image(tI*100)
    
    %% OPEN
    SE = strel('disk',3);
    oI = imerode(tI,SE);
    oI = imdilate(oI,SE);
    image((tI + oI)* 100);
    

    %%
    tProps = regionprops(oI,'Area','PixelIdxList');
    Areas = [tProps.Area];
    goodArea = find((Areas>=cbArea(1)) & (Areas<=cbArea(2)));
    
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
    
    CB(:,:,i) = colI(:,:,2)-colI(:,:,1);
    
end

