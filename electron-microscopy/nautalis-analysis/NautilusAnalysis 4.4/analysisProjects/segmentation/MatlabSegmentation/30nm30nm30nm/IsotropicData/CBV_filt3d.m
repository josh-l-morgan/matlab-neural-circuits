


colormap(gray(256))
SPN = GetMyDir;
TPN = [SPN(1:end-1) '_filtCBV3d_2' '\'];
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

colI = zeros(800,800,length(iNam),3,'uint8');

for i = 1:length(iNam)
   colI(:,:,i,:) = imread([SPN iNam{i}]); 
end
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%%Run for cell bodies

I = colI(:,:,:,2)>0;
G = colI(:,:,:,3);
%% get object properties

labI = bwlabeln(I,6);
iProps = regionprops(labI,G,'BoundingBox','PixelIdxList',...
    'PixelValues','MeanIntensity')

boxes = cat(1,iProps.BoundingBox);
widths= boxes(:,[4 5 6]);
minW = min(widths,[],2);
maxW = max(widths,[],2);
gVal = [iProps.MeanIntensity];

%% select properties

CBmin = 5;
CBmax = 40;
CBpass = find((minW>=CBmin) & (maxW<=CBmax));

Vmin = 4;
Vmax = 40;
Vpass = find((minW>=Vmin)& (maxW >= Vmax));

%% write labels
cbI = I * 0;
for i = 1:length(CBpass)
   cbI(iProps(CBpass(i)).PixelIdxList) = gVal(CBpass(i)); 
end
image(cbI(:,:,300)*1000)

vI = I * 0;
for i = 1:length(Vpass)
   vI(iProps(Vpass(i)).PixelIdxList) = gVal(CBpass(i)); 
end
image(vI(:,:,300)*1000)

%%

colI(:,:,:,2) = cbI*1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%%Run again for vessles

I = colI(:,:,:,1)>0;
G = colI(:,:,:,3);
%% get object properties

labI = bwlabeln(I,6);
iProps = regionprops(labI,G,'BoundingBox','PixelIdxList',...
    'PixelValues','MeanIntensity')

boxes = cat(1,iProps.BoundingBox);
widths= boxes(:,[4 5 6]);
minW = min(widths,[],2);
maxW = max(widths,[],2);
gVal = [iProps.MeanIntensity];

%% select properties

CBmin = 5;
CBmax = 40;
CBpass = find((minW>=CBmin) & (maxW<=CBmax));

Vmin = 4;
Vmax = 40;
Vpass = find((minW>=Vmin)& (maxW >= Vmax));

%% write labels
cbI = I * 0;
for i = 1:length(CBpass)
   cbI(iProps(CBpass(i)).PixelIdxList) = gVal(CBpass(i)); 
end
image(cbI(:,:,300)*1000)

vI = I * 0;
for i = 1:length(Vpass)
   vI(iProps(Vpass(i)).PixelIdxList) = gVal(CBpass(i)); 
end
image(vI(:,:,300)*1000)

%%

colI(:,:,:,1) = vI*1000;

%%
for i = 1:size(colI,3)
    
   image(uint8(squeeze(colI(:,:,i,:)))) ,pause(.1)
end


%%
for i = 1:size(colI,3)
    
   imwrite(uint8(squeeze(colI(:,:,i,:))),[TPN iNam{i}]);
 
end








