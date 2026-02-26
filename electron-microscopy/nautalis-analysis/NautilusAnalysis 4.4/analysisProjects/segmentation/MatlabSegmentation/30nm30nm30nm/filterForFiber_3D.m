

SPN = 'D:\LGNs1\testAlign\s8Sample_filt4_2DFiber\';
TPN = 'D:\LGNs1\testAlign\s8Sample_filt4_3DFiber\';
mkdir(TPN)

cbVol = [5 3000000];
ballness = .3;

dSPN = dir([SPN '*.tif']);
inams = {dSPN.name};

imageInfo = imfinfo([SPN inams{1}]);

W = imageInfo.Width;
H = imageInfo.Height;
D = length(inams);

I = zeros(H,W,D,'uint8');
for i = 1:D
   I(:,:,i) = imread([SPN inams{i}])>0;    
end

% SE = strel('ball',5,1);
% oI = imopen(double(I),SE);
oI = double(I);

image(mean(oI,3)*1656)
pause(.1)

%% 
labCB = bwconncomp(oI>.01,6);
tProps = regionprops(labCB,'Area','PixelIdxList');
Areas = [tProps.Area];
goodA = find((Areas>cbVol(1)) & (Areas<cbVol(2)));

CBsize = I*0;
for i = 1:length(goodA)
    CBsize(tProps(goodA(i)).PixelIdxList)= tProps(goodA(i)).Area;% CBsize(tProps(goodA(i)).PixelIdxList)+100;
end

image(CBsize(:,:,3))
labCB = bwlabeln(CBsize,6);

%%  Analyze morphology
%{
labCB = bwconncomp(CBsize,6);

ratCB = double(CB)*0;
for i  = 1:labCB.NumObjects
   [y x z] = ind2sub(size(I),labCB.PixelIdxList{i});
   y = y-mean(y);
   x = x-mean(x);
   z = z-mean(z);
   dists= sqrt(y.^2+x.^2+z.^2);
   rad  = mean(dists);
   eqVol = 4/3 * pi * rad^3;
   ratVol(i) = length(y)/eqVol;
   ratCB(labCB.PixelIdxList{i})  = ratVol(i);
end

ratCB = ratCB<ballness;

image(ratCB(:,:,3)*1000000)

C = zeros(H,W,3,'uint8');
C(:,:,1) = ratCB(:,:,3)*1000;
C(:,:,2) = CBsize(:,:,3) * 10000;

image(C)
labCB = bwlabeln(ratCB,6);
%}

%%

for i = 1:size(labCB,3)
    newName = sprintf('%05.0f.tif',i)
   %imwrite(uint16(labCB(:,:,i)),[TPN newName],'bitdepth',16)
       imwrite(uint8(labCB(:,:,i)),[TPN newName])

end


























