
clear all
res1 = 1;
res2 = 4;
id1 = 18;
id2 = 1;

SPN1 = 'D:\LGN\P32_highRes\SynapseSegChas\'
SPN2 = 'D:\LGN\P32_highRes\SynapseSegQuarterRes\'
TPN1 = 'D:\LGN\P32_highRes\SynapseSegQuarterResUpsample\'
TPNmissed1 = 'D:\LGN\P32_highRes\SynapseSegQuarterResUpsampleMissed1\'
TPNmissed2 = 'D:\LGN\P32_highRes\SynapseSegQuarterResUpsampleMissed2\'
if ~exist(TPN1,'dir'), mkdir(TPN1), end
if ~exist(TPNmissed1,'dir'), mkdir(TPNmissed1), end
if ~exist(TPNmissed2,'dir'), mkdir(TPNmissed2), end

SPN1dir = dir([SPN1 '*.png']);
SPN2dir = dir([SPN2 '*.png']);

inam1 = {SPN1dir.name};
inam2 = {SPN2dir.name};

I = imread([SPN1 inam1{1}]); 
seg1 = zeros(size(I,1),size(I,2),length(inam1),'logical');
for i = 1:length(inam1)
   I = imread([SPN1 inam1{i}]); 
   seg1(:,:,i) = I == id1; 
end
seg1 = seg1(1:1980,1:1980,1:100);


I = imread([SPN2 inam2{1}]); 
seg2 = zeros(size(I,1)*res2,size(I,2)*res2,length(inam2),'logical');
for i = 1:length(inam2)
   I = imread([SPN2 inam2{i}]); 
   I = imresize(I,4);
   seg2(:,:,i) = I == id2; 
end
seg2 = seg2(1:1980,1:1980,1:100);

overlap = seg1 & seg2;


sumCol = cat(3,sum(seg1,3),sum(seg2,3),sum(overlap,3));

image(uint8(sumCol*50))
% 
% for i = 1:size(overlap,3)
%     col(:,:,1) = seg1(:,:,i);
%     col(:,:,2) = seg2(:,:,i);
%     col(:,:,3) = overlap(:,:,i);
%     image(uint8(col*200))
%     pause
% end

%% segment segmentations

subplot(2,1,1)
se = strel('sphere',5);
seg1Di = imdilate(seg1,se);
image(squeeze(sum(seg1Di,2)))
[lab1 lab1num] = bwlabeln(seg1Di,26);
colormap colorcube(256)
image(mod(squeeze(max(lab1,[],3)),255)+1)

props1 = regionprops(lab1,'PixelIdxList');

subplot(2,1,2)
se = strel('sphere',5);
seg2Di = imdilate(seg2,se);
image(squeeze(sum(seg2Di,2)))
[lab2 lab2num] = bwlabeln(seg2Di,26);
colormap colorcube(256)
image(mod(squeeze(max(lab2,[],3)),255)+1)

props2 = regionprops(lab2,'PixelIdxList');

%% count overlap
minPix = 10000;
minPercent = 3;

sumOv1 = zeros(lab1num,1);
size1 = zeros(lab1num,1);
missed1 = overlap * 0;
for i = 1:lab1num
    %sprintf('%d of %d',i,lab1num)
   ovVal = seg2(props1(i).PixelIdxList);
   size1(i) = length(props1(i).PixelIdxList);
   sumOv1(i) = sum(ovVal);
   if ~sum(ovVal)
       missed1(props1(i).PixelIdxList) = 1;
   end
   
end
%sum(sumOv1>minPix)
lab1num

found1 = mean(((sumOv1./size1)*100)>minPercent)
%found1 = mean(sumOv1>minPix)



sumOv2 = zeros(lab2num,1);
size2 = zeros(lab2num,1);

missed2 = overlap * 0;
for i = 1:lab2num
    %sprintf('%d of %d',i,lab2num)
   ovVal = seg1(props2(i).PixelIdxList);
   size2(i) = length(props2(i).PixelIdxList);

   sumOv2(i) = sum(ovVal);
   if ~sum(ovVal)
       missed2(props2(i).PixelIdxList) = 1;
   end
   
end
%sum(sumOv2>minPix)
lab2num

found2 = mean(((sumOv2./size2)*100)>minPercent)
%found2 = mean(sumOv2>minPix)


seg2R = permute(seg2,[2 1 3]);
sumOvR = zeros(lab1num,1);
missed1 = overlap * 0;
for i = 1:lab1num
    %sprintf('%d of %d',i,lab1num)
   ovVal = seg2R(props1(i).PixelIdxList);
   sumOvR(i) = sum(ovVal);
   if ~sum(ovVal)
       missed1(props1(i).PixelIdxList) = 1;
   end
   
end
%sum(sumOv1>minPix)

foundR = mean(((sumOvR./size1)*100)>minPercent)
%foundR = mean(sumOv1>minPix)

colSeg = cat(3,sum(seg1,3),sum(seg2,3),sum(seg2R,3));
image(uint8(colSeg*10));

return
%% Write seg file
if 0
for i = 1:size(seg2,3)
   I =  seg2(:,:,i);
   I(1:100,100:200) = 1;
   I = I';
   I = flipud(I);
   I = fliplr(I);
   imwrite(uint8(I),[TPN1 inam2{i}]); 
   
   I =  missed1(:,:,i);
      I(1:100,100:200) = 1;
% 
%    I = I';
%    I = flipud(I);
%    I = fliplr(I);
   imwrite(uint8(I),[TPNmissed1 inam2{i}]); 
   
    I =  missed2(:,:,i);
   I(1:100,100:200) = 1;
% 
%     I = I';
%    I = flipud(I);
%    I = fliplr(I);
   imwrite(uint8(I),[TPNmissed2 inam2{i}]); 
end
end











