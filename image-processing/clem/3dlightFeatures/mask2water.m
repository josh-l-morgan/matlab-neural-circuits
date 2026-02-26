function[wI cI] = mask2water(I);


image(sum(I,3)*10)

voxelSize = [1 1 1];
[ys xs zs] = size(I);

cI = -I;
%% find buotns by watershed
%'watershed'

%smooth I
gKern = gaus3d([7 7 7],1,[1 1 1]);
cI = fastCon(I,gKern);
cI = cI>.5;
image(sum(cI,3)*10)



%find distanct to surrounding pixels
pI = bwperim(~I);
pI(:,:,1) = 0; pI(:,:,zs) = 0;
pI(:,1,:) = 0; pI(:,xs,:) = 0;
pI(1,:,:) = 0; pI(ys,:,:) = 0;
dI =bwdistsc(~I & ~pI,voxelSize);
image(sum(pI,3)*10);
image(sum(dI,3));

%%smooth distance transform
gKern = gaus3d([15 15 15],1,[1 1 1]);
cI = fastCon(dI,gKern);
image(sum(cI,3))

%%watershed
dI2 = max(cI(:))-cI;
dI2(~I) = -Inf;%min(dI2(:));\
%cI = imhmin(cI,.1);
dI2(pI) = max(dI2(:)); %paint surround highest
wI = watershed(dI2,26);
image(sum(-dI2,3))
image(max(wI,[],3))

%%Make background 1
wHist = hist(double(wI(:)),0:1:max(double(wI(:))));
wHist = wHist(2:end);
numWI = max(wI(:));
bg = find(wHist == max(wHist));
if bg ~=1
    wI(wI==1) = numWI +1;
    wI(wI==bg) = 1;
    wI(wI == (numWI+1)) = bg;
end


% for i = 1:size(wI,3)
%     imshow(label2rgb(wI(:,:,i),'jet','w')),pause
% end
