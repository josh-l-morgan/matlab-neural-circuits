

%% find membranes from segmentation data

SPN = GetMyDir
diThresh = 2;


S = readVast(SPN);
P = bwperim(S,6);

pause
%S = single(isoEM(S));

disp('distance transform')
[D idx] = bwdist(S);
disp('getting dilation')
dil = S(idx);
dil(D>diThresh) = 0;


%%
clear S D idx
dif1 = dil(1:end-1,:,:)-dil(2:end,:,:);
ids1 = find(abs(dif1)>0);
[y x z] = ind2sub(size(dif1),ids1);
ids1 = sub2ind(size(dil),y, x, z);
y = y+1;
ids1b = sub2ind(size(dil),y, x, z);

dif1 = dil(:,1:end-1,:)-dil(:,2:end,:);
ids2 = find(abs(dif1)>0);
[y x z] = ind2sub(size(dif1),ids2);
ids2 = sub2ind(size(dil),y, x, z);
y = y+1;
ids2b = sub2ind(size(dil),y, x, z);

dif1 = dil(:,:,1:end-1)-dil(:,:,2:end);
ids3 = find(abs(dif1)>0);
[y x z] = ind2sub(size(dif1),ids3);
ids3 = sub2ind(size(dil),y, x, z);
y = y+1;
ids3b = sub2ind(size(dil),y, x, z);

inds = unique(cat(1,ids1, ids1b, ids2, ids2b, ids3, ids3b));

M = zeros(size(dil),'uint8');
M(inds) = 1;


%%
for i = 1:size(S,3)
image(dil(:,:,i)*10000)
pause
end

%%
TPN = [SPN(1:end-1) '_mem\'];

if ~exist(TPN,'dir')
    mkdir(TPN)
end

for i = 1:size(M,3)
   iNam = sprintf('mem_%05.0f.tif',i);
   imwrite(M(:,:,i),[TPN iNam],'Compression','none')
end
   


%{

%% Get difs
for i = 1:3




dif1 = diff(dil);
dif2 = diff(dil,1);
dif3 = diff(dil,2);





diThresh = 10;


P = bwperim(S,6);


surfD = D<3;
sumSurf = sum(surfD,3);
image(sumSurf*256/max(sumSurf(:)))


%%
for i = 1:size(M,3)
image(M(:,:,i)*10000)
pause
end


%}
