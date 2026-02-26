

%% find membranes from segmentation data

SPN = GetMyDir
diThresh = 10;


S = readVast(SPN);
S = single(S);

disp('distance transform')
[D idx] = bwdist(~S);
disp('getting dilation')
ero = D>diThresh;


%%
for i = 1:size(ero,3)
image(ero(:,:,i)*10000)
pause
end

%%
TPN = [SPN(1:end-1) '_notMem\'];

if ~exist(TPN,'dir')
    mkdir(TPN)
end

for i = 1:size(ero,3)
   iNam = sprintf('notMem_%05.0f.tif',i);
   imwrite(ero(:,:,i),[TPN iNam],'Compression','none')
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
