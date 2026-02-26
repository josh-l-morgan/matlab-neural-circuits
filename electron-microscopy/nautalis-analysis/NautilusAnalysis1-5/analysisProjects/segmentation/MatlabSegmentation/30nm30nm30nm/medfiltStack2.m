

SPN = 'D:\LGNs1\testS8segmentation\s8Sample\'
FPN = 'D:\LGNs1\testS8segmentation\s8Sample_filt4test\'

downSamp = 4;


if ~exist(FPN,'dir'),mkdir(FPN);end

%Get  tifs
dSPN = dir([SPN '*.tif']);
inams = {dSPN.name};
inams = sort(inams');
D = fix(length(inams)/downSamp);

%% build reshape matrix
%%For every downSamp number of planes, a new matrix is created where cubes
%%are linear. 

imageInfo = imfinfo([SPN inams{1}]);
W = imageInfo.Width;
H = imageInfo.Height;
rawI = ones(H,W,downSamp);
%[y x z] = find(rawI);

nys = fix(H/downSamp);
nxs = fix(W/downSamp);
nzs = downSamp^3;

reshapeMat = zeros(nys, nxs, nzs);

for z = 1:size(rawI,3)
z
[y x] = find(rawI(:,:,z)==1);

newY = fix((y-1)/downSamp)+1;
newX = fix((x-1)/downSamp)+1;
%newZ = mod(y-1,downSamp)*downSamp +  mod(x-1,downSamp)+1 + (z-1)* downSamp^2;
newZ = mod(y-1,downSamp)*downSamp +  mod(x-1,downSamp) * downSamp^2 + (z);

sourceInd = sub2ind(size(rawI),y,x,ones(length(y),1)*z);
targetInd = sub2ind(size(reshapeMat),newY,newX,newZ);

reshapeMat(targetInd) = sourceInd;

end
clear y x z newY newX newZ sourceInd targetInd

%% Read data

for nz = 1:D;
    disp(sprintf('running %d of %d',nz,D))
    %Get Data
    tic
        disp('read data')

    parfor stepZ = 1:downSamp
        getZ = (nz-1)*downSamp + stepZ;
        rawI(:,:,stepZ) = imread([SPN inams{getZ}]);
    end
    toc
    %Reshape
    tic
    'reshape data'
    I = rawI(reshapeMat);
    toc
    
    %%first median
    %medI1 = zeros(nys,nxs,nzs/downSamp);
    tic
    disp('non linear filter')
    parfor m = 1:downSamp^2
        startMed = (m-1)*downSamp;
        Isamp = I(:,:,startMed+1:startMed+downSamp);
        %medI1(:,:,m) = median(Isamp,3);
        maxI1(:,:,m)=max(Isamp,[],3);
    end
    toc
    
    tic
    disp('linear filter')
    %medI2 = median(medI1,3);
    meanI = mean(maxI1,3);
    toc
    
    tic
    disp('writing')
    %filtI(:,:,nz) = meanI;
    newName = sprintf('filt%05.0f.tif',nz);
    imwrite(uint8(meanI),[FPN newName])
    toc
    
end
        
%medI= median(I,3);


% 
% for nz = 1:size(filtI,3)
%      newName = sprintf('filt%05.0f.tif',nz);
%     imwrite(filtI(:,:,nz),[FPN newName])
% end
% 
% 

% %%
% while 1
% for i = 1:size(filtI,3)
%     image(filtI(1500:1600,1500:1600,i))
%     pause
% end
% end








