

SPN = 'D:\LGNs1\testAlign\s8Sample\'

downSamp = 4;


%Get  tifs
dSPN = dir([SPN '*.tif']);
inams = {dSPN.name};
inams = sort(inams');

nzs = fix(length(inams)/downSamp);


%% build reshape matrix
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
newZ = mod(y-1,downSamp)*downSamp +  mod(x-1,downSamp)+1 + (z-1)* downSamp^2;

sourceInd = sub2ind(size(rawI),y,x,ones(length(y),1)*z);
targetInd = sub2ind(size(I),newY,newX,newZ);

reshapeMat(targetInd) = sourceInd;

end
clear y x z newY newX newZ sourceInd targetInd

%% Read data

for nz = 1:nzs;
    %Get Data
    for stepZ = 1:downSamp
        getZ = (nz-1)*downSamp + stepZ;
        rawI(:,:,stepZ) = imread([SPN inams{getZ}]);
    end
    %Reshape
    I = rawI(reshapeMat);
end
        
medI= median(I,3);







