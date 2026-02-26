

SPN = 'E:\IxQ_KarlsRetinaVG3_2019\Volumes\HighResSCM\Merge\';

load([SPN 'obI.mat']);
load([SPN 'dsObj.mat']);
load([SPN 'vastSubs.mat']);


%% connect to VAST

clear global vdata
if ~exist('vdata','var')
    vasttools
end

global vdata

if ~vdata.state.isconnected
    vdata.vast.connect('127.0.0.1',22081,1000)
end

%%
laytype = 2;
segName = 'mergeSeg';
vdata.vast.addnewlayer(laytype,segName);

%% reformat vastSubs to list [y x z i]

numVox = 0;
lastVox = 0;
voxList = zeros(10000,4);
for i = 1:length(vastSubs)
    sub = vastSubs{i};
    newNum = size(sub,1);
    numVox = numVox + newNum;
    if numVox > size(voxList,1);
       voxList(numVox+newNum *100,:)=0;
    end 
    voxList(lastVox+1:numVox,1:3) = sub;
    voxList(lastVox+1:numVox,4) = i;
    lastVox = numVox;     
end

voxList = voxList(1:numVox,:);
ySize = max(voxList(:,1));
xSize = max(voxList(:,2));
voxXY = sub2ind([ySize xSize],voxList(:,1), voxList(:,2));
voxID = voxList(:,4);
voxZ = voxList(:,3);

%% Write colors

cS = obI.colStruc;
    [id, res] = vdata.vast.addsegment(0,0,'mergeBlank');

for i = 1:length(cS.names)
    [id, res] = vdata.vast.addsegment(i,0,cS.names{i});
end


%% Write sections

zSize = max(voxList(:,3));

for i = 300:400%min(voxZ) : zSize
    isSec = find(voxZ == i);
    if length(isSec)
        X = voxList(isSec,2);
        Y = voxList(isSec,1);
        
        minx = min(X); miny = min(Y); minz = i;
        maxx = max(X); maxy = max(Y); maxz = i;
        
        X = X - minx+1;
        Y = Y - miny+1;
        I = zeros(max(Y), max(X));
        voxInd = sub2ind(size(I), Y, X);
        I(voxInd) = voxID(isSec)+1;
        
        res = vdata.vast.setsegimageraw(obI.em.mipLevel, minx,...
            maxx, miny, maxy, minz, maxz, I);
    end
    
end














