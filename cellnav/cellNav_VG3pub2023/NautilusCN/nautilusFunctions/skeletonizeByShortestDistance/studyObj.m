clear all

SPN = 'C:\Users\joshm\Documents\myWork\LGNs1\jlmHomeSeg\export1Mat\'
TPN = [SPN 'cell1\'];

%% Load object data
load([SPN 'objs.mat'])

%% Pick an object
for i = 1:length(o)
    objSize(i,:) = size(o{i});
end

objNum = find(objSize(:,1) == max(objSize(:,1)));
allVox.name = 'All object voxels';
allVox.subs = o{objNum};
numVox = size(allVox.subs,1);
clear o

%% Find center
midObj = median(allVox.subs,1);


%% turn object into voxel connectivity matrix
tic
disp('Turn object into connectivity matrix')
[allVox.conMat allVox.conDat] = obj2con(allVox.subs);
toc

%% Get obj surface
 [surfVox ] = subs2surf(allVox);
 numSurf = size(surfVox.subs,1);
 [moveObj maxSubs] = squeezeSubs(surfVox.subs);
 surf2cent = sqrt((surfVox.subs(:,1) - midObj(1)).^2 + (surfVox.subs(:,2) - midObj(2)).^2 + ...
     (surfVox.subs(:,3) - midObj(3)).^2);

 
%% Find Seed
dist2mid = sqrt((surfVox.subs(:,1)-midObj(1)).^2 + (surfVox.subs(:,2)-midObj(2)).^2 + ...
    (surfVox.subs(:,3)-midObj(3)).^2);

minDist = min(dist2mid);
seed = find(dist2mid==minDist,1);

%% find shortest path from seed to surface voxels

tic 
disp('Find shortest path from seed to all surfaces')
[seedPath] = conMat2allShortest(surfVox,seed);
toc

showMaxProp(moveObj,mod(seedPath.segID,256));

%% Distance to surface
% reps = 2;
% vProp  = (sum(conMat(:,1:6)>0,2)>5) * size(conMat,1) * 2; %make surface voxel distance 0;
% [passSurfVox passSurfProp passSurfPred] = passDist(conMat,conDat,vProp,reps);

%% Find tips by spreading ditances to seed


tic
reps = 20;
[longTip] = passMaxAndDist(surfVox,seedPath.dist,reps);
longTipness = hist(longTip.vox,1:numSurf);
toc

% reps = 15;
% [maxTipPath] = passMaxAndDist(surfVox,tipness,reps);
% tipness2 = hist(maxTipPath.vox,1:numSurf);


%% Analyze tip shape
voxEdges =  connectRegions(surfVox, longTip.vox);
tipShape = getTipShape(surfVox,voxEdges);

tipProp = (longTipness>100) & (tipShape.edge2cent > 5) & (seedPath.dist>100);
tips = find(tipProp);
colormap jet(256)
showMaxProp(moveObj,tipProp(longTip.vox)*100  + (longTipness>0)*50);



%% Shorten paths
[surfSkel] = shortenPath(surfVox,seedPath,tips);
showMaxProp(moveObj,surfSkel.prop *100  + (longTipness>0)*50);


%% Simplify skel
[skel] = simpleSkel(surfVox,surfSkel,10);
[skel] = bridgeGaps(skel,seedPath);
start = showSubs(surfVox.subs(skel.edges(:,1),:));

showSkel(skel);

%Consolidate and turn to arbor

%% Group voxels to skel
reps = 10;
[surfClose2node] = passMaxAndDist(surfVox,skel.prop,reps);
 [maxIm objVol] = showSubs(surfVox.subs,mod(surfClose2node.vox,256));
image(squeeze(max(objVol,[],1)))


% skelInVol = zeros(size(allVox.subs,1),1);
% skelInVol(surfVox.surf2all) = (skel.prop);
% [volClose2node] = passMaxAndDist(allVox,skelInVol,reps);
% showSubs(allVox.subs,mod(volClose2node.prop,256));



%% Show property
[shortSkel objVol] = showMaxProp(moveObj,skel.prop);
[fullSkel] = showMaxProp(moveObj,surfSkel.prop);

[maxIm objVol] = showMaxProp(moveObj,ones(1,size(moveObj,1)));
sumObj = sum(objVol,3);
sumSurf = sumObj*10 + ((sumObj>0)*20);

col(:,:,3) = sumSurf;
col(:,:,2) = shortSkel*1000;
col(:,:,1) = fullSkel * 1000;

image(uint8(col)), pause(.01)

