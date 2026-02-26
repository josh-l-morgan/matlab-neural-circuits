%function[] = connectSubs(objectSubs)
%%
longTipReps = 20; %connectivity length of tips
minTipEccentricity = 100;
interNodeSpacing = 10; % for skeleton simplification


%%
allVox.subs = objectSubs;
allVox.name = 'All object voxels';
numVox = size(allVox.subs,1);

%% Find Dimensions
midObj = median(allVox.subs,1);
allVox.mid = midObj;

minAV = min(allVox.subs,[],1);
offsetSubs = minAV-1;
maxAV = max(allVox.subs,[],1);
imageDims = maxAV-offsetSubs + 1;

minMax{1} = offsetSubs;
minMax{2} = imageDims;

allVox.minMax = minMax;



%% turn object into voxel connectivity matrix

disp('Turn object into connectivity matrix')
[allVox.conMat allVox.conDat] = obj2con(allVox.subs);


%% Get obj surface
disp('Get object surface')
 [surfVox ] = subs2surf(allVox);
 numSurf = size(surfVox.subs,1);
 
for i = 1:3
   moveObj(:,i) = surfVox.subs(:,i)- minMax{1}(i); 
end
 surf2cent = sqrt((surfVox.subs(:,1) - midObj(1)).^2 + (surfVox.subs(:,2) - midObj(2)).^2 + ...
     (surfVox.subs(:,3) - midObj(3)).^2);
 
 surfVox.dist2cent = surf2cent;
 surfVox.minMax = minMax;
 
 
%% Find Seed
dist2mid = sqrt((surfVox.subs(:,1)-midObj(1)).^2 + (surfVox.subs(:,2)-midObj(2)).^2 + ...
    (surfVox.subs(:,3)-midObj(3)).^2);

minDist = min(dist2mid);
seed = find(dist2mid==minDist,1);
surfVox.firstSeed = seed;

%% find shortest path from seed to surface voxels

disp('Find shortest path from seed to all surfaces')
[seedPath] = conMat2allShortest(surfVox,seed);

surfVox.seedPath = seedPath;
showMaxProp(moveObj,mod(seedPath.segID,256));
showMaxProp(moveObj,mod(seedPath.segID,256));


%% Find tips by spreading distances to seed

tic
disp('Find long tips')
[longTip] = countMaxAndDist(surfVox,seedPath.dist,longTipReps);
longTipness = hist(longTip.vox,1:numSurf);
longTipness = max(longTip.countVox(3:end,:),[],1);
toc

surfVox.longTip = longTip;

%% Analyze tip shape
voxEdges =  connectRegions(surfVox, longTip.vox);
tipShape = getTipShape(surfVox,voxEdges);

% tipShape = checkTipness(surfVox,longTip);
% [maxTipIm] = showMaxProp(moveObj,tipShape.edge2cent);
% [maxTipIm2] = showMaxProp(moveObj,tipShape.edge2cent(longTip.vox));
% image(maxTipIm+maxTipIm2)
% SE = strel('disk',10);
% dilateTipIm = imdilate(maxTipIm,SE);
% image(dilateTipIm)


tipProp = (longTipness>longTipReps*5) & (tipShape.edge2cent > (longTipReps/10)) & (surfVox.dist2cent'> minTipEccentricity);
showMaxProp(moveObj,tipProp(longTip.vox)*100  + (longTipness>0)*50);
tips = find(tipProp>0);

surfVox.longTip.voxEdges = voxEdges;
surfVox.longTip.tipShape = tipShape;
surfVox.longTip.tips = tips;
surfVox.longTip.tipProp = tipProp;
surfVox.longTip.longTipReps = longTipReps;
surfVox.longTip.minTipEccentricity = minTipEccentricity;

%% Shorten paths
disp('consolidating paths')
[surfSkel] = shortenPath(surfVox);
[imax ovol] = showMaxProp(moveObj,surfSkel.prop *100  + (longTipness>0)*50);
surfVox.surfSkel = surfSkel;

%% Simplify skel
disp('reducing number of skeleton nodes')
[skel] = simpleSkel(surfVox,interNodeSpacing);
[bridges] = bridgeGaps(skel,seedPath);
skel.bridge = bridges;

skel.minMax = surfVox.minMax;

skelIm = showSkel(skel,1);

%Consolidate and turn to arbor

%% Group voxels to skel
reps = 10;
[surfClose2node] = passMaxAndDist(surfVox,skel.prop,reps);
surfVox.surfClose2node = surfClose2node;

%  [maxIm objVol] = showSubs(surfVox.subs,mod(surfClose2node.vox,256));
% image(squeeze(max(objVol,[],1)))

%% Create cell images
clear skelIm
for i = 1:3
    skelIm{i} = showSkel(skel,i);
end
sumSurf = showSurf(surfVox);

for i = 1:3
    tempSkel = skelIm{i};
    tempSkel(:,:,3) = sumSurf{i} * 20;
%     filled = repmat(sum(tempSkel,3),[ 1 1 3]);
%     tempSurf = repmat(sumSurf{i},[1 1 3])*10;
%     tempSkel(filled==0) = tempSurf(filled == 0);
    image(uint8(tempSkel)),pause(1)
    combCell{i} = uint8(tempSkel);
end

image(combCell{1})

cellStruct.allVox = allVox;
cellStruct.surfVox = surfVox;
cellStruct.skel = skel;
cellStruct.sideViews = combCell;
cellStruct.date = date;

%% Unused analysis
%{

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


%% Distance to surface
% reps = 2;
% vProp  = (sum(conMat(:,1:6)>0,2)>5) * size(conMat,1) * 2; %make surface voxel distance 0;
% [passSurfVox passSurfProp passSurfPred] = passDist(conMat,conDat,vProp,reps);



%% Find tip to tipness
% reps = 15;
% [maxTipPath] = passMaxAndDist(surfVox,tipness,reps);
% tipness2 = hist(maxTipPath.vox,1:numSurf);


%% Find long tips
tic
disp('Find long tips')
[longTip] = passMaxAndDist(surfVox,seedPath.dist,longTipReps);
longTipness = hist(longTip.vox,1:numSurf);
toc



% skelInVol = zeros(size(allVox.subs,1),1);
% skelInVol(surfVox.surf2all) = (skel.prop);
% [volClose2node] = passMaxAndDist(allVox,skelInVol,reps);
% showSubs(allVox.subs,mod(volClose2node.prop,256));








%}


