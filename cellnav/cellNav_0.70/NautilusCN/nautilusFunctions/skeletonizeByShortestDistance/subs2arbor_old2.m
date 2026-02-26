function[cellStruct] = subs2arbor(objectSubs,seed)
%%

%objectSubs = objectSubs(3000:3500,:);
if 0 % clip subs
    useSub = find([(objectSubs(:,1)<1000) & (objectSubs(:,2)<1600) & ...
        (objectSubs(:,2)>500)] );
    objectSubs = objectSubs(useSub,:);
end

showSteps = 0;

longTipReps = 10; %connectivity length of tips
minTipEccentricity = 3;
interNodeSpacing = 2; % for skeleton simplification
usePrevious = 0;
useSurface = 0;

if ~exist('seed','var')
    seedSub = [];
else
    seedSub = seed;
end

%%


%% Clean object subs
objMax = max(objectSubs,[],1);
objInds = sub2ind(objMax,objectSubs(:,1),objectSubs(:,2),objectSubs(:,3));
objInds = unique(objInds);
[oy ox oz] = ind2sub(objMax,objInds);
objectSubs = cat(2,oy,ox,oz);
clear oy ox oz


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
allVox.voxFV = renderCon(allVox.subs,allVox.conMat); 
renderFV(allVox.voxFV)
tic
allVox = bridgeConMat(allVox,seed);
toc
%% Get obj surface

disp('Get object surface')
if useSurface
    [   surfVox ] = subs2surf(allVox);
else
    surfVox = allVox;
end

numSurf = size(surfVox.subs,1);
surfVox.minMax = allVox.minMax;

moveObj = surfVox.subs*0;
for i = 1:3
    moveObj(:,i) = surfVox.subs(:,i)- surfVox.minMax{1}(i);
end

surf2cent = sqrt((surfVox.subs(:,1) - midObj(1)).^2 + (surfVox.subs(:,2) - midObj(2)).^2 + ...
    (surfVox.subs(:,3) - midObj(3)).^2);

surfVox.dist2cent = surf2cent;


%[maxObj objVol] = showMaxProp(moveObj,ones(size(surfVox.subs,1),1)*1000,minMax);
%[labField labNum] = bwlabeln(objVol);


%% Find Seed
if isempty('seedSub') % get new seed in middle if none exists
    seedSub = midObj;
end

dist2mid = sqrt((surfVox.subs(:,1)-seedSub(1)).^2 + (surfVox.subs(:,2)-seedSub(2)).^2 + ...
    (surfVox.subs(:,3)-seedSub(3)).^2);

minDist = min(dist2mid);
seed = find(dist2mid==minDist,1);
surfVox.firstSeed = seed;


%% find shortest path from seed to surface voxels

[surfPath] = conMatRad2allShortestSurface(surfVox);
surfVox.path2surf = surfPath;
edgy  = 1./(surfVox.path2surf.d2sMat+1);
surfVox.path2surf.penalty = sum(edgy,2)*10;

%% Find path through middle to seed
clf
[seedPath] = conMat2allShortestMiddle3(surfVox,seed);
clf
surfVox.seedPath = seedPath;
%showLongPred3D(surfVox.subs,seedPath);

%% get distance to longer paths
[pathL] = drawPathDist(seedPath.pred,seedPath.dist,seedPath.predLength);  %% organize all paths according furthest from seed
[sortPaths pIDX]  = sort(pathL.lengths,'descend');
showTips = pIDX(1:min(100,length(pathL.tips)/100));
if showSteps, showPathL3D3(surfVox.subs,pathL,[],showTips), end

%prop = 1./(pathL.voxLengths + rand(1,length(pathL.voxLengths))/100); % penalize short paths


%% find nearest long path by max spread
prop = 1./(pathL.voxLengths/100); % penalize short paths

propPath = path2prop(surfVox,prop,2) %find nearest bigger property
[sortPaths pIDX]  = sort(propPath.lengths,'descend');
showTips = pIDX(1:min(100,length(propPath.tips)/100));
if showSteps, showPathL3D3(surfVox.subs,propPath,[],showTips), end


%% Get path morphology
% prop = surfVox.path2surf.rads.atMaxAxRad;
prop = surfVox.path2surf.rads.meanAll;
minTip = 5; %minimum length
ratThresh = 0.3; % rad to length

sProp = spreadPredProp(propPath,prop,1);
mProp = maxPredProp(propPath,sProp,10);
if showSteps, showPathL3D3(surfVox.subs,propPath,mProp,showTips), end

minRadRatProp = mProp*0;
minRadRat = zeros(length(propPath.tips),1);
getLength = minRadRat;
for t = 1:length(propPath.tips)
    isBranch =  find(propPath.owner == propPath.tips(t));
    radRat = mProp(isBranch)./(propPath.lengths(t));
    radRat(mProp(isBranch)>propPath.lengths(t)) = 1;
    minRadRat(t) = min(radRat);
    minRadRatProp(isBranch) = min(radRat);
    getLength(t) = sum(propPath.predLength(isBranch));
end

showTips = find((minRadRat< ratThresh) & (getLength>=minTip));
%showPathL3D3(surfVox.subs,propPath,minRadRatProp,showTips)
%     [sortPaths pIDX]  = sort(minRadRat,'descend');
%     showTips = pIDX(1:length(propPath.tips)/100);
%     showTips = pIDX(sortPaths<.3);
if showSteps, showPathL3D3(surfVox.subs,propPath,minRadRatProp,showTips), end

[surfSkel] = drawPathBones(propPath.pred,seedPath.dist,propPath.predLength,propPath.tips(showTips))
surfVox.surfSkel = surfSkel;

clf
if showSteps, showPathL3D3(surfVox.subs,surfSkel,prop,[]), end


%% Group to nodes remove non independent spurrs

[skel] = smoothBones(surfVox,1);
prop = ones(1,size(surfVox.subs,1));
prop(skel.node2surf) = 0;
lookupNode = prop*0;
lookupNode(skel.node2surf) = 1:length(skel.node2surf);
path2path = conMat2allShortestProp(surfVox,prop,1000);
vox2node = lookupNode(path2path.vox);
%showSkelChunks3D(surfVox.subs,skel,vox2node)

indi = checkNodeIndependence(skel,vox2node,surfVox.conMat);
useTips = skel.node2surf([skel.bones(find(indi.branches>2)).tip]);
[surfSkel] = drawPathBones(propPath.pred,seedPath.dist,propPath.predLength,useTips)


[skel] = smoothBones(surfVox,1);
prop = ones(1,size(surfVox.subs,1));
prop(skel.node2surf) = 0;
lookupNode = prop*0;
lookupNode(skel.node2surf) = 1:length(skel.node2surf);
path2path = conMat2allShortestProp(surfVox,prop,1000);
vox2node = lookupNode(path2path.vox);

surfSkel.vox2node = vox2node
surfVox.surfSkel = surfSkel;
clf
if showSteps, showPathL3D3(surfVox.subs,surfSkel,prop,[]), end



%% clean up
%[skel] = bones2skel(surfVox,skel);
skel.nodeRad = surfVox.path2surf.rads.atMaxAxRad(skel.node2surf);
[skel.bridges] = bridgeGaps(skel,surfVox.seedPath);
surfVox.skel = skel;
clf
if showSteps, showBones3D3(surfVox.subs,skel,skel.nodeRad), end

%[arbor] = group2nodes(surfVox); %????????????
arbor = makeArbor(surfVox);
%bridgeArbor(arbor)
%showArbor3D(arbor);
%showRads3D(arbor); 
%showRadSurf(pos,edge,rad,nodeCol)

%showArborChunks3D(arbor);
%showArborVox(arbor)



%%


cellStruct.allVox = allVox;
cellStruct.surfVox = surfVox;
cellStruct.skel = skel;
cellStruct.date = date;
cellStruct.arbor = arbor;


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


 %% Simplify skel
    % disp('reducing number of skeleton nodes')
    % [skel] = simpleSkel(surfVox,interNodeSpacing);
    % [bridges] = bridgeGaps(skel,seedPath);
    % skel.bridge = bridges;
    %
    % skel.minMax = surfVox.minMax;
    %
    % skelIm = showSkel(skel,1);
    
    %Consolidate and turn to arbor


 if 0
        %%Squeeze path
        [squeezedPathLN] = squeezePaths(surfVox)
        nodeCount = squeezedPathLN.owned;
        showPred3D(surfVox.subs,squeezedPathLN.pred,squeezedPathLN.seed)
        path = squeezedPathLN;
        prop = path.lengths;
        [sortPaths pIDX]  = sort(prop,'descend');
        clf
        for i = 1:40
            t = pIDX(i);
            disp(sprintf('pcRatio = %.2f, length = %.2f, cNum = %d',...
                prop(t),path.lengths(t),childNum(t)))
            showPred3D(surfVox.subs,path.pred,path.bases(t),path.tips(t))
            hold on
            pause(.01)
        end
        
        %showLongPred3D(surfVox.subs,squeezedPathLN);
        
        
        %%find local max
        [manyNode] = countMaxAndDist(surfVox,nodeCount+rand(size(nodeCount))/1000,5);
        showMaxProp(moveObj,nodeCount(manyNode.vox));
        
        minNode = 2;
        maxNodes = unique(manyNode.vox);
        maxNodeProp = nodeCount*0;
        maxNodeProp(maxNodes) = 1;
        tipProp = maxNodeProp & (nodeCount>minNode);
        tips = find(tipProp);
        
        showMaxProp(moveObj,maxNodeProp*1000+10);
        showMaxProp(moveObj,tipProp*1000+10);
        
        increment = 0;
        
          [surfSkel] = drawPathBones(squeezedPathLN.pred,seedPath.dist,increment,tips)
        clf
        showPath3D(surfVox,surfSkel)
        surfVox.surfSkel = surfSkel;
        
    end


%}


