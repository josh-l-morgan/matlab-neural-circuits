function[cellStruct] = subs2arbor(objectSubs,seed)
%%
longTipReps = 10; %connectivity length of tips
minTipEccentricity = 5;
interNodeSpacing = 3; % for skeleton simplification
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
    
    
    [maxObj objVol] = showMaxProp(moveObj,ones(size(surfVox.subs),1)*1000,minMax);
    [labField labNum] = bwlabeln(objVol);
    
    
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
    
    disp('Find shortest path from seed to all surfaces')
    [seedPath] = conMat2allShortest2(surfVox,seed);
    
    surfVox.seedPath = seedPath;
    showMaxProp(moveObj,mod(seedPath.segID,256));
    showMaxProp(moveObj,mod(seedPath.segID,256));
    %%
    
   [squeezedPathLN] = squeezePaths(surfVox)
   nodeCount = squeezedPathLN.owned;
    showMaxProp(moveObj,nodeCount*10+10);
    %%
    %find local max
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
    showMaxProp(moveObj,mod(surfSkel.owner*777,256)+10);
    
    surfVox.surfSkel = surfSkel;
    
    
    %% Simplify bones
    [skel] = simpleBones(surfVox,interNodeSpacing);
    %[skel] = bones2skel(surfVox,skel);
    
    [skel.bridges] = bridgeGaps(skel,surfVox.seedPath);
    %skel.bridges = [];
    surfVox.skel = skel;
    showBones(skel);
    
    
    %% Group to nodes
    [arbor] = group2nodes(surfVox);
    %bridgeArbor(arbor)
    %showMaxProp(surfVox.skel.node2subs,near.surfQual+1 * 100);
    showMaxProp(arbor.nodes.node2subs);
    showArbor(arbor,3);

    
    
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


%%


% %% Group voxels to skel
% reps = 20;
% nodeProp = zeros(size(surfVox.surfSkel.owner));
% nodeProp(skel.node2surf) = 1;
% [surfClose2node] = passMaxAndDist(surfVox,nodeProp,reps);
% surfVox.surfClose2node = surfClose2node;
% 
% 
% [skel.surf] = studySkelSurface(surfVox,skel);
% col1 = showBonesMid(skel);
% col2 = showBones(skel);
% image(uint8(col1 + col2))
% 
% [maxIm objVol] = showSubs(surfVox.subs,mod(surfClose2node.vox,256));
%  showMaxProp(moveObj,mod(surfClose2node.vox*777,256)+10);

% image(squeeze(max(objVol,[],1)))

%% Create cell images
clear skelIm
for i = 1:3
    skelIm{i} = showArbor(arbor,i);
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


cellStruct.allVox = allVox;
cellStruct.surfVox = surfVox;
cellStruct.skel = skel;
cellStruct.sideViews = combCell;
cellStruct.date = date;
cellStruct.arbor = arbor;

%% Show
showBones(skel,2);
showArbor(arbor,2);
image(combCell{1})

%showBones3D(skel,2);

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


