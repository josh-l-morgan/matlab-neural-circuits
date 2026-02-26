function[cellStruct] = subs2bones(objectSubs,seedSub)
%%
longTipReps = 10; %connectivity length of tips
minTipEccentricity = 5;
interNodeSpacing = 10; % for skeleton simplification
usePrevious = 0;
useSurface = 0;

if ~exist('seed','var')
    seedSub = [];
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
    if isempty(seedSub) % get new seed in middle if none exists
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
    
    
    
    %% Find tips by spreading distances to seed
    
    
    disp('Find long tips')
    [longTip] = countMaxAndDist(surfVox,seedPath.dist,longTipReps);
    longTipness = hist(longTip.vox,1:numSurf);
    longTipness = max(longTip.countVox(3:end,:),[],1);
    tipLengths = findFurthest(longTip);
    
    maxTip = showMaxProp(moveObj,longTipness(longTip.vox));
    maxTip = showMaxProp(moveObj,tipLengths(longTip.vox));
    
    
    surfVox.longTip = longTip;
    
    %% Analyze tip shape
    disp('Analyze tip shape')
    % voxEdges =  connectRegions(surfVox, longTip.vox);
    % tipShape = getTipShape(surfVox,voxEdges);
    
    % tipShape = checkTipness(surfVox,longTip);
    % [maxTipIm] = showMaxProp(moveObj,tipShape.edge2cent);
    % [maxTipIm2] = showMaxProp(moveObj,tipShape.edge2cent(longTip.vox));
    % image(maxTipIm+maxTipIm2)
    % SE = strel('disk',10);
    % dilateTipIm = imdilate(maxTipIm,SE);
    % image(dilateTipIm)
    
    
    %tipProp = (longTipness>(longTipReps*5)) & (tipShape.edge2cent > (longTipReps/10)) & (surfVox.dist2cent'> minTipEccentricity);
    %tipProp = longTipness>(longTipReps*25);
    %showMaxProp(moveObj,tipProp(longTip.vox)*100  + (longTipness>0)*20);
    
    tipProp = tipLengths>(longTipReps*2);
    %showMaxProp(moveObj,tipProp*400  + (longTipness>0)*20);
    showMaxProp(moveObj,tipProp*1000  + (tipLengths>0)*50);

    tips = find(tipProp>0);
    
    %surfVox.longTip.voxEdges = voxEdges;
    %surfVox.longTip.tipShape = tipShape;
    surfVox.longTip.tips = tips;
    surfVox.longTip.tipProp = tipProp;
    surfVox.longTip.longTipReps = longTipReps;
    surfVox.longTip.minTipEccentricity = minTipEccentricity;
    
    %% Shorten paths
    disp('consolidating paths')
    [surfSkel] = shortenPath(surfVox);
    surfVox.surfSkel = surfSkel;
    
    colormap colorcube(256)
    showMaxProp(moveObj,mod(surfSkel.prop,256));
    showMaxProp(moveObj,surfSkel.prop);

    [imax ovol] = showMaxProp(moveObj,surfSkel.prop *100  + (longTipness>0)*50);
    
    
    %% Remove spurs
    
    
    reps = 10;
    [path2surfSkel] = passMaxAndDist(surfVox,surfSkel.prop>0,reps);
    surfVox.surfSkel.path2surfSkel  = path2surfSkel;
    
    parentBone = surfSkel.prop(path2surfSkel.vox);
    %showMaxProp(moveObj,mod(parentBone,256));
    
    bones = surfSkel.bones;
    voxMeanDist = zeros(1,length(parentBone));
    boneThresh = 3;
    for i = 1:length(bones)
        bonePos = find(parentBone == bones(i).tip);
        dist2bone = path2surfSkel.dist(bonePos);
        bones(i).meanDist = mean(dist2bone);
        bones(i).voxNum = sum(dist2bone>0);
        bones(i).use = boneThresh <= (bones(i).length/bones(i).meanDist);
        vox2bone(bonePos) = i;
    end
    
    useBones = [bones.use];
    surfVox.surfSkel.bones = bones;
    
  

%% Simplify bones
[skel] = simpleBones(surfVox,interNodeSpacing);
%[skel] = bones2skel(surfVox,skel);
[skel.bridges] = bridgeGaps(skel,surfVox.seedPath);
surfVox.skel = skel;
showBones(skel);



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

%% Group voxels to skel
reps = 20;
nodeProp = zeros(size(surfVox.surfSkel.prop));
nodeProp(skel.node2surf) = 1;
[surfClose2node] = passMaxAndDist(surfVox,nodeProp,reps);
surfVox.surfClose2node = surfClose2node;


[skel.surf] = studySkelSurface(surfVox,skel);
col1 = showBonesMid(skel);
col2 = showBones(skel);
image(uint8(col1 + col2))

[maxIm objVol] = showSubs(surfVox.subs,mod(surfClose2node.vox,256));
% image(squeeze(max(objVol,[],1)))

%% Create cell images
clear skelIm
for i = 1:3
    skelIm{i} = showBones(skel,i);
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

%% Show
image(combCell{1})

showBones(skel,2);
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


