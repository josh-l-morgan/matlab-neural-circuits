function[cellStruct] = subs2arbor(allVox)
%%

showSteps = 1;

longTipReps = 10; %connectivity length of tips
minTipEccentricity = 3;
interNodeSpacing = 2; % for skeleton simplification
usePrevious = 0;
useSurface = 0;



%% Get obj surface


numSurf = size(allVox.subs,1);
allVox.minMax = allVox.minMax;
seedSub = double(allVox.seedSub);
seed = allVox.seedInd;

moveObj = allVox.subs*0;
for i = 1:3
    moveObj(:,i) = allVox.subs(:,i)- allVox.minMax{1}(i);
end

midObj = allVox.mid;
surf2cent = sqrt((allVox.subs(:,1) - midObj(1)).^2 + (allVox.subs(:,2) - midObj(2)).^2 + ...
    (allVox.subs(:,3) - midObj(3)).^2);

allVox.dist2cent = surf2cent;


%[maxObj objVol] = showMaxProp(moveObj,ones(size(allVox.subs,1),1)*1000,minMax);
%[labField labNum] = bwlabeln(objVol);


%% Find Seed
if isempty('seedSub') % get new seed in middle if none exists
    seedSub = midObj;
end

dist2mid = sqrt((allVox.subs(:,1)-seedSub(1)).^2 + (allVox.subs(:,2)-seedSub(2)).^2 + ...
    (allVox.subs(:,3)-seedSub(3)).^2);

allVox.firstSeed = seed;


%% find shortest path from seed to surface voxels

[surfPath] = conMatRad2allShortestSurface(allVox);
allVox.path2surf = surfPath;
edgy  = 1./(allVox.path2surf.d2sMat+1);
allVox.path2surf.penalty = sum(edgy,2)*100;

%% Find path through middle to seed
clf
[seedPath] = conMat2allShortestMiddle3(allVox,seed);
clf
allVox.seedPath = seedPath;
if showSteps, showLongPred3D(allVox.subs,seedPath); end

%% get distance to longer paths
[pathL] = drawPathDist(seedPath.pred,seedPath.dist,seedPath.predLength);  %% organize all paths according furthest from seed
[sortPaths pIDX]  = sort(pathL.lengths,'descend');
showTips = pIDX(1:round(min(100,length(pathL.tips)/100)));
if showSteps, showPathL3D3(allVox.subs,pathL,pathL.voxLengths,showTips), end

%prop = 1./(pathL.voxLengths + rand(1,length(pathL.voxLengths))/100); % penalize short paths


% %% find nearest long path by max spread
% spreadL = [5 3  1];
% prop = 1./(pathL.voxLengths/100); % penalize short paths
% propPath = path2prop(allVox,prop,spreadL) %find nearest bigger property
% [sortPaths pIDX]  = sort(propPath.lengths,'descend');
% showTips = pIDX(1:min(100,length(propPath.tips)/10));
% if showSteps, showPathL3D3(allVox.subs,propPath,prop*1000,showTips), end
% 
% [tipLength voxLength] = path2seed(propPath);
% [pathPL] = drawPathTip2Seed(propPath.pred,tipLength,propPath.predLength,propPath.tips);  %% organize all paths according furthest from seed
% 
% %[pathPL] = drawPathDist(propPath.pred,voxLength,propPath.predLength,pathL.tips);  %% organize all paths according furthest from seed
% [sortPaths pIDX]  = sort(pathPL.lengths,'descend');
% showTips = pIDX(1:min(100,length(pathPL.tips)/10));
% if showSteps, showPathL3D3(allVox.subs,pathPL,pathPL.voxLengths,showTips), end
% 
% %drawPathDistTest(propPath.pred,pathL.voxLengths,propPath.predLength,propPath.tips,allVox.subs);
% 
% %sortPaths = pathL;
% 
% for r = 2:length(spreadL) %repeat
%       
%     prop = 1./(pathPL.voxLengths/100); % penalize short paths
%     propPath = path2prop(allVox,prop,spreadL(r)); %find nearest bigger property
%     [tipLength voxLength] = path2seed(propPath);
% 
%     [sortPaths pIDX]  = sort(tipLength,'descend');
%     showTips = pIDX(1:min(100,length(tipLength)/10));
%     if showSteps, showPathL3D3(allVox.subs,propPath,voxLength,showTips), end
%     [pathPL] = drawPathTip2Seed(propPath.pred,tipLength,propPath.predLength,propPath.tips);  %% organize all paths according furthest from seed
%     [sortPaths pIDX]  = sort(pathPL.lengths,'descend');
%     showTips = pIDX(1:min(100,length(pathPL.tips)/10));
%     if showSteps, showPathL3D3(allVox.subs,pathPL,pathPL.voxLengths,showTips), end
%   
%     
% end
%% find nearest long path by max spread
spreadL = [5 3  1];
prop = pathL.voxLengths; % penalize short paths
propPath = path2prop(allVox,-prop,spreadL); %find nearest bigger property
newTips = setdiff(1:length(propPath.prop),propPath.pred);
adjustDist = -propPath.prop - (propPath.path2path.dist);
adjustDist = 100 - (propPath.path2path.dist);

% stable = find(propPath.path2path.dist==0);
% scatter3(allVox.subs(stable,1),allVox.subs(stable,2),allVox.subs(stable,3),'.')

%[tipLength voxLength] = path2seed(propPath);
%[pathPL] = drawPathTip2Seed(propPath.pred,prop,propPath.predLength,propPath.tips);  %% organize all paths according furthest from seed

[pathPL] = drawPathDist(propPath.pred,adjustDist,propPath.predLength,newTips);  %% organize all paths according furthest from seed
[sortPaths pIDX]  = sort(adjustDist(pathPL.tips),'descend');
showTips = pIDX(1:min(100,length(adjustDist)/10));
if showSteps, showPathL3D3(allVox.subs,pathPL,adjustDist,showTips), end

%drawPathDistTest(propPath.pred,pathL.voxLengths,propPath.predLength,propPath.tips,allVox.subs);

%sortPaths = pathL;


%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Error fix tip number?
disp('error, fix tip number?')
% for r = 2:length(spreadL) %repeat
%       
%     prop = 1./(pathPL.voxLengths/100); % penalize short paths
%     propPath = path2prop(allVox,prop,spreadL(r)); %find nearest bigger property
%     [tipLength voxLength] = path2seed(propPath);
% 
%     [sortPaths pIDX]  = sort(tipLength,'descend');
%     showTips = pIDX(1:min(100,length(tipLength)/10));
%     if showSteps, showPathL3D3(allVox.subs,propPath,voxLength,showTips), end
%     [pathPL] = drawPathTip2Seed(propPath.pred,tipLength,propPath.predLength,propPath.tips);  %% organize all paths according furthest from seed
%     [sortPaths pIDX]  = sort(pathPL.lengths,'descend');
%     showTips = pIDX(1:min(100,length(pathPL.tips)/10));
%     if showSteps, showPathL3D3(allVox.subs,pathPL,pathPL.voxLengths,showTips), end
%   
%     
% end

%% Get path morphology
% prop = allVox.path2surf.rads.atMaxAxRad;
%pathPL = propPath;
prop = allVox.path2surf.rads.meanAll;
minTip = 5; %minimum length
ratThresh = 0.3; % rad to length

sProp = spreadPredProp(pathPL,prop,1);
mProp = maxPredProp(pathPL,sProp,10);
if showSteps, showPathL3D3(allVox.subs,pathPL,mProp,showTips), end

minRadRatProp = mProp*0;
minRadRat = zeros(length(pathPL.tips),1);
getLength = minRadRat;
for t = 1:length(pathPL.tips)
    isBranch =  find(pathPL.owner == pathPL.tips(t));
    radRat = mProp(isBranch)./(pathPL.lengths(t));
    radRat(mProp(isBranch)>pathPL.lengths(t)) = 1;
    minRadRat(t) = min(radRat);
    minRadRatProp(isBranch) = min(radRat);
    getLength(t) = sum(pathPL.predLength(isBranch));
end

useTips = find((minRadRat< ratThresh) & (getLength>=minTip));

%%Check Soma
useTips = useTips(~allVox.isSoma(pathPL.tips(useTips)));

%     [sortPaths pIDX]  = sort(minRadRat,'descend');
%     showTips = pIDX(1:length(propPath.tips)/100);
%     showTips = pIDX(sortPaths<.3);
if showSteps, showPathL3D3(allVox.subs,pathPL,minRadRatProp,useTips), end

[surfSkel] = drawPathBones(pathPL.pred,seedPath.dist,pathPL.predLength,pathPL.tips(useTips))
allVox.surfSkel = surfSkel;

clf
if showSteps, showPathL3D3(allVox.subs,surfSkel,prop,[]), end


%% Group to nodes remove non independent spurrs

[skel] = smoothBones(allVox,1);
prop = ones(1,size(allVox.subs,1));
prop(skel.node2surf) = 0;
lookupNode = prop*0;
lookupNode(skel.node2surf) = 1:length(skel.node2surf);
path2path = conMat2allShortestProp(allVox,prop,1000);
vox2node = lookupNode(path2path.vox);
%showSkelChunks3D(allVox.subs,skel,vox2node)

indi = checkNodeIndependence(skel,vox2node,allVox.conMat);
useTips = skel.node2surf([skel.bones(find(indi.branches>2)).tip]);
[surfSkel] = drawPathBones(propPath.pred,seedPath.dist,propPath.predLength,useTips)


[skel] = smoothBones(allVox,1);
prop = ones(1,size(allVox.subs,1));
prop(skel.node2surf) = 0;
lookupNode = prop*0;
lookupNode(skel.node2surf) = 1:length(skel.node2surf);
path2path = conMat2allShortestProp(allVox,prop,1000);
vox2node = lookupNode(path2path.vox);

surfSkel.vox2node = vox2node
allVox.surfSkel = surfSkel;
clf
if showSteps, showPathL3D3(allVox.subs,surfSkel,prop,[]), end



%% clean up
%[skel] = bones2skel(allVox,skel);
skel.nodeRad = allVox.path2surf.rads.atMaxAxRad(skel.node2surf);
[skel.bridges] = bridgeGaps(skel,allVox.seedPath);
allVox.skel = skel;
clf
if showSteps, showBones3D3(allVox.subs,skel,skel.nodeRad), end

%[arbor] = group2nodes(allVox); %????????????
arbor = makeArbor(allVox);
%bridgeArbor(arbor)
%showArbor3D(arbor);
%showRads3D(arbor); 
%showRadSurf(pos,edge,rad,nodeCol)

%showArborChunks3D(arbor);
%showArborVox(arbor)



%%


cellStruct.allVox = allVox;
cellStruct.allVox = allVox;
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
% [maxTipPath] = passMaxAndDist(allVox,tipness,reps);
% tipness2 = hist(maxTipPath.vox,1:numSurf);


%% Find long tips
tic
disp('Find long tips')
[longTip] = passMaxAndDist(allVox,seedPath.dist,longTipReps);
longTipness = hist(longTip.vox,1:numSurf);
toc



% skelInVol = zeros(size(allVox.subs,1),1);
% skelInVol(allVox.surf2all) = (skel.prop);
% [volClose2node] = passMaxAndDist(allVox,skelInVol,reps);
% showSubs(allVox.subs,mod(volClose2node.prop,256));


 %% Simplify skel
    % disp('reducing number of skeleton nodes')
    % [skel] = simpleSkel(allVox,interNodeSpacing);
    % [bridges] = bridgeGaps(skel,seedPath);
    % skel.bridge = bridges;
    %
    % skel.minMax = allVox.minMax;
    %
    % skelIm = showSkel(skel,1);
    
    %Consolidate and turn to arbor


 if 0
        %%Squeeze path
        [squeezedPathLN] = squeezePaths(allVox)
        nodeCount = squeezedPathLN.owned;
        showPred3D(allVox.subs,squeezedPathLN.pred,squeezedPathLN.seed)
        path = squeezedPathLN;
        prop = path.lengths;
        [sortPaths pIDX]  = sort(prop,'descend');
        clf
        for i = 1:40
            t = pIDX(i);
            disp(sprintf('pcRatio = %.2f, length = %.2f, cNum = %d',...
                prop(t),path.lengths(t),childNum(t)))
            showPred3D(allVox.subs,path.pred,path.bases(t),path.tips(t))
            hold on
            pause(.01)
        end
        
        %showLongPred3D(allVox.subs,squeezedPathLN);
        
        
        %%find local max
        [manyNode] = countMaxAndDist(allVox,nodeCount+rand(size(nodeCount))/1000,5);
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
        showPath3D(allVox,surfSkel)
        allVox.surfSkel = surfSkel;
        
    end


%}


