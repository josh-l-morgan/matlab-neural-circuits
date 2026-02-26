function[cellStruct] = subs2arbor(allVox)
%%Skeletonize cell using a shortest path algorithm on a list of xyz
%%possitions

%%2020 + 02 + 24 strategy
%%find path to seed.  Sort into paths with path lengths. Make permanant tip list.  Find Parent child relationship by
%%spreading to longer paths. Find new shortest paths from childrent to
%%parents. Store new shortest paths as bones.

showSteps = 0;

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


%% Redefine paths according to tips, starting with longest tip to seed path
[pathL] = drawPathDist(seedPath.pred,seedPath.dist,seedPath.predLength);  %% organize all paths according furthest from seed
[sortPaths pIDX]  = sort(pathL.lengths,'descend');
showTips = pIDX(1:round(min(100,length(pathL.tips)/100)));
if showSteps, showPathL3D3(allVox.subs,pathL,pathL.voxLengths,showTips), end

if 0 % attempt to fix parrallel paths

    tips1 = pathL.tips; %remember first tips
    lengths1 = pathL.lengths;
    owner1 = pathL.owner;
    tipParent1 = pathL.owner(pathL.bases);

    %% repeat spread.  Figure out parent child relationship of tips
    if 1
        spreadL = [2];%[2 5 2 1 1]; %[ 5 3 1];

        propPathI = pathL;
        [sortPaths pIDX]  = sort(propPathI.lengths,'descend');
        showTips = pIDX(1:30);
        showPathL3D3(allVox.subs,propPathI,propPathI.voxLengths,showTips)
        pause(.01)
        for r = 1:length(spreadL)
            size(propPathI.tips)
            propPathI = setDist2Seed(propPathI);
            allVox.seedPath = propPathI;
            prop = propPathI.voxLengths; % penalize short paths
            penalty = allVox.path2surf.penalty';
            penalty = penalty/max(penalty(:)) * .00001;
            %prop = prop + penalty;
            propPathI = path2prop(allVox,-prop,spreadL(r)); %find nearest bigger property
            propPathI = setDist2Seed(propPathI);
            [sortPaths pIDX]  = sort(propPathI.lengths,'descend');
            showTips = pIDX(1:30);
            %showPathL3D3(allVox.subs,propPathI,propPathI.voxLengths,showTips)


            isTip = propPathI.owner * 0;
            isTip(propPathI.tips) = 1;
            tip2tip = isTip(propPathI.bases);
            [sortPaths pIDX]  = sort(propPathI.lengths,'descend');
            showPathL3D3(allVox.subs,propPathI,isTip,pIDX(1:5))

            pause(.1)
        end

        isTip = propPathI.owner * 0;
        isTip(propPathI.tips) = 1;
        tip2tip = isTip(propPathI.bases);
        [sortPaths pIDX]  = sort(propPathI.lengths,'descend');
        showTips = pIDX(1:30);
        showPathL3D3(allVox.subs,propPathI,propPathI.voxLengths,pIDX(1:30))
        propPath = propPathI;
    end



    %% Remember tip parent relationships
    owner2 = propPathI.owner; %new owners of paths
    tips2 = owner2(tips1); %get current owner of path of that tips1 are in
    one2two = tips1 *0; %make list to find index of path that original tip is in
    for i = 1:length(tips2)
        one2two(i) = find(propPathI.tips==tips2(i));
    end
    bases2 = propPathI.bases(one2two);
    tipParent2 = tips1 * 0;

    if 1
        tipParent2 = owner2(bases2); %Find parent tip vox (old path) of each tip

        %Pick longest tip1 from tips in tip parent2
        tipParent = tipParent2 * 0;
        for i = 1:length(tipParent2)
            pidx = find(tips2 == tipParent2(i));
            if isempty(pidx)
                %%folow the path until a valid tip is found
                v = bases2(i); %start voxel
                for p = 1:length(tipParent) %follow pred
                    v = propPathI.pred(v);
                    testTip = propPathI.owner(v);
                    if sum(tips2==testTip) %is the owner of the new voxel a new tip mapped to an old tip?
                        pidx = find(tips2==testTip);
                        break
                    end
                end
            end

            if length(pidx>0)
                pL = lengths1(pidx);
                pidx = pidx(find(pL == max(pL),1));
            end
            tipParent(i) = tips1(pidx);
        end
    else
        tipParent = owner1(bases2); %Find parent tip vox (old path) of each tip
        %%Will generate some bad streatches by using original path owners

    end



    for i = 1:10
        t = round(rand *length(tips1));
        prop = propPathI.voxLengths * 0;
        %prop = prop + rand(size(prop))*100;
        i1 = one2two(t)
        i2 = one2two(find(tips1 == tipParent(t)))
        prop(propPathI.owner == propPathI.tips(i1)) = 50;
        prop(propPathI.owner ==  propPathI.tips(i2)) = 166;
        showPathL3D3(allVox.subs,propPathI,prop,[i1 i2])
        pause(.2)
    end

    pidx = find(tipParent==tips1');
    showPathL3D3(allVox.subs,propPathI,isTip,pidx)




    %% Find shortest path from tip to base path
    %%Now solve for tips1
    %%Doesnt work

    %%find main parent
    tipParent(tipParent == tips1') = 0;
    source = tips1(find(tipParent==0));

    %%Make starting Paths
    clear shortPath
    for i = 1:length(source)
        targ = find(tips1==source(i));
        shortPath(targ).tip = source(i);
        shortPath(targ).vox = find(pathL.owner == source(i));
        shortPath(targ).pred = pathL.pred(shortPath(i).vox);
        shortPath(targ).parentTip = 0;
        shortPath(targ).childrenTip = tips1(find(tipParent==source(i)));
    end

    checkTips = [];
    for r = 1: length(tipParent)
        checkTips = [checkTips source];
        disp(sprintf('%d of %d',length(checkTips),length(tips1)))
        nextSource = [];
        for s = 1:length(source)
            prop = prop * 0;
            targ = find(tips1==source(s));
            prop(shortPath(targ).vox) = 1;
            childrenIdx = find(tipParent==source(s));
            childrenVox = tips1(childrenIdx);
            nextSource = [nextSource childrenVox];


            if length(shortPath(targ).vox) > 0
                if ~isempty(childrenVox)
                    propPathI = path2prop(allVox,-prop,length(prop)); %find nearest bigger property
                    %                 pidx = [targ; childrenIdx];
                    %                 pidx = pidx(1:min(30,length(pidx)));
                    %                 showPathL3D3(allVox.subs,propPathI,prop,pidx )
                    %                 pause
                    shortOwner = propPathI.owner;
                    shortOwner(prop>0) = 0;
                    for i = 1:length(childrenVox)
                        targ = find(tips1==childrenVox(i));
                        newTip = propPathI.owner(childrenVox(i));
                        if newTip~=childrenVox(i) %if old tip is now on a path
                            vox = childrenVox(i);
                            for p = 1:length(prop)
                                next = propPathI.pred(vox(end));
                                if shortOwner(next) == newTip
                                    vox = [vox next];
                                else
                                    break
                                end
                            end
                            shortPath(targ).vox = vox;
                        else
                            shortPath(targ).vox = find(shortOwner == childrenVox(i));
                        end
                        shortPath(targ).tip = childrenVox(i);
                        shortPath(targ).pred = propPathI.pred(shortPath(targ).vox);
                        shortPath(targ).parentTip = source(s);
                        shortPath(targ).childrenTip = tips1(find(tipParent==childrenVox(i)));
                    end
                end
            end
        end
        source = nextSource;
        if isempty(source)
            break
        end
    end
    subs = allVox.subs;


    clf
    %scatter3(subs(:,1),subs(:,2),subs(:,3)'.','k')
    hold on
    set(gca,'clipping','off')
    axis 'equal'
    for i = 1:length(shortPath)
        vox1 = shortPath(i).vox;
        vox2 = shortPath(i).pred;
        vox2(vox2==0) = vox1(vox2==0);
        if length(vox1)>4
            %scatter(subs(vox,1),subs(vox,2),'.','r')
            plot3([subs(vox1,1) subs(vox2,1)]',[subs(vox1,2) subs(vox2,2)]',...
                [subs(vox1,3) subs(vox2,3)]')
            pause(.01)
        end

    end




    if 0
        col = jet(100);
        clf
        hold on
        isTip = [shortPath.tip]
        for i = 1:length(shortPath)

            vox1 = shortPath(i).vox;
            vox2 = subs(shortPath(i).pred,:);
            prop = propPathI.owner * 0;
            prop(vox) = 1;
            showPathL3D3(allVox.subs,propPathI,prop,shortPath(i).tip)



            pause

        end
    end

    propS = pathL.owner * 0;
    %propS(pathL.owner == tips1(find(tipParent==0))) = 1;
    propS(shortPath(1).vox) = 1;
    childrenIdx = find(tipParent==shortPath(1).tip);
    childrenVox = tips1(childrenIdx);
    for t = 1:length(childrenVox)
        tip = childrenVox(t)
        prop = propPathI.owner * 0;
        prop(shortOwner == tip) = 1;
        %         showPathL3D3(allVox.subs,propPathI,prop,tip)
        %

        if sum(prop)>2
            clf
            scatter(subs(:,1),subs(:,2),'.','k')
            hold on
            scatter(subs(propS>0,1),subs(propS>0,2),'o','b')
            scatter(subs(prop>0,1),subs(prop>0,2),'+','r')
            hold off
            pause
        end
    end




    %% Record a voxel list for the path of every tip to its parent

    for i = 1:length(tips2)





    end


    if 0
        'hi'
        if ~exist('propPath','var')
            propPathI = pathL;
        end
        [sortPaths pIDX] = sort(propPathI.lengths,'descend');
        tipOwner = propPathI.owner(propPathI.bases);
        tips = propPathI.tips;
        for i = 1:length(pIDX)
            if propPathI.owned(pIDX(i)) % if there is a child

                prop = propPathI.owner == propPathI.tips(pIDX(i));
                [sortPaths pIDX] = sort(propPathI.lengths,'descend');
                showPathL3D3(allVox.subs,propPathI,prop,unique([pIDX(1:30) pIDX(i)]))
                allVox.seedPath = propPathI;
                propPathT = path2prop(allVox,-prop,length(prop)); %find nearest bigger property
                childrenIdx = find(tipOwner == tips(pIDX(i))); % get old tips attached to current base
                pause(.10)
                %showPathL3D3(allVox.subs,propPathI,prop,childrenIdx(1:10))
                %%update children
                for c = 1:length(childrenIdx)
                    cTip = tips(childrenIdx(c)); % get old tip ID
                    cOwner = propPathT.owner(cTip); % get new tip for path of old tip

                    addPath = find(propPathT.owner == cOwner); %find all new voxels that are part of the identified child path
                    newTipIdx = find(propPathT.tips == cOwner);






                    %             prop2 = propPathT.owner == cOwner;
                    %             showPathL3D3(allVox.subs,propPathT,prop2,newTipIdx)

                    %             propPathI.pred(addPath) = propPathT.pred(addPath);
                    %             propPathI.owned(addPath) = propPathT.owned(addPath);
                    %             propPathI.owner(addPath) = propPathT.owner(addPath);
                    %             propPathI.pathDist(addPath) = propPathT.pathDist(addPath);
                    %             propPathI.predLength(addPath) = propPathT.predLength(addPath);
                    %             propPathI.dist(addPath) = propPathT.dist(addPath);
                    %             propPathI.voxLengths(addPath) = propPathT.voxLengths(addPath);
                    %             oldTipIdx = find(propPathI.tips == cTip);
                    %             propPathI.tips(oldTipIdx) = propPathT.tips(newTipIdx);
                    %             propPathI.bases(oldTipIdx) = propPathT.bases(newTipIdx);
                    %             propPathI.lengths(oldTipIdx) = propPathT.lengths(newTipIdx);
                    %prop(addPath) = 30;
                end

                [sortPaths pIDX] = sort(propPathI.lengths,'descend');

                showPathL3D3(allVox.subs,propPathI,prop,unique([pIDX(1:30) pIDX(i)]))

                pause(.01)


            end
        end
    end
end


%% find nearest long path by max spread
if 1
    spreadL = [4];%[2 5 2 1 1]; %[ 5 3 1];
    prop = pathL.voxLengths; % penalize short paths
    propPath = path2prop(allVox,-prop,spreadL); %find nearest bigger property
    [sortPaths pIDX]  = sort(propPath.lengths,'descend');
    showTips = pIDX(1:round(min(100,length(propPath.tips)/100)));
    showPathL3D3(allVox.subs,propPath,propPath.voxLengths,showTips)
end


newTips = setdiff(1:length(propPath.prop),propPath.pred);
adjustDist = -propPath.prop - (propPath.path2path.dist);
adjustDist = 100 - (propPath.path2path.dist);
adjustDist = propPath.voxLengths + propPath.path2path.dist;
% stable = find(propPath.path2path.dist==0);
% scatter3(allVox.subs(stable,1),allVox.subs(stable,2),allVox.subs(stable,3),'.')

%[tipLength voxLength] = path2seed(propPath);
%[pathPL] = drawPathTip2Seed(propPath.pred,prop,propPath.predLength,propPath.tips);  %% organize all paths according furthest from seed

[pathPL] = drawPathDist(propPath.pred,adjustDist,propPath.predLength,newTips);  %% organize all paths according furthest from seed
[sortPaths pIDX]  = sort(adjustDist(newTips),'descend');
showTips = pIDX(1:min(200,length(adjustDist)/10));
if showSteps, showPathL3D3(allVox.subs,pathPL,pathPL.voxLengths,showTips), end

%drawPathDistTest(propPath.pred,pathL.voxLengths,propPath.predLength,propPath.tips,allVox.subs);

%sortPaths = pathL;


%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Needs improvement
disp('needs improvement')
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


