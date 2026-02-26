%function[skel] = shortenAllPaths(surfVox)

colormap gray(255)
moveObj = surfVox.subs;
seedPath = surfVox.seedPath;

voxList = ones(size(surfVox.subs,1),1);
usedPred = unique(seedPath.pred);
usedPred = usedPred(usedPred>0);
voxList(usedPred) = 0;
tips = find(voxList);

numSurf = size(surfVox.conMat,1);


showMaxProp(moveObj,voxList*100+100);
pause(.01)


%% for each vox on path from tip, enter distance of longest source tip


pathDist = zeros(numSurf,1);
for t = 1:length(tips)
    tip = tips(t);
    tipDist = seedPath.dist(tip);
    for r = 1:numSurf*2
        tipDist = tipDist + 1.01;
        if tip<1,tip;break,end
        pathDist(tip) = max(pathDist(tip),tipDist);
        tip =  seedPath.pred(tip);
    end
end

showMaxProp(moveObj,pathDist+10);
pause(.01)
%% find path to voxel that is path to greatest distnce tip
reps = 10;
%[path2maxDist] = passMaxAndDist(surfVox,pathDist,reps);
[path2maxDist] = passMaxAndDist(surfVox,pathDist,reps);

newPred = path2maxDist.pred;
notChanged = sum(newPred<1)
newPred(notChanged) =  seedPath.pred(notChanged);
showMaxProp(moveObj,(path2maxDist.dist==0)*100);

showMaxProp(moveObj,mod(path2maxDist.prop,256)+10);

% 
% bases = path2maxDist.vox; %find new base for each tip
% useOld = path2maxDist.pred<1;
% %bases(useOld) = seedPath.vox(useOld);
% maxPred = path2maxDist.pred; % combine pred lists
% maxPred(useOld) = seedPath.pred(useOld); %keep preds with no new solutions
% minDist = path2maxDist.dist;
% minDist(useOld) = seedPath.dist(useOld);

%% for each vox on path from tip, enter distance of longest source tip


newPathDist = zeros(numSurf,1);
for t = 1:length(tips)
    tip = tips(t);
    tipDist = seedPath.dist(tip);
    for r = 1:numSurf*2
        if tip<1,tip;break,end
        newPathDist(tip) = max(newPathDist(tip),tipDist);
        tip =  newPred(tip);
    end
end

showMaxProp(moveObj,newPathDist+10);
pause(.01)



%% Create bones
skelTip = zeros(1,numSurf);  % keep track of which tip each voxel belongs to
numNodes = zeros(1,numSurf);
showActiveNode = numNodes;
% 
% boneLengths = path2maxDist.dist(tips);
% boneLengths(boneLengths == 0) = pathDist(tips(boneLengths == 0));
% %[sortLengths bLidx] = sort(boneLengths,'descend');
[sortLengths bLidx] = sort(seedPath.dist(tips),'descend');
sortTips = tips(bLidx);

colormap gray(256)
% colormap jet(256)
for t = 1:length(sortTips)
    tip = sortTips(t); 
    bones(t).tip = tip; %sorted space
%     boneLength = path2maxDist.dist(tip);
%     if ~boneLength
%         boneLength = newPathDist(tip);
%     end
%     bones(t).length = boneLength;
    nodes = [];
    for r = 1:numSurf*2
        if tip<1, %% if you reach the seed
            bones(t).base = nodes(end);
            bones(t).parent = 0;
            break,
        end
        
        if skelTip(tip)
            if newPathDist(skelTip(tip)) >= newPathDist(sortTips(t)) %if new path found longer path
            %if pathDist(skelTip(tip)) >= pathDist(sortTips(t)) %if new path found longer path

                bones(t).base = tip;
                bones(t).parent = skelTip(tip);
                break
            else
                nodes = [nodes tip];
                skelTip(tip) = sortTips(t); % keep track of which tip each voxel belongs to
                disp('hit smaller process, parent child and nodes not valid')
            end
        else %ir empty space
            nodes = [nodes tip];
            skelTip(tip) = sortTips(t); % keep track of which tip each voxel belongs to
        end
        
        tip =  newPred(tip); %use new pred list
       %tip = seedPath.pred(tip);
        
%         if tip>0
%           %  showActiveNode = showActiveNode * .5;
%             showActiveNode(tip) = 300;
%          end
        
        
    end
    bones(t).nodes = nodes;
    bones(t).numNodes = length(nodes);
%      if length(nodes)>10
%                 showMaxProp(moveObj,showActiveNode+(showActiveNode*0)+20);
%                 %showMaxProp(moveObj,mod(showActiveNode,256)+(showActiveNode*0)+20);
%                 pause
%         end
%                     showActiveNode(showActiveNode>0) = 100;

    
end

numNodes = hist(skelTip,[1:1:numSurf]);

skel.tips = sortTips;
skel.skelTip = skelTip; %original index space
skel.numNodes = numNodes;
skel.pred = newPred;
skel.allBones = bones;
skel.bones = bones;
skel.sub2tipIndex = bLidx;

%skel.tip2subIndex(bLidx) = 1:length(bLidx);
    showMaxProp(moveObj,numNodes+20);
pause(.01)



%%