%load skeletons if necessary
skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\SMs\';
vgcCidList=[2 3 4 5 13 14];
loadSkel=0;
if ~exist('allSkels')
    allSkels=cell(1,6);
    for i=1:length(vgcCidList)
        curCid=vgcCidList(i);
        skelFileName=['sm_cid' + string(curCid)+'.mat'];
        curSkel=load([skelDir+skelFileName]);
        allSkels{i}=curSkel;
    end
end

if ~exist('curTis')
    global glob tis
    curTis=tis;
end

%% Trying to visualize sections of branch
vobj=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\dsObj.mat');
vtis=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\tis.mat');
vobj=vobj.dsObj;
vtis=vtis.tis;
figure();
volD=8;
test=getCidVox(14,1,vobj,vtis);
testVox=test{1};
cleanVox=testVox(testVox(:,1)>0,:);
emptyMat=false(2000,2000,500);
inds=sub2ind(size(emptyMat),cleanVox(:,1),cleanVox(:,2),cleanVox(:,3));
emptyMat(inds)=true;
testSkel=bwskel(emptyMat);
skelVox=[];
[skelVox(:,1),skelVox(:,2),skelVox(:,3)]=ind2sub(size(testSkel),find(testSkel));
volPlot=scatter3(cleanVox(1:volD:end,1),cleanVox(1:volD:end,2),cleanVox(1:volD:end,3),20,'mo','filled');
hold on
alphaV = 0.1;
set(volPlot, 'MarkerEdgeAlpha', alphaV, 'MarkerFaceAlpha', alphaV);
skelPlot=scatter3(skelVox(:,1),skelVox(:,2),skelVox(:,3),3,'co','filled');
alphaS = 0.8;
set(skelPlot, 'MarkerEdgeAlpha', alphaS, 'MarkerFaceAlpha', alphaS);

%% Using prebaked skeletons
for i=2
curSkel=allSkels{i};
curSkelD=curSkel.sm.skel2skel.linDist;
curSynD=curSkel.sm.syn2Skel.syn2SkelDist;
curSyn=curSkel.sm.syn;
curSM=curSkel.sm;
curSyn.closestNode=curSkel.sm.syn2Skel.closest;
closeNodePos=curSkel.sm.arbor.nodes.pos(curSyn.closestNode,:);
synPos=curSyn.pos*10;


figure();
scatter3(closeNodePos(:,1),closeNodePos(:,2),closeNodePos(:,3));
hold on
scatter3(synPos(:,1),synPos(:,2),synPos(:,3));

ptA=[105.6303  135.5316   25.1500];
ptB=[125.1914  123.2726   17.9312];

boneSize=zeros(size(curSM.arbor.skel.bones,2),1);
for j=1:size(curSM.arbor.skel.bones,2)
    boneSize(j)=length(curSM.arbor.skel.bones(j).nodes);
end
% get the ordered ID list for all of the bones.
boneNumber=32;

[boneSrt boneSrtIdx]=sort(boneSize);
% draw the best 256 bones by order of length.
biggestBranches=boneSrtIdx(end-boneNumber+1:end);
f5=figure();
title('boneTest');
%scatter3(1000,1000,350);
hold on
colours=turbo(boneNumber);
for k=1:length(biggestBranches)
%bigBoneID=find(boneSize==max(boneSize));
bigBoneID=biggestBranches(k);
bigBoneNodeIDs=curSM.arbor.skel.bones(bigBoneID).nodes;
bigBoneNodePos=curSM.arbor.skel.bones(bigBoneID).nodePos;
scatter3(bigBoneNodePos(:,1),bigBoneNodePos(:,2),bigBoneNodePos(:,3),5,colours(k,:));
end


end



%% new test plot
tf=figure();
hold on
for i=4
curSkel=allSkels{i};
curSkelD=curSkel.sm.skel2skel.linDist;
curSynD=curSkel.sm.syn2Skel.syn2SkelDist;
curSyn=curSkel.sm.syn;
curSM=curSkel.sm;
curSyn.closestNode=curSkel.sm.syn2Skel.closest;
closeNodePos=curSkel.sm.arbor.nodes.pos(curSyn.closestNode,:);
synPos=curSyn.pos*10;

[boneSrt boneSrtIdx]=sort(boneSize);

brhDat=curSkel.sm.arbor.branches;
parents=cell2mat({brhDat(:).parent});
bases=cell2mat({brhDat(:).base});
tips=cell2mat({brhDat(:).tip});
brhIDs=[6,8];
for j=1:length(brhIDs)
bID=brhIDs(j);
curBrhDat=curSkel.sm.arbor.branches(bID);
scatter3(curBrhDat.nodePos(:,1),curBrhDat.nodePos(:,2),curBrhDat.nodePos(:,3));
pause();
%figure();
%subBrh=ismember(tips,curBrhDat.nodes);
%sum(subBrh)
%tg=graph(curSkel.sm.arbor.edges(:,1),curSkel.sm.arbor.edges(:,2));
%tp = plot(tg,'Layout','force');
end


end

%%
parPos=curSkel.sm.arbor.nodes.pos(parents,:);

%% test code

% figure();
% tg=graph(curSkel.sm.arbor.edges(1:1000,1),curSkel.sm.arbor.edges(1:1000,2));
% tp = plot(tg,'Layout','force');

allEdges=curSkel.sm.arbor.edges;
uniqueNodes=unique(allEdges(:));
nodeCounts=histcounts(allEdges,uniqueNodes);
tipIDs=find(nodeCounts==1);
forkIDs=find(nodeCounts==3);

firstNodeID=find(curSkel.sm.arbor.nodes.pos(:,3)==min(curSkel.sm.arbor.nodes.pos(:,3)));
rootNode=firstNodeID(1);

brhStruct=struct();
for k=1:length(tipIDs)
    brhStruct(k).nodeList=[];
    previousEdgeID=0;
    curTip=tipIDs(k);
    if curTip~=rootNode
    curNode=curTip;
    while ~ismember(curNode,forkIDs)
        curEdgeID=find(any(allEdges == curNode, 2));
        curEdgeID=curEdgeID(curEdgeID~=previousEdgeID);
        curEdgeNodes=allEdges(curEdgeID,:);
        nextNode=curEdgeNodes(curEdgeNodes~=curNode);
        brhStruct(k).nodeList=[brhStruct(k).nodeList; curNode];
        curNode=nextNode;
        previousEdgeID=curEdgeID;
    end
    end
end


