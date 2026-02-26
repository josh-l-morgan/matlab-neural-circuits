

editDir = [MPN 'skel\editSkel\']
editedDir = [MPN 'skel\editSkelDone\']

mkdir(editDir);

%%
for i = 1:length(axList)
    preNodes = axNodes{i};
    maxDist(i)
    axLength(i)
    if maxDist(i)>100
        scatter(preNodes(:,3),preNodes(:,2),'.','k')
    else
        scatter(preNodes(:,3),preNodes(:,2),'.','r')
        
    end
    
    roundNodes = round(preNodes)+1;
    roundNodes(:,1) = roundNodes(:,1) * -1;
    nodeImage = zeros(max(roundNodes(:,2),[],1),max(roundNodes(:,3),[],1));
    nodeInd = sub2ind(size(nodeImage),roundNodes(:,2),roundNodes(:,3));
    nodeImage(nodeInd) = 1;
    colorNode = nodeImage*0;
    colorNode(:,:,2) = nodeImage*1000;
    colorNode(:,:,3) = nodeImage*1000;
    colorNode = uint8(colorNode);
    
    editSkelStruc.axName = axList(i);
    editSkelStruc.preNodes = preNodes;
    editSkelStruc.roundNodes = roundNodes;
    editSkelStruc.nodeImage = nodeImage;
    editSkelStruc.colorNode = colorNode;
    
    image(colorNode)
    strucName = sprintf('editSkelStruc_%d.mat',axList(i));
    imageName = sprintf('newEditSkelImage_%d.png',axList(i));
    
    imwrite(colorNode,[editDir imageName])
    save([editDir strucName],'editSkelStruc')
    
    pause(.01)
    
end

plot(maxDist,'k')
hold on
plot(minDist,'r')
hold off

%% Test completeness
missingAx = [];
for i = 1:length(axList)
    
    strucName = sprintf('editSkelStruc_%d.mat',axList(i));
    imageName = sprintf('editSkelImage_%d.png',axList(i));
    if ~exist([editedDir imageName],'file')
                missingAx = [missingAx axList(i)];
    end
end
missingAx

%%  Read
tertAxLength = zeros(length(axList),1);
axLength = zeros(length(axList),1);


editedDir
minLength = 3;
for i = 1:length(axList)
    
    strucName = sprintf('editSkelStruc_%d.mat',axList(i));
    imageName = sprintf('editSkelImage_%d.png',axList(i));
    
    load([editedDir strucName]);
   
    I = imread([editedDir imageName]);
    redChan = I(:,:,1);
    roundNodes = editSkelStruc.roundNodes;
    redChan(max(roundNodes(:,2))+3,max(roundNodes(:,3))+3) = 0;
    nodeInd = sub2ind(size(redChan),roundNodes(:,2),roundNodes(:,3));
    cutOut = redChan(nodeInd)>0
    
%     preNodes = editSkelStruc.preNodes;
%     scatter(preNodes(~cutOut,3),preNodes(~cutOut,2),'.','k')
%     hold on
%     scatter(preNodes(cutOut,3),preNodes(cutOut,2),'.','r')
%     hold off
%     pause(.01)
    
 %   tertSkel{i} = preNodes(~cutOut,:);
    
    
    
    %cellLength = zeros(length(cellList),1);
    
    
    skel = axSkel(i).skel;
    nodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keepNodes = [];
    for b = 1:length(skel.bones);
        bone = skel.bones(b);
        L = 0;
        tL = 0;
        scaleMids = scaleSubs(bone.edgeMids,voxelScale);
        scaleMids = round(scaleMids)+1;
        scaleMids(:,1) = scaleMids(:,1) * -1;
        
        for e = 1:size(bone.edges,1);
           
                edgeNodes = nodeSubs(bone.edges(e,:),:);
                edgeDif = diff(edgeNodes,1);
                edgeLength = sqrt(sum(edgeDif.^2));
                
                L = L + edgeLength;
                
                try
                    if ~redChan(scaleMids(e,2),scaleMids(e,3));
                        tL = tL + edgeLength
                    end
                catch err
                    err
                end
        end
        disp(L)
        boneLengths(b) = L;
        tertLengths(b) = tL;
        if L>minLength
            keepNodes = cat(2,keepNodes,bone.nodes);
            axLength(i) = axLength(i) + L;
            tertAxLength(i) = tertAxLength(i) + tL;
        end
    end
    keepNodes = unique(keepNodes);
    
    
    allNodeSubs = scaleSubs(skel.node2subs,voxelScale);
    keptNodeSubs = allNodeSubs(keepNodes,:);
    scatter(allNodeSubs(:,3),allNodeSubs(:,2),'.','k')
    hold on
    scatter(keptNodeSubs(:,3),keptNodeSubs(:,2),'.','r')
    hold off
    pause(.01)
    
    axNodes{i} = keptNodeSubs;
    
    
end


%% Get tert length

targ108 = find(cellList == 108);
pre108 = find(con(:,targ108)>0);
tertLengthDist = tertAxLength(pre108);
sum(tertLengthDist)

    
    
    
    
    