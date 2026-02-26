function[newPos,nPred] = regularEdges(pos,pred,gap,showSteps)


if ~exist('gap','var')
    gap = 0.1;
end

nodes = 1:size(pos,1);
nCount = hist(pred,nodes);
startNode = find(pred<0);

if ~exist('showSteps','var')
showSteps = 1;
end


if showSteps
    clf
    hold on
    [sortPred idx] = sort(pred,'ascend');
    for iX = 2:length(idx)
        i = idx(iX);
        
        plot([pos(i,1) pos(pred(i),1)],[pos(i,2) pos(pred(i),2)],'b')
    end
    drawnow
end


%%
clear bone boneType
bone{1}(1) = startNode;
freshBones = 1;
recNode = [];
boneParent(1) = -1;
for i = 1:length(nodes)
    if ~sum(freshBones)
        break
    end
    runBone = find(freshBones);
    for br = 1:length(runBone) %% run each bone
        b = runBone(br);
        lastNode =  bone{b}(end);
        
        for n = 1:length(nodes) %check everything downstream
            nextNodes = find(pred==lastNode);
            recNode = cat(1,recNode,nextNodes(:)); %track which nodes have been used
            %disp(sprintf('percent bonified %0.1f',length(recNode)/length(nodes)*100));
            if length(nextNodes)==1 %keep going if internode
                bone{b} = cat(1, bone{b},nextNodes);
                lastNode = nextNodes;
            elseif length(nextNodes)<1 %Stop if at tip
                boneType(b) = 1;
                freshBones(b) = 0;
                break
            elseif length(nextNodes)>1 %stop if at branch point
                boneType(b) = 0;
                numBone = length(bone);
                freshBones(b) = 0;
                for nB = 1:length(nextNodes)
                    bone{numBone+nB}=[bone{b}(end); nextNodes(nB)];
                    boneParent(numBone+nB) = b;
                    freshBones(numBone+nB) = 1;;
                end
                break
            end
            
            if 0
                clf
                scatter(pos(:,1),pos(:,2),'.','k')
                hold on
                scatter(pos(recNode,1),pos(recNode,2),'r')
                scatter(pos(nextNodes,1),pos(nextNodes,2),'g')
                hold off
                pause(.01)
            end
            
        end
        
        
    end % run Bone
    
end


if 0%showSteps
    clf
    hold on
    for b = 1:length(bone)
        bn = bone{b};
        plot([pos(bn(1:end-1),1) pos(bn(2:end),1)]',...
            [pos(bn(1:end-1),2) pos(bn(2:end),2)]','r','linewidth',5)
        bP = boneParent(b);
        if bP>0
            n1 = bone{bP}(end);
            n2 = bn(1);
            plot([pos(n1,1) pos(n2,1)], [pos(n1,2) pos(n2,2)],'b')
        end
    end
    drawnow
end


%% Change gaps within bones
clear newBone
for b = 1:length(bone)
    
    nb = bone{b};
    newMidPts = [];
    if length(nb)>1
        p1 = pos(nb(1),:);
        
        for n = 2:length(nb)-1
            p2 = pos(nb(n+1),:);
            
            dif1 = p2(1)-p1(1);
            dif2 = p2(2)-p1(2);
            dif3 = p2(3)-p1(3);
            dist = sqrt(dif1^2+dif2^2+dif3^2);
            
            
            if dist>=gap %cut up the length once you excede the gap
                steps = [0:gap:dist]/dist;
                difP = [dif1 dif2 dif3];
                pN = zeros(length(steps)-1,3);
                for s = 2:length(steps)
                    difS = difP * steps(s);
                    pN(s-1,:) = p1+difS;
                    
                    if 0
                        clf
                        scatter3(pos(nb,1),pos(nb,2),pos(nb,3),'k','.')
                        hold on
                        scatter3(p1(1),p1(2),p1(3),'g')
                        scatter3(p2(1),p2(2),p2(3),'b')
                        scatter3(pN(1:s-1,1),pN(1:s-1,2),pN(1:s-1,3),'r')
                        if ~isempty(newMidPts)
                            scatter3(newMidPts(:,1),newMidPts(:,2),newMidPts(:,3),'.','m')
                        end
                        pause
                    end
                    
                end
                p1 = pN(end,:);
                newMidPts = cat(1,newMidPts,pN);
                
                if 0
                    clf
                    scatter3(pos(nb,1),pos(nb,2),pos(nb,3),'k','.')
                    hold on
                    scatter3(p1(1),p1(2),p1(3),'g')
                    scatter3(p2(1),p2(2),p2(3),'b')
                    scatter3(pN(1:s-1,1),pN(1:s-1,2),pN(1:s-1,3),'r')
                    if ~isempty(newMidPts)
                        scatter3(newMidPts(:,1),newMidPts(:,2),newMidPts(:,3),'.','m')
                    end
                    pause
                end
                
            end
            
        end
        
        newBone{b} = cat(1,pos(nb(1),:),newMidPts,pos(nb(end),:));
        
        
    else
        newBone{b} = pos(nb(1),:);
    end
end


if 0%showSteps
    clf
    hold on

    for b = 1:length(newBone)
        bn = newBone{b};
                 scatter3(bn(:,1),bn(:,2)+.1000,bn(:,3),'.','k')
   
        plot3([bn(1:end-1,1) bn(2:end,1)],[bn(1:end-1,2) bn(2:end,2)],...
            [bn(1:end-1,3) bn(2:end,3)],'r')
        
    end
    drawnow
end




%% Rebuild skeleton with new bones and boneParent



[bpSort idx] = sort(boneParent,'ascend');


%%make new root node

root =  newBone{idx(1)}(1,:);
newPos = newBone{idx(1)}(1,:);
nPred = [-1];
countNodes = 1;
lastNodeID = zeros(length(bone),1);
lastNodeID(1) = 1;

for bX = 1:length(idx);
    b = idx(bX);
    
    oldNum = size(newPos,1); %old number of nodes
    bPos = newBone{b}(2:end,:);
    newPos = cat(1,newPos,bPos);
    
    %%set preds
    num = size(bPos,1);
    bPred = oldNum:oldNum+num-1;
    if ~isempty(bPred)
        if boneParent(b)>0
            bPred(1) = lastNodeID(boneParent(b));
        end
    end
    nPred = cat(1,nPred,bPred');
    
    lastNodeID(b) = size(newPos,1); %enter ending node ID for that bone
    
end


if showSteps

scatter(root(:,1),root(:,2),130,'m','p')
    hold on
    [sortPred idx] = sort(nPred,'ascend');
    for iX = 1:length(idx)
        i = idx(iX);
        scatter(newPos(i,1),newPos(i,2),'.','k')
        if nPred(i)>0
        plot([newPos(i,1) newPos(nPred(i),1)],[newPos(i,2) newPos(nPred(i),2)],'r')
        end
        
    end
  drawnow
  pause(1)
end



