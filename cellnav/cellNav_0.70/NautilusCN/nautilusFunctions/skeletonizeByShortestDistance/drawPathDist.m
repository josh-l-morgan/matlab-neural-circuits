function[pathLengths] = drawPathDist(pred,dists,increment,tips)


%%Draws path distance and owner based on seed paths. second input is voxel values used to sort tips.
%By default, distances
%%are incrimented slightly so they can be climbed towards seed
%%by default, all voxels are used as tips.
showProgress = 1;


%%
% if ~exist('increment','var')
%     increment = 0.000001;
% else
%     increment = 0;
% end
if ~exist('increment','var') | isempty(increment)
    increment = ones(length(pred),1);
end

numSurf = length(dists);

if ~exist('tips','var')
    tips = [1:numSurf];
    tips = setdiff(tips,pred);
end


%% Sort tips
[sortDists distIDX] = sort( dists(tips),'descend');
sortTips = tips(distIDX);


pathDist = zeros(numSurf,1);
pathOwner = pathDist;
pathLengthOwned = pathDist;

for t = 1:length(sortTips)
    sourceTip = sortTips(t);
    tip = sourceTip;
    tipDist(t) = 0;
    bases(t) = tip;
    for r = 1:numSurf*2
                  
        
        
        if tip<1,
            tip;break,
        end
            bases(t) = tip;

        if pathOwner(tip)
            break
        else
            pathDist(tip) = tipDist(t);
            pathOwner(tip) = sourceTip;
        end
        tipDist(t) = tipDist(t) + increment(tip);
        tip =  pred(tip);
       

    end
%         if showProgress
%             
%         end

end



%% recalculate base length







%% Count owned
uOwner = unique(pathOwner);
uOwner = uOwner(uOwner>0);  %zeros came up?
countOwner = hist(pathOwner,uOwner);
pathLengthOwned(uOwner) = countOwner;

pathLengths.tips = sortTips;
pathLengths.bases = bases;
pathLengths.pred = pred;
pathLengths.lengths = tipDist;%dists(sortTips)-dists(bases)+1;
pathLengths.owned = pathLengthOwned;
pathLengths.owner = pathOwner;
pathLengths.pathDist = pathDist;

lookUpTip(sortTips) = [1:length(sortTips)];
pathLengths.voxLengths = pathOwner'*0;
pathLengths.voxLengths(pathOwner>0) = pathLengths.lengths(lookUpTip(pathOwner(pathOwner>0)));
pathLengths.predLength = increment;

% %% Show
% showMaxProp(moveObj,pathDist+10);
% pause(.1)
% showMaxProp(moveObj,pathLengthOwned+10);
% % pause(.1)
% showPath = 1;
% if showPath
%     [sortPaths pIDX]  = sort(pathLengths.lengths,'descend');
%     for i = 1:length(pIDX)
%         useNodes = find(pathLengths.owner == pathLengths.tips(pIDX(i)));
%         clf
%         showPred3D(surfVox.subs,seedPath.pred,seedPath.seed,pathLengths.tips(pIDX(i)))
%         hold on
%         subs = surfVox.subs;
%         subs = subs - repmat(min(subs,[],1),[size(subs,1) 1])+1;
%         %scatter3(subs(useNodes,1),subs(useNodes,2),subs(useNodes,3),10,'g','o','filled')
%         pause(.1)
%     end
% end
% 








