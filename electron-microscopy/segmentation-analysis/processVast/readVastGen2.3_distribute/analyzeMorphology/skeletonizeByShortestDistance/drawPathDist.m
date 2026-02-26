function[pathLengths] = drawPathDist(pred,dists,increment,tips)


%%Draws path distance and owner based on seed paths. By default, distances
%%are incrimented slightly so they can be climbed towards seed
%%by default, all voxels are used as tips.


%%
% if ~exist('increment','var')
%     increment = 0.000001;
% else
%     increment = 0;
% end
if ~exist('increment','var')
    increment = 0;
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
    tipDist = dists(tip);
    for r = 1:numSurf*2
        tipDist = tipDist + increment;
        if tip<1,tip;break,end
        
        if pathDist(tip)
            break
        else
            pathDist(tip) = tipDist;
            pathOwner(tip) = sourceTip;
        end
        
        tip =  pred(tip);
    end
%             showMaxProp(moveObj,pathOwner*7777+30);
%             pause(.01)

end

%% Count owned
uOwner = unique(pathOwner);
uOwner = uOwner(uOwner>0);  %zeros came up?
countOwner = hist(pathOwner,uOwner);
pathLengthOwned(uOwner) = countOwner;

pathLengths.owned = pathLengthOwned;
pathLengths.owner = pathOwner;
pathLengths.pathDist = pathDist;

% %% Show
% showMaxProp(moveObj,pathDist+10);
% pause(.1)
% showMaxProp(moveObj,pathLengthOwned+10);
% pause(.1)





