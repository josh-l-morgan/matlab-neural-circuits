function[skel] = drawPathBones(pred,dists,increment,tips)


%%Draws path distance and owner based on seed paths. By default, distances
%%are incrimented slightly so they can be climbed towards seed
%%by default, all voxels are used as tips.


%%

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
    bones(t).tip = sourceTip;
    nodes = [];
    
    
    for r = 1:numSurf*2
        tipDist = tipDist + increment;
        if tip<1, tip; 
                    bones(t).base = nodes(end);
                    bones(t).parent = 0;

                    break, 
        end %stop at seed
                    
                  
        if pathDist(tip)>=tipDist %stop at longer/parent branch
           
            bones(t).base = tip;
            bones(t).parent = pathOwner(tip);

            break
        else
            pathDist(tip) = tipDist;
            pathOwner(tip) = sourceTip;
            nodes = [nodes tip];
        end
        
        tip =  pred(tip);
    end
    bones(t).nodes = nodes;  
    bones(t).use = 1;


end


%% Count owned
uOwner = unique(pathOwner);
uOwner = uOwner(uOwner>0);  %zeros came up?
if length(uOwner) <2
    countOwner = sum(pathOwner>0);
else
    countOwner = hist(pathOwner,uOwner);
end
pathLengthOwned(uOwner) = countOwner;

pathLengths.owned = pathLengthOwned;
pathLengths.owner = pathOwner;
pathLengths.pathDist = pathDist;

skel.tips = sortTips;
skel.owner = pathOwner;
skel.pathDist = pathDist;
skel.owned = pathLengthOwned;
skel.pred = pred;
skel.bones = bones;


