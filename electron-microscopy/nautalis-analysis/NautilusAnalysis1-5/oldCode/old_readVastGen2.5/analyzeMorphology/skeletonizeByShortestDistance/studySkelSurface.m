function[near] = studySkelSurface(surfVox,skel);

surfClose2node = surfVox.surfClose2node;

for i = 1:length(skel.node2surf)
    nearVox = find(surfClose2node.vox == skel.node2surf(i) );
    nearNum(i) = length(nearVox);
    nearSubs = surfVox.subs(nearVox,:);
    nearMid(i,:) = mean(nearSubs,1);
    nearDist = sqrt((nearSubs(:,1) - nearMid(i,1)).^2 + (nearSubs(:,2)-nearMid(i,2)).^2 + ...
        (nearSubs(:,3)-nearMid(i,3)).^2);
    nearMeanDist(i) = mean(nearDist);
    
end

near.num = nearNum;
near.mid = nearMid;
near.dist = nearMeanDist;