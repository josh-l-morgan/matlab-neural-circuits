

%arbor1, arbor2

mid1 = cat(1,arbor1.branches.edgeMids);
mid2 = cat(1,arbor2.branches.edgeMids);

dists = zeros(size(mid1,1),1);

tic
for i = 1:length(dists)
    dists = sqrt((mid2(:,1)-mid1(i,1)).^2 + ...
        (mid2(:,2)-mid1(i,2)).^2 + (mid2(:,3)-mid1(i,3)).^2);
    minClosest(i) = min(dists);
    closestMid2(i) = find(dists == minClosest(i),1);
end


closestDist = min(minClosest);
bestMid1 = find(minClosest == closestDist);
bestMid2 = closestMid2(bestMid1);

bestEdge = [ mid1(bestMid1,:) ; mid2(bestMid2,:)];
sprintf('%.0f ',bestEdge)

toc




