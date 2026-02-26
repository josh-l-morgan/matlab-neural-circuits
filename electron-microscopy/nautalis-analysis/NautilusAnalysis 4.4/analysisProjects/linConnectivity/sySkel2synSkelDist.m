function[sm] = sySkel2synSkelDist(sm)


edges = sm.skelEdges;
nodes = sm.syn2skel.closestSkel;
num = length(nodes);

pos = sm.skelPos;
lengths = sqrt((pos(edges(:,1),1)-pos(edges(:,2),1)).^2 + ...
    (pos(edges(:,1),2)-pos(edges(:,2),2)).^2 + ...
    (pos(edges(:,1),3)-pos(edges(:,2),3)).^2);
eucDist = zeros(num,max(edges(:)));
linDist = zeros(num,max(edges(:)));
parfor y = 1:num
    y
    pp = node2nodeDist(edges,lengths,nodes(y));
    eucDist(y,:) = sqrt((pos(:,1)-pos(nodes(y),1)).^2 + (pos(:,2)-pos(nodes(y),2)).^2 + ...
        (pos(:,3)-pos(nodes(y),3)).^2);
    linDist(y,:) = pp.dists;
    
end


linDist2 = max(linDist,eucDist);

plot(1:1000,1:1000)
hold on
scatter(eucDist(:),linDist(:),'.');
hold off

sm.syn2skel.syn2skelLinDist = linDist2;
sm.syn2skel.syn2skelEucDist = eucDist;

sm.syn2skel.syn2synLinDist = linDist2(:,sm.syn2skel.closestSkel);
sm.syn2skel.syn2synEucDist = eucDist(:,sm.syn2skel.closestSkel);