function[edgeLength] = getLengths(edges,pos);


if size(pos,2) == 3;
    edgeLength = sqrt((pos(edges(:,1),1)-pos(edges(:,2),1)).^2 + ...
        (pos(edges(:,1),2)-pos(edges(:,2),2)).^2 + ...
        (pos(edges(:,1),3)-pos(edges(:,2),3)).^2);
elseif size(pos,1) == 2
    edgeLength = sqrt((pos(edges(:,1),1)-pos(edges(:,2),1)).^2 + ...
        (pos(edges(:,1),2)-pos(edges(:,2),2)).^2);
end