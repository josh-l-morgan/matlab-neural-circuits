function[p] = mapSkelPropIntoPatch(p,pos,col,alph)


pts = p.Vertices; 


if 0
    hold on
    scatter3(pos(:,1),pos(:,2),pos(:,3),1,'r')
    hold on
    scatter3(pts(:,1),pts(:,2),pts(:,3),1,'g')
    hold off
end

% create vertices(row) by skeleton pos (column) matrix
% dif = cat(3,pos(:,1)'-pts(:,1),pos(:,2)'-pts(:,2),pos(:,3)'-pts(:,3));
% dist = sqrt(sum(dif.^2,3)); 

dists = zeros(size(pts,1),size(pos,1));
'mapping vertices to positions'
for d = 1:size(pos,1)
    sprintf('%d of %d',d,size(pos,1))
    dists(:,d) = sqrt((pts(:,1)-pos(d,1)).^2 + (pts(:,2)-pos(d,2)).^2 +...
        (pts(:,3)-pos(d,3)).^2); 
end
    
minDist = min(dists,[],2);
closest = minDist * 0;
for i = 1:length(minDist)
    closest(i) = find(dists(i,:) == minDist(i),1);
end

A = alph(closest)';
C = col(closest,:);
p.FaceColor = 'interp';
p.FaceVertexCData = C;
p.FaceVertexAlphaData = A;
p.FaceAlpha = 'interp';


    
    