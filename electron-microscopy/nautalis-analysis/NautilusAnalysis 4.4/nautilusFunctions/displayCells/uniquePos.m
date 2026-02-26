function[uPos] = posOr(pos1,pos2,minDist)
%% combine two lists of x,y,z coordinates and eliminate redundant positions within minDist

if ~exist(minDist,'var')
    minDist = 0;
end

use1 = ones(size(pos1,1),1);
for i = 1:size(pos1,1)
   dists = sqrt((pos2(:,1)-pos1(i,1)).^2 + (pos2(:,2)-pos1(i,2)).^2  + ...
    (pos2(:,3)-pos1(i,3)).^2 );
    nearest = min(dists);
    if nearest <= minDist
        use1(i) = 0;
    end 
end

uPos = cat(1,pos1(use1>0,:),pos2);