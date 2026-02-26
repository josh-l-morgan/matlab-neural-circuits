function[lengthOverlap] = threshOverlap(arbor1,arbor2,lookDist);

%%find length of overlap (at lookDist) of arbor1 and arbor2


if (length(arbor1)>0) & (length(arbor2)>0)
    
    mid1 = cat(1,arbor1.branches.edgeMids);
    mid2 = cat(1,arbor2.branches.edgeMids);
    
    midLength1 = cat(2,arbor1.branches.edgeLengths);
    midLength2 = cat(2,arbor2.branches.edgeLengths);
    
    overlaps = zeros(length(mid1),lookDist);

    parfor i = 1:length(mid1)
        dists = sqrt((mid2(:,1)-mid1(i,1)).^2 + ...
            (mid2(:,2)-mid1(i,2)).^2 + (mid2(:,3)-mid1(i,3)).^2);
        for o = 1:lookDist
            closeEnough = dists<=o;
            overlaps(i,o) = midLength1(i)/lookDist * sum(midLength2(closeEnough));
        end
    end
    
else
    disp('Error: No branches found for one or more arbor')
end

if exist('overlaps','var')
lengthOverlap = sum(overlaps,1);
else
   lengthOverlap = zeros(1,lookDist); 
end

