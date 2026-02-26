function[sameGrid] = countGrids(mid1,mid2,gridSize);

%%find length of overlap (at lookDist) of arbor1 and arbor2

    sameGrid.gridSize = gridSize;


if (length(mid1)>0) & (length(mid2)>0)
   
    mid1 = ceil(mid1/gridSize);
    mid2 = ceil(mid2/gridSize);
    
    maxMid = [max([mid1(:,1); mid2(:,1)])  max([mid1(:,2); mid2(:,2)])  max([mid1(:,3); mid2(:,3)])]+gridSize; 
    
    inds1 = sub2ind(maxMid,mid1(:,1),mid1(:,2),mid1(:,3));
    inds2 = sub2ind(maxMid,mid2(:,1),mid2(:,2),mid2(:,3));
    
    gridOverlapInd = intersect(inds1,inds2);
    gridCount = length(gridOverlapInd);
    [gY gX gZ] = ind2sub(maxMid,gridOverlapInd);
    
    sameGrid.gridCount = gridCount;
    sameGrid.overlapSubs = [gY gX gZ];
    
    
else
    disp('Error: No branches found for one or more arbor')
    sameGrid.gridCount = 0;
    sameGrid.overlapSubs = [];
    
    
    
end

