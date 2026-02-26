function[furthest] = findFurthest(path2max)


%%

vNum = length(path2max.vox);
uPred = unique(path2max.vox);

furthest = zeros(vNum,1);

for i = 1:length(uPred)
    
    targ = uPred(i);
    
    findPred = find(path2max.vox==targ);
    
    if ~(isempty(findPred))
        foundDist = path2max.dist(findPred);
        furthest(targ) = max(foundDist)-min(foundDist);
    else
        furthest(targ) = 0;
    end
        
end
