function[subs ia ic] = uniqueSubs(subs)

if ~isempty(subs)
    minSub = min(subs(:));
    subs = subs - minSub + 1;
    
    maxSub = max(subs,[],1);
    inds = sub2ind(maxSub,subs(:,1),subs(:,2),subs(:,3));
    [uInd ia ic] = unique(inds);
    subs = subs(ia,:);
    subs = subs + minSub - 1;
else
    subs = subs;
    ia = [];
    ic = [];
end



