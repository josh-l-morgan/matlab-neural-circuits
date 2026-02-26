function[dsSubs] = downSampSubs(subs,dSamp);


sub = ceil(double(subs)/dSamp);
sub(sub<1) = 1; %%!!! Problem getting subs should never get zeros
maxSub = max(sub,[],1);
if ~isempty(sub)
    inds = sub2ind(maxSub,sub(:,1),sub(:,2),sub(:,3));
    uInds = unique(inds);
    if length(uInds>1)
        hInds = hist(inds,uInds);
    else
        hInds = length(inds);
    end
    [y x z] = ind2sub(maxSub,uInds);
    dsSubs = cat(2,uint16(y),uint16(x),uint16(z));
else
    dsSubs = [];
end