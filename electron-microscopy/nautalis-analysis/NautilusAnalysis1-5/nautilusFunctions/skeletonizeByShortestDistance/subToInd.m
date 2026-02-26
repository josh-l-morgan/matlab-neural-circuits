function[inds] = subToInd(siz,subs);

if length(siz) == 3
inds = sub2ind(siz,subs(:,1),subs(:,2),subs(:,3));
else
    inds = sub2ind(siz,subs(:,1),subs(:,2));

end