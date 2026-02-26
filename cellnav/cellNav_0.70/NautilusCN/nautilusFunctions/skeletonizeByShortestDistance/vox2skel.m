function[vox2node] = vox2skel(subs,skel)


if 0 % nearest by euclidian distance
nodePos = skel.nodePos;
dists = sqrt((subs(:,1)-nodePos(:,1)').^2 + (subs(:,2)-nodePos(:,2)').^2 + (subs(:,3)-nodePos(:,3)').^2);
ismin = dists == repmat(min(dists,[],2),[1 size(dists,2)]);
[yIDX xIDX] = find(ismin);
vox2node = zeros(size(subs,1),1);
vox2node(yIDX) = xIDX;

else
    

    
end




































