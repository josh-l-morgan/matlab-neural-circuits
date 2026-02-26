 function[con2,conDat] = obj2con(objSubs)
%%transform voxel list into n x 26 matrix of connectivity. 
%%also returns conDat that lists the distances and relative positions of
%%the 26 x entries of con2

%% Turn objSubs in to more compact space. define virtDim as a space around the object buffered with one voxel

vNum = size(objSubs,1);
minSub = min(objSubs,[],1);
for i = 1:3
    objSubs(:,i) = objSubs(:,i) - minSub(i) + 2;
end


maxSub = max(objSubs,[],1)+10;
virtDim = maxSub+10;
maxInd = prod(virtDim);




%% requires a real space index to connectivity list index lookup table.
objInds = sub2ind(virtDim,objSubs(:,1),objSubs(:,2),objSubs(:,3)); % linear index of objSubs
objInds = unique(objInds);
objSubs = zeros(length(objInds),3);
[objSubs(:,1) objSubs(:,2) objSubs(:,3)] = ind2sub(virtDim,objInds);

%% Make shifts
shiftMat = ones(3,3,3);
shiftMat(2,2,2) = 0;
shiftInd = find(shiftMat);

[y x z] = ind2sub(size(shiftMat),shiftInd);
shift = [y x z]-2;
shiftDist = sqrt(shift(:,1).^2 + shift(:,2).^2 + shift(:,3).^2);
[sortDist sortOrder] = sort(shiftDist,'ascend');
shift = shift(sortOrder,:);

conDat.shift = shift;
conDat.dists = sortDist;




%% Each voxel will have 26 potential locations to register in the table. 
con1 = zeros(length(objInds),26); %Create matrix for connections


parfor i = 1:size(shift,1) %fill in each column of con1 with new subs
   con1(:,i) = sub2ind(virtDim,objSubs(:,1)+shift(i,1),objSubs(:,2)+shift(i,2),objSubs(:,3)+shift(i,3)); 
end


% % sparse matrix
% lookupInd = sparse(objInds,ones(1,vNum),1:vNum,max(con1(:)),1);
% con2 = full(lookupInd(con1));
% con2 = lookupInd(con1);
% clear con1
% con2 = full(con2);
% %}

lookupInd = zeros(max(con1(:)),1,'single');
lookupInd(objInds) = 1:vNum;
con2 = lookupInd(con1);


disp(sprintf('reciprical failures = %d',testReciprocal(con2)))

