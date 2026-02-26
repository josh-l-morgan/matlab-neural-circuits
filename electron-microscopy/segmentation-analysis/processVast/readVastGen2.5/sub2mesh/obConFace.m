 function[fv] = obConFace(objSubs)
%%transform voxel list into n x 26 matrix of connectivity. 
%%also returns conDat that lists the distances and relative positions of
%%the 26 x entries of con2

%% make unique
obSiz = max(objSubs,[],1);
obInd = sub2ind(obSiz,objSubs(:,1),objSubs(:,2),objSubs(:,3));
uObInd = unique(obInd);
[y x z] = ind2sub(obSiz,uObInd);
objSubs = [y x z];


%% Turn objSubs in to more compact space. define virtDim as a space around the object buffered with one voxel

vNum = size(objSubs,1);
minSub = min(objSubs,[],1);
for i = 1:3
    objSubs(:,i) = objSubs(:,i) - minSub(i) + 5;
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
% shiftMat = ones(3,3,3);
% shiftMat(2,2,2) = 0;
% shiftInd = find(shiftMat);
% 
% [y x z] = ind2sub(size(shiftMat),shiftInd);
% shift = [y x z]-2;

shift = [ -1 0 0 ; 1 0 0 ; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1];

shiftDist = sqrt(shift(:,1).^2 + shift(:,2).^2 + shift(:,3).^2);
%[sortDist sortOrder] = sort(shiftDist,'ascend');
sortDist = shiftDist;
%shift = shift(sortOrder,:);

conDat.shift = shift;
conDat.dists = sortDist;


%% Face shift
aSquare = [0 -.5 -.5; 0 -.5 .5; 0 .5 .5; 0 .5 -.5] ;


for i = 1:size(shift,1);
    rotSquare = circshift(aSquare',ceil(i/2)-1)';
    newVert = rotSquare + repmat(shift(i,:)/2,[4,1]);
    shiftVerts(:,:,i) = newVert;
end

shiftFace = [1 2 3; 3 4 1];


%% Each voxel will have 26 potential locations to register in the table. 
con1 = zeros(length(objInds),size(shift,1)); %Create matrix for connections


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

lookupInd = zeros(max(con1(:)),1);
lookupInd(objInds) = 1:vNum;
con2 = lookupInd(con1);


%% Find faces

isFace = con2==0;
[faceID whichFace] = find(isFace);
% inFace = objInds(faceID);
% outFace = con1(isFace);

allVerts = zeros(length(faceID)*4,3);
allFaces = zeros(length(faceID)*2,3);
numFaces = 0; numVerts = 0;
for i = 1:length(faceID);
    getSub = objSubs(faceID(i),:);
    getVerts = shiftVerts(:,:,whichFace(i));
    allVerts((i-1)*4+1:(i-1)*4+4,:) = repmat(getSub,[4,1])+getVerts;
    allFaces(i*2-1,:) = shiftFace(1,:)+(i-1)*4;
    allFaces(i*2,:) = shiftFace(2,:)+(i-1)*4;
end

% 
% fv.faces = allFaces;
% fv.vertices = allVerts;
% 
% 
% clf
% 
% p = patch(fv);
% view(30,-15);
% axis vis3d;
% colormap copper
% set(p,'FaceColor','red','EdgeColor','none');
% daspect([1,1,1])
% view(3); axis tight
% camlight 
% lighting gouraud
        



%%

allVerts = (allVerts + .5);
vertInd = round(sub2ind(virtDim,allVerts(:,1),allVerts(:,2),allVerts(:,3)));
uVert = unique(vertInd);

lookupVert = zeros(max(vertInd),1);
lookupVert(round(vertInd)) = 1:length(vertInd);

face2Ind = vertInd(allFaces);
uniqueFace = lookupVert(face2Ind);

numVert = sum(lookupVert>0);

verts = zeros(numVert,3);
pickVerts = lookupVert(lookupVert>0);
lookupNewVerts = zeros(size(allVerts,1),1);
lookupNewVerts(pickVerts) = 1:length(pickVerts);

verts = allVerts(pickVerts,:)-.5;
faces = lookupNewVerts(uniqueFace);


fv.faces = faces;
fv.vertices = verts;

% 
% clf
% 
% p = patch(fv);
% view(30,-15);
% axis vis3d;
% colormap copper
% set(p,'FaceColor','red','EdgeColor','none');
% daspect([1,1,1])
% view(3); axis tight
% camlight 
% lighting flat
%         
% 
% 
% 
% 
% 
% 
% 
