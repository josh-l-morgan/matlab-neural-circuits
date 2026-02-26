
%% make unique
obSiz = max(obSub,[],1);
obInd = sub2ind(obSiz,obSub(:,1),obSub(:,2),obSub(:,3));
uObInd = unique(obInd);
[y x z] = ind2sub(obSiz,uObInd);
obSub = [y x z];

%% Get surface voxels
conMat = obj2con(obSub);
isSurf = sum(conMat(:,1:6)>0,2)<6;
numSurf = sum(isSurf);

objSurf = obSub(isSurf,:);
conSurf = conMat(isSurf,:);

lookupSurf = zeros(length(isSurf)+1,1);
surf2all = find(isSurf);
lookupSurf(surf2all+1) = 1:numSurf;

conSurf = lookupSurf(conSurf+1);
all2surf = lookupSurf(2:end);

%% Get faces 
conMat = obj2con(obSub);
isFace = conMat(:,1:6) == 0;

isSurf = sum(conMat(:,1:6)>0,2)<6;
numSurf = sum(isSurf);

objSurf = obSub(isSurf,:);
conSurf = conMat(isSurf,:);

lookupSurf = zeros(length(isSurf)+1,1);
surf2all = find(isSurf);
lookupSurf(surf2all+1) = 1:numSurf;

conSurf = lookupSurf(conSurf+1);
all2surf = lookupSurf(2:end);


%%

verticies = objSurf;
edges = conSurf;



