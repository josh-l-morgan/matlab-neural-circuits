function[surfVox] = subs2surf(allVox);


isSurf = sum(allVox.conMat(:,1:6)>0,2)<6;
numSurf = sum(isSurf);


objSurf = allVox.subs(isSurf,:);
conSurf = allVox.conMat(isSurf,:);

lookupSurf = zeros(length(isSurf)+1,1);
surf2all = find(isSurf);
lookupSurf(surf2all+1) = 1:numSurf;

conSurf = lookupSurf(conSurf+1);

all2surf = lookupSurf(2:end);

surfVox.name = 'surface voxels';
surfVox.source = allVox.name;
surfVox.subs = objSurf;
surfVox.conMat = conSurf;
surfVox.conDat = allVox.conDat;
surfVox.all2surf = all2surf;
surfVox.surf2all = surf2all;




