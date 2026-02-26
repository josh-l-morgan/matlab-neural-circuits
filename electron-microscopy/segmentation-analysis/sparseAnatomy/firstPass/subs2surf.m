function[objSurf conSurf] = subs2surf(objSubs, conMat);


surfVox = sum(conMat(:,1:6)>0,2)<6;
numSurf = sum(surfVox);


objSurf = objSubs(surfVox,:);
conSurf = conMat(surfVox,:);

lookupSurf = zeros(length(surfVox)+1,1);
lookupSurf(find(surfVox)+1) = 1:numSurf;

conSurf = lookupSurf(conSurf+1);

%testSurf = sum(conSurf(:,1:6)>0,2)<6;

