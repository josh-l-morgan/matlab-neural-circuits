function[conInd] = lookupCon()



shiftMat = ones(3,3,3);
shiftMat(2,2,2) = 0;
shiftInd = find(shiftMat);
[y x z] = ind2sub(size(shiftMat),shiftInd);
vec = [y x z] -2;
op = y*0;
vecNum = size(vec,1);

vecVol = zeros(3,3,3);
vecVol(shiftInd) = 1:26;


lookUpVec(:,1) = vecVol(:);
pVecVol = permute(vecVol,[2 1 3]);
lookUpVec(:,2) = pVecVol(:);
pVecVol = permute(vecVol,[3 2 1]);
lookUpVec(:,3) = pVecVol(:);

showVol = zeros(3,3,3);

% 0 tilt
axVol = zeros(3,3,3);
axVol(1,2,2) = 1; axVol(3,2,2) = 2; 
planeVol = zeros(3,3,3);
planeVol(1,:,:) = 1; planeVol(2,:,:) = 2; planeVol(3,:,:) = 3; planeVol(2,2,2) = 0;
edgeVol = zeros(3,3,3);
edgeVol(:,1,1) = 1; edgeVol(:,1,3) = 2; edgeVol(:,3,1) = 3; edgeVol(:,3,3) = 4;
c = 0;
clear ax mid
for i = 1:3   
    c = c+1;
    tilt0.axI(:,i) = lookUpVec(find(axVol>0),i);
    tilt0.frontPlane(:,i) = lookUpVec(find(planeVol == 1),i);
    tilt0.midPlane(:,i) = lookUpVec(find(planeVol == 2),i);
    tilt0.backPlane(:,i) = lookUpVec(find(planeVol == 3),i);
    tilt0.edges(:,i) = lookUpVec(find(edgeVol > 0),i);
    ax(c,:) = tilt0.axI(:,i);
    mid{c} = tilt0.midPlane(:,i);
    
%       showVox = showVox * 0;
%     showVox(shiftInd(tilt0.axI(:,i))) = 2000;
%     showVox(shiftInd(tilt0.midPlane(:,i))) = 50;
%     image(sum(showVox,3))
%     pause
    
    
end

% 1 tilt
axVol = zeros(3,3,3);
axVol(2,1,1) = 1; axVol(2,3,3) = 2;
planeVol = zeros(3,3,3);
%planeVol(1,3,:) = 1; planeVol(3,1,:) = 1; planeVol(2,2,1) = 1; planeVol(2,2,3) = 1;
planeVol(:,1,3) = 1; planeVol(:,3,1) = 1; planeVol(1,2,2) = 1; planeVol(3,2,2) = 1;

showVox = zeros(3,3,3);
for i = 1:3
    i
    c = c+1;
    tilt1.axI(:,i) = lookUpVec(find(axVol>0),i);
    tilt1.midPlane(:,i) = lookUpVec(find(planeVol == 1),i);
    ax(c,:) = tilt1.axI(:,i);
    mid{c} = tilt1.midPlane(:,i);
%     showVox = showVox * 0;
%     showVox(shiftInd(tilt1.axI(:,i))) = 100;
%     showVox(shiftInd(tilt1.midPlane(:,i))) = 50;
%     image(sum(showVox,3))
%     pause
end


% 2 tilt


lookUpVec = vecVol(:);
pVecVol = flip(vecVol,1);
lookUpVec(:,2) = pVecVol(:);
pVecVol = flip(vecVol,2);
lookUpVec(:,3) = pVecVol(:);
pVecVol = flip(flip(vecVol,1),2);
lookUpVec(:,4) = pVecVol(:);


axVol = zeros(3,3,3);
axVol(3,1,1) = 1; axVol(1,3,3) = 2;
equatorVol = zeros(3,3,3);
equatorV = [1 1 3;1 3 1; 3 1 1; 1 3 3; 3 1 3; 3 3 1; 1 2 3; 3 2 1; 2 1 3; 2 3 1;1 3 2; 3 1 2];
equatorVol(sub2ind([3 3 3],equatorV(:,1),equatorV(:,2),equatorV(:,3))) = 1;

for i = 1:4
        c = c+1;
        tilt2.axI(:,i) = lookUpVec(find(axVol>0),i);
        tilt2.midPlane(:,i) = lookUpVec(find(equatorVol == 1),i);
        ax(c,:) = tilt2.axI(:,i);
        mid{c} = tilt2.midPlane(:,i);
%         showVox = showVox * 0;
%         showVox(shiftInd(tilt2.axI(:,i))) = 200;
%         %showVox(shiftInd(tilt2.midPlane(:,i))) = 50;
%         subplot(2,1,1),image(sum(showVox,3))
%         subplot(2,1,2), image(squeeze(sum(showVox,2)))
%         pause
end

conInd.ind = shiftInd;
conInd.lookUpVec = lookUpVec;
conInd.ax = ax;
conInd.mid = mid;
conInd.tilt0 = tilt0;
conInd.tilt1 = tilt1;
conInd.tilt2 = tilt2;










