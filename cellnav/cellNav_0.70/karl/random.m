%fxPts=zeros(6,3);
%mvPts=zeros(6,3);

geoT=fitgeotrans(mvPts(:,[1 2]),fxPts(:,[1 2]),'affine');