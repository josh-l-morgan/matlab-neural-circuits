function[equivilantDiameter] = vol2diam(vols);


equivilantDiameter = (vols*3/4/pi).^(1/3) * 2;
% 
%         butVol = (butVols*3/4/pi).^(1/3)*butSize.voxLength * 2
%         butVol = ((butVols*butSize.voxVol)*3/4/pi).^(1/3) * 2
