function center=findwafer(img,radius,expectedcenter,maxtrans,startfigure)
%A function which attempts to find the position of a wafer in image img
%The expected radius and center of the wafer are given
%A translatory search is made in x/y using +-maxtrans as range
%returns the coordinates of the center of the wafer in img as [y,x]
%By Daniel Berger for MIT-BCS Seung, February 2010

%use for example center=findwafer(cregistered,505,[575 575],100)

%compute derivative of image
% dxi=(img(:,2:end)-img(:,1:end-1));
% dyi=(img(2:end,:)-img(1:end-1,:));
% %interpolate to get derivative at pixel positions
% dxic=[dxi(:,1) dxi(:,2:end)+dxi(:,1:end-1) dxi(:,end)];
% dyic=[dyi(1,:); dyi(2:end,:)+dyi(1:end-1,:); dyi(end,:)];
% 
% di=sqrt(dxic.*dxic+dyic.*dyic);
di=imageabsderivative(img);
figure(startfigure);
imagesc(di); axis equal;

%Generate circle ring mask
cmask=zeros(size(img,1)-2*maxtrans,size(img,2)-2*maxtrans);
midy=(size(cmask,1)+1)/2;
midx=(size(cmask,2)+1)/2;
sqrad=radius*radius;
tolerance=((radius+1)*(radius+1))-(radius*radius); %(about 2 pixel tolerance)
for y=1:1:size(cmask,1)
  dy=midy-y;
  for x=1:1:size(cmask,2)
    dx=midx-x;
    sqdist=dx*dx+dy*dy; %compute square of distance to center
    if abs(sqrad-sqdist)<tolerance
      cmask(y,x)=1;
    end;
  end;
end;

figure(startfigure+1);
imagesc(cmask); axis equal;

%Do cross-correlation between gradient image and ring mask
transmap=convn_fast(di,flipdims(cmask),'valid');
figure(startfigure+2);
imagesc(transmap);
title('convn_fast translation map');

tmidy=(size(transmap,1)+1)/2; tmidx=(size(transmap,2)+1)/2;
[A,B]=max(transmap);
[C,D]=max(A);
transx=D-tmidx;  %X-translation of second image leading to maximal correlation
transy=B(D)-tmidy; %y-translation of second image leading to maximal correlation
transcorr=C;

center=[(size(img,1)+1)/2+transy (size(img,2)+1)/2+transx];