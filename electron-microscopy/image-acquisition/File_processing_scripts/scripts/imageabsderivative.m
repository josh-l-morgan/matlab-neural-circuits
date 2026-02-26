function di=imageabsderivative(img)
%A function which computes an absolute image derivative
%By Daniel Berger for MIT-BCS Seung, February 2010

%compute derivative of image
dxi=(img(:,2:end)-img(:,1:end-1));
dyi=(img(2:end,:)-img(1:end-1,:));
%interpolate to get derivative at pixel positions
dxic=[dxi(:,1) dxi(:,2:end)+dxi(:,1:end-1) dxi(:,end)];
dyic=[dyi(1,:); dyi(2:end,:)+dyi(1:end-1,:); dyi(end,:)];

di=sqrt(dxic.*dxic+dyic.*dyic);
