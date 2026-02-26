function outputimage=rotateimage(inputimage, radangle, targetxsize, targetysize, method)
%rotates the inputimage by angle radangle and pastes the rotated image 
%centrally into rotimage. method is the method for interp2, for example
%'linear'.
%by Daniel Berger, March 21 2009 (for MIT-BCS Seung)

sinrotang=sin(radangle);
cosrotang=cos(radangle);

xmin=round((size(inputimage,2)-targetxsize)/2); %in image coordinates
xmax=xmin+targetxsize-1;
ymin=round((size(inputimage,1)-targetysize)/2);
ymax=ymin+targetysize-1;

%generate coordinate grid
XI=ones(targetysize,1)*(xmin:xmax);
YI=(ymin:ymax)'*ones(1,targetxsize);
midx=size(inputimage,2)/2;
midy=size(inputimage,1)/2;
  
%rotate coordinate grid
XI=XI-midx; YI=YI-midy;
XIr=cosrotang*XI-sinrotang*YI;
YIr=sinrotang*XI+cosrotang*YI;
XIr=XIr+midx;
YIr=YIr+midy;
        
%compute rotated slice
outputimage=interp2(inputimage,XIr,YIr,method);