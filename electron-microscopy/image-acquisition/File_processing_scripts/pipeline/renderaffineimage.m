function outputimage=renderaffineimage(inputimage, A, minx,miny,maxx,maxy, targetxsize, targetysize, method)
%Renders inputimage affinely transformed by A into outputimage
%method is the method for interp2, for example 'linear'.
%by Daniel Berger, June 4 2009 (for MIT-BCS Seung)

%global sourceimage; %this is the source image
%global source; %this is the target image (called source for further processing)

%xmin=-round((targetxsize-size(sourceimage,2))/2);
%xmax=xmin+targetxsize-1;
%ymin=-round((targetysize-size(sourceimage,1))/2);
%ymax=ymin+targetysize-1;

% minx=-100; maxx=200; targetxsize=500;
% (maxx-minx)=300
% -100+[0..499]*300/499

xcoords=minx+[0:targetxsize-1]*(maxx-minx)/(targetxsize-1);
ycoords=miny+[0:targetysize-1]*(maxy-miny)/(targetysize-1);

%generate coordinate grid
XI=ones(targetysize,1)*xcoords;
YI=ycoords'*ones(1,targetxsize);
% XI=ones(targetysize,1)*(minx:maxx);
% YI=(miny:maxy)'*ones(1,targetxsize);
midx=size(inputimage,2)/2;
midy=size(inputimage,1)/2;
  
%transform coordinate grid
XI=XI-midx; YI=YI-midy;
%XIr=cosrotang*XI-sinrotang*YI;
%YIr=sinrotang*XI+cosrotang*YI;
%XIr=XIr+midx;
%YIr=YIr+midy;
XIa=XI; YIa=YI; %initialize to same size for increased speed
for y=1:1:size(XI,1)
  for x=1:1:size(XI,2)
    XIa(y,x)=XI(y,x)*A(1,1)+YI(y,x)*A(1,2)+A(1,3);
    YIa(y,x)=XI(y,x)*A(2,1)+YI(y,x)*A(2,2)+A(2,3);
  end;
end;

XIa=XIa+midx; YIa=YIa+midy;
        
%compute rotated slice
outputimage=interp2(inputimage,XIa,YIa,method);


% ----------------------------
% 
% sinrotang=sin(radangle);
% cosrotang=cos(radangle);
% 
% xmin=round((size(inputimage,2)-targetxsize)/2); %in image coordinates
% xmax=xmin+targetxsize-1;
% ymin=round((size(inputimage,1)-targetysize)/2);
% ymax=ymin+targetysize-1;
% 
% %generate coordinate grid
% XI=ones(targetysize,1)*(xmin:xmax);
% YI=(ymin:ymax)'*ones(1,targetxsize);
% midx=size(inputimage,2)/2;
% midy=size(inputimage,1)/2;
%   
% %rotate coordinate grid
% XI=XI-midx; YI=YI-midy;
% XIr=cosrotang*XI-sinrotang*YI;
% YIr=sinrotang*XI+cosrotang*YI;
% XIr=XIr+midx;
% YIr=YIr+midy;
%         
% %compute rotated slice
% outputimage=interp2(inputimage,XIr,YIr,method);