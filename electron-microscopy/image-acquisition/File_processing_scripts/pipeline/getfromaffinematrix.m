function [scalex, scaley, shearx, rotang, transx, transy]=getfromaffinematrix(A)
%retrieves scaling, shearing, rotation and translation from an affine
%transformation matrix A (which has translation values in the right column)
%by Daniel Berger for MIT-BCS Seung, April 19 2009

%We assume that sheary=0
scalex=sqrt(A(1,1)*A(1,1)+A(2,1)*A(2,1));
rotang=atan2(A(2,1)/scalex,A(1,1)/scalex);

R=[[cos(-rotang) -sin(-rotang)];[sin(-rotang) cos(-rotang)]];
%rotate back shearx and scaley
v=R*[A(1,2) A(2,2)]';
shearx=v(1);
scaley=v(2);

%rotate back translation
v=R*[A(1,3) A(2,3)]';
transx=v(1);
transy=v(2);



% transxpre=A(1,3);
% transypre=A(2,3);
% 
% %compute rotation angle by using atan2 between original and transformed
% %vector
% v=[0 1 0]';
% av=A*v;
% rotang=atan2(av(1),av(2));
% 
% %remove rotation from matrix
% R=[[cos(rotang) -sin(rotang) 0]; [sin(rotang) cos(rotang) 0]; [0 0 1]];
% A2=R*A;
% 
% scalex=A2(1,1);
% shearx=A2(1,2);
% scaley=A2(2,2);
% 
% transxpost=A2(1,3);
% transypost=A2(2,3);