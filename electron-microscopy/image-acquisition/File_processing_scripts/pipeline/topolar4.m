function B=topolar4(A,dist,ang,method)
%Transforms matrix A to polar coordinates and stores result in matrix B.
%Assumes that A is a square matrix with an odd number of rows and columns.
%dist is an array of distance values and ang an array of angles to use.
%'method' is the method used by interp2; for example 'linear' or 'cubic'
%By Daniel Berger for MIT-BCS Seung, April 20 2009

B=[]; 
sa=size(A);
if (size(sa,1)~=1) || (size(sa,2)~=2) || sa(1)~=sa(2)
  return; %only 2D square matrices allowed.
end;
isodd=floor(sa(1)/2)~=sa(1)/2;
if ~isodd
  return; %only odd-row-column-matrices allowed
end;

%make column vectors
if size(dist,2)>size(dist,1)
  dist=dist';
end;
if size(ang,2)>size(ang,1)
  ang=ang';
end;

m=(sa(1)+1)/2; %midpoint coordinate

%compute coordinates for interp2
% xc=zeros(sa(1),sa(1));
%xc=zeros(nrd,nrp);
xc=zeros(size(dist,1),size(ang,1));
yc=xc;

% for p=1:nrp
%   for d=1:nrd
%     phi=(p-1)*(2*pi/nrp);
%     dist=(d-1)/nrd*(sa(1)/2);
%     xc(d,p)=cos(phi)*dist;
%     yc(d,p)=sin(phi)*dist;
%   end;
% end;

for p=1:1:size(ang,1)
  for d=1:1:size(dist,1)
    xc(d,p)=cos(ang(p))*dist(d);
    yc(d,p)=sin(ang(p))*dist(d);
  end;
end;

%figure(10);
%plot(xc+m,yc+m,'*');

B=interp2(A,xc+m,yc+m,method);