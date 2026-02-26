function d=getavgtransdist(xcoord,ycoord,matrix1,matrix2)
%This function computes the average distance between the target points of the points specified in xcoord,ycoord when transformed by matrix1 and matrix2
%It is assumed that xcoord and ycoord are vectors of the same length, and that matrix1 and matrix2 are 3x3 matrices (affine)
%Example usage:
%d=getavgtransdist([1 param.scaledsize(2) param.scaledsize(2) 1],[1 1 param.scaledsize(1) param.scaledsize(1)],crmatrix1,crmatrix2);
%
%By Daniel Berger for MIT-BCS Seung, June 17th

d=-1;
if (size(xcoord,1)~=size(ycoord,1))||(size(xcoord,2)~=size(ycoord,2))
  disp('ERROR: Input point coordinate numbers do not match.');  
  return;
end;

if (size(matrix1,1)~=3)||(size(matrix1,2)~=3)||(size(matrix2,1)~=3)||(size(matrix2,2)~=3)
  disp('ERROR: Matrices are not 3x3');
  return;
end;

nrofpoints=max(size(xcoord));

d=0;
for point=1:1:nrofpoints
  vec=[xcoord(point) ycoord(point) 1]';
  tvec1=matrix1*vec;
  tvec2=matrix2*vec;
  
  d=d+ sqrt((tvec2(1)-tvec1(1))*(tvec2(1)-tvec1(1))+(tvec2(2)-tvec1(2))*(tvec2(2)-tvec1(2)));
end;

d=d/nrofpoints;