function outmat=makerigidmatrix(inmat)
%This function constrains a 3x3 matrix inmat to be a rigid transformation 
%and returns the result in outmat. Assumes that the rotation is in the upper left corner and
%the translation in the right column.
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

%estimate rotation angle
cphi=(inmat(1,1)+inmat(2,2))/2; 
sphi=(-inmat(1,2)+inmat(2,1))/2;
phi=atan2(sphi,cphi);

tx=inmat(1,3); ty=inmat(2,3);

outmat=[[cos(phi) -sin(phi) tx]; [sin(phi) cos(phi) ty]; [0 0 1]];