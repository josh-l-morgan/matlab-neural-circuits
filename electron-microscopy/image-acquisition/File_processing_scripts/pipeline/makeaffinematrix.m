function A=makeaffinematrix(scalex, scaley, shearx, sheary, rotang, transx, transy)
%makes an affine transformation matrix from the given scale, shear,
%rotation and translation values
%if you want a uniquely retrievable matrix, give sheary=0
%by Daniel Berger for MIT-BCS Seung, April 19 2009

A=[[scalex shearx transx];[sheary scaley transy];[0 0 1]];
A=[[cos(rotang) -sin(rotang) 0];[sin(rotang) cos(rotang) 0];[0 0 1]] * A;
