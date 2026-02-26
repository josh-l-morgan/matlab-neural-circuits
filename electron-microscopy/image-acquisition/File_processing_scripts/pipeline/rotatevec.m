function [xr,yr]=rotatevec(x,y,radangle)
%A function to rotate the vector (x,y) by angle radangle (in radians)
%By Daniel Berger, March 21 2009 (for MIT-BCS Seung)

sina=sin(radangle);
cosa=cos(radangle);
xr=cosa*x-sina*y;
yr=sina*x+cosa*y;