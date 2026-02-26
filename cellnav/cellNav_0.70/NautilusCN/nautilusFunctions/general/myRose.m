function [tout,rout] = myRose(nn,zz)
%ROSE   Angle histogram plot.
%   ROSE(THETA) plots the angle histogram for the angles in THETA.  
%   The angles in the vector THETA must be specified in radians.
%
%   ROSE(THETA,N) where N is a scalar, uses N equally spaced bins 
%   from 0 to 2*PI.  The default value for N is 20.
%
%   ROSE(THETA,X) where X is a vector, draws the histogram using the
%   bins specified in X.
%
%   ROSE(AX,...) plots into AX instead of GCA.
%
%   H = ROSE(...) returns a vector of line handles.
%
%   [T,R] = ROSE(...) returns the vectors T and R such that 
%   POLAR(T,R) is the histogram.  No plot is drawn.
%
%   See also HIST, POLAR, COMPASS.

%   Clay M. Thompson 7-9-91
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 5

%Modified by Josh Morgan


% Form radius values for histogram triangle
if min(size(nn))==1, % Vector
  nn = nn(:); 
end
[m,n] = size(nn);
mm = 4*m;
r = zeros(mm,n);
r(2:4:mm,:) = nn;
r(3:4:mm,:) = nn;

% Form theta values for histogram triangle from triangle centers (xx)
%zz = edges;

t = zeros(mm,1);
t(2:4:mm) = zz(1:m);
t(3:4:mm) = zz(2:m+1);
   h = polar(t,r,'-r');
% 
% cax = [];
% if nargout<2
%   if ~isempty(cax)
%     h = polar(cax,t,r);
%   else
%   end
%   
%   % Register handles with m-code generator
%   if ~isempty(h)
%      mcoderegister('Handles',h,'Target',h(1),'Name','rose');
%   end
%   
%   if nargout==1, tout = h; end
%   return
% end
% 
% if min(size(nn))==1,
%   tout = t'; rout = r';
% else
%   tout = t; rout = r;
% end


