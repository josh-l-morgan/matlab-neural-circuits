function [zhat] = gaussian2d (beta, xy)
%[zhat] = Gaussian2D (beta, xy)
%This function contains the equation for fitting the 2D Gaussian. I takes
%as input a row vector of initial parameters (beta) and a matrix of x and y
%coordinates (xy) for the fitting. The output is a vector of z values
%(zhat) which get passed to nlinfit in mastergauss
%See also: GAUSSPREPARATION, MASTERGAUSS

%% ASSIGNING VARIABLES TO BETA
amplitude = beta(1);
xOffset = beta(2);
yOffset = beta(3);
theta = beta(4);
sigmaX = beta(5);
sigmaY = beta(6);
background = beta(7);

%% THE EQUATION
zhat = ( amplitude * exp(-0.5 *( ( ( (xy(:,1)-xOffset)*cos(theta) - (xy(:,2)-yOffset)*sin(theta) ) / sigmaX ).^2....
    + ( ( (xy(:,2)-yOffset)*cos(theta) + (xy(:,1)-xOffset)*sin(theta) ) / sigmaY ).^2 ) ) ) + background;
