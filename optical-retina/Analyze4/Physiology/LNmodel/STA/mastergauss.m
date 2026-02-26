% This script assembles the different functions to fit a 2D
% Gaussian to STA data. 
%See also: GAUSSIAN2D, GAUSSPREPARATION

%% GAUSSPREPARATION
[z beta0 xy] = gausspreparation (STA);

%% NLINFIT
gaussOptions = statset ('MaxIter', 10000);
beta = nlinfit(xy, z, @gaussian2d, beta0, gaussOptions);

%% FINAL OUT
[zhat] = gaussian2d (beta, xy);

% subplot(2,1,1)
% zPlot = reshape(z, 60, 80);
% surf(zPlot)
% shading flat
% colormap jet
% axis ([25 55 20 40])
% 
% subplot(2,1,2)
zhatPlot = reshape(zhat, 60, 80);
a = zhatPlot(23:43, 25:55);
surf(a)
shading flat
colormap gray
caxis([-0.5  1.5])
axis tight