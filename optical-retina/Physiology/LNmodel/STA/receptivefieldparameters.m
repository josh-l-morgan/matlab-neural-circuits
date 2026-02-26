% This script calculates spatial and temporal receptive field parameters
% from STA matrices. 
% spatial parameters: sigma, area, eccentricity
% temporal parameters: kernel, peak
% these parameters are saved in the nested structure
% spikes.space.sigma/area/eccentricity and spikes.time.kernel/peak

%See also: GAUSSPREPARATION, GAUSSIAN2D, RECEPTIVEFIELDTIME

%% LOADING FILES AND GETTING NUMBER OF CHANNELS
uiopen('load')
cd('C:\Program Files\MATLAB\R2006b\work')
s = whos('spikes');
nChannels = s.size(2);
pixelSize = 53;

%% COMPUTING SPATIAL AND TEMPORAL RECEPTIVE FIELD PARAMETERS
for i = 1:nChannels
    STA =  spikes(i).STA;
    
    %first-pixel anomaly correction
    [pixel timeWindow] = size(STA);
    meanGray = 0.5;
    STA(1,:) = meanGray*ones(1,timeWindow);

    %spatial receptive field parameters
    [z beta0 xy] = gausspreparation (STA);
    gaussOptions = statset ('MaxIter', 10000);
    beta = nlinfit(xy, z, @gaussian2d, beta0, gaussOptions);
    sigmaX = abs(beta(5));
    sigmaY = abs(beta(6));
    spikes(i).space.sigma = pixelSize * sqrt(sigmaX * sigmaY);
    spikes(i).space.area = pi * pixelSize^2 * (sigmaX * sigmaY);
    if sigmaX > sigmaY
        spikes(i).space.eccentricity = sqrt(1- sigmaY^2 / sigmaX^2);
    elseif sigmaY > sigmaX
        spikes(i).space.eccentricity = sqrt(1- sigmaX^2 / sigmaY^2);
    else
        spikes(i).space.eccentricity = 0;
    end

    %temporal receptive field parameters
    [kernel peak] = receptivefieldtime(STA);
    spikes(i).time.kernel = kernel;
    spikes(i).time.peak = peak;
end

[filename, pathname]=uiputfile('*.mat', 'save file as');
save([pathname filename], 'spikes');

