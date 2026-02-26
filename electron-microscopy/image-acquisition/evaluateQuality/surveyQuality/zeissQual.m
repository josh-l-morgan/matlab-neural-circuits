
% resolution = RRE(I,threshold)
%
% Input:
%   I: [n x n]
%   threshold
%	
% Output:
%   resolution
%
% Principle (opt):
%   Use: resolution = RRE(I,threshold)
%
% Uses (self written and third-party functions):
%
% Comments on current status (opt):
% 
% ToDo (opt):
% 
% History:
%   2013 CZ Microscopy GmbH, confidential

function[zqual] = zeissQual(I);

    %% Get SNR
     tic
     [~,~,~,tmp] = dwt2(I,'db3','mode','sym');
     sigma = median(abs(tmp(:)))/.6745;
     snrResult = mean(mean(I))/sigma;
     'snr time';
     toc
     tic
    
     %% Get Resolution
     tic
    threshold = .00001; 
        
        % Image windowing in real space before calculating FFT
        [dimY,dimX] = size(I);
        hI = uint8(double(I).*(hann(dimY)*hann(dimX)'));

        % Calculate FFT
        powerfft = fftshift(fft2(hI));
        powerfft = powerfft.*conj(powerfft);
                
        % Noise reduction
        se = [1 6 15 20 15 6 1]'*[1 6 15 20 15 6 1];
        se = se/sum(se(:));
        powerfft = imfilter(powerfft,se,'circular','same');        

        % Normalize FFT
        powerfft = powerfft/max(powerfft(:));

        %% Thresholding
        mp = mode(log(powerfft(:)));
        for t = 1:10;    
            threshold = mp + (mp*.1)*t;
        BW = powerfft>threshold;
        BW = imfill(BW,'holes');
        
        % Calculate resolution
        parameters = regionprops(BW,'MinorAxisLength','Area');
        if length(parameters) > 1
            [~,index] = max(cat(1,parameters.Area));
            parameters = parameters(index);
        end

        resolution = 1/parameters.MinorAxisLength
        memRes(t) = resolution;
        end
        plot(memRes)
       %%
               plot(powerfft(round(size(powerfft,1)/2),:))

        plot(log(powerfft(round(size(powerfft,1)/2),:)))
        'resolution time'
    toc
    
    
    zqual.resolution = resolution;
    zqual.snr = snrResult;