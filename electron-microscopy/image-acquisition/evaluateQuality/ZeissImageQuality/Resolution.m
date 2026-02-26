
% resolution = RRE(image,threshold)
%
% Input:
%   image: [n x n]
%   threshold
%	
% Output:
%   resolution
%
% Principle (opt):
%   Use: resolution = RRE(image,threshold)
%
% Uses (self written and third-party functions):
%
% Comments on current status (opt):
% 
% ToDo (opt):
% 
% History:
%   2013 CZ Microscopy GmbH, confidential

    function resolution = Resolution(image,threshold)

    if ~isnan(threshold)
        
        % Image windowing in real space before calculating FFT
        [dimY,dimX] = size(image);
        image = uint8(double(image).*(hann(dimY)*hann(dimX)'));

        % Calculate FFT
        powerfft = fftshift(fft2(image));
        powerfft = powerfft.*conj(powerfft);
                
        % Noise reduction
        se = [1 6 15 20 15 6 1]'*[1 6 15 20 15 6 1];
        se = se/sum(se(:));
        powerfft = imfilter(powerfft,se,'circular','same');        

        % Normalize FFT
        powerfft = powerfft/max(powerfft(:));

        % Thresholding
       
        BW = powerfft>threshold;
        BW = imfill(BW,'holes');
        
        % Calculate resolution
        parameters = regionprops(BW,'MinorAxisLength','Area');
        if length(parameters) > 1
            [~,index] = max(cat(1,parameters.Area));
            parameters = parameters(index);
        end

        resolution = 1/parameters.MinorAxisLength;
        
    else
        disp('No threshold.')        
    end
    
end
