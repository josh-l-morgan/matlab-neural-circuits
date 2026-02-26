% result = SNR(image)
%
% Input:
%   image
%	
% Output:
%   result = SNR
%
% Principle (opt):
%   Use: result = SNR(image)
%
% Uses (self written and third-party functions):
%
% Comments on current status (opt):
% 
% ToDo (opt):
% 
% History:
%   2013 CZ Microscopy GmbH, confidential

function result = SNR(image)

    [~,~,~,tmp] = dwt2(image,'db3','mode','sym');
    
    sigma = median(abs(tmp(:)))/.6745;

    result = mean(mean(image))/sigma;

end