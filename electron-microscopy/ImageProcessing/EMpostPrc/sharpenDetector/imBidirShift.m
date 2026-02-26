function outImage = imBidirShift(inImage, nShiftPix)
%imBidirShift (calcium): offset every other line; changes bidir scan alignment
%   outImage = imBidirShift(inImage, nShiftPix)
%
%   Shift alternate lines of an image to correct for wrong phase/delay in
%   a bidirectional scan
%
%$Id: imBidirShift.m 92 2008-03-14 21:39:15Z histed $

[nRows nCols nFrames nPlanes] = size(inImage);

assert(length(nShiftPix) == 1 || length(nShiftPix)==nFrames);


outImage = inImage;
if length(nShiftPix) > 1
   assert(length(nShiftPix)==nFrames);
   % need to do separately for each frame
   for iF = 1:nFrames
       tNShiftPix = nShiftPix(iF);
       outImage(:,:,iF) = subShiftByConst(outImage(:,:,iF), tNShiftPix);
   end
else
   outImage = subShiftByConst(outImage, nShiftPix);
end



%%%%%%%%%%%%%%%%

function outImage = subShiftByConst(outImage, nShiftPix)
if nShiftPix > 0
   outImage(2:2:end, 1:(end-nShiftPix), :) ...
       = outImage(2:2:end, (1+nShiftPix):end, :);
elseif nShiftPix < 0
   nS = abs(nShiftPix);
   outImage(1:2:end, 1:(end-nS), :) ...
       = outImage(1:2:end, (1+nS):end, :);
elseif nShiftPix == 0
   % pass
end