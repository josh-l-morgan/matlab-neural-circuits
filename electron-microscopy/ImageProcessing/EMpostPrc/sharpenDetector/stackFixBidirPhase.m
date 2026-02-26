function [outStack, lagP, lagSmoothed, lagSmRound, figH] ...
   = stackFixBidirPhase(stack, nSmoothFrames, doPlot, nAvgFrames, debug)
%STACKFIXBIDIRPHASE (calcium): automatically detect, remove bidir scan offsets
%
%
%   doPlot: boolean; true to plot each frame's lag and smoothed version
%       used to correct shifts.
%
%$Id: stackFixBidirPhase.m 383 2008-11-03 17:32:58Z histed $

if nargin < 2 || isempty(nSmoothFrames), nSmoothFrames = 10; end
if nargin < 3 || isempty(doPlot), doPlot = false; end
if nargin < 4 || isempty(nAvgFrames), nAvgFrames = 1; end
if nargin < 5 || isempty(debug), debug = false; end

% t
%nAvgFrames = 10;
%nSmoothFrames = nSmoothFrames/10;


[nRows nCols nFrames] = size(stack);

%% find optimal shift for each frame
fprintf(1, '%s: list of lags: ', mfilename);

% subset of frames
if nAvgFrames > nFrames, nAvgFrames = nFrames; end
if nFrames == 1
   frameList = 1;
else
   frameList = 1:nAvgFrames:(nFrames-nAvgFrames);
end
nFrToDo = length(frameList);

lagP = repmat(NaN, 1, nFrToDo);
for iF = 1:length(frameList)
   %rFrame = double(stack(:,:,iF));
   tFrN = frameList(iF);
   if nAvgFrames > 1
       frNs = tFrN:(tFrN+nAvgFrames);
       tFrame = mean(stack(:,:,frNs),3);
   else
       tFrame = stack(:,:,iF);
   end

   maxLag = 15;
   lagNs = -maxLag:1:maxLag;
   nLags = length(lagNs);
   centerN = max(floor(nCols/7), maxLag+1);  % chop off this many points at either end
                              % of frame

   if iF == 1
       lagsToDo = lagNs;
       % for first frame, do all lags
   else
       % use only a subset of lags centered on last lag
       lastLag = lagP(iF-1);
       lagsToDo = lastLag-3:1:lastLag+3;  % set small
   end

   [tLagP isEdge] = subComputeDiff(tFrame, lagsToDo, centerN);

   if isEdge
       lagsToDo = lagNs;  % do all on this repeat
       [tLagP isEdge] = subComputeDiff(tFrame, lagsToDo, centerN);

       if isEdge
           tLagP = lagP(iF-1);
           if debug
               warning('algorithm instability in phase calculation, skipping this frame');
           else
               % not debug mode
               % this should never happen for good data; make it a crit error
               error('algorithm instability: After two repetitions, still at edge');
           end
       end
   end

   lagP(iF) = tLagP;

   fprintf(1, '%d ', lagP(iF));
end

%% compute a low-pass filtered average
% first smooth
smLagP = smooth(lagP, nSmoothFrames, 'lowess')';
lagAdjust = round(smLagP);

if any(abs(diff(lagAdjust)) > 1)
   % smooth again: bigger span
   smLagP2 = smooth(lagP, nSmoothFrames*3, 'lowess')';
   lagAdjust2 = round(smLagP2);
   if any(abs(diff(lagAdjust2)) > 1)
       wMsg = 'after 2x: lagP changes by more than 1 unit: algorithm or SNR?';
       if debug
           warning(wMsg);
       else
           error(wMsg);
       end

   end
   smLagP = smLagP2;
   lagAdjust = lagAdjust2;
end

if nFrames > 1
   % interp back to all frames
   yi = interp1(frameList,smLagP,1:nFrames);
   % fill in NaN at end
   yi(end-nAvgFrames+1:end) = yi(end-nAvgFrames);
   lagAdjust = round(yi);
else
   lagAdjust = lagP;
end




%% plot it if desired
if doPlot
   figH=figure('Tag', 'stackFixBidirPhase');
   hold on;
   lH(1) = plot(frameList, lagP, 'r');
   lH(2) = plot(frameList, smLagP, 'k');
   lH(3) = plot(lagAdjust,'b', ...
                'LineWidth', 3);

   mLag = mean(lagAdjust);
   fixLim = 5*[-1 1]+mLag;
   maxLim = [min(min(lagP), fixLim(1)) max(max(lagP), fixLim(2))];
   if range(maxLim>0)
       set(gca, 'YLim', maxLim);
   end
   xlabel('Frames');
   ylabel('Pixel shift');
   legH = legend(lH, { 'computed best shift', ...
                       'smoothed best', ...
                       'rounded to int; used' });

end

%% create new stack
outStack = imBidirShift(stack, lagAdjust);
% save outputs
lagSmoothed = smLagP;
lagSmRound = lagAdjust;


% debug stats
fprintf(1, '\nstd around smoothed mean: %g\n', ...
       chop(std(lagP - smLagP), 5));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tLagP isEdge meanDiff] = subComputeDiff(frame, lagsToDo, centerN)

[nRows nCols] = size(frame);


nLToDo = length(lagsToDo);
meanDiff = repmat(0, [1, nLToDo]);
lastL = floor(nRows/2)*2;  % allow odd number of rows
for iL = 1:nLToDo
   tL = lagsToDo(iL);

   %tM1 = rFrame(1:2:end, (centerN):(nCols-centerN));
   %tM2 = rFrame(2:2:end, (centerN+tL):(nCols-centerN+tL));
   % old difference metrics
   %tDiff = abs(tM2-tM1).*(tM2 + tM1);
   %tDiff = (tM2-tM1).^2./(tM2.^2 + tM1^2 + 1);% 0.59 std around smoothed(lag,10)
   %tDiff = abs(tM2-tM1);                    % 0.34 std+faster, use this

   % 4.2s to do 100 fr
   % I think this is faster not because of the cast itself but
   % because indexing into a double matrix is slower than indexing
   % into a uint8
   tDiff = abs( double(frame(1:2:lastL,(centerN):(nCols-centerN))) - ...
                double(frame(2:2:lastL,(centerN+tL):(nCols-centerN+tL))) );

   % multiply by mean line intensity for weighting
   lineMean = mean(frame(1:2:end),2);

   % 7.2s
   %  tDiff = abs( rFrame(1:2:end,(centerN):(nCols-centerN)) - ...
   %               rFrame(2:2:end,(centerN+tL):(nCols-centerN+tL)) );

   meanDiff(iL) = mean(mean(tDiff,2).*lineMean);
end

% find min, giving minimum fuzz for this overlap
[minV minN] = min(meanDiff);
tLagP = lagsToDo(minN);

isEdge = false;
if (minN == 1 || minN == nLToDo)
   isEdge = true;
end


