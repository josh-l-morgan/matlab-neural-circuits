function STA = firstpixelcorrection (STA);

timeWindow = size(STA);
meanGray = 0.5;
STA(1,:) = meanGray*ones(1,timeWindow(2));