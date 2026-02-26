function [kernel peak] = receptivefieldtime(STA)

meanGray = 0.5;
peakTime = find (std(STA) == max(std(STA)));
peakFrame = mean(STA(:,peakTime-1:peakTime+1),2);

if ( mean(STA(:))- min(STA(:))  <  max(STA(:))- mean(STA(:)) ) %for On cells
    center = find ( peakFrame >= ( mean(peakFrame) + 4*std(peakFrame) ) );
    kernel = mean(STA(center,:)) - meanGray;
    peak = 500 - find(kernel == max(kernel));
    
elseif ( mean(STA(:))- min(STA(:))  >  max(STA(:))- mean(STA(:)) ) %for Off cells
    center = find ( peakFrame <= ( mean(peakFrame) - 4*std(peakFrame) ) );
    kernel = mean(STA(center,:)) - meanGray;
    peak = 500 - find(kernel == min(kernel));

else
    disp(ERROR)
end

