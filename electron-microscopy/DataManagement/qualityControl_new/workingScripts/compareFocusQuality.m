TPN = GetMyDir
%qual = checkQuality(TPN);
 
%get quality and focus
if ~exist([TPN 'qual.mat'],'file')
    qual = checkQuality(TPN);
else
    load([TPN 'qual.mat'])
end

focusMos = grabFocusWB(TPN);

%% display

subplot(2,2,1)
image(fitH(focusMos))

subplot(2,2,2)
qualMos = qual.mos.qualMos;
image(fitH(qualMos))

subplot(2,2,3)
mosUM = (focusMos - median(focusMos(:))) * 10^6;
mesh(flipud(mosUM))
axis([1 7 1 10 median(mosUM(:))-2 median(mosUM(:))+2])


subplot(2,2,4) 
mesh(flipud(qual.mos.qualMos))

internalSpread = median(abs(mosUM(:)-mean(mosUM(:))));
%
%
% 
% newmos = mos-mean(mos(:));
% internalSpread = median(abs(mos(:)-mean(mos(:))));
% image(newmos/internalSpread * 50 + 128);
