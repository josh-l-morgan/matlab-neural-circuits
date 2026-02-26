


%% set variables
useVgc = [2 3 4 5 13 14];


for v = 1:length(useVgc)
    cid = useVgc(v);
    useSM(v) = 1;
    disp(sprintf('running sm %d, %d of %d',cid,v,length(useVgc)))
    
    targ = find(runCids==cid);
    sm = sms(targ).sm;%load([smDir fileName],'syn','syn2Skel','nep');

    syn = sm.syn; % get syn information
    synD = sm.syn2Skel.syn2SynDist;
    skelD = sm.syn2Skel.syn2SkelDist;



    isBip = find(syn.preClass == 7);
    isRGC = find(syn.postClass == 1);

    synBR = synD(isBip,:);
    synBR = synBR(:,isRGC);

    allBR{v} = synBR;
end

poolBR = [];
for v = 1:length(allBR)
    poolBR = cat(1,poolBR,allBR{v}(:));
end

hRange = [0:1:160];
hBR = hist(poolBR,hRange);

clf
bar(hRange,hBR);


W = exp(-poolBR/lc); % Apply length constant
hRange2 = [0:.01:1];
hW = hist(W,hRange2);
bar(hRange2,hW)

hWscaled = hW .* hRange2;

bar(hRange2,hWscaled/sum(hWscaled(:)))

cumW = cumsum(hWscaled);
plot(cumW)


scaleRange = exp(-hRange/lc); % Apply length constant

hBRscaled = hBR .* scaleRange; 

bar(hRange, hBRscaled)

cumHBRscaled = cumsum(hBRscaled);
cumHBRscaled = cumHBRscaled/max(cumHBRscaled);

clf, hold on
plot(hRange,cumHBRscaled,'k')
bar(hRange,hBR/max(hBR(:)),'facecolor','none','edgecolor','k');

hRange(min(find(cumHBRscaled>=0.5)))














