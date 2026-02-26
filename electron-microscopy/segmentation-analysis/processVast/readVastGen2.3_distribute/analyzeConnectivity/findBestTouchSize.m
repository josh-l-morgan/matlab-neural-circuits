

%%find best touch size based on appositions


MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\export+14+04+27_mat\'
TPN = [MPN 'skel\'];

load([TPN 'terSubs.mat'])
minTouchRange = 2.^[0:20];
maxTouchRange = [1:20];
synMat = terSubs.synMat;
histMat = terSubs.histMat;
for i = 1:length(minTouchRange)
    i
    for p = 1:length(maxTouchRange)
        minTouch = minTouchRange(i);
        maxTouch = minTouch * maxTouchRange(p);
        
        
        disp(sprintf('Testing minTouch %d and maxTouch %d',minTouch,maxTouch))
        
        biggestTouch = 0;
        collectTouches = [];
        for pre = 1:size(histMat,1)
            for post = 1:size(histMat,2)
                histVal = histMat{pre,post} ;
                if ~isempty(histVal)
                    biggestTouch = max(biggestTouch,max(histVal));
                end
                atLeast = histVal>minTouch;
                atMost = floor(histVal/maxTouch);
                countTouch = max(atLeast,atMost);
                touchMat(pre,post) = sum(countTouch);
                collectTouches = [collectTouches histVal];
            end
        end
        
        
        useVals = touchMat>=0;
        if sum(useVals)
        Xdat = touchMat(useVals);
        Xdat = Xdat-mean(Xdat);
        Xdat = Xdat/std(Xdat);
        Ydat = synMat(useVals);
        Ydat = Ydat - mean(Ydat);
        Ydat = Ydat/std(Ydat);
        [rho pc] = corr(Xdat,Ydat);
        
                [fit1 gof] = fit(Xdat,Ydat,'poly1');

        scatter(Xdat,Ydat,'.','r')
        hold on
        maxX = max(Xdat);
        X = [0 maxX];
        Y = X* fit1.p1 + fit1.p2;
        line(X ,Y);
        hold off
        pause(.01)
        
        
        fitDat.Xdat = Xdat;
        fitDat.Ydat = Ydat;
        fitDat.fit1 = fit1;
        fitDat.gof = gof;
        fitDat.rho = rho;
        fitDat.corrP = pc;
        
        trackCorr(i,p) = rho;
        trackSSE(i,p) = gof.sse;
        else
            
        
        trackCorr(i,p) = 0;
        trackSSE(i,p) = 0;
        end
        
        
    end
end

hist(collectTouches,[1:100:max(collectTouches)])
image(trackCorr*200/max(trackCorr(:)))
