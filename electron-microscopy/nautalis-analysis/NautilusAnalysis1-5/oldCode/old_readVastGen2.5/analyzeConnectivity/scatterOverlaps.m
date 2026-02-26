







%clear all
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\'
TPN = [MPN 'manyTerSubs\'];
reps = 100000;


downSamps = [ 1 1 1 1 1  2 4 4  8 8 8  16 32 64]
lookDists = [ 0 1 2 3 4  4 2 4  4 8 16 16 16 16]



downSamps = [ 1  1 4 8  64]
lookDists = [ 0  4 4  8  16]
lookMicrons =  .03 * 4  + lookDists .* downSamps * 8 * 0.03



for s = 1: length(lookDists)
    
    downSamp = downSamps(s);
    lookDist = lookDists(s);
    
    fileName = sprintf('terSubs_Ds%d_Ds%d_Look%d.mat',8,downSamp,lookDist)
    
    load([TPN fileName]);
    synMat = terSubs.synMat;
    touchMat = terSubs.touchMat;
    
    subplot(length(lookDists),1,length(lookDists)-s+1)
    
    useVals = touchMat>=0;
    Xdat = touchMat(useVals)/sum(touchMat(useVals))*100;
    Ydat = synMat(useVals)/sum(synMat(useVals))*100;
    [rho pc] = corr(Xdat,Ydat)
    
    hold on
    maxX = max(Xdat);
    X = [0 maxX];
    [fit1 gof] = fit(Xdat,Ydat,'poly1');
    Y = X* fit1.p1 + fit1.p2;
    line(X ,Y,'lineWidth',2);
    
       scatter(Xdat,Ydat,'.','r','MarkerSize',3)

    ylim([0 6])
    xlim([0 maxX+maxX*.1])
    hold off
    
    
    fitDat.Xdat = Xdat;
    fitDat.Ydat = Ydat;
    fitDat.fit1 = fit1;
    fitDat.gof = gof;
    fitDat.rho = rho;
    fitDat.corrP = pc;
    
    terSubs.fitDat = fitDat;
    
    
end

