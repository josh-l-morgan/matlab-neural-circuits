%% Plot synapses, prediction, and functional ROIs

global tis glob

if 0
    
    SPN =  [glob.datDir 'Analysis\Data\'];
    
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'NOI.mat']);
end

roiCid = ptDat(:,3);

runCids = [2 3 4 13];
numCid = length(runCids);
sp = ceil(sqrt(numCid));

lc = 10;

for i = 1:length(runCids);
    
    cid = runCids(i);
    v = find(NOI.cids==cid,1);
    nep = NOI.neps(v).nep;
    syn = NOI.syns(v).syn;
    preSign = NOI.preSigns(v).preSign;
    d = NOI.ds(v).d;
    W = exp(-d/lc); % Apply length constant
    
    
    
    
    
    pos = nep.pos;
    clf
    scatSkel = scatter3(pos(:,1),pos(:,2),pos(:,3),'.','clipping','off');
    ylim([min(pos(:,2)) max(pos(:,2))]);
    xlim([min(pos(:,1)) max(pos(:,1))]);
    zlim([min(pos(:,3)) max(pos(:,3))]);
      axis equal
        view([ 0 90])
    
    
    for p = 1:length(syn.pre)
        preSign = preSign*0;
        preSign(p) = 1;
        
        Won = W .* repmat(preSign==1,[1 size(W,2)]);
        sumOn = sum(Won,1);
        
        synPos = syn.pos;
        onSynPos = synPos(p,:);
        
        colormap(jet(100))
        sumOnN = sumOn * 50/mean(sumOn(:));
        
        
        onCol = ceil(sumOn * 99/max(sumOn))+1;
        
        %% roi
        useRois = find(roiCid == cid);
        roiPos = SOI.pos(useRois,:);
        roiPol = ROI.Polarity(useRois);
        roiSN = ROI.SNR(useRois);
        roiPolCol = roiPol * 50 + 50;
        roiPolCol(roiPolCol<1) = 1;
        roiPolCol(roiPolCol>100) = 1;
        
        %% show
        markerSize = 150;
        %figure
        
        %subplot(sp,sp,i)
        %subplot(1,4,i)
        
      
        scatSkel.CData = ceil(onCol);
        hold on
        try scatOn.delete; end
        scatOn = scatter3(onSynPos(:,1),onSynPos(:,2),onSynPos(:,3),markerSize,'k','o','linewidth',2,'clipping','off');
        
        hold off
        
        title(sprintf('cell %d',cid))
        drawnow
        colormap(jet(100))
        pause(.1)
    end
    
end

%% Convert um to em


X = 104.892
Y = 129.66
Z = 21.52


xEM = X / 0.004;
yEM = Y / 0.004;
zEM = Z / 0.04;

[yEM xEM zEM]













