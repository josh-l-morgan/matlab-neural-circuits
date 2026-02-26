%% Plot synapses, prediction, and functional ROIs

global tis glob

if 0
    SPN =  [glob.datDir 'Analysis\Data\preproc\'];
    
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'NOI.mat']);
    load([SPN 'MOI.mat']);
    load('MPN.mat')
    swcDir = [WPN 'swc\'];
end

cmap = jet(100);
lc = 10


roiCid = ptDat(:,3);
runCids = 2; %[2 3 4 13];%MOI.cids;
numCid = length(runCids);
sp = ceil(sqrt(numCid));


for i = 1:length(runCids);
    
    cid = runCids(i);
    load(sprintf('%snrn_cid%d.mat',swcDir,cid))
    
    v = find(NOI.cids==cid,1);
    nep = NOI.neps(v).nep;
    syn = NOI.syns(v).syn;
    preSign = NOI.preSigns(v).preSign;
    d = NOI.ds(v).d;
    
    W = exp(-d/lc); % Apply length constant
    
    Won = W .* repmat(preSign==2,[1 size(W,2)]);
    Woff = W .* repmat(preSign==1,[1 size(W,2)]);
    Winhib = W.* repmat(preSign==-1,[1 size(W,2)]);
    
    mv = find(MOI.cids == cid,1);
    sumOn = MOI.c(mv).mOn;
    sumOff = MOI.c(mv).mOff;
    sumInhib = MOI.c(mv).mUnk;
    
    synPos = syn.pos;
    onSynPos = synPos(preSign == 2,:);
    offSynPos = synPos(preSign == 1,:);
    inhibSynPos = synPos(preSign == -1,:);
    
    colormap(jet(100))
    sumOnN = sumOn * 50/mean(sumOn(:));
    sumOffN = sumOff * 50/mean(sumOff(:));
    sumInhibN = sumInhib * 50/mean(sumInhib(:));
    
    maxI = max([sumOn(:); sumOff(:)]);
    col2 = [sumOff(:)/max(sumOff(:))  sumOn(:) * 0 sumOn(:)/max(sumOn(:))];
    
    offBias = (sumOff-sumOn)./(sumOn+sumOff);
    
    pos = nep.pos;
    
    onCol = sumOn;
    onCol = onCol - mean(onCol(:));
    onCol = onCol/std(onCol(:));
    onCol = round(onCol * 25 + 50);
    onCol(onCol<1) = 1;
    onCol(onCol>100) = 100;
    onCol = cmap(onCol,:);
    
    offBiasCol = round(offBias * 50 + 50);
    offBiasCol(offBiasCol<1) = 1;
    offBiasCol(offBiasCol>100) = 1;
    offBiasCol = cmap(offBiasCol,:);
    
    %% roi
    
    useRois = find(roiCid == cid);
    roiPos = SOI.pos(useRois,:);
    roiPol = ROI.Polarity(useRois);
    roiSN = ROI.SNR(useRois);
    roiPolCol = round(roiPol * 50 + 50);
    roiPolCol(roiPolCol<1) = 1;
    roiPolCol(roiPolCol>100) = 1;
    roiPolCol = cmap(roiPolCol,:);
    
    %% show
    markerSize = 150;
    
    %subplot(sp,sp,i)
    %subplot(1,4,i)
    
    scatSkel = scatter3(pos(:,1),pos(:,2),pos(:,3),'.','clipping','off');
    axis equal
    view([ 0 90])
    scatSkel.CData = offBiasCol;
    hold on
    scatOn = scatter3(onSynPos(:,1),onSynPos(:,2),onSynPos(:,3),markerSize,'k','+','linewidth',2,'clipping','off');
    scatOff = scatter3(offSynPos(:,1),offSynPos(:,2),offSynPos(:,3),markerSize,'k','d','linewidth',2,'clipping','off');
    scatRoi = scatter3(roiPos(:,1),roiPos(:,2),roiPos(:,3),markerSize,'k','o','linewidth',2,'clipping','off');
    scatRoi.MarkerFaceColor = 'flat';
    scatRoi.CData = roiPolCol;
    scatRoi.SizeData = roiSN * 10 * markerSize;
    
    scatRoiE = scatter3(roiPos(:,1),roiPos(:,2),roiPos(:,3),markerSize,'k','o','linewidth',1,'clipping','off');
    scatRoiE.SizeData = roiSN * 10 * markerSize;
    hold off
    ylim([min(pos(:,2)) max(pos(:,2))]);
    xlim([min(pos(:,1)) max(pos(:,1))]);
    zlim([min(pos(:,3)) max(pos(:,3))]);
    title(sprintf('cell %d',cid))
    drawnow
    colormap(jet(100))
    pause(1)
    
    
end

%% Convert um to em


X = 104.892
Y = 129.66
Z = 21.52


xEM = X / 0.004;
yEM = Y / 0.004;
zEM = Z / 0.04;

[yEM xEM zEM]













