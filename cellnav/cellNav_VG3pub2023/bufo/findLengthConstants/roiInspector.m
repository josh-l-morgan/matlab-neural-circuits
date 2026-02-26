%% Compare each ROI to surrounding synapses
%% display cell with ROI highlighted and color coded by polarity.  
%% color code synapses by influence/distance, marker type shows synapse type

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
lc = 14


roiCid = ptDat(:,3);
runCids = [2 3 4 13];%MOI.cids;
numCid = length(runCids);
sp = ceil(sqrt(numCid));


for i = 1:length(runCids);
    
    cid = runCids(i);
    %load(sprintf('%snrn_cid%d.mat',swcDir,cid))
    
    
    v = find(NOI.cids==cid,1);
    nep = NOI.neps(v).nep;
    syn = NOI.syns(v).syn;
    preSign = NOI.preSigns(v).preSign;
    d = NOI.ds(v).d;
    
    W = exp(-d/lc); % Apply length constant
    Won = W .* repmat(preSign==2,[1 size(W,2)]);
    Woff = W .* repmat(preSign==1,[1 size(W,2)]);
    Winhib = W.* repmat(preSign==-1,[1 size(W,2)]);
    
    
    if 1
        sumOn = sum(Won,1);
        sumOff = sum(Woff,1);
        sumInhib = sum(Winhib,1);
    else
        mv = find(MOI.cids == cid,1);
        sumOn = MOI.c(mv).mOn;
        sumOff = MOI.c(mv).mOff;
        sumInhib = MOI.c(mv).mUnk;
    end
    
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
    closeNodes = SOI.closeNode(useRois);

    roiPolCol = round(roiPol * 50 + 50);
    roiPolCol(roiPolCol<1) = 1;
    roiPolCol(roiPolCol>100) = 1;
    roiPolCol = cmap(roiPolCol,:);
    
%     %% test ROI pos
%     scatSkel = scatter3(pos(:,1),pos(:,2),pos(:,3),3,'.','clipping','off');
%     hold on
%     scatSkel = scatter3(pos(SOI.closeNode(useRois),1),pos(SOI.closeNode(useRois),2),pos(SOI.closeNode(useRois),3),30,'x','clipping','off');
%     scatSkel = scatter3(roiPos(useRois,1),roiPos(useRois,2),roiPos(useRois,3),30,'o','clipping','off');
%     hold off
    
    %% show
    markerSize = 300;

    %%Run each Roi
    for r = 1:size(roiPos,1)
        
        subplot(1,1,1)
        closeNode = closeNodes(r); %skeleton node to use for ROI
        synDist = W(:,closeNode);
        synDistInd = round(synDist*100);
        synDistInd(synDistInd<1) = 1;
        synDistInd(synDistInd>100) = 100;
        synCol = cmap(synDistInd,:);
        %subplot(sp,sp,i)
        %subplot(1,4,i)
        cla
        scatSkel = scatter3(pos(:,1),pos(:,2),pos(:,3),3,'.','clipping','off');
        axis equal
        view([ 0 0])
        scatSkel.CData = [0,0,0];%offBiasCol;
        hold on

        scatSkel = scatter3(pos(closeNode,1),pos(closeNode,2),pos(closeNode,3),300,'x','clipping','off');       
        
        
        scatOn = scatter3(onSynPos(:,1),onSynPos(:,2),onSynPos(:,3),markerSize,[0 0 1],'o','linewidth',2,'clipping','off');
        %scatOn.CData = synCol(preSign == 2,:);
        scatOn.SizeData = synDist(preSign == 2)*markerSize;
        scatOn.MarkerFaceColor = 'flat';        
        scatOnE = scatter3(onSynPos(:,1),onSynPos(:,2),onSynPos(:,3),markerSize,'k','p','linewidth',.5,'clipping','off');
        scatOnE.SizeData = synDist(preSign == 2)*markerSize;
        
        scatOff = scatter3(offSynPos(:,1),offSynPos(:,2),offSynPos(:,3),markerSize,[1 0 0],'o','linewidth',2,'clipping','off');
        %scatOff.CData = synCol(preSign == 1,:);
        scatOff.SizeData = synDist(preSign == 1)*markerSize;
        scatOff.MarkerFaceColor = 'flat';
        scatOffE = scatter3(offSynPos(:,1),offSynPos(:,2),offSynPos(:,3),markerSize,'k','o','linewidth',.5,'clipping','off');
        scatOffE.SizeData = synDist(preSign == 1)*markerSize;
        
        scatRoi = scatter3(roiPos(r,1),roiPos(r,2),roiPos(r,3),markerSize * 20,'k','s','linewidth',2,'clipping','off');
        scatRoi.MarkerFaceColor = 'flat';
        scatRoi.CData = offBiasCol(closeNode,:);    
        scatRoi.MarkerFaceAlpha = .2;    
       
        scatRoi2 = scatter3(roiPos(r,1),roiPos(r,2),roiPos(r,3),markerSize * 5,'k','s','linewidth',2,'clipping','off');
        scatRoi2.MarkerFaceColor = 'flat';
        scatRoi2.CData = roiPolCol(r,:);
        scatRoi2.MarkerFaceAlpha = .2;
        
        
     %   scatRoi.SizeData = roiSN * 10 * markerSize;
        
%         scatRoiE = scatter3(roiPos(r,1),roiPos(r,2),roiPos(r,3),markerSize,'k','+','linewidth',1,'clipping','off');
%         scatRoiE.SizeData = roiSN * 10 * markerSize;
        hold off
        ylim([min(pos(:,2)) max(pos(:,2))]);
        xlim([min(pos(:,1)) max(pos(:,1))]);
        zlim([min(pos(:,3)) max(pos(:,3))]);
        title(sprintf('cell %d, roi %d, OffPol %0.2f, Pred %0.2f, SNR %0.2f',...
            cid,useRois(r),roiPol(r),offBias(closeNode),roiSN(r)*100))
        colormap(jet(100))
        view([0 90])
        drawnow
        
        pause
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













