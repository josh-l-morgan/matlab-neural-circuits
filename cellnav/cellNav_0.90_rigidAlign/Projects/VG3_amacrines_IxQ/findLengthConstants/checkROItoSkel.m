global tis glob

if 0
    
    SPN =  [glob.datDir 'Analysis\Data\preproc\'];
    
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'NOI.mat']);
end

f = figure

cids = NOI.cids;
cids = cids(1:length(NOI.neps));

onCid = ptDat(:,3);
roiPos = ptDat(:,[7 6 8]);
roiPos(:,[1 2]) = roiPos(:,[1 2]) * 0.004;
roiPos(:,3) = roiPos(:,3) * .04;
roiPlane = ptDat(:,1);
roiPlaneID = unique(roiPlane);
closeNode = SOI.closeNode;

%% Get ROI close positions
closePos = zeros(length(onCid),3);
useROI = zeros(length(onCid),1);
for i = 1:length(onCid)
    cidTarg = find(cids == onCid(i),1);
    if ~isempty(cidTarg)
        nep = NOI.neps(cidTarg).nep;
        if ~isempty(nep)
            useROI(i) = 1;
            closePos(i,:) = nep.pos(closeNode(i),:);
        else
            useROI(i) = 0;
        end
    end
end




%%
clf
hold on
nepCol = hsv(length(cids)) * 0;


for n  = 1:length(NOI.neps)
    
    nep = NOI.neps(n).nep;
    if ~isempty(nep)
    scatter3(nep.pos(:,1),nep.pos(:,2),nep.pos(:,3),10,'o','markerfacecolor',...
        nepCol(n,:),'markeredgealpha',.02,'markerfacealpha',.1)
    else
        disp('nep is epmpty')
    end
end


planeCol = hsv(length(roiPlaneID));
planeCol = planeCol(randperm(size(planeCol,1)),:);
for p = 1:length(roiPlaneID)
    
    isPlane = find((roiPlane == roiPlaneID(p)) & (useROI));
    scatter3(closePos(isPlane,1),closePos(isPlane,2),closePos(isPlane,3),...
        150,'o','markerfacecolor',planeCol(p,:),'markeredgecolor','k',...
        'markeredgealpha',1,'markerfacealpha',.4)
    
    isPlane = find((roiPlane == roiPlaneID(p)) & (useROI+1));
    scatter3(roiPos(isPlane,1),roiPos(isPlane,2),roiPos(isPlane,3),...
        40,'d','markerfacecolor',planeCol(p,:),'markeredgecolor','k',...
        'markeredgealpha',1,'markerfacealpha',.4)

    
end

























