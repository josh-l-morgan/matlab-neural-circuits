function[] = makeGOI()

%% group ROIs according to proximity 
%%Use itterative algorithm that pools rois les than X from one another and
%%create 'heavier' averaged roi possition
%%M

global tis glob
SPN =  [glob.datDir 'Analysis\Data\preproc\'];

if 1
    SPN =  [glob.datDir 'Analysis\Data\preproc\'];
    
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'NOI.mat']);
    load([SPN 'MOI.mat']);
    load('MPN.mat')
    swcDir = [WPN 'swc\'];
    load([SPN 'RawOI.mat']);

end

poolDist = 1;
frames = ptDat(:,1);

roiCid = ptDat(:,3);
runCids = unique(roiCid);%MOI.cids;
numCid = length(runCids);

posCal = ptDat(:,[7 6 8]);
posCal(:,1) = posCal(:,1) * 0.004;
posCal(:,2) = posCal(:,2) * 0.004;
posCal(:,3) = posCal(:,3) * .04;

smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];

clear GOI
GOI.cids = runCids; % create structure containing grouped rois
GOI.poolDist = poolDist; 

for i = 1:length(runCids);
    cid = runCids(i);
    
    disp(sprintf('grouping rois for cell %d. %d of %d',cid,i,length(runCids)))
    
    %load(sprintf('%snrn_cid%d.mat',swcDir,cid))
    
    %%Get distances between nodes
  
    useSM(i) = 1;
    if exist('sms','var')
    sm = sms(i).sm;
    else
        fileName = sprintf('sm_cid%d.mat',cid);
        load([smDir fileName]);
    end
    dAll = sm.skel2skel.linDist;
    pos = sm.nep.pos;
       
    %%Get Rois
    useRois = find(roiCid == cid);
    roiPos = SOI.pos(useRois,:);
    roiPol = ROI.Polarity(useRois);
    roiSN = ROI.SNR(useRois);
    closeNodes = SOI.closeNode(useRois);
    
    
    %%Build grouped rois
    gPos = posCal(useRois,:);
    gClose = closeNodes;
    clear gRoiID
    for r = 1:length(closeNodes)
        gRoiID{r} = useRois(r);
    end
    gWeight = ones(length(closeNodes),1);
    
    %%Show what is going on
    scatter(pos(:,1),pos(:,2),'.','k')
    hold on
    scatter(posCal(useRois,1),posCal(useRois,2),100,'g','linewidth',3)
    scatter(gPos(:,1),gPos(:,2),200,'r','linewidth',3)
    pause(.1)
    hold off
    
    %%Run itterative grouping of rois on one cell
    for r = 1: length(closeNodes)
        
        %%Get distance between nodes
        d = dAll(gClose,gClose);
        d(d==0) = inf;
        minD = min(d(:));
        if minD <= poolDist
            [y x] = find(d==minD,1);
            
            pos1 = gPos(y,:) * gWeight(y);
            pos2 = gPos(x,:) * gWeight(x);
            newWeight = (gWeight(y) + gWeight(x));
            newPos = (pos1 + pos2)/newWeight;
            newFrame = frames(x);
            
            newIDs = cat(2,gRoiID{y},gRoiID{x});
            
            nDist = sqrt((pos(:,1)-newPos(1)).^2 + (pos(:,2)-newPos(2)).^2 + (pos(:,3)-newPos(3)).^2);
            closest = find(nDist == min(nDist),1);
            
            %%Make new group
            gPos(y,:) = newPos;
            gClose(y) = closest;
            gRoiID{y} = newIDs;
            gWeight(y) = newWeight;
            gFrame{y} = newFrame;
            
            %Cut old group
            useG = setdiff([1:length(gClose)],x);
            gPos = gPos(useG,:);
            gClose = gClose(useG);
            gRoiID = gRoiID(useG);
            gWeight = gWeight(useG);
            
            scatter(pos(:,1),pos(:,2),'.','k')
            hold on
            scatter(posCal(useRois,1),posCal(useRois,2),100,'g','linewidth',3)
            scatter(gPos(:,1),gPos(:,2),200,'r','linewidth',3)
            pause(.1)
            hold off
            
            
        else % leave loop if no roi nodes are close together
            break
        end
    end
    
    GOI.c(i).roiID = gRoiID;
    GOI.c(i).closeNode = gClose;
    GOI.c(i).pos = gPos;
    GOI.c(i).weight = gWeight;
    
end

save([SPN 'GOI.mat'],'GOI');

%% Combine physiologies

gNum = 0;
GOI.roiRaw = ROI;
for c = 1:length(GOI.c) %each cid/cell
   for r = 1:length(GOI.c(c).roiID)
      ids = GOI.c(c).roiID{r};
      GOI.c(c).Polarity(r) = mean(ROI.Polarity(ids));
      GOI.c(c).Polarity1(r) = mean(ROI.Polarity1(ids));
      GOI.c(c).Polarity2(r) = mean(ROI.Polarity2(ids));
            GOI.c(c).Polarity3(r) = mean(ROI.Polarity3(ids));
      GOI.c(c).Polarity4(r) = mean(ROI.Polarity4(ids));

      GOI.c(c).SNR(r) = sum(ROI.SNR(ids));
      GOI.c(c).polVar(r) = var(ROI.Polarity(ids));
      GOI.c(c).frames{r} = frames(ids);
      GOI.c(c).frame(r) = mode(frames(ids));
      GOI.c(c).varPol(r) = mean(ROI.varPol(ids));
      GOI.c(c).qual(r) = mean(ROI.qual(ids));

      gNum = gNum + 1;
      
      GOI.Polarity(gNum,1) = GOI.c(c).Polarity(r);
      GOI.Polarity1(gNum,1) = GOI.c(c).Polarity1(r);
      GOI.Polarity2(gNum,1) = GOI.c(c).Polarity2(r);
      GOI.Polarity3(gNum,1) = GOI.c(c).Polarity3(r);
      GOI.Polarity4(gNum,1) = GOI.c(c).Polarity4(r);
      GOI.varPol(gNum,1) = GOI.c(c).varPol(r);
      GOI.qual(gNum,1) = GOI.c(c).qual(r);
      GOI.roiCids(gNum,1) = GOI.cids(c);
      GOI.SNR(gNum,1) = GOI.c(c).SNR(r);
      GOI.closeNode(gNum,1) = GOI.c(c).closeNode(r);
      GOI.pos(gNum,:) = GOI.c(c).pos(r,:);
      GOI.weight(gNum,1) = GOI.c(c).weight(r);
      GOI.roiID{gNum,1} = GOI.c(c).roiID{r};
      GOI.frame(gNum,1) = GOI.c(c).frame(r);
      GOI.frames{gNum,1} = GOI.c(c).frames{r};


   end    
end


for y = 1:size(RawOI.PixResp,1)
    for x = 1:size(RawOI.PixResp,2)
        useRois = RawOI.useRois{y,x};
        GOI.PixResp{y,x} = [];
        GOI.useRois{y,x} = [];
        GOI.rawPixResp{y,x} = [];
        gf = 0;
        for g = 1:length(GOI.roiID)
            groupRoi = GOI.roiID{g};
            [a hit] = intersect(useRois,groupRoi);
            if length(hit)
                gf = gf+1;
                getResp = RawOI.PixResp{y,x}(:,:,:,hit);
                GOI.PixResp{y,x}(:,:,:,gf) = mean(getResp,4);
                GOI.useRois{y,x}(gf) = g;
                GOI.rawPixResp{y,x}(gf,:) = mean(RawOI.rawPixResp{y,x}(hit,:),1);

            end
        end
    end
end





save([SPN 'GOI.mat'],'GOI');











