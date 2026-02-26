clear all
SPN = 'D:\LGNs1\Segmentation\VAST\S8\S8_ds16\CellBodies\matOut_AutoCB\';
load([SPN 'cbObj.mat'])



%% load data
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\'
TPN = [MPN 'cbVSnetwork\'];
load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])


downSamp = 2;


%% load paths
mainDir = 'C:\Users\joshm\Documents\MATLAB\jlm_Code\EM\SegmentationAnalysis\processVast\readVastGen2.1\'
%mainDir = pwd;
slash = regexp(mainDir,'\');
mainDir = mainDir(1:slash(end));
path(path,mainDir)
path(path,sprintf('%sanalyzeMorphology\\skeletonizeByShortestDistance',mainDir))
path(path,sprintf('%sanalyzeMorphology',mainDir))
path(path,sprintf('%sviewCells',mainDir))
path(path,sprintf('%sanalyzeColors',mainDir))
path(path,sprintf('%sanalyzeConnectivity',mainDir))
path(path,sprintf('%svast2ob',mainDir))




%% Define variables
allSubs = cat(1,dsObj.subs);
maxSub = max(allSubs,[],1)+10;
clear allSubs



%% Get syn data and select synapses
targCell = 108;
synapses = obI.nameProps.edges;
preWithTarg = unique(synapses(synapses(:,1)==targCell,2));
preWithTarg = preWithTarg(preWithTarg>0);
synWithPre = [];
for i = 1:length(preWithTarg)
    synWithPre = cat(1,synWithPre, find((synapses(:,2)==preWithTarg(i)) & ...
        (synapses(:,1)>0)));
end

%postList = ([108  129	109	117	162	131	116	137	130	135	106]);


synPre = synapses(synWithPre,2);
synPost = synapses(synWithPre,1);
synObj = synapses(synWithPre,3);
%
% synPos = ceil(obI.colStruc.anchors(synObj,:)/8);
% synPos(:,1:2) = synPos(:,1:2); %%% Error in tif to sub assignment



%% find closest point and total length within 2x dist of closest point
%
% preList = unique(synPre);
% preList = preList(preList>0);
postList = unique(synPost);
postList = postList(postList>0);
postList = postList(postList<1000);
preList = preWithTarg;
%bigList = [ 1001 1006 1009 1012 1014 1021 1023 1025 1027 1028 1029 1030 1031 1032 1033 1034 1036 1037 1041 1050 1051 1053 1054 1055 1056 1058]
%preList = intersect(preList,bigList);


goodCheck = synObj*0 + 1;
showPlots = 0;
synMat = zeros(length(preList),length(postList));


% for s = 1:length(lookDists)
%     tic
%     disp(sprintf('running %d of %d',s,length(lookDists)))
%     clear cb2net
%     downSamp = downSamps(s);
%     lookDist = lookDists(s);
%     fileName = sprintf('terSubs_Ds%d_Ds%d_Look%d.mat',8,downSamp,lookDist)
%     
%     if ~exist([TPN fileName])
        
        cb2net.targCell  = targCell;
        cb2net.downSamp1 = 8;
        cb2net.downSamp2 = downSamp;
        cb2net.preList = preList;
        cb2net.postList = postList;
        
            
        
        %%
        cellList = unique([postList(:) ;preList(:)]);
        cb2net.cellList = cellList;
        dsMaxSub = cbObj.ImageSize;
        for cellNum = 1: length(cellList)
            mid2 = getCellSubs(obI,dsObj,cellList(cellNum));
            if cellList(cellNum) == 106
%             mid2(:,1:2) = mid2(:,1:2) + 1024/8;  %% adjustment for previous error !!!!!!!!
%              mid2(:,1:2) = mid2(:,1:2) + 300;  %% adjustment for previous error !!!!!!!!

            end
            dsSubs = double(downSampSubs(mid2,downSamp));
            
            dsSubs(dsSubs(:,1)>dsMaxSub(1),1) = dsMaxSub(1);
            dsSubs(dsSubs(:,2)>dsMaxSub(2),2) = dsMaxSub(2);
            dsSubs(dsSubs(:,3)>dsMaxSub(3),3) = dsMaxSub(3);
%             
            cb2net.cells(cellNum).dsSubs = dsSubs;
            
        end
        
        %%
        viewVol = zeros(cbObj.ImageSize,'uint8');
        
      
            
            for post = 1: length(postList)
                
                dsSubs = cb2net.cells(find(cb2net.cellList == postList(post))).dsSubs;
                dsInd = sub2ind(cbObj.ImageSize,dsSubs(:,1),dsSubs(:,2),dsSubs(:,3));
                
                subplot(2,1,1)
                cbVals = cbObj.vol(dsInd);
                cbVals = cbVals(cbVals>0);
                uCB = unique(cbVals);
                showOverlap = 0;
                if length(uCB)== 0
                  'missing CB'
                  postCB(post) = 0;
                  showOverlap = 1;
                elseif length(uCB) == 1
                        postCB(post) = uCB;
                else
                    histCB = hist(cbVals,uCB); 
                    postCB(post) = uCB(find(histCB == max(histCB),1));
                end
                
                
                if showOverlap
                    %%temp view combo
                    hist(cbVals)

                    viewVol = viewVol * 0;
                    viewVol(dsInd) = 1;
                    rangeZ = 10;
                    middleZ = mode(dsSubs(:,3));

                    sampCell = double(viewVol(:,:,max(1,middleZ-rangeZ):min(dsMaxSub(3),middleZ + rangeZ)));
                    sampCB =  double(cbObj.vol(:,:,max(1,middleZ-rangeZ):min(dsMaxSub(3),middleZ + rangeZ))>0);

                    colSamp(:,:,1) = sum(sampCell,3)*30;
                    colSamp(:,:,2) = sum(sampCB,3)*10;
                    colSamp(:,:,3) = sum(sampCB,3) * 10;

                    subplot(2,1,2)
                    image(uint8(colSamp))
                    postCB(post)
                    pause(.12)
                end
                %}
               
               % synMat(pre,post) = sum(synPre==preList(pre) & synPost == postList(post));
                
            end
            %synMat(pre,post) = length(uniqueInd);  % how many unique synapse positions are there
            %}
     
            
            %% Get syn mat
            for pre = 1:length(preList)
                for post = 1:length(postList)
                    synMat(pre,post) = sum(synPre==preList(pre) & synPost == postList(post));
                end
            end
            
            
            %% assign cell to cb
            cbID = zeros(1,cbObj.NumObjects);
            cbID(postCB(postCB>0)) = postList(postCB>0);
            
            
            
            cbSynMat = zeros(length(preList),length(cbID));
            for i = 1:length(cbID)
                
                if cbID(i)>0
                    targ = find(postList == cbID(i));
                    cbSynMat(:,i) = synMat(:,targ);
                end
            end
            
            %% Find distances of CB to seed cell
            voxelSize = 0.480;
            targCB = find(cbID == targCell);
            cbProps = regionprops(cbObj,'Centroid');
            cbCenters = cat(1,cbProps.Centroid);
            cbDists = sqrt((cbCenters(:,1)-cbCenters(targCB,1)).^2 + ...
                (cbCenters(:,2)-cbCenters(targCB,2)).^2 + ...
                (cbCenters(:,3)-cbCenters(targCB,3)).^2) * voxelSize;
            
            %% find shared Syn
            targInputs = cbSynMat(:,targCB);
            cbSharedSyn = cbSynMat*0;
            for i = 1:length(cbID);
               cbSharedSyn(:,i) = min(cbSynMat(:,i),targInputs); 
            end
            
            cbSyn = sum(cbSharedSyn,1);
            cbBridge = sum(cbSharedSyn>0,1);
            
            %%
            
              inNetwork = sum(cbSynMat,1)>0;
            
            useCB = inNetwork &  (cbID ~= targCell);
            subplot(2,2,1)
            scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB),1),'.')
                        xlim([0 200])

            subplot(2,2,3)
            scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB)>0,1),'.')
                        xlim([0 200])

            
            useCB = (cbID ~= targCell);
            subplot(2,2,2)
            scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB),1),'.')
                        xlim([0 200])

            subplot(2,2,4)
            scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB)>0,1),'.')
                        xlim([0 200])

        %%
           useCB = (cbID ~= targCell);
            subplot(2,1,1)
            scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB),1),'.')
                        xlim([0 200])
                        ylabel('Shared synapse number')
                        xlabel('Distance ove cell body to origin cell')

            subplot(2,1,2)
            scatter(cbDists(useCB),sum(cbSharedSyn(:,useCB)>0,1),'.')
                        xlim([0 200])
                        xlabel('Distance ove cell body to origin cell')
                        ylabel('Shared axon number')
                        
                        
            
            %% Display distance to synapse relationship
            
            inNetwork = sum(cbSynMat,1)>0;
            
            useCB = inNetwork &  (cbID ~= targCell);
            subplot(2,2,1)
            scatter(cbDists(useCB),sum(cbSynMat(:,useCB),1),'.')
                        xlim([0 200])

            subplot(2,2,3)
            scatter(cbDists(useCB),sum(cbSynMat(:,useCB)>0,1),'.')
                        xlim([0 200])

            
            useCB = (cbID ~= targCell);
            subplot(2,2,2)
            scatter(cbDists(useCB),sum(cbSynMat(:,useCB),1),'.')
                        xlim([0 200])

            subplot(2,2,4)
            scatter(cbDists(useCB),sum(cbSynMat(:,useCB)>0,1),'.')
            
                        xlim([0 200])
            %% Threshold look area
            
            
            maxLook = 50;
            closeCB = find((cbDists'<maxLook) & (cbID ~= targCell));
            
            closeSyn = sum(cbSharedSyn(:,closeCB),1);
            closeBridge = sum(cbSynMat(:,closeCB)>0,1);
            
            
            subplot(2,1,1)
            hist(closeSyn,[0:1:20])
            subplot(2,1,2)
            hist(closeBridge,[0:1:8])
            
            
            %% Draw volume
            subplot(1,1,1)
            closeVol = zeros(cbObj.ImageSize);
            conVol = closeVol;
            targVol = closeVol;
            for i = 1:length(closeCB)
                
                conVol(cbObj.PixelIdxList{closeCB(i)}) = cbSyn(closeCB(i));
                closeVol(cbObj.PixelIdxList{closeCB(i)}) = 1;

            end
            
            targVol(cbObj.PixelIdxList{targCB}) = 1;
            
            viewDim = 1;
            maxConVol = squeeze(max(conVol,[],viewDim));
            maxCloseVol = squeeze(sum(closeVol,viewDim));
            sumTargVol = squeeze(sum(targVol,viewDim));
            maxCloseVol = (maxCloseVol>0) * 20 + maxCloseVol * 2;
            maxCloseVol(maxCloseVol>150) = 150;
            
            colComb = maxCloseVol +maxConVol*10 + sumTargVol*30;
            colComb(:,:,2) = maxCloseVol +maxConVol*100;
            colComb(:,:,3) = maxCloseVol;
            
            image(uint8(colComb))
            
            %%{
            fileName = sprintf('conSphere_look%d_Dim%d.png',maxLook,viewDim)
            imwrite(uint8(colComb),[TPN fileName])
            %}
            
            
            %%  Reference sphere
            
            sumAllCB = squeeze(sum(cbObj.vol>0,viewDim));
            sumAllCB = (sumAllCB>0) * 40 + sumAllCB * 2;
            sumAllCB(sumAllCB>150) = 150;
            
            colComb =  maxCloseVol*2 +maxConVol*10 + sumTargVol*30;
            colComb(:,:,2) = maxCloseVol * 2 +maxConVol*100;
            colComb(:,:,3) = maxCloseVol  + sumAllCB;
            
            image(uint8(colComb))
            
            
            fileName = sprintf('conSphere_Reference_Dim%d.png',maxLook,viewDim)
            imwrite(uint8(colComb),[TPN fileName])
            
            
            
            
            
            
            
            
            
            
            
            
            
            targCell
        %%
        cb2net.synMat = synMat;
%         cb2net.touchMat = touchMat;
        % cb2net.histMat = histMat;
        % cb2net.minTouch = minTouch;
        
        
        
        
        
        %%
        
        useVals = touchMat>=0;
        Xdat = touchMat(useVals);
        Ydat = synMat(useVals);
        [rho pc] = corr(Xdat,Ydat)
        
        scatter(Xdat,Ydat,'.','r')
        hold on
        maxX = max(Xdat);
        X = [0 maxX];
        [fit1 gof] = fit(Xdat,Ydat,'poly1');
        Y = X* fit1.p1 + fit1.p2;
        line(X ,Y);
        hold off
        
        
        pause(2)
        
        fitDat.Xdat = Xdat;
        fitDat.Ydat = Ydat;
        fitDat.fit1 = fit1;
        fitDat.gof = gof;
        fitDat.rho = rho;
        fitDat.corrP = pc;
        
        cb2net.fitDat = fitDat;
        
        
        %% Show syn vs length
        
        colMat = synMat*256/max(synMat(:));
        colMat(:,:,2) = touchMat * 256/max(touchMat(:));
        colMat(:,:,3) = touchMat * 0;
        
        %image(uint8(colMat))
        
        %%
        
        
        toc
        