
%% load data
clear all
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\'
TPN = [MPN 'manyTerSubs\'];
load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])


downSamps = [ 1 1 1 1 1 2 2 4 4 8 8 8 8  16 32 64]
lookDists = [ 0 1 2 3 4 2 4 2 4 2 4 8 16 16 16 16]




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

postList = ([108  129	109	117	162	131	116	137	130	135	106]);
postList = ([129	109	117	162	131	116	137	130	135	106]);


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
% postList = unique(synPost);
% postList = postList(postList>0);
preList = preWithTarg;
bigList = [ 1001 1006 1009 1012 1014 1021 1023 1025 1027 1028 1029 1030 1031 1032 1033 1034 1036 1037 1041 1050 1051 1053 1054 1055 1056 1058]
preList = intersect(preList,bigList);


goodCheck = synObj*0 + 1;
showPlots = 0;
synMat = zeros(length(preList),length(postList));


for s = 1:length(lookDists)
    tic
    disp(sprintf('running %d of %d',s,length(lookDists)))
    clear terSubs
    downSamp = downSamps(s);
    lookDist = lookDists(s);
    fileName = sprintf('terSubs_Ds%d_Ds%d_Look%d.mat',8,downSamp,lookDist)
    
    if ~exist([TPN fileName])
        
        terSubs.targCell  = targCell;
        terSubs.downSamp1 = 8;
        terSubs.downSamp2 = downSamp;
        terSubs.lookDist = lookDist;
        
        
        terSubs.preList = preList;
        terSubs.postList = postList;
        
        %% make ball filter
        ballSize = [lookDist*2+1  lookDist*2+1 lookDist*2+1];
        ball = strel('ball',lookDist,lookDist);
        ball = ones(ballSize);
        ballInd = find(ball);
        [bY bX bZ] = ind2sub(ballSize,ballInd);
        dists = sqrt((bY - lookDist - 1).^2 + (bX - lookDist - 1).^2 + (bZ - lookDist - 1).^2);
        ball(dists>lookDist) = 0;
        %ball(dists<lookDist-1.8) = 0;
        colormap gray(256)
        
        image(sum(ball,3)*10)
        
        %%
        cellList = unique([preList(:); postList(:)]);
        terSubs.cellList = cellList;
        terSubs.lookDist = lookDist;
        terSubs.downSamp = downSamp;
        dsMaxSub = maxSub/downSamp+1;
        for cellNum = 1: length(cellList)
            mid2 = getCellSubs(obI,dsObj,cellList(cellNum));
            mid2(:,1:2) = mid2(:,1:2) + 1024/8;
            dsSubs = double(downSampSubs(mid2,downSamp));
            
            %%lesmid
            %     maxMid = max(mid2,[],1);
            %     midInd = sub2ind(maxMid,mid2(:,1),mid2(:,2),mid2(:,3));
            %     midInd = unique(midInd);
            %     mid3 = zeros(length(midInd),3);
            %     [mid3(:,1) mid3(:,2) mid3(:,3)] = ind2sub(maxMid,midInd);
            conMat = obj2con(dsSubs);
            getSurf = sum(conMat>0,2)<23;
            mid4 = dsSubs(getSurf,:);
            
            tic
            diSub2 = dilateSub(mid4,ball);
            pause(.01)
            toc
            
            
            diSub2(diSub2(:,1)>dsMaxSub(1),1) = dsMaxSub(1);
            diSub2(diSub2(:,2)>dsMaxSub(2),2) = dsMaxSub(2);
            diSub2(diSub2(:,3)>dsMaxSub(3),3) = dsMaxSub(3);
            
            terSubs.cells(cellNum).terSubs = diSub2;
            
        end
        
        %%
        viewVol = zeros(dsMaxSub,'uint8');
        
        for pre = 1 : length(preList)
            
            %disp(sprintf('running pre %d of %d',pre,length(preList)))
            
            diSub1 = terSubs.cells(find(terSubs.cellList == preList(pre))).terSubs;
            diInd1 = sub2ind(dsMaxSub,diSub1(:,1),diSub1(:,2),diSub1(:,3));
            
            for post = 1: length(postList)
                
                diSub2 = terSubs.cells(find(terSubs.cellList == postList(post))).terSubs;
                diInd2 = sub2ind(dsMaxSub,diSub2(:,1),diSub2(:,2),diSub2(:,3));
                
                sharedVox = intersect(diInd1,diInd2);
                
                %         viewVol = viewVol * 0;
                %         viewVol(diInd2) = 1;
                %         viewVol(diInd1) = 2;
                %         viewVol(diInd2) = 1;
                %         viewVol(sharedVox) = 5;
                %         image(squeeze(max(viewVol,[],1))*50);
                %         pause(.01)
                
                touchMat(pre,post) = length(sharedVox);
                synMat(pre,post) = sum(synPre==preList(pre) & synPost == postList(post));
                
            end
            %synMat(pre,post) = length(uniqueInd);  % how many unique synapse positions are there
            %}
        end
        
        %%
        terSubs.synMat = synMat;
        terSubs.touchMat = touchMat;
        % terSubs.histMat = histMat;
        % terSubs.minTouch = minTouch;
        
        
        
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
        
        terSubs.fitDat = fitDat;
        
        
        %% Show syn vs length
        
        colMat = synMat*256/max(synMat(:));
        colMat(:,:,2) = touchMat * 256/max(touchMat(:));
        colMat(:,:,3) = touchMat * 0;
        
        %image(uint8(colMat))
        
        %%
        
        save([TPN fileName],'terSubs')
        toc
    end
end
