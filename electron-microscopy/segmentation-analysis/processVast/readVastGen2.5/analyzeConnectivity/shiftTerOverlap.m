
%% load data
clear all
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\'
TPN = [MPN 'manyTerSubs\'];
load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])

downSamp = 1;
lookDist = 0;
shiftUm = [-50:1:50];
shiftDim = 1;


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

%postList = ([ 108 129	109	117	162	131	116	137	130	135	106]);


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

tic
clear shiftTer

fileName = sprintf('terSubs_Ds%d_Ds%d_Look%d.mat',8,downSamp,lookDist)

%%fileName = sprintf('terSubs_Ds%d_Ds%d_Look%d.mat',8,downSamp,lookDist)
load([TPN fileName])
%}

shiftTer = terSubs;

preList = shiftTer.preList;
postList = shiftTer.postList;
dsMaxSub = maxSub/downSamp+1;
synMat = zeros(length(preList),length(postList));


%%

viewVol = zeros(dsMaxSub,'uint8');
shiftVec = round(shiftUm/( .03 * 8 * downSamp));
shiftTer.shiftVec = shiftVec;
shiftTer.shiftUm = shiftUm;
shiftTer.shiftDim = shiftDim;

for s = 1:length(shiftVec);
    disp(sprintf('running shift %d of %d',s,length(shiftVec)))
    
    for pre = 1 : length(preList)
        
        %disp(sprintf('running pre %d of %d',pre,length(preList)))
        
        diSub1 = shiftTer.cells(find(shiftTer.cellList == preList(pre))).terSubs;
        diSub1(:,shiftTer.shiftDim) = diSub1(:,shiftTer.shiftDim)+ shiftVec(s);
        useSub = (diSub1(:,1)>0) & (diSub1(:,1)<=dsMaxSub(1)) & ...
            (diSub1(:,2)>0) & (diSub1(:,2)<=dsMaxSub(2)) & ...
            (diSub1(:,3)>0) & (diSub1(:,3)<=dsMaxSub(3));
        diSub1 = diSub1(useSub,:);
        
        diInd1 = sub2ind(dsMaxSub,diSub1(:,1),diSub1(:,2),diSub1(:,3));
        
        for post = 1: length(postList)
            
            diSub2 = shiftTer.cells(find(shiftTer.cellList == postList(post))).terSubs;
            %diSub2(:,shiftTer.shiftDim) = diSub2(:,shiftTer.shiftDim)+ shiftVec(s);
            diInd2 = sub2ind(dsMaxSub,diSub2(:,1),diSub2(:,2),diSub2(:,3));
            
            sharedVox = intersect(diInd1,diInd2);
            
            %                         viewVol = viewVol * 0;
            %                         viewVol(diInd2) = 1;
            %                         viewVol(diInd1) = 2;
            %                         viewVol(diInd2) = 1;
            %                         viewVol(sharedVox) = 5;
            %                         image(squeeze(max(viewVol,[],1))*50);
            %                         pause(.01)
            %
            shiftMat(pre,post,s) = length(sharedVox);
            synMat(pre,post) = sum(synPre==preList(pre) & synPost == postList(post));
            
        end
        %synMat(pre,post) = length(uniqueInd);  % how many unique synapse positions are there
        %}
    end
    
end

%%
shiftTer.synMat = synMat;
shiftTer.shiftMat = shiftMat;
% shiftTer.histMat = histMat;
% shiftTer.minTouch = minTouch;



%%
%{
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
        
        shiftTer.fitDat = fitDat;
        
        
        %% Show syn vs length
        
        colMat = synMat*256/max(synMat(:));
        colMat(:,:,2) = touchMat * 256/max(touchMat(:));
        colMat(:,:,3) = touchMat * 0;
        
        %image(uint8(colMat))
        
        %%
        
        %save([TPN fileName],'shiftTer')
        toc
%}
writeName = sprintf('shiftTer_b_Ds%d_Ds%d_Look%d_Dim%d.mat',8,downSamp,lookDist,shiftDim)

save([TPN writeName],'shiftTer');




