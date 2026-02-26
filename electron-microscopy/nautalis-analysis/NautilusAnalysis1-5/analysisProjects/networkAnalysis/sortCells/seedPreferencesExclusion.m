function[conPref] = seedPreferences(seedCells, useList);

axList = useList.preList;
cellList = useList.postList;
%cellList = [cellList setdiff(seedCells,cellList)];
con = useList.con;



%% find cell prefs without exclusions
covMat= cov(con);

for s = 1:length(seedCells)
    targ = find(cellList==seedCells(s));
    targMat = repmat(con(:,targ),[1 size(con,2)]);
    ax2seed(s,:) = con(:,targ)';
    
    covarMat = con.*targMat;
    geometricMeanMat = sqrt(covarMat);
    sharedAxMat  = (con>0) & (targMat>0);
    sharedSynMat  = con .* (targMat>0);
    zcovMat = con.* targMat;
    
    seedCov(s,:) = covMat(targ,:);
    
    geoMean(s,:) = sum(geometricMeanMat,1);
    sharedAx(s,:) = sum(sharedAxMat,1);
    sharedSyn(s,:) = sum(sharedSynMat,1);
    zcov(s,:) = mean(zcovMat,1);
    
    geoMeanNorm(s,:) = geoMean(s,:)./max(geoMean(s,:));
    sharedAxNorm(s,:) = sharedAx(s,:)./max(sharedAx(s,:));
    sharedSynNorm(s,:) = sharedSyn(s,:)./max(sharedSyn(s,:));
    zcovNorm(s,:) = zcov(s,:)./max(zcov(s,:));
    %
    %     scatter(sharedAx(s,:),sharedSyn(s,:),'r');
    %     hold on
    %     scatter(sharedAx(s,:),geoMean(s,:),'b');
    %     hold off
    %
    
    %image(geoMean*10)
end



%bar(cellPref')

%% Calculate what the preference of each post would be, absent a particular pre

postPrefMinusPre = zeros(size(con,1),size(con,2),2);
for i = 1: length(axList)
    usePre = setdiff(1:length(axList),i);
    sampCon = con(usePre,:);
    for s = 1:length(seedCells)
        targ = find(cellList==seedCells(s));
        simMat = sqrt(sampCon .* repmat(sampCon(:,targ),[1 size(con,2)]));
        postPrefMinusPre(i,:,s) = sum(simMat,1);
    end
end
% 
% subplot(2,2,1)
% image(postPrefMinusPre(:,:,1));
% subplot(2,2,2)
% image(postPrefMinusPre(:,:,2));

%% find ax prefrence given exclusion of post cell
clear sumMat
prePrefFromPost = zeros(size(con,1),size(con,2),2);
%normCellPref = cellPref./repmat(sum(cellPref,1),[2 1]);
for i = 1:length(axList)
    for p = 1:length(cellList)
        usePost = setdiff(1:length(cellList),p);
        sampPref = squeeze(postPrefMinusPre(i,usePost,:));
        sampCon = con(i,usePost);
        for s = 1:length(seedCells)
            simMat = sqrt(sampCon .* sampPref(:,s)');
            prePrefFromPost(i,p,s) = sum(simMat);
        end
    end
end

% subplot(2,2,3)
% image(prePrefFromPost(:,:,1));
% subplot(2,2,4)
% image(prePrefFromPost(:,:,2));

%% struct

conPref.seedCells = seedCells;
conPref.axList = axList;
conPref.cellList = cellList;
conPref.con = con;
conPref.ax2seed = ax2seed;

conPref.covMat = covMat;

conPref.seedCov = seedCov
conPref.zcov = zcov;
conPref.geoMean = geoMean;
conPref.sharedAx =sharedAx;
conPref.sharedSyn = sharedSyn;

conPref.geoMeanNorm = geoMeanNorm;
conPref.sharedAxNorm = sharedAxNorm;
conPref.sharedSynNorm = sharedSynNorm;
conPref.zcovNorm = zcovNorm;


conPref.cellPrefNoExclusion = geoMean;
conPref.cellPrefAxExcluded =  postPrefMinusPre;
conPref.axPrefCellExcluded = prePrefFromPost;
