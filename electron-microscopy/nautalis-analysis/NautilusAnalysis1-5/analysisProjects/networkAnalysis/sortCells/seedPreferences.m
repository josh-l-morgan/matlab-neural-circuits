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



%% struct

conPref.seedCells = seedCells;
conPref.axList = axList;
conPref.cellList = cellList;
conPref.con = con;
conPref.ax2seed = ax2seed;

conPref.covMat = covMat;

conPref.seedCov = seedCov;
conPref.zcov = zcov;
conPref.geoMean = geoMean;
conPref.sharedAx =sharedAx;
conPref.sharedSyn = sharedSyn;

conPref.geoMeanNorm = geoMeanNorm;
conPref.sharedAxNorm = sharedAxNorm;
conPref.sharedSynNorm = sharedSynNorm;
conPref.zcovNorm = zcovNorm;

conPref.cellPrefNoExclusion = geoMean;
