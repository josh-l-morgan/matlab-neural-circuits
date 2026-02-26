

clear all
MPN = GetMyDir;
load([MPN 'obI.mat']);
allSeeds = [ 108 201 109 903 907];
%allSeeds = [108 201 109]
%seedList = [ 108  201 109 ];

useList = obI2cellList_seedInput_RGC_TCR(obI,allSeeds);

%%
clf
scatCol={'r','g','b','b','b'};
%scatCol={[1 0 0],[0 1 0], [ 0 0 1],[ 0 0 1],[ 0 0 1],[0 .7 .7],[.7 0 .7]};
scatSize = [ 10 10 20 20 20 20 20 20 ];
clear uL histAxSyn histAxEdge histAxMean histPostSyn histPostEdge histPostMean

for i = 1:length(allSeeds)
    seedList = allSeeds(i);
    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    targ = find(useList.postList == useList.seedList(1))
    con = useList.con;
    withSeed = con(:,targ);
    [sortWithSeed idx] = sort(withSeed,'descend');
    preList = useList.preList(idx);
    con = con(idx,:);
    sumPost = sum(con,1);
    [sortSumPost idx] = sort(sumPost,'descend')
    con = con(:,idx);
    postList = useList.postList(idx);
    %subplot(length(allSeeds),1,i)
    %image(con*10)
    usePre = sum(con,2)>4;
    usePost = sum(con,1)>4;
    useCon = (repmat(usePre,[1 size(con,2)]) & repmat(usePost,[size(con,1) 1]));
   

    clear allGAT allGNT allMean
    targ = postList == seedList;
    for a = 1:size(con,1)
        
        repCon = repmat(con(a,:),[size(con,1) 1]);
        geoMean = sqrt(con.*repCon);
        allMean(a,:,:) = geoMean;  
        
        geoAtTarg = geoMean(:,find(targ));
        sumGeo = sum(geoMean,2);
        countGeo = sum(geoMean>0,2)-1;
        geoNotTarg = sumGeo-geoAtTarg;
        geoNatTarg = geoNotTarg./countGeo;
        
        allGAT(:,a) = geoAtTarg;
        allGNT(:,a) = geoNotTarg;
        

    end
    
    fullGeoTarg = allMean(:,:,targ);
    fullGeoNotTarg = allMean(:,:,setdiff(1:size(allMean,3),targ));
    repFullGeoTarg = repmat(fullGeoTarg,[1 1 size(fullGeoNotTarg,3)]);
    useG = ~diag(ones(size(con,1),1),0);
    
    repUseG = repmat(useG,[1 1 size(fullGeoNotTarg,3)]);
    
    filtCon = con;
    axSyn = sum(con,2);
    filtCon = con(axSyn>0,:);
    synBin = [0:3:30];
    edgeBin = [0:3:30];
    meanBin = [0:1.5:15];
    
    axSyn = sum(filtCon,2);
    axEdge  = sum(filtCon>0,2);
    axMean = axSyn./axEdge;
   
    subplot(2,1,1)
    sg = scatter(axSyn,axEdge,scatCol{i})
    %set(sg,'markerFacerColor',scatCol{i});
    hold on
    
    postSyn = sum(filtCon,1);
    postEdge = sum(filtCon>0,1);
    postMean = postSyn./postEdge;
    
    subplot(2,1,2)
    scatter(postSyn,postEdge,scatCol{i});
    hold on
    
    
    
    
    
    histAxSyn(:,i) = hist(axSyn,synBin);
    histAxEdge(:,i) = hist(axEdge,edgeBin);
    histAxMean(:,i) = hist(axMean,meanBin);
    
    histPostSyn(:,i) = hist(postSyn,synBin);
    histPostEdge(:,i) = hist(postEdge,edgeBin);
    histPostMean(:,i) = hist(postMean,meanBin);
    
    
end
groupHist = {[1] [2] [3 4 5]};
for g = 1:length(groupHist)
    histAxSy(:,g) = sum(histAxSyn(:,groupHist{g}),2);
    histAxEdge(:,g) = sum(histAxEdge(:,groupHist{g}),2);
    histAxMean(:,g) = sum(histAxMean(:,groupHist{g}),2);
    histPostSyn(:,g) = sum(histPostSyn(:,groupHist{g}),2);
    histPostEdge(:,g) = sum(histPostEdge(:,groupHist{g}),2);
    histPostMean(:,g) = sum(histPostMean(:,groupHist{g}),2);

end


%{

subplot(3,2,1) 
bar(synBin,histAxSyn(:,1:length(groupHist)))

subplot(3,2,2)
bar(synBin,histPostSyn(:,1:length(groupHist)))

subplot(3,2,3)
bar(edgeBin,histAxEdge(:,1:length(groupHist)))

subplot(3,2,4)
bar(edgeBin,histPostEdge(:,1:length(groupHist)))

subplot(3,2,5)
bar(meanBin,histAxMean(:,1:length(groupHist)))

subplot(3,2,6)
bar(meanBin,histPostMean(:,1:length(groupHist)))
%}



% xlim([0 26])
% ylim([0 26])
hold off











