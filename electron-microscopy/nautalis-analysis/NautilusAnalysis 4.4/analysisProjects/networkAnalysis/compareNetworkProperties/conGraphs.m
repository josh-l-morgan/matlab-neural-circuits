

clear all
MPN = GetMyDir;
load([MPN 'obI.mat']);
allSeeds = [ 108 201 109 903 907];
%allSeeds = [108 201 109]
%seedList = [ 108  201 109 ];

useList = obI2cellList_seedInput_RGC_TCR(obI,allSeeds);

%%
clf
scatCol={'r','g','b','c','m'};
scatCol={[1 0 0],[0 1 0], [ 0 0 1],[ 0 0 1],[ 0 0 1],[0 .7 .7],[.7 0 .7]};
scatSize = [ 10 10 20 20 20 20 20 20 ];
clear uL
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
    
    sg = scatter(repFullGeoTarg(repUseG),fullGeoNotTarg(repUseG),scatSize(i),'o','filled')
    
    
%     subplot(length(allSeeds),3,(i-1)*3+1)
%     image(allGAT*2);
%     subplot(length(allSeeds),3,(i-1)*3+2)
%     image(allGNT*2);
%     subplot(length(allSeeds),3,(i-1)*3+3)
%     scatter(allGNT(:),allGAT(:),'.')
%     
    %sg = scatter(allGAT(useG),allGNT(useG),scatSize(i),'o','filled')
                    set(sg,'MarkerEdgeColor',scatCol{i},'MarkerFaceColor',scatCol{i});

    hold on
    
    seedPreferences(seedList,useList);
    
end
xlim([0 26])
ylim([0 26])
hold off











