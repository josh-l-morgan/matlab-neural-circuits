
MPN = GetMyDir

%load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])

seedList = [108 201 109];
useList = obI2cellList_all(obI);
%seedPref = seedPreferences(seedList,useList);
allEdges = obI.nameProps.edges(:,[2 1]);
conTo = makeConTo(obI,seedList);


%%
clf
allPres = [];

plotPos = 1:(length(seedList)*4);

for s = 1:length(seedList)
    preList = conTo(s).rgcList
    allPres = cat(2,allPres,preList);
    %preList = preList(randperm(length(preList)));
    allPost = [];
    postLength = zeros(size(preList));
    for p = 1:length(preList)
        isPre = allEdges(:,1) == preList(p);
        isPost = allEdges(isPre,2);
        allPost = cat(1,allPost,isPost);
        uPost = unique(allPost);
        postLength(p) = length(uPost);
    end
    
    reps = 100
    randLength = zeros(reps,length(preList));
    
    for r = 1:reps
        preList = preList(randperm(length(preList)));
        allPost = [];
        for p = 1:length(preList)
            isPre = allEdges(:,1) == preList(p);
            isPost = allEdges(isPre,2);
            allPost = cat(1,allPost,isPost);
            uPost = unique(allPost);
            randLength(r,p) = length(uPost);
        end
    end
    
    
    usePos = (s-1)*4+1
    plotPos = setdiff(plotPos,usePos);
    %subplot(length(seedList),4,usePos)
    subplot(2,2,s)
    plot(randLength','c')
    hold on
    plot(postLength,'k','LineWidth',2);
    text(1,100,sprintf('seed cell %d',seedList(s)))
    xlim([0 40]);
    ylim([0 150]);
    
    lengthChanges{s} = postLength;
    
    
end



    allPost = [];
    postLength = zeros(size(allPres));
    for p = 1:length(allPres)
        isPre = allEdges(:,1) == allPres(p);
        isPost = allEdges(isPre,2);
        allPost = cat(1,allPost,isPost);
        uPost = unique(allPost);
        postLength(p) = length(uPost);
    end
    
    reps = 100
    randLength = zeros(reps,length(allPres));
    
    for r = 1:reps
        allPres = allPres(randperm(length(allPres)));
        allPost = [];
        for p = 1:length(allPres)
            isPre = allEdges(:,1) == allPres(p);
            isPost = allEdges(isPre,2);
            allPost = cat(1,allPost,isPost);
            uPost = unique(allPost);
            randLength(r,p) = length(uPost);
        end
    end
    
    
    
    %subplot(length(seedList),4,plotPos)
    subplot(2,2,length(seedList)+1);
    plot(randLength','c')
    hold on
    plot(postLength,'k','LineWidth',2);
    text(1,190,sprintf('all axons connected to the five seed cells'))
    xlim([0 100]);
    ylim([0 200]);
    
    lengthChanges{s} = postLength;
    









