function[res] = scrambleSyn(preList,postList,synapses);

    manySyn = 5;
s = 1;

  %seed = conTo(s).targ;
   % preList = conTo(s).rgcList;
   % postList = setdiff(conTo(s).tcrList,seed);
   % synapses = conTo(s).syn;
    
    con = zeros(length(preList),length(postList));
    useSyn = zeros(size(synapses,1),1);
    syn = [];
    allSyn = [];
    for y = 1:length(preList)
        for x = 1:length(postList)
            goodSyn = synapses(:,2)==postList(x) & (synapses(:,1) == preList(y));
            useSyn = useSyn+ goodSyn;
            syn = cat(1,syn,repmat([y x],sum(goodSyn),1));
            allSyn = cat(1,allSyn,[y x]);
            con(y,x) = sum(goodSyn);
        end
    end
    
    
    %% measure real
    
    
    allID = allSyn(:,1) * 1000 + allSyn(:,2);
    preid = syn(:,1) * 1000;
    idsyn = preid + syn(:,2);
    histsyn = hist(idsyn,allID);
    fullCVRMSE = cvrmse(histsyn);
    
    
    
    sortPreCon = con*0;
    for y = 1:size(con,1)
        sortPreCon(y,:) = sort(con(y,:),'descend');
    end
    
    bar(mean(sortPreCon,1))
    bar(sort(con(:),'descend'))
    hist(con(:),[0:1:max(con(:))]);
    
    conList = ones(size(con));
    conList = conList(:);
    
    
    idsyn = syn(:,1) * 1000 + syn(:,2);
    usyn = unique(idsyn);
    histsyn = hist(idsyn,usyn);
    sortHist = sort(histsyn,'descend');
    
    res.realMany(s) = sum(histsyn>=manySyn)
    res.realCVRMSE(s) = fullCVRMSE;
    
    %%
    
    numSyn = size(syn,1);
    reps = 10000;
    allHist = zeros(reps,length(allID));
    sortAllHist = allHist;
    constrainPost = 1;
    for r = 1:reps
        if constrainPost
            newpost = syn(randperm(numSyn),2);
        else
            newpost = ceil(rand(numSyn,1)*size(con,2));
        end
        newid = preid + newpost;
        newhist = hist(newid,allID);
        randCVRMSE(r) = cvrmse(newhist);
        allHist(r,:) = newhist;
        sortAllHist(r,:) = sort(newhist,'descend');
    end
    
    %%
    meanRandHist = median(sortAllHist,1);
    doubleSort = sort(sortAllHist);
    minHist = doubleSort(round(reps* 0.025),:);
    maxHist = doubleSort(round(reps* 0.975),:);
    
    
    freqMany = sum(allHist>=manySyn,2);
    res.randMany(s) = mean(freqMany);
    res.randManySTD(s) = std(freqMany);
    res.randManyP(s) = sum(freqMany>=res.realMany(s))/reps;
    
    res.meanRandCVRMSE(s) = mean(randCVRMSE);
    res.cvrmseP(s) = sum(randCVRMSE>=fullCVRMSE)/reps;
    
    res.con{s} = con;
    res.reps(s) = reps;
    
    subplot(2,1,1)
    hist(randCVRMSE);
    hold on
    scatter(fullCVRMSE,1,'r')
    xlim([0 fullCVRMSE + 1])
    hold off
    
    subplot(2,1,2)
    plot(minHist,'b')
    hold on
    plot(maxHist,'b')
    plot(meanRandHist,'k')
    plot(sortHist,'r')
    xlim([0 numSyn])
    hold off
    
%     
%     synBin = (0:30);
%     hist(sortHist,synBin,'k')
%     hold on
%     
%     hold off
    
    
   pause(.01)


