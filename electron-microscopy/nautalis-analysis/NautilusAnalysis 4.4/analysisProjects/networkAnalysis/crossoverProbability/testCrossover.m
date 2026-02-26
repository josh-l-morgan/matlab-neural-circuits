
%%Test the probability that crossover synapses would be concentrated on few
%%axons
clear all
load('MPN.mat')
load([MPN 'obI.mat'])
allEdges = obI.nameProps.edges(:,[2 1]);
seedList = [201 108 903 907];
useCells = obI2cellList_seedInput_RGC_TCR(obI,seedList);
conTo = makeConTo(obI,seedList);

edges = unique(allEdges,'rows');


preA = setdiff(conTo(1).rgcList,conTo(2).rgcList);
preB = setdiff(conTo(2).rgcList,conTo(1).rgcList);

postA = [];
for i = 1:length(preA);
    postT = postTo(allEdges,preA(i));
    postA = cat(1, postA, postT(:,1));
end
    postA = intersect(postA,conTo(1).tcrList);
    postA = setdiff(postA,seedList);
    
    
postB = [];
for i = 1:length(preB);
    postT = postTo(allEdges,preB(i));
    postB = cat(1, postB, postT(:,1));
end
    postB = intersect(postB,conTo(2).tcrList);
    postB = setdiff(postB,seedList);

    
crossed = intersect(postA,postB);

preT = preB;
postT = unique([postA postB]);

conA = [];
for y = 1:length(preA)
    for x = 1:length(postT)
        conA(y,x) = sum((allEdges(:,1) == preA(y)) & (allEdges(:,2) == postT(x)));
    end
end


conB = [];
for y = 1:length(preB)
    for x = 1:length(postT)
        conB(y,x) = sum((allEdges(:,1) == preB(y)) & (allEdges(:,2) == postT(x)));
    end
end

tot = sum(conB>0,2);
%conB = conB(tot>3,:);

subplot(2,1,1)
  image(conA * 20)
  subplot(2,1,2)
  image(conB*20)
  
  %%
  
 crossT = sum(conA>0,1);
 crossT = repmat(crossT,[size(conB,1) 1]);
 image(crossT * 10);
 [y x] = find(conB);
 realEdge = [y x];
 
 
 realDist = sum(crossT .* (conB>0),2);
 realTot = sum((conB>0),2);
 realNorm = realDist./realTot;
 
 image((crossT .* (conB>0)) * 10)
 realTest =  var(realDist);
 
 sortRealDist = sort(realDist);
 bar(sortRealDist)
 
 
 %%
 reps = 10000;
 randDists = zeros(reps,size(conB,1));
 sortRandDists = randDists;
 for r = 1:reps
    pick = randperm(size(realEdge,1)); 
    newEdge = [realEdge(pick,1) realEdge(:,2)];
    newCon = conB*0;
    newCon(sub2ind(size(newCon),newEdge(:,1),newEdge(:,2))) = 1;
    newDist = sum(crossT .* (newCon>0),2);
    normDist = newDist./realTot;
    randTest(r) = var(newDist);
    randMedian(r) = median(newDist);
    %image(newCon*20), pause(.1)
   % image((crossT .* newCon)*10),pause(.1)
    randDists(r,:) = newDist;
    sortRandDists(r,:) = sort(newDist);
    %bar(sortRandDists(r,:)),pause(.1)
 end
 
 meanRand = mean(sortRandDists,1);
 
 %%
 subplot(2,1,1)
 bar([meanRand' sortRealDist])
 subplot(2,1,2)
 histRand = hist(randTest);
 bar(histRand)
 hold on
 scatter(realTest,0,'r','filled')
 hold off
 
 P = sum(randTest>=realTest)/length(randTest)
 
