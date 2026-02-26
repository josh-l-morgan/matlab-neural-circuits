%{
Modularity is sum of over all pairs of vertices 
(i,j) falling in the same group

jm Modularity definition to use is mean(Ein/E - Pin/P)

%}

load('MPN.mat')
load([MPN 'obI.mat'])
seedCells = [108 201];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedCells);
conPref = seedPreferences(seedCells,useList)

reps = 1000;

%%

rawCon = useList.con;
[a isSeed] = intersect(useList.postList,seedCells);
notSeed = setdiff(1:length(useList.postList),isSeed);
seedCon = rawCon(:,isSeed);
con = rawCon(:,notSeed);


%% get links
links = zeros(size(con,1),size(con,1));
ingroup = links;
for y = 1:size(con,1)
    for x = 1:size(con,1)
        links(y,x) = sum((con(y,:)>0) & (con(x,:)>0));
        if y == x
           % links(y,x) = 0;
        end
        ingroup(y,x) = (seedCon(y,1)>0) == (seedCon(x,1)>0);
    end
end


image(links*20)
image(ingroup*30)



adj = links;
modules = {[find(seedCon(:,1)>0)] [find(seedCon(:,2)>0)]};
Q=modularity_metric(modules,adj)






%%

[pre post] = find(links>0);
edgeNum = length(post);

edge2Seed1 = rawCon(sub2ind(size(rawCon),pre,pre*0+isSeed(1)));
edge2Seed2 = rawCon(sub2ind(size(rawCon),pre,pre*0+isSeed(2)));
edge1 = find(edge2Seed1>0);
edge2 = find(edge2Seed2>0);

%%

    post1 = post(edge1);
    post2 = post(edge2);
    realSharedPost = intersect(post1,post2);
    realShared = length(realSharedPost);
    

postNum = length(post);
P = links * 0;
[ys xs] = size(links);
for r = 1:reps
    
    newPost = post(randperm(ys));
    post1 = newPost(edge1);
    post2 = newPost(edge2);
    sharedPost = intersect(post1,post2);
    
    randShared(r) = length(sharedPost);
    
    addInd = sub2ind([ys xs],pre,newPost);
    
end

%%

maxShared = max([randShared realShared])
histBin = [0 : maxShared/20 : maxShared]
histRand = hist(randShared,histBin)
subplot(1,1,1)



bar(histBin,histRand/max(histRand(:)),'BarWidth',1)
hold on
scatter(realShared,0,'r','LineWidth',5)
xlabel('LinkedNumber')
ylabel('Montecarlo distribution')
hold off
xlim([0 maxShared])


sortRandShared = sort(randShared,'ascend');
rand99 = [sortRandShared(round(reps*.005)) sortRandShared(end - round(reps*.005))] 
rand95 = [sortRandShared(round(reps*.025)) sortRandShared(end - round(reps*.025))] 
realShared
meanRand = mean(sortRandShared);
P = sum(sortRandShared<=realShared)/reps




