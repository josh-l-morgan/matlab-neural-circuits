function[listCells cellCols] = seedMarkov(MPN)


%MPN = GetMyDir;
load([MPN 'obI.mat']);

seedList = [108 201 109 903 907];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
con = useList.con;
%%

reps = 1000000000
%%
preCount = zeros(size(con,1),1);
transPre= {};

for i = 1:size(con,1)
    transSub = [];
    for p = 1:size(con,2)
        if con(i,p)
            %pause
        
        transSub = [transSub ones(1,con(i,p))*p];
        end
    end
    transPre{i} = transSub;
    preCount(i) = length(transSub);
end


postCount = zeros(size(con,1),1);
transPost = {};
for p = 1:size(con,2)
    transSub = [];
    for i = 1:size(con,1)
        if con(i,p)
            %pause
        
        transSub = [transSub ones(1,con(i,p))*i];
        end
    end
    transPost{p} = transSub;
    postCount(p) = length(transSub);
end

trans = cat(1,transPre',transPost');
pCount = cat(1,preCount,postCount);

%% 

num = length(transPost);
P = ceil(rand*num);
recPass = con*0;
memLength = 1000
mem = zeros(memLength,1);
newMem = mem;


wasHit = zeros(1,num);
lastHit = zeros(1,num);
countHits = zeros(num);
trackHit = zeros(num);

numHit = zeros(num);
sumLength = zeros(num,1);

for r = 1:reps
    I = transPost{P}(ceil(rand*postCount(P)));
    P = transPre{I}(ceil(rand*preCount(I)));
    
    wasHit(P) = 1;
    lastHit = lastHit + 1;
    countHits(P,wasHit>0) = countHits(P,wasHit>0) + 1; %record how many measurements are mode for each cell
    trackHit(P,wasHit>0) = trackHit(P,wasHit>0) + lastHit(wasHit>0); %record last time each P was passed
    lastHit(P) = 0; %reset hit memory
    
    
    if ~mod(r,10000)
        
        image(trackHit./countHits*.2),pause(.01)
        minCount = min(countHits(:));
        disp(sprintf('min count is %d after %d jumps',minCount,r))
        if minCount>1000,break,end
    end
    
    
%     showCon = con;
%     showCon(I,:) = showCon(I,:) +5;
%     showCon(:,P) = showCon(:,P) + 5;
%     image(showCon*10),pause(.01)
%     recPass(I,P) = recPass(I,P) +1;
end



meanPath = trackHit./countHits;
image(meanPath*.1)

%%

clear toSeed fromSeed
for s = 1:length(seedList)
    targ = (useList.postList == seedList(s));
    toSeed(s,:) = meanPath(targ,:);
    fromSeed(s,:) = meanPath(:,targ)';
end

seedWalk(1,:) = toSeed(1,:);
seedWalk(2,:) = toSeed(2,:);
seedWalk(3,:) = min(toSeed(3:end,:),[],1);

listCells = useList.postList;
cellCols = seedWalk'
for i = 1:3;
    cellCols(:,i) = cellCols(:,i)/max(cellCols(:,i));
end




















