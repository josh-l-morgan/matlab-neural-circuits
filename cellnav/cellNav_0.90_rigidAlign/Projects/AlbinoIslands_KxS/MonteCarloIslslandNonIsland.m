
reps = 100;

cidPrimNames = {'cid' 'number of primaries with RGC inputs' 'number of primaries'};
allCidPrim = [1001 7  3;
    1003 4 5;
    1004 5 6;
    1006 8 8;
    1016 7 7;
    1015 4 6;
    1002 4 5;
    1005 2 6;
    1007 4 5;
    1008 0 5;
    1011 3 5;
    1013 7 8;
    1014 5 7;
    1009 2 4;
    1019 4 4;
    1017 2 3;];
[useCid idx] = setdiff(allCidPrim(:,1),[1006 1013]);

cidPrim = allCidPrim(idx,:);
primConFrac = sum(cidPrim(:,2))/sum(cidPrim(:,3));


cNum = size(cidPrim,1);
pNum = max(cidPrim(:,3));
prefVec = zeros(1,cNum) + 0.5; %Preference of each cell for island 

prefMat = repmat(prefVec,[pNum 1]);
usePrim = repmat([1:pNum]',[1 cNum]) <= repmat(cidPrim(:,3)',[pNum 1]);



for r = 1:reps

pick =  rand(pNum,cNum) <= prefMat; %Primaries that chose island
makesCon = rand(pNum,cNum) <= primConFrac; % Choose primaries to use 
hitsIsland = sum(pick & makesCon & usePrim,1);
hitsNon = sum(~pick & makesCon & usePrim,1);


end







