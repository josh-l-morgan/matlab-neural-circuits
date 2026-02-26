


if 0
    clear all
    
    obMovDir = 'D:\LGNs1\Analysis\movies\subLin125_TargFac\'
    if ~exist(obMovDir,'dir'),mkdir(obMovDir),end
    % else
    %     'directory already exists'
    %     return
    % end
    
    load('MPN.mat')
    load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])
    
    sm = addDatToSynMat(obI)
    sm = getTopoEucDistBetweenSyn(sm);
    sm = getTopoEucDistBetweenSkelAndSyn(sm);
    sm = labelShaftSkel(sm);
    sm = labelSubTypes(sm);
    
end

%% Types {'targ ax shaft body'}
 if 1
nodeNum = length(sm.skelNodes);

tcList = sm.post(sm.postClass == 2);


%% Distances between synapses
dists = zeros(nodeNum);

for i = 1:size(sm.skelPos,1)
   
    dists(i,:) = sqrt((sm.skelPos(:,1)-sm.skelPos(i,1)).^2 + ...
        (sm.skelPos(:,2)-sm.skelPos(i,2)).^2 + ... 
        (sm.skelPos(:,3)-sm.skelPos(i,3)).^2);
        
end

%% Process type

nodeProcType = sm.isTarg + sm.isAx * 2 + sm.isShaft * 3 + sm.isBody * 4;
procTypeMatch = zeros(nodeNum);
goodCheck = zeros(nodeNum);
for y = 1:nodeNum
    for x = 1:nodeNum
        procTypeMatch(y,x) = nodeProcType(y) * 10 + nodeProcType(x);
        goodCheck(y,x) = y>x;
        
    end
end

%% PrePost match

hasSyn = zeros(nodeNum,1);
hasSyn(sm.syn2skel.closestSkel) = 1;
orSyn = hasSyn | hasSyn';
andSyn = hasSyn & hasSyn';

[y x] = find(andSyn);
samePre = zeros(nodeNum);
samePost = zeros(nodeNum);
difPre = zeros(nodeNum);
difPost = zeros(nodeNum);

for i = 1:length(y)
    
    ySyn = sm.syn2skel.syn(find(sm.syn2skel.closestSkel == y(i)),:);
    xSyn = sm.syn2skel.syn(find(sm.syn2skel.closestSkel == x(i)),:);
    if y(i) ~= x(i)
        for s = 1:size(ySyn,1)
            if (ySyn(s,1) ~= 125) & (ySyn(s,1) ~= 0)
                samePre(y(i),x(i)) = samePre(y(i),x(i)) + sum(xSyn(:,1)==ySyn(s,1));
                difPre(y(i),x(i)) = difPre(y(i),x(i)) + sum((xSyn(:,1)~=ySyn(s,1)) & ...
                    (xSyn(:,1) ~= 0) & (xSyn(:,1) ~= 125));
            end
            if (ySyn(s,2) ~= 125) & (ySyn(s,2) ~= 0) & sum(tcList == ySyn(s,2))
                samePost(y(i),x(i)) = samePost(y(i),x(i)) + sum(xSyn(:,2)==ySyn(s,2));
                difPost(y(i),x(i)) = difPost(y(i),x(i)) + sum((xSyn(:,2)~=ySyn(s,2)) & ...
                    (xSyn(:,2) ~= 0) & (xSyn(:,2) ~= 125));
            end
        end
    end
    
end

%% Count synapses on each process type





%% Ask your questions

bin = 10;

bins = [0:bin:400];
clear ddPostMatch ddPostAll ssPostMatch ssPostAll ttPostMatch ttPostAll
clear aaPostMatch aaPostAll stPostMatch stPostAll adPostMatch adPostAll
clear ssPostDif ttPostDif aaPostDif stPostDif ddPostDif adPostDif


for i = 1:length(bins)-1
    

    subplot(1,1,1)
    scatter(sm.skelPos(sm.isShaft,2),sm.skelPos(sm.isShaft,1),'.','b')
    hold on
    scatter(sm.skelPos(sm.isTarg,2),sm.skelPos(sm.isTarg,1),'.','g')
    scatter(sm.skelPos(sm.isAx,2),sm.skelPos(sm.isAx,1),'.','r')


   useDist =  (dists >= bins(i)) & (dists < bins(i+1)) ;
    
   ssPostMatch(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 33) & samePost));
   ssPostAll(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 33)));
   ssPostDif(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 33) & difPost));

   ttPostMatch(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 11) & samePost));
   ttPostAll(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 11)));
   ttPostDif(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 11) & difPost));
   
   aaPostMatch(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 22) & samePost));
   aaPostAll(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 22)));
   aaPostDif(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 22) & difPost));
   
   stPostMatch(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 31) & samePost));
   stPostAll(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 31)));
   stPostDif(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 31) & difPost));
   
   ddPostMatch(i) = sum(sum( useDist &  goodCheck & ...
       ((procTypeMatch == 11) | (procTypeMatch == 33) | (procTypeMatch == 13) | (procTypeMatch == 31))...
        & samePost));
   ddPostAll(i) = sum(sum(useDist &  goodCheck & ...
       ((procTypeMatch == 11) | (procTypeMatch == 33) | (procTypeMatch == 13) | (procTypeMatch == 31))));
   ddPostDif(i) = sum(sum( useDist &  goodCheck & ...
       ((procTypeMatch == 11) | (procTypeMatch == 33) | (procTypeMatch == 13) | (procTypeMatch == 31))...
        & difPost));
    
   adPostMatch(i) = sum(sum( useDist &  goodCheck & ...
       ((procTypeMatch == 12) | (procTypeMatch == 21) | (procTypeMatch == 23) | (procTypeMatch == 32))...
        & samePost));
   adPostAll(i) = sum(sum( useDist &  goodCheck & ...
       ((procTypeMatch == 12) | (procTypeMatch == 21) | (procTypeMatch == 23) | (procTypeMatch == 32))));
   adPostDif(i) = sum(sum( useDist &  goodCheck & ...
       ((procTypeMatch == 12) | (procTypeMatch == 21) | (procTypeMatch == 23) | (procTypeMatch == 32))...
        & difPost));   
   
      testMatch = useDist &  goodCheck & ...
       ((procTypeMatch == 12) | (procTypeMatch == 21) | (procTypeMatch == 23) | (procTypeMatch == 32));
   
   
    [y x] = find(testMatch);
    posY = sm.skelPos(y,:);
    posX = sm.skelPos(x,:);
   
    for p = 1:size(posY,1)
    plot([posY(p,2) posX(p,2)],[posY(p,1) posX(p,1)]);
    end
        hold off
    pause(.01)
end



%%count synapses on axon


subplot(4,1,1)
%plot(bins(2:end),adPostMatch./adPostAll,'r')
hold on
%plot(bins(2:end),ddPostMatch./ddPostAll,'b')
plot(bins(2:end),ssPostMatch./ssPostAll,'k')
%plot(bins(2:end),aaPostMatch./aaPostAll,'c')
plot(bins(2:end),stPostMatch./stPostAll,'m')
plot(bins(2:end),ttPostMatch./ttPostAll,'g')

hold off

subplot(4,1,2)

%plot(bins(2:end),adPostAll,'r')
hold on
%plot(bins(2:end),ddPostAll,'b')
plot(bins(2:end),ssPostAll,'k')

%plot(bins(2:end),aaPostAll,'c')
plot(bins(2:end),stPostAll,'m')
plot(bins(2:end),ttPostAll,'g')

hold off

subplot(4,1,3)

%plot(bins(2:end),adPostMatch,'r')
hold on
plot(bins(2:end),ssPostMatch,'k')
%plot(bins(2:end),ddPostMatch,'b')
%plot(bins(2:end),aaPostMatch,'c')
plot(bins(2:end),stPostMatch,'m')
plot(bins(2:end),ttPostMatch,'g')



subplot(4,1,4)
%plot(bins(2:end),adPostMatch./(adPostMatch + adPostDif),'r')
hold on
%plot(bins(2:end),ddPostMatch./(ddPostMatch + ddPostDif),'b')
plot(bins(2:end),ssPostMatch./(ssPostMatch + ssPostDif),'k')
%plot(bins(2:end),aaPostMatch./(aaPostMatch + aaPostDif),'c')
plot(bins(2:end),stPostMatch./(stPostMatch + stPostDif),'m')
plot(bins(2:end),ttPostMatch./(ttPostMatch + ttPostDif),'g')

hold off





% 
% 
% plot(bins(2:end),adPostMatch)
% hold on
% plot(bins(2:end),ddPostMatch)
% plot(bins(2:end),aaPostMatch)
% hold off

end
%% Monte carlo that uses the dend probabiltiy of a synapse match at a node 
%% to predict how often a axon/dend match should be expected

subplot(3,1,1)

reps = 10000
testEffect = 2

plot(bins(2:end),ssPostMatch./(ssPostMatch + ssPostDif),'b')
hold on
plot(bins(2:end),ttPostMatch./(ttPostMatch + ttPostDif),'g')
hold off

conProb = stPostMatch./stPostAll;
conOpportunities = stPostAll;
matchRat = stPostMatch./stPostAll;
testStat = sum(stPostMatch);


randHit = zeros(length(bins)-1,reps);
randHitHigh = zeros(length(bins)-1,reps);
randHitLow = zeros(length(bins)-1,reps);
for r = 1:reps
    for i = 1:length(bins)-1
        randHit(i,r) = sum(rand(conOpportunities(i),1)<=conProb(i));
    end    
    
    for i = 1:length(bins)-1
        randHitHigh(i,r) = sum(rand(conOpportunities(i),1)<=(conProb(i)* testEffect));
    end  
    
    for i = 1:length(bins)-1
        randHitLow(i,r) = sum(rand(conOpportunities(i),1)<=(conProb(i)* 1/testEffect));
    end  
    
    
    plot(bins(2:end),randHit(:,r),'k')
    hold on

end
 
subplot(1,1,1)

plot(bins(2:end),stPostMatch./stPostAll,'r')

hold off

rankRes  = sort(randHit,2)
%randStat = 

ci95 = [rankRes(:,round(reps*.025)) rankRes(:,round(reps * .975))]
ci50 = [rankRes(:,round(reps*.25)) rankRes(:,round(reps * .75))]

mean(rankRes,2)
plot(bins(2:end),ci95(:,1),'k');
hold on
plot(bins(2:end),ci50(:,1),'b');
plot(bins(2:end),ci50(:,2),'b');

plot(bins(2:end),ci95(:,2),'k');
plot(bins(2:end),matchRat,'r')
hold off

randResSum = sum(randHit,1);
randResHighSum = sum(randHitHigh,1);
randResLowSum = sum(randHitLow,1);


histRange = [0:.1:20];
randRes = hist(randResSum,histRange);
bar(histRange,randRes)
hold on
scatter(testStat,0,50,'r','o','filled')
hold off

rankRes = sort(randResSum,'ascend')
ci95 = [rankRes(:,round(reps*.025)) rankRes(:,round(reps * .975))]
ci50 = [rankRes(:,round(reps*.25)) rankRes(:,round(reps * .75))]
mean(rankRes)

P = mean(randRes>=testStat)


%%Find Power
sig = 0.05;
pHigh = zeros(reps,1);
pLow = zeros(reps,1);
for i = 1:length(randResSum)
    pHigh(i) = mean(randResSum>=randResHighSum(i)) ;
    pLow(i) = mean(randResSum<=randResLowSum(i)) ;
end
mean(randResHighSum)
mean(randResLowSum)
Bhigh = mean(pHigh<=sig)
Blow = mean(pLow<=sig)



