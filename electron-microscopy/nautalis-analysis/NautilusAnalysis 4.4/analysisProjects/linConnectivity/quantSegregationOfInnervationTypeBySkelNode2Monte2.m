


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
    load('D:\LGNs1\Analysis\sm.mat')
    
end

%% Types {'targ ax shaft body'}


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

isTC = zeros(max(sm.syn2skel.syn(:)+1),1);
isTC(tcList+1) = 1;
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
                    (xSyn(:,2) ~= 0) & (xSyn(:,2) ~= 125)& (isTC(xSyn(:,2)+1)));
            end
        end
    else % if the same node
        
        for s = 1:size(ySyn,1)
            if (ySyn(s,1) ~= 125) & (ySyn(s,1) ~= 0)
                samePre(y(i),x(i)) = samePre(y(i),x(i)) + sum(xSyn(:,1)==ySyn(s,1))-1;
                difPre(y(i),x(i)) = difPre(y(i),x(i)) + sum((xSyn(:,1)~=ySyn(s,1)) & ...
                    (xSyn(:,1) ~= 0) & (xSyn(:,1) ~= 125));
            end
            if (ySyn(s,2) ~= 125) & (ySyn(s,2) ~= 0) & sum(tcList == ySyn(s,2))
                samePost(y(i),x(i)) = samePost(y(i),x(i)) + sum(xSyn(:,2)==ySyn(s,2))-1;
                difPost(y(i),x(i)) = difPost(y(i),x(i)) + sum((xSyn(:,2)~=ySyn(s,2)) & ...
                    (xSyn(:,2) ~= 0) & (xSyn(:,2) ~= 125) & (isTC(xSyn(:,2)+1)));
            end
            
        end
        samePre(y(i),x(i)) = samePre(y(i),x(i))/2;
        difPre(y(i),x(i)) = difPre(y(i),x(i))/2;
        samePost(y(i),x(i)) = samePost(y(i),x(i))/2;
        difPost(y(i),x(i)) =  difPost(y(i),x(i))/2;
        %           if sum(isTC(xSyn(:,2)+1))>1
        
        %             ySyn,xSyn,isTC(xSyn(:,2)+1)
        %             samePost(y(i),x(i)),difPost(y(i),x(i)),
        %             pause
        %           end
        
    end
    
end






%% Ask your questions

bin = 10;

bins = [0:bin:10];
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
    
    %       testMatch = useDist &  goodCheck & ...
    %        ((procTypeMatch == 12) | (procTypeMatch == 21) | (procTypeMatch == 23) | (procTypeMatch == 32));
    %
    %
    %     [y x] = find(testMatch);
    %     posY = sm.skelPos(y,:);
    %     posX = sm.skelPos(x,:);
    %
    %     for p = 1:size(posY,1)
    %     plot([posY(p,2) posX(p,2)],[posY(p,1) posX(p,1)]);
    %     end
    %         hold off
    %     pause(.01)
end



subplot(4,1,1)
plot(bins(2:end),adPostMatch./adPostAll,'r')
hold on
plot(bins(2:end),ddPostMatch./ddPostAll,'b')
plot(bins(2:end),ssPostMatch./ssPostAll,'k')
plot(bins(2:end),aaPostMatch./aaPostAll,'c')
plot(bins(2:end),stPostMatch./stPostAll,'m')
plot(bins(2:end),ttPostMatch./ttPostAll,'g')

hold off

subplot(4,1,2)

plot(bins(2:end),adPostAll,'r')
hold on
plot(bins(2:end),ddPostAll,'b')
plot(bins(2:end),ssPostAll,'k')

plot(bins(2:end),aaPostAll,'c')
plot(bins(2:end),stPostAll,'m')
plot(bins(2:end),ttPostAll,'g')

hold off

subplot(4,1,3)

plot(bins(2:end),adPostMatch,'r')
hold on
plot(bins(2:end),ssPostMatch,'k')
plot(bins(2:end),ddPostMatch,'b')
plot(bins(2:end),aaPostMatch,'c')
plot(bins(2:end),stPostMatch,'m')
plot(bins(2:end),ttPostMatch,'g')



subplot(4,1,4)
plot(bins(2:end),adPostMatch./(adPostMatch + adPostDif),'r')
hold on
plot(bins(2:end),ddPostMatch./(ddPostMatch + ddPostDif),'b')
plot(bins(2:end),ssPostMatch./(ssPostMatch + ssPostDif),'k')
plot(bins(2:end),aaPostMatch./(aaPostMatch + aaPostDif),'c')
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


%% Monte carlo that switches process match types

subplot(3,1,1)

plot(bins(2:end),ssPostMatch./(ssPostMatch + ssPostDif),'b')
hold on
plot(bins(2:end),ttPostMatch./(ttPostMatch + ttPostDif),'g')
hold off

ssConRats = ssPostMatch./(ssPostMatch + ssPostDif);
ttConRats  =  ttPostMatch./(ttPostMatch + ttPostDif);
ttConRat = mean(ttConRats(1))
ttConTot = mean(ttPostMatch(1) + ttPostDif(1))
ssConRat = mean(ssConRats(1))
ssConTot = mean(ssPostMatch(1) + ssPostDif(1))

testRat = ttConRat/ssConRat;




ssConRats = ssPostMatch./(ssPostAll);
ttConRats  =  ttPostMatch./(ttPostAll);
ttConRat = mean(ttConRats(1))
ttConTot = mean(ttPostAll(1))
ssConRat = mean(ssConRats(1))
ssConTot = mean(ssPostAll(1))

testRat2 = ttConRat/ssConRat;




reps = 10000
testRatRand = zeros(reps,1);
testRatRand2 = zeros(reps,2);
for r = 1:reps
    sprintf('running itteration %d of %d', r,reps)
    ind33 = find(procTypeMatch == 33);
    ind11 = find(procTypeMatch == 11);
    mixInd = [ind33; ind11];
    mixedInd = mixInd(randperm(length(mixInd)));
    
    procTypeMatchRand = procTypeMatch;
    procTypeMatchRand(mixedInd(1:length(ind33))) = 33;
    procTypeMatchRand(mixedInd(length(ind33)+1:end)) = 11;
    
    
clear ddPostMatchRand  ssPostMatchRand ssPostAllRand ssPostDifRand ttPostMatchRand ttPostAllRand ttPostDifRand
clear testRatRand testRatRand2
    
    
    for i = 1:length(bins)-1
        
        useDist =  (dists >= bins(i)) & (dists < bins(i+1)) ;
        
        ssPostMatchRand(i) = sum(sum( useDist &  goodCheck & (procTypeMatchRand == 33) & samePost));
        ssPostAllRand(i) = sum(sum( useDist &  goodCheck & (procTypeMatchRand == 33)));
        ssPostDifRand(i) = sum(sum( useDist &  goodCheck & (procTypeMatchRand == 33) & difPost));
        
        ttPostMatchRand(i) = sum(sum( useDist &  goodCheck & (procTypeMatchRand == 11) & samePost));
        ttPostAllRand(i) = sum(sum( useDist &  goodCheck & (procTypeMatchRand == 11)));
        ttPostDifRand(i) = sum(sum( useDist &  goodCheck & (procTypeMatchRand == 11) & difPost));
        
        
        
        
        
        %    aaPostMatch(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 22) & samePost));
        %    aaPostAll(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 22)));
        %    aaPostDif(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 22) & difPost));
        %
        %    stPostMatch(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 31) & samePost));
        %    stPostAll(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 31)));
        %    stPostDif(i) = sum(sum( useDist &  goodCheck & (procTypeMatch == 31) & difPost));
        %
        %    ddPostMatch(i) = sum(sum( useDist &  goodCheck & ...
        %        ((procTypeMatch == 11) | (procTypeMatch == 33) | (procTypeMatch == 13) | (procTypeMatch == 31))...
        %         & samePost));
        %    ddPostAll(i) = sum(sum(useDist &  goodCheck & ...
        %        ((procTypeMatch == 11) | (procTypeMatch == 33) | (pro rocTypeMatch == 13) | (procTypeMatch == 31))...
        %         & difPost));
        %
        %    adPostMatchRand(i) = sum(sum( useDist &  goodCheck & ...
        %        ((procTypeMatch == 12) | (procTypeMatch == 21) | (procTypeMatch == 23) | (procTypeMatchRand == 32))...
        %         & samePost));
        %    adPostAllRand(i) = sum(sum( useDist &  goodCheck & ...
        %        ((procTypeMatch == 12) | (procTypeMatch == 21) | (procTypeMatch == 23) | (procTypeMatchRand == 32))));
        %    adPostDifRand(i) = sum(sum( useDist &  goodCheck & ...
        %        ((procTypeMatch == 12) | (procTypeMatch == 21) | (procTypeMatch == 23) | (procTypeMatchRand == 32))...
        %         & difPost));
        %
        %       testMatch = useDist &  goodCheck & ...
        %        ((procTypeMatch == 12) | (procTypeMatch == 21) | (procTypeMatch == 23) | (procTypeMatch == 32));
        %
        %
        %     [y x] = find(testMatch);
        %     posY = sm.skelPos(y,:);
        %     posX = sm.skelPos(x,:);
        %
        %     for p = 1:size(posY,1)
        %     plot([posY(p,2) posX(p,2)],[posY(p,1) posX(p,1)]);
        %     end
        %         hold off
        %     pause(.01)
        
        
    end
    subplot(3,1,2)
    plot(bins(2:end),ssPostMatchRand./(ssPostMatchRand + ssPostDifRand),'b')
    hold on
    plot(bins(2:end),ttPostMatchRand./(ttPostMatchRand + ttPostDifRand),'g')
    hold off
    
    
    
    ssConRatsRand = ssPostMatchRand./(ssPostMatchRand + ssPostDifRand);
    ttConRatsRand  =  ttPostMatchRand./(ttPostMatchRand + ttPostDifRand);
    ttConRatRand = mean(ttConRatsRand(1))
    ttConTotRand = mean(ttPostMatchRand(1) + ttPostDifRand(1))
    ssConRatRand = mean(ssConRatsRand(1))
    ssConTotRand = mean(ssPostMatchRand(1) + ssPostDifRand(1))
    
    testRatRand(r) = ttConRatRand/ssConRatRand;
    
    
    ssConRatsRand = ssPostMatchRand./(ssPostAllRand);
    ttConRatsRand  =  ttPostMatchRand./(ttPostAllRand);
    ttConRatRand = mean(ttConRatsRand(1));
    ttConTotRand = mean(ttPostAllRand(1));
    ssConRatRand = mean(ssConRatsRand(1));
    ssConTotRand = mean(ssPostAllRand(1));
    
    testRatRand2(r) = ttConRatRand/ssConRatRand;
    
    
end

subplot(4,1,3)
histRange = [0:.1:4];
randRes = hist(testRatRand,histRange);
bar(histRange,randRes)
hold on
scatter(testRat,0,50,'r','o','filled')
hold off

rankRes = sort(testRatRand,'ascend');
ci95 = [rankRes(round(reps*.025)) rankRes(round(reps * .975))]
mean(rankRes)
P = sum(rankRes>=testRat)/reps


subplot(4,1,4)
histRange = [0:.1:4];
randRes2 = hist(testRatRand2,histRange);
bar(histRange,randRes2)
hold on
scatter(testRat2,0,50,'r','o','filled')
hold off

rankRes2 = sort(testRatRand2,'ascend');
ci95 = [rankRes2(round(reps*.025)) rankRes2(round(reps * .975))]
mean(rankRes2)
P = sum(rankRes2>=testRat2)/reps








