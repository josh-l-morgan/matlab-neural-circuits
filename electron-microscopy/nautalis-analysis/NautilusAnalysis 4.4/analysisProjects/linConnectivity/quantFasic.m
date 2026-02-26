


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


%%
fascicDist = 6;

distRange = 1:1:fascicDist;
distBin = .5;


targCell = 125;
filtPre = sm.pre == targCell;
filtPost = sm.postClass == 2;

testSyn = find(filtPre & filtPost);

voxSiz = 1;

targSub = double(sm.subs(sm.subType == 1,:))*obI.em.dsRes(1);
axSub = double(sm.subs(sm.subType == 2,:))*obI.em.dsRes(1);
shaftSub = double(sm.subs(sm.subType == 3,:))*obI.em.dsRes(1);
allSub = double(sm.subs)*obI.em.dsRes(1);

clear typeDist
closePlot = zeros(length(testSyn),length(distRange));

goodTest = testSyn * 0;
for i = 1:length(testSyn)
    
    sprintf('running %d of %d',i,length(testSyn))
    
    synPos = sm.pos(testSyn(i),:);
    postID = sm.post(testSyn(i),:);
    
    %Find process type
    typeDist{1} = min(sqrt((targSub(:,1)-synPos(1)).^2 + ...
        (targSub(:,2)-synPos(2)).^2 + (targSub(:,3)-synPos(3)).^2));
    typeDist{2}  = min(sqrt((axSub(:,1)-synPos(1)).^2 + ...
        (axSub(:,2)-synPos(2)).^2 + (axSub(:,3)-synPos(3)).^2));
    typeDist{3}  = min(sqrt((shaftSub(:,1)-synPos(1)).^2 + ...
        (shaftSub(:,2)-synPos(2)).^2 + (shaftSub(:,3)-synPos(3)).^2));
    
    minType = [typeDist{1} typeDist{2} typeDist{3}];
    
    subDist  = (sqrt((allSub(:,1)-synPos(1)).^2 + ...
        (allSub(:,2)-synPos(2)).^2 + (allSub(:,3)-synPos(3)).^2)); %distance of targ cell to synapse
    
    synType(i) = find(minType == min(minType),1);
    
    
    checkSub = find(subDist<=fascicDist); %Targ cell 125 subs
    
    subCell = names2Subs(obI,dsObj,postID);
    subTarg = double(subCell{1})*obI.em.dsRes(1);
    
    targSubDist  = (sqrt((subTarg(:,1)-synPos(1)).^2 + ...
        (subTarg(:,2)-synPos(2)).^2 + (subTarg(:,3)-synPos(3)).^2));
    subTarg = subTarg(targSubDist <=(fascicDist*3),:); %positions of nearby pieces of targeted cells
    closeSub = allSub(checkSub,:); %position of all nearby pieces of target cell
    
    
    if ~isempty(subTarg) & ~isempty(closeSub)
         goodTest(i) = 1;
        clear targDist
        for s = 1:length(checkSub)
            cSub = closeSub(s,:);
            targDist(s)  = min(sqrt((subTarg(:,1)-cSub(1)).^2 + ...
                (subTarg(:,2)-cSub(2)).^2 + (subTarg(:,3)-cSub(3)).^2));
        end
        
        closeDist = subDist(checkSub); %distance of nearby targ125 to synapse
        
        
        %%Analyze
        
        
        closest = distRange * 0;
        for d = 1:length(distRange) %run through distances to filter target cell positions
            
            isRange = (closeDist >= (distRange(d)-distBin)) & ...
                (closeDist < (distRange(d) + distBin));
            if sum(isRange)
                closest(d) = min(targDist(isRange')); %find closest position of targ to post cell
            else
                goodTest(i) = 0;
            end
        end
        
        closePlot(i,:) = closest;
        
        subplot(1,2,1)
        scatter3(subTarg(:,1),subTarg(:,2),subTarg(:,3),'.','b')
        daspect([1 1 1])
        hold on
        scatter3(closeSub(:,1),closeSub(:,2),closeSub(:,3),'.','r')
        scatter3(synPos(:,1),synPos(:,2),synPos(:,3),200,'s','filled','g')
        hold off
        subplot(1,2,2)
        hold off
        plot(distRange,closest,'color','r','linewidth',5)
        hold on
        scatter(closeDist,targDist,'.','b');
        
        ylim([0 fascicDist]);
        xlim([0 fascicDist]);
        pause(.01)
        
    else
        'missing positions'
    end
end

%%

scatter(sm.pos(testSyn(synType == 1),1),sm.pos(testSyn(synType == 1),2),'.','g')
hold on
scatter(sm.pos(testSyn(synType == 2),1),sm.pos(testSyn(synType == 2),2),'.','r')
scatter(sm.pos(testSyn(synType == 3),1),sm.pos(testSyn(synType == 3),2),'.','b')

hold off
%%

closeTarg =  closePlot((synType' == 1) & goodTest,:);
closeShaft =  closePlot((synType' == 3) & goodTest,:);



L = size(closeTarg,1);
sortClose = sort(closeTarg,1,'ascend');
targ05 = sortClose(round(L*.05),:);
targ50 = sortClose(round(L*.50),:);
targ95 = sortClose(round(L*.95),:);
targ025 = sortClose(round(L*.025),:);
targ975 = sortClose(round(L*.975),:);



L = size(closeShaft,1);
sortClose = sort(closeShaft,1,'ascend');
shaft05 = sortClose(round(L*.05),:);
shaft50 = sortClose(round(L*.50),:);
shaft95 = sortClose(round(L*.95),:);
shaft975 = sortClose(round(L*.975),:);
shaft025 = sortClose(round(L*.025),:);


hold off
plot(distRange,targ025,'g')
hold on
plot(distRange,targ50,'color','g','linewidth',5)
plot(distRange,targ975,'g')
plot(distRange,shaft025,'r')
plot(distRange,shaft50,'color','r','linewidth',5)
plot(distRange,shaft975,'r')

hold off

%% Test fasci
%definition is amount of neurite within 1um 5-15 um from synapse

closeEnough = .5;
sampTarg = closeTarg(:,(distRange>5) & (distRange <= 15))<=closeEnough;
sampShaft = closeShaft(:,(distRange>5) & (distRange <= 15))<=closeEnough;

mean(sampTarg(:))
mean(sampShaft(:))



%% p
subplot(1,1,1)
closestRange = [0:.5:5];
histShaft = hist(closeShaft(:,end),closestRange)/sum(histShaft);
histTarg = hist(closeTarg(:,end),closestRange)/sum(histTarg);
bar(closestRange,[histTarg' histShaft'])
hold on
bar(closestRange +.2,histShaft,'r')

hold off

realDif = shaft50(end) - targ50(end);
realFacDif = mean(closeTarg(:,end)<1) - mean(closeShaft(:,end)<1) ;
mean(closeTarg(:,end)<1)
mean(closeShaft(:,end)<1)


%%

reps = 100000;
for r = 1:reps
    choosen = randperm(length(synType));
    randSyn =  synType(choosen);
    closeTarg =  closePlot((synType(randSyn)' == 1) & goodTest(randSyn),:);
    closeShaft =  closePlot((synType(randSyn)' == 3) & goodTest(randSyn),:);
    testTarg(r) = median(closeTarg(:,end));
    testShaft(r) = median(closeShaft(:,end));
    randFacDif(r) = mean(closeTarg(:,end)<1) - mean(closeShaft(:,end)<1) ;

end
randDif = testShaft-testTarg;


P = sum(randDif>=realDif)/length(randDif)
Pfac = sum(randFacDif>=realFacDif)/length(randDif)

%
% hold off
% lineCol = [0 1 0; 1 0 0; 0 0 1; 0 0 0; 0 0 0; 0 0 0];
% for i = 1:length(synType)
%     plot(distRange,closePlot(i,:),'color',lineCol(synType(i),:),'linewidth',1)
%     hold on
% end
%











