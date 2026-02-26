
clear all

SPN = 'D:\LGNs1\Segmentation\VAST\ArtSub_1kCore\export\export_14+11+14\';
showViews = 0;



processString = 'tube'
startAspect = [4 4 30];

MPN = [SPN(1:end-1) '_mat\']
TPN = [MPN 'tube\'];
if ~exist(TPN,'dir'),mkdir(TPN);end
load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])

obNames = obI.nameProps.names;
isUse = zeros(length(obNames),1);
for i = 1:length(obNames)
    if sum(regexp(obNames{i},processString))
        if size(dsObj(i).subs,1)>100
            isUse(i) = 1;
        end
    end
end

%% Group tubes

allOb = find(isUse);
groupNum = 2;
allOb = allOb(randperm(length(allOb)));
groupSize = fix(length(allOb)/groupNum);
for g = 1:groupNum
groupOb{g} = allOb((g-1)*groupSize+1:g*groupSize);
end

for g = 1:groupNum
useOb = groupOb{g};
clear tubes
for i = 1:length(useOb)
    tubes(i).sub =  dsObj(useOb(i)).subs;
end




%%
xrange = [4];
yrange = [3.5:.1:5.5];
zrange = [25:.5:35]
testAs = [];
for z = 1:length(zrange)
    for y = 1:length(yrange)
        for x = 1:length(xrange)
            
            testAs(size(testAs,1)+1,:) = [yrange(y) xrange(x) zrange(z)];
        end
    end
end

%% check
clear recVar minVar bestVals
addDif = zeros(size(testAs,1),1);
allDif = zeros(size(testAs,1),length(tubes));
for i = 1:length(tubes)
    
    sub = double(tubes(i).sub);
    sub(:,1) = sub(:,1)-mean(sub(:,1));
    sub(:,2) = sub(:,2)-mean(sub(:,2));
    sub(:,3) = sub(:,3)-mean(sub(:,3));
    
    
    
    fixSub = [sub(:,1) * startAspect(1,1) sub(:,2) * startAspect(1,2) sub(:,3)* startAspect(1,3)];
    showSub = fixSub;
    mins = min(showSub,[],1);
    showSub = showSub - repmat(mins,[size(showSub,1) 1])+10;
    showSub = fix(showSub/30)+1;
    maxes = max(showSub,[],1)+10;
    subVol = zeros(maxes);
    subVol(sub2ind(maxes,showSub(:,1),showSub(:,2),showSub(:,3))) = 1;
    for d = 1:3
        sumView = squeeze(sum(subVol,d));
        subplot(2,6,d)
        imshow(sumView,[0 max(sumView(:))] );
    end
    
    [Coef,Score,latent,tsquare] = princomp(fixSub);
    showSub = Score;
    mins = min(showSub,[],1);
    showSub = showSub - repmat(mins,[size(showSub,1) 1])+10;
    showSub = fix(showSub/30)+1;
    maxes = max(showSub,[],1)+10;
    subVol = zeros(maxes);
    subVol(sub2ind(maxes,showSub(:,1),showSub(:,2),showSub(:,3))) = 1;
    for d = 1:3
        sumView = squeeze(sum(subVol,d));
        subplot(2,6,d+6)
        imshow(sumView,[0 max(sumView(:))] );
    end
    
    pause(.01)
    
    
    
    parfor r = 1:size(testAs,1)
        fixSub = [sub(:,1) * testAs(r,1) sub(:,2) * testAs(r,2) sub(:,3)* testAs(r,3)];
        fixSub = round(fixSub-min(fixSub(:))+1);
        fixSub = shrinkSub(fixSub,2);
        fixSub = [fixSub(:,1)-mean(fixSub(:,1)) ...
            fixSub(:,2)-mean(fixSub(:,2)) fixSub(:,3)-mean(fixSub(:,3)) ];
        
        [Coef,Score,latent,tsquare] = princomp(fixSub);
        recVar(r,:) = std(Score);
        if showViews
            showSub = Score;
            mins = min(showSub,[],1);
            showSub = showSub - repmat(mins,[size(showSub,1) 1])+10;
            showSub = fix(showSub/30)+1;
            maxes = max(showSub,[],1)+10;
            subVol = zeros(maxes);
            subVol(sub2ind(maxes,showSub(:,1),showSub(:,2),showSub(:,3))) = 1;
            for d = 1:3
                sumView = squeeze(sum(subVol,d));
                subplot(2,3,d+3)
                imshow(sumView,[0 max(sumView(:))] );
            end
            pause(.01)
        end
    end
    
    difVar = (recVar(:,2)-recVar(:,3))./recVar(:,3);
    minVar(i) = min(difVar);
    bestVals(i,:) = testAs(find(difVar==minVar(i),1),:);
    
    addDif = addDif + difVar/max(difVar);
    
    subplot(2,6,[4 5 6 10 11 12])
    plot(recVar(:,2)/max(recVar(:,2)),'b')
    hold on
    plot(recVar(:,3)/max(recVar(:,2)),'r')
    plot(difVar/max(difVar),'g')
    plot(addDif/i,'k')
    hold off
    pause(.1)
    
    allDif(:,i) = difVar;
end

%%

%scatter(bestVals(:,1),bestVals(:,3))

minAdd = find(addDif == min(addDif),1);
bestAdd = testAs(minAdd,:)

bestAspects(g,:) = bestAdd;


meanDif = mean(allDif,2);
medDif = median(allDif,2);

plot(medDif,'r')
hold on
plot(meanDif,'g')
hold off

testAs(find(meanDif == min(meanDif)),:)
testAs(find(medDif == min(medDif)),:)


end
bestAspects

