
if 0
    clear all
    [TFN TPN] = uigetfile
    load([TPN TFN])
else
    clear all
    TFN = 'sm_cid4.mat';
    TPN = 'G:\IxQ\Matlab\Analysis\SMs\';
    load([TPN TFN])
end



pos = sm.nep.pos;
edges = sm.nep.edges;
groupEdges = sm.nep.groupEdges;
bones = sm.nep.bones;
rad = sm.nep.nodeRad;
rad = sm.nep.meanNodeRad;

sm.nep.props

%% Pick seed


global globPt
%globalPt.figH = figure;

scatter3(pos(:,1),pos(:,2),pos(:,3),50,'o','k','filled');
globPt.Y = nan;
guiPt
for t = 1:1000
    pause(.1)
    if ~isnan(globPt.Y)
        set(globPt.dcm ,'DisplayStyle','datatip',...
            'SnapToDataVertex','off','Enable','off')
        delete(findall(globPt.figH,'Type','hggroup'));
        break
    end
end

dif1 = sum(abs(pos - repmat([globPt.Y globPt.X globPt.Z],[size(pos,1) 1])),2);
seed = find(dif1 == min(dif1));
%close(globPt.figH)

%% Turn edges into predlist
%
% nCnt = hist(edges(:),uN);
% tips = uN(nCnt==1);
%
% for t = 1:length(tips)
%     t
%     pred(:,t) = edges2pred(edges,tips(t));
% end
%

%% Pred list
pred = edges2pred(edges,seed);


pre = pred;
showN = seed;
for p = 1:length(pre)
    
    scatter3(pos(:,1),pos(:,2),pos(:,3),10,'o','k','filled');
    hold on
    
    pre(showN) = 0;
    
    hit = [];
    for s = 1:length(showN)
        hit = cat(1,hit,find(pre==showN(s)));
    end
    
    showN = hit;
    scatter3(pos(showN,1),pos(showN,2),pos(showN,3),150,'o','r','filled');
    hold off
    pause(.01)
    if isempty(showN)
        break
    end
end

%% define special nodes


uP = 1:length(pred);
pCnt = hist(pred,uP);
tips = uP(pCnt==0);
branches = setdiff(uP(pCnt==2),seed);

scatter3(pos(uP,1),pos(uP,2),pos(uP,3),10,'o','k','filled');
hold on
scatter3(pos(seed,1),pos(seed,2),pos(seed,3),50,'o','b','filled');
scatter3(pos(tips,1),pos(tips,2),pos(tips,3),50,'o','r','filled');
scatter3(pos(branches,1),pos(branches,2),pos(branches,3),50,'o','g','filled');
hold off


%% check branches
nBin = 5; % number of nodes to bin for radius measure
clear uRad uVar dRad
for b = 1:length(branches)
    
    branchN = branches(b);
    
    %%Get nodes downstream of branch
    downN = find(pred==branchN);
    for d = 1:length(downN)
        showN = downN(d);
        downListN{d} = showN;
%         scatter3(pos(uP,1),pos(uP,2),pos(uP,3),10,'o','k','filled');
%         hold on
        for n = 1:nBin
            hit = [];
            for s = 1:length(showN)
                hit = cat(1,hit,find(pre==showN(s)));
            end
            showN = hit;
            downListN{d} = cat(1,downListN{d},showN);
%             scatter3(pos(downListN{d},1),pos(downListN{d},2),pos(downListN{d},3),50,'o','r','filled');
%             pause(.01)
            if isempty(showN)
                break
            end
        end
        hold off
        pause(.1) 
        
        dRads = rad(downListN{d});
        dRad(b,d) = mean(dRads);
        
    end
    
    %%Get list of nodes upstream of branch
    upList = [];
    showN = branches(b);
    for n = 1:nBin
        showN = pred(showN);
        upList = cat(1,upList,showN);
    end
    
    %%Get rads
    upRads = rad(upList);
    uVar(b) = var(upRads);
    uRad(b) = mean(upRads);
    
end

scatter(uRad,dRad(:,1),'markeredgealpha',0,'markerfacecolor','k','markerfacealpha',.2)
hold on
scatter(uRad,dRad(:,2),'markeredgealpha',0,'markerfacecolor','k','markerfacealpha',.2)
scatter(dRad(:,1),dRad(:,2),'markeredgealpha',0,'markerfacecolor','r','markerfacealpha',.2)
scatter(dRad(:,2),dRad(:,1),'markeredgealpha',0,'markerfacecolor','r','markerfacealpha',.2)

line([0 1.2],[0 1.2])
ylim([0 1.5])
xlim([0 1.5])
hold off

%% compare 
useN = pos(branches,3) < 30; %select nodes to use by stratification

dRad32 = dRad(useN,1).^3/2 + dRad(useN,2).^3/2;
dRad32mean = (dRad(useN,1).^3/2 + dRad(useN,2).^3/2)/2;
uRad32 = uRad(useN).^3/2;
scatter(uRad32,dRad32,'markeredgealpha',0,'markerfacecolor','k','markerfacealpha',.2)
hold on
scatter(uRad32,dRad32mean,'markeredgealpha',0,'markerfacecolor','r','markerfacealpha',.2)
hold off
line([0 1.2],[0 1.2])
ylim([0 1])
xlim([0 1])


difSum = (dRad32-uRad32)/uRad32;
difMean = (dRad32mean-uRad32)/uRad32;

rng = -1:.1:2;
hSum = hist(difSum,rng);
hMean = hist(difMean,rng);

bar(rng,[hSum ;hMean]')


%% compare depth to radius

prop = pos(seed,3) - pos(:,3);
scatter(prop,rad,'markeredgealpha',0,'markerfacecolor',...
    'k','markerfacealpha',.2)

%% compare radius to topo distance

linDist = sm.skel2skel.linDist(seed,:);
eucDist = sm.skel2skel.eucDist(seed,:);
scatter(linDist,rad,'markeredgealpha',0,'markerfacecolor',...
    'k','markerfacealpha',.4)
hold on
% scatter(eucDist,rad,'markeredgealpha',0,'markerfacecolor',...
%     'k','markerfacealpha',.2)
rng = [0:2:60];
[N,edges,bin] = histcounts(linDist,rng)

bRad = 2;
rng = 0:1:100;
for i = 1:length(rng)
    isThere = (linDist>(rng(i)-bRad) & (linDist<=(rng(i)+bRad)));
    meanRad(i) = mean(rad(isThere));
    
end
plot(rng,meanRad,'linewidth',4,'color','k')
hold off


showN  = linDist<=10;
dotSize = rad.^2*80;
scatter3(pos(:,1),pos(:,2),pos(:,3),dotSize,'marker','o','markerfacecolor',[0 .7 0],'markeredgecolor','k');
hold on
scatter3(pos(showN,1),pos(showN,2),pos(showN,3),dotSize(showN),'marker','o','markerfacecolor',[.8 0 .8],'markeredgecolor','k');
hold off
grid off
axis off
axis square
set(gca,'clipping', 'off')
set(gcf,'color','w')
%%

pre = pred;
showN = seed;
for p = 1:length(pre)
    
    scatter3(pos(:,1),pos(:,2),pos(:,3),10,'o','k','filled');
    hold on
    
    pre(showN) = 0;
    
    hit = [];
    for s = 1:length(showN)
        hit = cat(1,hit,find(pre==showN(s)));
    end
    
    showN = hit;
    scatter3(pos(showN,1),pos(showN,2),pos(showN,3),150,'o','r','filled');
    hold off
    pause(.01)
    if isempty(showN)
        break
    end
end




big = rad>1.1;
prop = rad;
prop(prop>2) = 2;
prop = prop * 100;
cm = jet(ceil(max(prop)));
s1 = scatter3(pos(uP,1),pos(uP,2),pos(uP,3),10,cm(ceil(prop),:),'filled');













