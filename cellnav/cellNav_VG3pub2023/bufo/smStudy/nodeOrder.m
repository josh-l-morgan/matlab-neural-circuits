function[ord] = nodeOrder(sm,seed)


%% Get properties
pos = sm.nep.pos;
edges = sm.nep.edges;
groupEdges = sm.nep.groupEdges;
bones = sm.nep.bones;
rad = sm.nep.nodeRad;
rad = sm.nep.meanNodeRad;

%% Pick seed
global globPt

if ~exist('seed','var')
    
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
    
else
    
end

ord.seed = seed;

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

%% Make pred list.  node previous to each node
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

ord.pred = pred;

%% define special nodes

ord.uP = 1:length(pred);
ord.pCnt = hist(pred,ord.uP);
ord.tips = uP(ord.pCnt==0);
ord.branches = setdiff(uP(ord.pCnt==2),seed);

scatter3(pos(ord.uP,1),pos(ord.uP,2),pos(ord.uP,3),10,'o','k','filled');
hold on
scatter3(pos(seed,1),pos(seed,2),pos(seed,3),50,'o','b','filled');
scatter3(pos(ord.tips,1),pos(ord.tips,2),pos(ord.tips,3),50,'o','r','filled');
scatter3(pos(ord.branches,1),pos(ord.branches,2),pos(ord.branches,3),50,'o','g','filled');
hold off



%% check branches
nBin = 5; % number of nodes to bin for radius measure
clear uRad uVar dRad
for b = 1:length(branches)
    hold off
%     scatter3(pos(:,1),pos(:,2),pos(:,3),10,'o','k','filled');
%     hold on
    branchN = branches(b);
    
    %%Get nodes downstream of branch
    downN = find(pred==branchN);
    for d = 1:length(downN)
        showN = downN(d);
        downListN{d} = showN;
        for n = 1:nBin
            hit = [];
            for s = 1:length(showN)
                hit = cat(1,hit,find(pred==showN(s)));
            end
            showN = hit;
            downListN{d} = cat(1,downListN{d},showN);
            
            if isempty(showN)
                break
            end
        end
        pause(.1)
        
        dRads = rad(downListN{d});
        dRad(b,d) = mean(dRads);
        %scatter3(pos(downListN{d},1),pos(downListN{d},2),pos(downListN{d},3),50,'o','r','filled');
    end
    
    %%Get list of nodes upstream of branch
    upList = [];
    showN = branches(b);
    for n = 1:nBin
        showN = pred(showN);
        upList = cat(1,upList,showN);
    end
    upList = cat(1,upList,showN);
    %scatter3(pos(upList,1),pos(upList,2),pos(upList,3),50,'o','g','filled');
    %%Get rads
    upRads = rad(upList);
    uVar(b) = var(upRads);
    uRad(b) = mean(upRads);
    
    ord.branch(b).uVar = uVar(b);
    ord.branch(b).uRad = uRad(b);
    ord.branch(b).downList = downListN;
    ord.branch(b).upList = upList;
    %pause(.02)
end

%% Plot split
hold off
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


%%  get order
% 
% minLength = 5; %minimum length
% 
% tipLength = [];
% for t = 1:length(ord.tips);
%     last = ord.tips(t);
%     L = 0;
%     for i = 1:length(pCnt);
%         next = pred(last);
%         L = L + sqrt((pos(last,1)-pos(next,1)).^2 + (pos(last,2)-pos(next,2)).^2 + ...
%             (pos(last,3)-pos(next,3)).^2);
%         if pCnt(next) == 1;
%             last = next;
%         else
%             tipParent(t) = next;
%             break
%         end
%     end
%     tipLength(t) = L;
% end




%% Simple branch order 

pre = pred;
showN = seed;
ordN = pred * 0;
ordN(showN) = 0;
for p = 1:length(pre)
    
    scatter3(pos(:,1),pos(:,2),pos(:,3),10,'o','k','filled');
    hold on
    
    pre(showN) = 0;
    
    hit = [];
    for s = 1:length(showN)
        h = find(pre==showN(s));
        if length(h) == 1
            ordN(h) = ordN(showN(s));
        elseif length(h) > 1
            ordN(h) = ordN(showN(s))+1;
        end
        hit = cat(1,hit,h);
    end
    
    showN = hit;
    scatter3(pos(showN,1),pos(showN,2),pos(showN,3),150,'o','r','filled');
    hold off
    pause(.01)
    if isempty(showN)
        break
    end
end
% 
% 
% prop = ordN+1;
% cm = jet(ceil(max(prop)));
% s1 = scatter3(pos(uP,1),pos(uP,2),pos(uP,3),20,cm(ceil(prop),:),'filled');
% 
% ord.ordN = ordN;











