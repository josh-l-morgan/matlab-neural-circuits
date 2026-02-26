%%Define Island and nonIsland
%%First run loadDataFromKxS


fig = figure;
ax = subplot(1,1,1,'parent',fig);

%% Find object ID for planes defining exlusion zone and island
islandName = 'Island top plane';
exclusionName = 'Exclusion top plane';

oNames = obI.nameProps.names;
for i = 1:length(oNames)
    isReg = regexp(oNames{i},islandName);
    if ~isempty(isReg)
       islandPlaneObID = i;
    end
    isReg = regexp(oNames{i},exclusionName);
    if ~isempty(isReg)
       exclusionPlaneObID = i;
    end
end

islandSubs = dsObj(islandPlaneObID).subs *  obI.em.dsRes;
exclusionSubs = dsObj(exclusionPlaneObID).subs *  obI.em.dsRes;

%% Get cell synapses


fromRGC = tis.syn.preClass == 1;
notFromRGC = ~fromRGC;
sPosAll = tis.syn.pos;

cla(ax)
ax.NextPlot = 'add';
scatIsland = scatter3(ax,islandSubs(:,2),islandSubs(:,1),islandSubs(:,3),5,'markeredgecolor','none','markerfacecolor','flat');
scatIsland.CData = [0 1 0];
scatIsland.MarkerFaceAlpha = .1;
scatExlusion = scatter3(ax,exclusionSubs(:,2),exclusionSubs(:,1),exclusionSubs(:,3),5,'markeredgecolor','none','markerfacecolor','flat');
scatExlusion.CData = [1 0 1];
scatExlusion.MarkerFaceAlpha = .1;
scatter3(ax,sPosAll(fromRGC,2),sPosAll(fromRGC,1),sPosAll(fromRGC,3),8,'markeredgecolor','none','markerfacecolor','b');
%scatter3(ax,sPosAll(notFromRGC,2),sPosAll(notFromRGC,1),sPosAll(notFromRGC,3),3,'markeredgecolor','none','markerfacecolor','k');


for s = 1 : length(sms)
    sm = sms(s).sm;
    pos = sm.nep.pos;
    sPos = sm.syn.pos;
    preType = sm.syn.preClass;
    isRGC = preType== 1;
    notRGC = ~isRGC;
   
    title(ax,sprintf('cell cid%d',sm.cid))
    s1 = scatter3(pos(:,2),pos(:,1),pos(:,3),5,'markeredgecolor','none','markerfacecolor','k')
    s2 = scatter3(sPos(isRGC,2),sPos(isRGC,1),sPos(isRGC,3),20,'markeredgecolor','none','markerfacecolor','g')
    s3 = scatter3(sPos(notRGC,2),sPos(notRGC,1),sPos(notRGC,3),10,'markeredgecolor','none','markerfacecolor','r')

    pause

    s1.delete
    s2.delete
    s3.delete

end















