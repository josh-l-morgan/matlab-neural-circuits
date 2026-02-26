%%Show distribution of synapses on TCs
%%First make SMs then run loadDataFromKxS.mat

global tis
load('MPN.mat');
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])


fig = figure;
ax = subplot(1,1,1,'parent',fig)

for o = 1:length(tis.syn.obID)
    nam = obI.nameProps.names{tis.syn.obID(o)};
    findRGC = regexp(lower(nam),'rgc');
    if ~isempty(findRGC)
        tis.syn.preClass(o) = 1;
    end
end


for s = 1 : length(sms)

    cla(ax)
    sm = sms(s).sm;
    pos = sm.nep.pos;
    scatter3(ax,pos(:,1),pos(:,2),pos(:,3),'.','k')
    
    synToCid = find(tis.syn.post==sm.cid);
    sPos = tis.syn.pos(synToCid,:);
    preType = tis.syn.preClass(synToCid);
    isRGC = preType== 1;
    notRGC = ~isRGC;
    hold on
    scatter3(ax,sPos(isRGC,1),sPos(isRGC,2),sPos(isRGC,3),'o','lineWidth',1,'MarkerEdgeColor','g')
    scatter3(ax,sPos(notRGC,1),sPos(notRGC,2),sPos(notRGC,3),'o','lineWidth',1,'MarkerEdgeColor','r')
    title(ax,sprintf('cell cid%d',sm.cid))
  


    %%Mark cell body
    cidID = find(tis.cids==sm.cid,1);
    anchor = tis.cells.anchors(cidID,:);
    aPos = double(anchor([2 1 3])) .* obI.em.res/1000;
    scatter3(ax,aPos(1),aPos(2),aPos(3),100,'MarkerFaceColor','b')
    scatter3(ax,aPos(1),aPos(2),aPos(3),400,'MarkerFaceColor','none','MarkerEdgeColor','b','lineWidth',2)


  pause

end
close(fig)




