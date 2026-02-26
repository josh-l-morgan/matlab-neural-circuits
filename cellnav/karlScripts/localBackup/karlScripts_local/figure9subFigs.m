%% loading
highTis=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\tis.mat');
highTis=highTis.tis;
highFV='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
curTis=highTis;
%% lists n such
allVGCids=type2cid({'amc'},{'vgc'},curTis);
allVGCids=allVGCids{1};
alPre=cid2type(curTis.syn.edges(:,2),curTis);
alPost=cid2type(curTis.syn.edges(:,1),curTis);
sbCids=unique(curTis.syn.edges(find(ismember(curTis.syn.edges(:,1),allVGCids)&alPre{1}'==7),2));
allbCids=unique(curTis.syn.edges(alPre{1}'==7,2));
nsbCids=setdiff(allbCids,sbCids);
%% parameters
cidList=[2002 3051 3119];
%vin=yellow; sb=cyan; nsb=mag;
synPlotCols=[[0 1 1];[1 1 0];[1 0 1]];
markSize=40;
barlocs=[[155,155];[165,150];[180,155]];
barSize=10;
barWidth=2;
%% loop
figs={};
patches={};
synPlots={};
for cidIt=1:length(cidList)
    curCid=cidList(cidIt);
    %figs{cidIt}=figure();
    curMorphFig=compareMorph([],curCid,highFV);
    figs{cidIt}=curMorphFig;
    
end
    
for figIt=1:length(figs)
    curFig=figs{figIt};
    curCid=cidList(figIt);
    figure(curFig.figNam);
    view(90,180);
    zlim=curFig.axes.ZLim;
    xlim=curFig.axes.XLim;
    ylim=curFig.axes.YLim;
    %put scale bar in Y plane
    scaleBar=[[barlocs(figIt,1) barlocs(figIt,1)];[barlocs(figIt,2) barlocs(figIt,2)+barSize];[mean(zlim) mean(zlim)]];
    %scaleBar=[[mean(zlim) mean(zlim)];[ylim(2)-15 ylim(2)-5];[ylim(2)-8 ylim(2)-8]];
    sBar=plot3(scaleBar(3,:),scaleBar(2,:),scaleBar(1,:),'LineWidth',barWidth,'Color','white'); %'Color','cyan',
    curFig.texts.Visible=0;
    curFig.figNam.Color=[0 0 0];
    curFig.patches.FaceColor=[1 1 1];
    vIn=find(curTis.syn.edges(:,1)==curCid&ismember(curTis.syn.edges(:,2),allVGCids));
    sbIn=find(curTis.syn.edges(:,1)==curCid&ismember(curTis.syn.edges(:,2),sbCids));
    nsbIn=find(curTis.syn.edges(:,1)==curCid&ismember(curTis.syn.edges(:,2),nsbCids));
    synPlotInds={vIn,sbIn,nsbIn};
    
    for i=1:3
        curInds=synPlotInds{i};
        curCol=synPlotCols(i,:);
        scatter3(curTis.syn.pos(curInds,3),curTis.syn.pos(curInds,1),curTis.syn.pos(curInds,2), ...
            markSize,curCol,'o','filled');
    end
    
    curFig.figNam.Position = [100+(500*(figIt-1)) 100 600 800];
    curFig.axes.Position = [0 0 1 1];
    filNam2=['C:\work\figs\fig9\' num2str(curCid) '_prnt.png'];
    filNam=['C:\work\figs\fig9\' num2str(curCid) '_prnt.eps'];
    filNam3=['C:\work\figs\fig9\' num2str(curCid) '_ex.eps'];
    %hgexport(curFig.figNam,filNam);
    %hgexport(curFig.figNam,filNam3);
    %exportgraphics(curFig.figNam,filNam3,'ContentType','vector')
    set(curFig.figNam, 'InvertHardcopy', 'off')
    print(curFig.figNam,filNam2,'-r1200','-dpng');
    print(curFig.figNam,filNam,'-r1200','-depsc');
    
end