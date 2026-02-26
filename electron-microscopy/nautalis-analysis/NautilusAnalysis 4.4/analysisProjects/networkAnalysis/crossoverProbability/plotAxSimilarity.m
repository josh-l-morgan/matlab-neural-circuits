%function[] = seedPreferences(seedCells, useList);

axList = useList.preList;
cellList = useList.postList;
con = useList.con;


plotConDir = [MPN 'connectivity\compareToSeed\'];
if ~exist(plotConDir,'dir'),mkdir(plotConDir),end

%%

seedCol = {'r','b','g','c','m','k'}
maxCon = max(con(:));

for i = 1: length(cellList)
        hold off

    
    
    if ~(sum(seedCells == cellList(i)))
        
    for s = 1: length(seedCells)
            targ = find(cellList==seedCells(s));
    
        seedAx = con(:,targ);
        cellAx = con(:,i);
        scatter(jitter(seedAx,.2),jitter(cellAx,.2),seedCol{s},'o')
        title(sprintf('synapses each axon forms with cell %d compared to seed cell',...
            cellList(i)))
        text(maxCon* .8, maxCon * .8+s,num2str(seedCells(s)),'edgecolor',seedCol{s})
        ylabel('thalamocortical test cell')
        xlabel('seed cell')
        hold on
    end
    end
    ylim([-0.5 maxCon]);
    xlim([-0.5 maxCon]);
pause(.1)

imageName = sprintf('%scompare2seed_%03.0f',plotConDir,cellList(i));
saveas(gcf,imageName,'png')
end
    hold off

