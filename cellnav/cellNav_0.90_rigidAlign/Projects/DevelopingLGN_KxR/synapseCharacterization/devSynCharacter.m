
 global  glob tis
%% Clear dat
if 0
    dat = {};
    %%Paste KxR cell Registry -> SynapseTyping
s = parseExcelSynProps(dat)

end



fs = figure;
ax = subplot(1,1,1,'parent',fs);
ax.NextPlot = 'add';
ax.Clipping = 'off';
fs.Color = [0 0 0];
ax.Color = [0 0 0];

allPost = setdiff(unique(s.postCid),0);
allProp = setdiff(unique(s.postCid(s.hasSynProp==1)),0);




showCells = 3209;




clear ps
for i = 1:length(showCells)

    cla(ax)
    cid = showCells(i);
    
    %% Draw cell
    fileName = sprintf('%s%d.mat',glob.useFvDir,showCells(i));
    fv = loadFV(fileName);
    try
        delete(ps.cell(i))
    end
    volName = glob.vol.activeName;
    if exist([glob.dir.Volumes volName '\volTransform.mat'],'file');
        load([glob.dir.Volumes volName '\volTransform.mat']);
    else
        volTransform = [];
    end
    ps.cell(i) = renderFVnav(fv,[1 1 1],.1,cid,volTransform,ax);

    %% Pick synapses
    mSize = 60;
    isCid = s.postCid == cid;
    allSyn = isCid & s.isPoint;
    isRGC = allSyn & (s.synType == 1);
    hasGlia = allSyn & (s.hasGlia == 1);
    hasSpine = allSyn & (s.hasSpine == 1);
    isLarge = allSyn & (s.large == 1);
    
   
    scatter3(s.pos(allSyn,1),s.pos(allSyn,2),s.pos(allSyn,3),ceil(mSize^(.7)),'markerfacecolor',[1 0 0],'markeredgecolor','none');
    scatter3(s.pos(isRGC,1),s.pos(isRGC,2),s.pos(isRGC,3),mSize,'markerfacecolor',[0 1 0],'markeredgecolor','none');
    scatter3(s.pos(hasGlia,1),s.pos(hasGlia,2),s.pos(hasGlia,3),mSize^1.5,'markerfacecolor','none','markeredgecolor',[1 1 0]);
    scatter3(s.pos(hasSpine,1),s.pos(hasSpine,2),s.pos(hasSpine,3),mSize^1.5,'+','markerfacecolor','none','markeredgecolor',[0 1 1]);
    scatter3(s.pos(isLarge,1),s.pos(isLarge,2),s.pos(isLarge,3),mSize^1.3,'s','markerfacecolor','none','markeredgecolor',[0 1 0]);


    disp(sprintf('cid%d',cid))
    drawnow
    pause


end


