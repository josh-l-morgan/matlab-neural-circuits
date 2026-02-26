

MPN = GetMyDir;
load([MPN 'dsObj.mat']);
load([MPN 'obI.mat']);


%%

[cellIDs cellProps] = getList_prinCompRats;
[cellIDs2 cellProps2] = getList_pcLatRats;
scatter(cellProps(:,1),cellProps(:,2))
hold on
scatter(cellProps2(:,1),cellProps2(:,2),'r')
hold off
ylim([0 1]);
xlim([0 1]);





%%


col = [1 1 1];

fsize = max(cat(1,dsObj.subs),[],1);

viewProps.dim = 1;
viewProps.maxScaleFactor = .1;
viewProps.sumScaleFactor = 3;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.viewWindow = double([1 1 1; fsize]);


%% Display Cells

dsamp = 8;
spread = 4000;
clf
xlim([0 spread]);
ylim([0 spread]);
hold on
for i = 1:length(cellIDs)
    
    viewProps.cellId = {cellIDs(i)};
    
    Iraw = showCellsAndMore(viewProps);
    [ind] = find(Iraw>0);
    [y x z] = ind2sub(size(Iraw),ind);
%     I = Iraw(min(y):max(y),min(x):max(x),:);
%     I = permute(I,[2 1 3]);
%     image(uint8(I*1000))
    
    y = ceil(y/dsamp);
    x = ceil(x/dsamp);
    
    ind = sub2ind([fsize(1) fsize(2)],y,x);
    ind = unique(ind);
    [y x] = ind2sub([fsize(1) fsize(2)],ind);
    
    cx = y + cellProps(i,1) * spread;
    cy = x + cellProps(i,2) * spread;
    scatter(cx,cy,5,'.')
    
    %     cellProps(i,:)
    %     cellProps2(i,:)
    pause(.01)
    
end
hold off




