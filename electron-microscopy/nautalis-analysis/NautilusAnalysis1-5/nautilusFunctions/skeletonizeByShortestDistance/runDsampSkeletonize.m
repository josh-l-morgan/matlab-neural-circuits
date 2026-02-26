MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\export+14+04+27_mat\'

TPN = [MPN 'skel\']
if ~exist(TPN,'dir'),mkdir(TPN);end
TPNview = [TPN 'view\'];
TPNmat = [TPN 'mat\'];
if ~exist(TPNview,'dir'),mkdir(TPNview);end
if ~exist(TPNmat,'dir'),mkdir(TPNmat);end




load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])

cellList = obI.cell.name;


cellList = 1033;
for i = 1:length(cellList)
    disp(sprintf('skeletonizing cell %d of %d',i,length(cellList)))
    cellTarg = cellList(i);
   
    objectSubs = getCellSubs(obI,dsObj,cellTarg);
    %connectSubs(objectSubs);
    
    if size(objectSubs,1)>100
    startTime = clock;
    
    pass = 1;
    try
    cellStruct = subs2bones(objectSubs);
    catch err
        err
        pass = 0;
    end
    if pass
    cellView = cellStruct.sideViews{1};
    skel = cellStruct.skel;
    
    skelFile = sprintf('%s%d.m',TPNmat,cellTarg);
    viewFile = sprintf('%s%d.png',TPNview,cellTarg);
    save(skelFile,'skel')
    imwrite(cellView,viewFile)
    end
   
    startTime
    stopTime = clock
     end
end








