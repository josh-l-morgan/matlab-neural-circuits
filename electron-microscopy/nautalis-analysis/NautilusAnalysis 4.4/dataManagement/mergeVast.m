
clear all
SPN = 'E:\IxQ_KarlsRetinaVG3_2019\Merge\'; %GetMyDir;

%SPN = GetMyDir;
MPN = [SPN(1:end-1) '_mat\']
if ~exist(MPN,'dir'),mkdir(MPN),end
         

%% parse names
subNames = dir([SPN '*vastSubs*.mat']);
subNames = {subNames.name}
mainFile = 1;
for i = 1:length(subNames);
    if strcmp(lower(subNames{i}),'vastsubs_joshm1.mat')
        mainFile = i;
        break
    end
end

vastNames = {subNames{mainFile} subNames{setdiff([1:length(subNames)],mainFile)}}
fuse.SPN = SPN;
for i = 1:length(vastNames)
    subTag = vastNames{i}(9:end-4);
    fuse.files(i).subTag = subTag;
    fuse.files(i).vastSub = vastNames{i};
    fuse.files(i).obI = ['obI' subTag '.mat'];
    
end

%% merge subs
lastInd = 0;
allSubs = {};
clear allOb
allOb.fuse = fuse;
obSource = [];
for i = 1:length(fuse.files)
    load([SPN fuse.files(i).vastSub]);
    load([SPN fuse.files(i).obI]);
    
    numOb = length(obI.colStruc.names);
    
    obSource(length(obSource)+1:length(obSource)+numOb) = i;
    %     for n = 1:numOb
    %     allOb.colStruc.names{startInd+n-1} = obI.colStruc.names{n};
    %     allOb.colStruc.ids{startInd+n-1,:} = obI.colStruc.ids{n}+startInd-1;
    %     allOb.colStruc.anchors(startInd+n-1,:) = obI.colStruc.anchors(n,:);
    %     allOb.colStruc.boundBox(startInd+n-1,:) = obI.colStruc.boundBox(n,:);
    %     allOb.colStruc.col1(startInd+n-1,:) = obI.colStruc.col1(n,:);
    %     allOb.colStruc.col2(startInd+n-1,:) = obI.colStruc.col2(n,:);
    %     end
    
    if i == 1
        allOb.colStruc.names = obI.colStruc.names;
        allOb.colStruc.ids = obI.colStruc.ids;
        allOb.colStruc.anchors = obI.colStruc.anchors;
        allOb.colStruc.boundBox =  obI.colStruc.boundBox;
        allOb.colStruc.col1 = obI.colStruc.col1;
        allOb.colStruc.col2 = obI.colStruc.col2;
        
        allSubs =vastSubs(1:numOb);
    else
        
        allOb.colStruc.names = cat(2,allOb.colStruc.names,obI.colStruc.names)
        
        ids = obI.colStruc.ids;
        for n = 1:size(ids,1)
           ids{n,1} = ids{n,1} + lastInd; 
        end
        allOb.colStruc.ids = cat(1,allOb.colStruc.ids,ids);
        allOb.colStruc.anchors = cat(1,allOb.colStruc.anchors,obI.colStruc.anchors)
        allOb.colStruc.boundBox = cat(1,   allOb.colStruc.boundBox, obI.colStruc.boundBox)
        allOb.colStruc.col1 = cat(1,allOb.colStruc.col1,obI.colStruc.col1);
        allOb.colStruc.col2 = cat(1,allOb.colStruc.col2,obI.colStruc.col2);
        
        allSubs = cat(2,allSubs,vastSubs(1:numOb));
        
    end
    
    
    lastInd = lastInd + numOb;
end

allOb.fuse.obSource = obSource;
obI = allOb;
vastSubs = allSubs;

%%

obI.nameProps = getNameProps(obI.colStruc.names);

cellIDs = obI.nameProps.cellNum;
allNames = unique(cellIDs);
c = 0;
for i = 1:length(allNames)
    obIDs = find((cellIDs == allNames(i)) & obI.nameProps.ofID);
    
    if ~isempty(obIDs)
        c = c+1;
        obI.cell.name(c) = allNames(i);
        obI.cell.obIDs{c} = obIDs;
        mainCell = obI.nameProps.cell(obIDs);
        if sum(mainCell) & (i>9)
            mainObID = obIDs(find(mainCell>0,1,'first'));
            obI.cell.isCell(c) = 1;
        else
            
            mainObID = obIDs(1);
            obI.cell.isCell(c) = 0;
        end
        obI.cell.anchors(c,:) = obI.colStruc.anchors(mainObID,:);
        obI.cell.mainObID(c) = mainObID;
    end
end


save([MPN 'vastSubs.mat'],'vastSubs','-v7.3')
save([MPN 'obI.mat'],'obI')

%% Down sample vastSubs


res = [6 4 30]; 
vRes = [4 4 1] .* [8 8 1].* res;
dsRes = [200 200 200];
dsDim = dsRes./vRes;

% useSaved = 0;
% if useSaved & exist([MPN 'dsObj.mat'],'file')
%     load([MPN 'dsObj.mat'])
% else
dsObj = downSampObj(MPN, dsDim);
% end


obI.em.res = res; %or [4.6 4 30]?
obI.em.vRes =vRes;
obI.em.dsRes =dsRes/1000;
save([MPN 'obI.mat'],'obI')

return

%% make golgi library
cellList = obI.cell.name;
%cellList = [10	129	162	170]

colMap = hsv(256);
for d = 1:3
    Dim = d;
    fsize = max(cat(1,dsObj.subs),[],1);
    golgiDir = sprintf('%sgolgi%d\\',MPN,Dim);
    if ~exist(golgiDir,'dir'),mkdir(golgiDir),end
    
    parfor i = 1:length(cellList)
    
    %i = find(cellList == 107)
    
    showCellNames = cellList(i);
    
    fileName = sprintf('golgi_%05.0f_dim%d.png',showCellNames,Dim);
    if 1%~exist([golgiDir fileName],'file')
        col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
        col = col(randperm(size(col,1)),:);
        col = [1 1 1];
        I_partCell = showCellSum(obI,dsObj,showCellNames,col,Dim,fsize);
        I_partCell = 256-I_partCell*3;
        image(uint8(I_partCell))
        disp(showCellNames)
        fileName = sprintf('golgi_%05.0f_dim%d.png',showCellNames,Dim);
        imwrite(uint8(I_partCell),[golgiDir fileName])
        pause(.01)
    end
    
    
    end
end

%% make OBJ
cellList = obI.cell.name;
cellList = getList_tracedCells;
renderOb = 1;

objDir = [MPN 'obFiles\']
if ~exist(objDir,'dir'),mkdir(objDir),end

cellSubs = names2Subs(obI,dsObj,cellList);
downSamp = 4;
for i = 1:length(cellSubs)
    sub = cellSubs{i};
    obName = cellList(i);
    if ~ischar(obName),obName = num2str(obName);end
    if ~isempty(sub)
    smallSub = shrinkSub(sub,downSamp);
   
    fv = subVolFV(smallSub,[],renderOb);
    fileName = sprintf('%sdSamp%d_%s.obj',objDir,downSamp,obName);
    vertface2obj(fv.vertices,fv.faces,fileName,obName);
   
    cellDat(i).subs = sub;
    cellDat(i).fv = fv;
    end
    
end




