

SPN = GetMyDir;
MPN = [SPN(1:end-1) '_mat\']
if ~exist(MPN,'dir'),mkdir(MPN),end

%% load paths
path(path,sprintf('%s\\analyzeMorphology\\skeletonizeByShortestDistance',pwd))
path(path,sprintf('%s\\viewCells',pwd))
path(path,sprintf('%s\\analyzeColors',pwd))
path(path,sprintf('%s\\analyzeConnectivity',pwd))
path(path,sprintf('%s\\vast2ob',pwd))
path(path,sprintf('%s\\sub2mesh',pwd))
path(path,sprintf('%s',pwd))



%% parse names
subNames = dir([SPN '*vastSubs*.mat']);
subNames = {subNames.name}
for i = 1:length(subNames);
    if sum(regexp(subNames{i},'jm'))
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
startInd = 1;
allSubs = {};
clear allOb
allOb.fuse = fuse;
obSource = [];
for i = 1:length(fuse.files)
    load([SPN fuse.files(i).vastSub]);
    load([SPN fuse.files(i).obI]);
        numOb = length(obI.colStruc.names);
    obSource(length(obSource)+1:length(obSource)+numOb) = i;
    for n = 1:numOb
    allOb.colStruc.names{startInd+n-1} = obI.colStruc.names{n};
    allOb.colStruc.ids{startInd+n-1,:} = obI.colStruc.ids{n}+startInd-1;
    allOb.colStruc.anchors(startInd+n-1,:) = obI.colStruc.anchors(n,:);
    allOb.colStruc.boundBox(startInd+n-1,:) = obI.colStruc.boundBox(n,:);
    allOb.colStruc.col1(startInd+n-1,:) = obI.colStruc.col1(n,:);
    allOb.colStruc.col2(startInd+n-1,:) = obI.colStruc.col2(n,:);
    end
        
  allSubs = cat(2,allSubs,vastSubs);
  startInd = startInd + numOb;
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
useSaved = 0;
if useSaved & exist([MPN 'dsObj.mat'],'file')
    load([MPN 'dsObj.mat'])
else
   dsDim = [1 1 4];
   dsObj = downSampObj(MPN, dsDim);
end


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
        showCellNames = cellList(i);

        fileName = sprintf('golgi_%05.0f_dim%d.png',showCellNames,Dim);
 if ~exist([golgiDir fileName],'file')
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



