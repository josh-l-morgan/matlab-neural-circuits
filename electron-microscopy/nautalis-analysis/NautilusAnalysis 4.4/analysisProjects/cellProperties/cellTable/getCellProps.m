function[cellList cellProp] = getAttributes(obI);

%MPN = GetMyDir

%load([MPN 'dsObj.mat'])
%load([MPN 'obI.mat'])


%% get cell list
  useList = obI2cellList_all(obI);
  cellList = useList.postList;
  cellNum = length(cellList);
  
%% Cell types

cellClass = zeros(cellNum,1);

cellClass(obI.nameProps.rgc(obI.cell.mainID)) = 1;
cellClass(obI.nameProps.tcr(obI.cell.mainID)) = 2;
cellClass(obI.nameProps.lin(obI.cell.mainID)) = 3;
  

%% get connectivity

  allEdges = obI.nameProps.edges(:,[2 1]);

  seedList = [108 201 109 903 907 ];
  conList = obI2cellList_seedInput(obI,seedList);
  seedPref = seedPreferences(seedList,conList);
  conTo = makeConTo(obI,seedList);
  
  seedLink = zeros(length(cellList),length(seedList));
  for i = 1:length(seedPref.cellList)
    targ = find(cellList == seedPref.cellList(i));
    if ~isempty(targ)
        seedLink(targ,:) = seedPref.sharedAx(:,i);
    end   
  end
  
  preSeed = zeros(length(cellList),length(seedList));
  for s = 1:length(seedList)
      targ = find(cellList == seedList(s));
      preSeed(:,s) = useList.con(:,targ);
  end
  
  cell2seed = preSeed;
  cell2seed(cellClass == 2,:) = seedLink(cellClass == 2,:);
  
%% get anchors
hold on

dSamp = [.016 0.016 0.030];

anchors = obI.colStruc.anchors;
anchors(:,1) = anchors(:,1).*dSamp(1);
anchors(:,2) = anchors(:,2).*dSamp(2);
anchors(:,3) = anchors(:,3).*dSamp(3);
anchors = round(anchors);
anchors(anchors<1) = 1;


synapses = obI.nameProps.edges;
synAnchors = anchors(synapses(:,3),:);

edges = synapses(:,1:2);




%%  Get possitions

ancPos = zeros(length(cellList),3);
for i = 1:length(cellList)
    targ = find(obI.cell.name == cellList(i));
    mainID = obI.cell.mainID(targ);
    ancPos(i,:) = anchors(mainID,:);
end
%scatter(ancPos(:,2),ancPos(:,1),30,'filled')

center = ancPos(cellList == 108,:);

conPos = ancPos;
for i = 1:length(cellList);
   isPre = synapses(:,2) == cellList(i);
   isPost = synapses(:,1) == cellList(i);
   isEither = isPre+ isPost;
   if sum(isEither)
       conPos(i,:) = mean(synAnchors(find(isEither),:),1);
   end
end
%scatter(conPos(:,2),conPos(:,1),30,'filled')


cellPos = ancPos;
cellPos(cellClass == 1,:) = conPos(cellClass == 1,:);



dist2center = sqrt((cellPos(:,1)-center(1)).^2 + ...
    (cellPos(:,2)-center(2)).^2 + ...
    (cellPos(:,3)-center(3)).^2);







%% cell Color

cellCol = zeros(length(cellList),3);

colParam = (cellPos(:,1) - center(1)+75)/150;
colParam = (cellPos(:,3) - center(3)+75)/150;

cellCol(:,1) = colParam;
cellCol(:,3) = 1-colParam;
%cellCol(find(cellClass == 2),2) = 1;


%cellCol(cellClass == 1,:) = .3;

cellCol(cellCol<0) = 0;
cellCol(cellCol>1) = 1;



cellAtt = cat(2,cellList',cellCol,cellClass,cellPos)


[cellList' cellClass,cellCol]

%% Make table


for i = 1:length(seedList)
    seedListCell{i} = num2str(seedList(i));
end
header = cat(2,{'cellList' 'dist2center'}, seedListCell)
rawTable  = cat(2,cellList', dist2center, cell2seed);




sortParam = dist2center;
sortParam(cellClass == 1) = sortParam(cellClass == 1) * 1000;
sortParam(cellClass == 3) = sortParam(cellClass == 3) * 10000000;
sortParam(cellClass == 0) = sortParam(cellClass == 0) * 100000000000000;

[sortParam cellOrder] = sort(sortParam,'ascend');
cellList(cellOrder(1:50))

sortTable = rawTable(cellOrder,:);

cellTable = mat2cell(sortTable,[ones(1,size(sortTable,1))],[ones(1,size(sortTable,2))])
fullTable = cat(1,header,cellTable);

for y = 1:size(fullTable,1)
    for x = 1:size(fullTable,2)
        if fullTable{y,x} == 0
            fullTable{y,x} = ' ';
        end
    end
end


cellProp = cellPos(:,1);




