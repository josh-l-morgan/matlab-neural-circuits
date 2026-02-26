

%MPN = GetMyDir

%load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])

  seedList = [108 201 109 903 907];
  useList = obI2cellList_all(obI);
  seedPref = seedPreferences(seedList,useList);
  allEdges = obI.nameProps.edges(:,[2 1]);
  conTo = makeConTo(obI,seedList);
  
  cellList = useList.postList;
  
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


%%
ancPos = zeros(length(cellList),3);
for i = 1:length(cellList)
    targ = find(obI.cell.name == cellList(i));
    mainID = obI.cell.mainID(targ);
    ancPos(i,:) = anchors(mainID,:);
end
scatter(ancPos(:,2),ancPos(:,1),30,'filled')


conPos = ancPos;
for i = 1:length(cellList);
   isPre = synapses(:,2) == cellList(i);
   isPost = synapses(:,1) == cellList(i);
   isEither = isPre+ isPost;
   if sum(isEither)
       conPos(i,:) = mean(synAnchors(find(isEither),:),1);
   end
end
scatter(conPos(:,2),conPos(:,1),30,'filled')

%% cell Color


cellCol = conPos;
cellCol(:,1) = cellCol(:,1)/max(cellCol(:,1));
cellCol(:,2) = cellCol(:,2)/max(cellCol(:,2))*0;
cellCol(:,3) = cellCol(:,3)/max(cellCol(:,3))*0;



cellCol(cellCol<0) = 0;
cellCol(cellCol>1) = 1;
















