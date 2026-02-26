

load([MPN 'obI.mat'])

%%
conTo(1).targ = 108;
conTo(2).targ = 201;
conTo(3).targ = 109;


synapses = obI.nameProps.edges;

tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);
rgcList = obI.nameProps.cellNum(obI.nameProps.rgc);
linList = obI.nameProps.cellNum(obI.nameProps.lin);


for t = 1:length(conTo)
  
preWithTarg = unique(synapses(synapses(:,1)==conTo(t).targ,2));
preWithTarg = preWithTarg(preWithTarg>0);
rgcWithTarg = intersect(preWithTarg,rgcList)
linWithTarg = intersect(preWithTarg,linList)


synWithPre = [];
for i = 1:length(rgcWithTarg)
    synWithPre = cat(1,synWithPre, find((synapses(:,2)==rgcWithTarg(i)) & ...
        (synapses(:,1)>0)));
end
    
synPre = synapses(synWithPre,2);
synPost = synapses(synWithPre,1);
synObj = synapses(synWithPre,3);
postList = unique(synPost(synPost>0));

conTo(t).preList = preWithTarg;
conTo(t).rgcList = rgcWithTarg;
conTo(t).postList = postList;
conTo(t).tcrList = intersect(postList,tcrList);
conTo(t).postLIN = intersect(postList,linList);
conTo(t).preLIN = linWithTarg

conTo(t).syn = [synPre synPost];
conTo(t).synObj = synObj;

end




%%
colMap = hsv(255);
colMap = cat(1,[0 0 0],colMap)
cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end

cellList = obI.cell.name


showCellNames = num2cell(conTo(1).tcrList);
showCellNames = num2cell(conTo(2).tcrList);
showCellNames = num2cell(cellList((cellList>1900)& (cellList<3000)));

testPre = cellList((cellList>1900)& (cellList<3000));


for p = 1:length(testPre)

showCellNames = num2cell(testPre(p));
    
col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
col = col(randperm(size(col,1)),:);


%{
%% Glom A split cell types
showCellNames = num2cell(intersect(cellList,[201 203 204 2000:2099]));
col = zeros(length(showCellNames),3);
col(4:end,2) = 1;
col(1,:) = [.5 0 1]*1.3;
col(2,:) = [0 .5 1]*1.3;

col(3,1) = 5;
%}



fsize = double(max(cat(1,dsObj.subs),[],1));
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [minVal; fsize];
%}

%{
%%two pops
pop1 = [conTo(1).tcrList conTo(2).tcrList];
pop2 = [conTo(1).postLIN conTo(2).postLIN];
showCellNames = num2cell([pop1 pop2]);
col = zeros(length(pop1)+length(pop2),3);
col(1:length(pop1),1) = 1;
col(1:length(pop1),3) = 1;

col(length(pop1)+1:end,2) = 1;
%}

%col = colMap(ceil((1:length(dsObj))*256/length(dsObj)),:);
% 

downSamp = [1 1 1];

mid = [12492, 21454, 5875] %massive input on 30001 225
mid = [14815, 18606, 2581] %glomG
mid = [11422, 19706, 4689]; % glomB
mid = [13496, 19858, 3743] %glomA

mid = round([mid(1)/downSamp(2) mid(2)/downSamp(1) mid(3)/downSamp(3)]);
rad = [120 200 120*4];

viewWindow = [mid(2) - rad(1)  mid(1) - rad(2)  mid(3) - rad(3) ; ...
    mid(2) + rad(1)  mid(1) + rad(2)  mid(3) + rad(3) ];

viewProps.viewWindow = viewWindow;

%}




viewProps.maxScaleFactor = 1;
viewProps.sumScaleFactor = 1;
viewProps.viewWindow = double([1 1 viewWindow(1,3); fsize]);
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.dim = 1;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;


I_topSum = showCellsAndMore(viewProps);
image(uint8(I_topSum * .4))

image(uint8(backGroundAx * .2 + I_topSum*6))
pause(.01)
%{
imwrite(uint8(I_topSum ),[cellPicDir 'niceGlomA2.png'])
imwrite(uint8(backGroundAx * .4 ),[cellPicDir 'niceGlomA2.png'])

%}

%%Count synapses

colPos = obI.colStruc.anchors;
colName = obI.colStruc.names;
synPos = colPos(synapses(:,3),:);
synName = colName(synapses(:,3));

vW = viewWindow;
synInWindow = (synPos(:,1)/downSamp(1) > vW(1,1)) &  (synPos(:,1)/downSamp(1)< vW(2,1)) & ...
    (synPos(:,2)/downSamp(2) > vW(1,2)) &  (synPos(:,2)/downSamp(2) < vW(2,2)) & ...
    (synPos(:,3)/downSamp(3) > vW(1,3)) &  (synPos(:,3)/downSamp(3) < vW(2,3));

preCell = testPre(p);
postCells = unique(synapses(synapses(:,2)==preCell,1));
disp(preCell)


if ~isempty(postCells)
clear synCount
for tc = 1:length(postCells)
    postCell = postCells(tc);
synCount(tc) = sum((synapses(:,1)==postCell) & (synapses(:,2) == preCell) );

%synCount(tc) = sum((synapses(:,1)==postCell) & (synapses(:,2) == preCell) & synInWindow);
end
%disp(sprintf('%d synapses between %d and %d',synCount,preCell,postCell))



showRes = num2cell([postCells synCount'])
end

if sum(I_topSum(:))
   'paused' 
pause
end

end


