

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
showCellNames = num2cell(cellList);
showCellNames = num2cell(cellList);
showCellNames = num2cell([108 201]);
showCellNames = cat(2,num2cell([108 201]), {'boutons'}, num2cell([1000 : 1200]));
showCellNames = num2cell([2032 2033 2034 2035 201 266]);
showCellNames = num2cell([2001:2038]);
showCellNames = num2cell([2032 2033 2034 2035 201 266]);
showCellNames = num2cell([2000:2040]);
showCellNames = cat(2,num2cell([2037:2040]),{'myelin'});
showCellNames = num2cell([2000:2031]);
showCellNames = num2cell([1000:1060]);
showCellNames = num2cell(conTo(1).preList);
showCellNames = num2cell(cellList);


col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
col = col(randperm(size(col,1)),:);
%col = col * 0 + .3;
% 
% col = ones(length(showCellNames),3)/5;
% col(1,:) = [1 0 0];
% col(2,:) = [0 1 0];
% col(3,:) = [0 0 1];




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
pop1 = [conTo(1).tcrList conTo(1).tcrList];
pop2 = [conTo(1).postLIN conTo(2).tcrList];
showCellNames = num2cell([pop1 pop2]);
col = zeros(length(pop1)+length(pop2),3);
col(1:length(pop1),1) = 1;
col(1:length(pop1),3) = 1;

col(length(pop1)+1:end,2) = 1;
%}

%{
%% three axon populations
pop1 = [2032 2033 2034 2035];
pop2 = [2001:2031];
pop3 = [1001:1100];
pop4 = [201];
showCellNames = num2cell([pop1 pop2 pop3 pop4]);
col = zeros(length(pop1)+length(pop2)+ length(pop3) + length(pop4),3);
col(1:length(pop1),1) = 3;
col(length(pop1)+length(pop2)+1:end,3) = 1;
col(length(pop1)+1:length(pop1)+length(pop2),2) = 3;
col(end,:) = [3 3 3];
%}


% 
% pop1 = [conTo(1).tcrList];
% pop2 = [200:];
% pop3 = [1001:1100];
% pop4 = [201];
% showCellNames = num2cell([pop1 pop2 pop3 pop4]);
% col = zeros(length(pop1)+length(pop2)+ length(pop3) + length(pop4),3);
% col(1:length(pop1),1) = 3;
% col(length(pop1)+length(pop2)+1:end,3) = 1;
% col(length(pop1)+1:length(pop1)+length(pop2),2) = 3;
% col(end,:) = [3 3 3];
% 
% 
% 
% 
% 
% 
%col = colMap(ceil((1:length(dsObj))*256/length(dsObj)),:);
% 





mid = [11422, 19706, 4689]; % glomB
mid = [14815, 18606, 2581] %glomG
mid = [13496, 19858, 3743] %glomA
mid = [12492, 21454, 5875] %massive input on 30001 225
mid = [15529, 8045, 3656]

mid = round([mid(1)/8 mid(2)/8 mid(3)/4]);
rad = [100 100 40];
viewProps.viewWindow = [mid(2) - rad(1)  mid(1) - rad(2)  mid(3) - rad(3) ; ...
    mid(2) + rad(1)  mid(1) + rad(2)  mid(3) + rad(3) ];
%}


% viewProps.viewWindow = [2300  1400  600 ; ...
%     2400  1600  1200 ];

%Target 2350 1500 700, Y 18800, x 1200, z = 2800

viewProps.maxScaleFactor = 0;
viewProps.sumScaleFactor = 1;
%viewProps.viewWindow = double([1 1 1; fsize]);
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.dim = 3;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;


I_topSum = showCellsAndMore(viewProps);
image(uint8(I_topSum*3))
%{
imwrite(uint8(I_topSum),[cellPicDir 'axo108view2gray.png'])
%}




