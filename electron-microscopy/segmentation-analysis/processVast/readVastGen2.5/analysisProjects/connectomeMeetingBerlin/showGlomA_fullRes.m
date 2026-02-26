

load([MPN 'obI.mat'])

%%
conTo(1).targ = 108;
conTo(2).targ = 201;
conTo(3).targ = 109;


synapses = obI.nameProps.edges;



for t = 1:length(conTo)
  
preWithTarg = unique(synapses(synapses(:,1)==conTo(t).targ,2));
preWithTarg = preWithTarg(preWithTarg>0);
synWithPre = [];
for i = 1:length(preWithTarg)
    synWithPre = cat(1,synWithPre, find((synapses(:,2)==preWithTarg(i)) & ...
        (synapses(:,1)>0)));
end
    
synPre = synapses(synWithPre,2);
synPost = synapses(synWithPre,1);
synObj = synapses(synWithPre,3);
postList = unique(synPost(synPost>0));

conTo(t).preList = preWithTarg;
conTo(t).postList = postList;
conTo(t).syn = [synPre synPost];
conTo(t).synObj = synObj;

end




%%

cellPicDir = [MPN '\cellPic\'];
if ~exist(cellPicDir,'dir'), mkdir(cellPicDir), end

cellList = obI.cell.name;

showCellNames = num2cell([conTo(2).postList]);
showCellNames = num2cell([201 203]); %Glom a post
showCellNames = num2cell([201 237]); %Glom B post
showCellNames = num2cell([201 232]); %Glom g post
showCellNames = num2cell([conTo(2).preList]);


showCellNames = num2cell(intersect(cellList,[201 203 2000:2099]));


%%Axons plus dim TC
showCellNames = num2cell(intersect(cellList,[201 203  2000:2099]));
% 
% col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
% col = col(randperm(size(col,1)),:);
% 
% col(1,:) = [1 .7 1] *.2;
% col(2,:) = [.7 1 1]*.2;
% col = col /2;



pop1 = intersect(200:300,obI.nameProps.cellNum(obI.nameProps.lin));
pop2 = [201 203];
showCellNames = num2cell([pop1 pop2]);
%}

pop1 = intersect(200:300,obI.nameProps.cellNum(obI.nameProps.tcr));
pop2 = [201 203];
showCellNames = num2cell([pop1 pop2]);


pop1 = intersect(2000:3000,obI.nameProps.cellNum(obI.nameProps.rgc));
pop2 = [201 203 207 255];
showCellNames = num2cell([pop1 pop2]);

%}

col = colMap(ceil((1:length(pop1))*256/length(pop1)),:);
col = col(randperm(size(col,1)),:)/2;
col = cat(1,col,[1 .7 1; .7 1 1; 1 1 .7; 1 1 .7]/8);

viewProps.maxScaleFactor = 1;
viewProps.sumScaleFactor = .2;
%}



%{
%% Glom A split cell types
showCellNames = num2cell(intersect(cellList,[201 203 204 2000:2099]));
col = zeros(length(showCellNames),3);
col(4:end,2) = 1;
col(1,:) = [.5 0 1]*1.3;
col(2,:) = [0 .5 1]*1.3;

col(3,1) =5;
col = col/4;


viewProps.maxScaleFactor = 0;
viewProps.maxSumFactor = 1;
%}



fsize = double(max(cat(1,dsObj.subs),[],1));
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [minVal; fsize];
%}

%{
%%two pops
pop1 = conTo(1).postList';
pop2 = conTo(2).postList';
showCellNames = num2cell([pop1 pop2]);
col = zeros(length(pop1)+length(pop2),3);
col(1:length(pop1),1) = 1;
col(1:length(pop1),3) = 1;

col(length(pop1)+1:end,2) = 1;
%}

%col = colMap(ceil((1:length(dsObj))*256/length(dsObj)),:);

% col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
% col = col(randperm(size(col,1)),:);

viewProps.viewWindow = double([1 1 1; fsize]);
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.dim = 3;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;




%{
mid = [11422, 19706, 4689]; % glomB
mid = [14815, 18606, 2581] %glomG
mid = [13396, 20000, 3743] %glomA

%Y = 19040 20960 X = 11796 , 14996, z = 3263, 4223

mid = round([mid(1)/8 mid(2)/8 mid(3)/4]);
rad = [120 200 120];
viewProps.viewWindow = [mid(2) - rad(1)  mid(1) - rad(2)  mid(3) - rad(3) ; ...
    mid(2) + rad(1)  mid(1) + rad(2)  mid(3) + rad(3) ];
%}


I_topSum = showCellsAndMore(viewProps);
image(uint8(I_topSum*1))
%{
imwrite(uint8(I_topSum),[cellPicDir 'glomA_theFuckingShit.png'])
%}




