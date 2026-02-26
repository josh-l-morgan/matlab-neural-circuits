load('MPN.mat')

synDir = [MPN 'synPos3\'];
mkdir(synDir);

load([MPN 'obI.mat'])
synapses = obI.nameProps.edges;
allPostCells = unique(synapses(:,1));

showCell = [10 5001 5003 ];

%%Reference lists
axnet1 =  [1000:1500];
axnet2 = [2000:2031 2038:2300];
axnet3 = [2032:2037];
oneAx = 2040;

prePop{2} = [5003];
postPop{2} = [10];

prePop{3} = [5001];
postPop{3} =  [10];

prePop{1} = [];
postPop{1} =  [];




viewProps.dim = 1;

showCellNames = cat(2,num2cell(showCell), {'myelin 5003'});
% 
% col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
% col = col(randperm(size(col,1)),:);
% 
% col = col + .5;

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
rgcWithTarg = intersect(preWithTarg,rgcList);
linWithTarg = intersect(preWithTarg,linList);


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

cellList = obI.cell.name;


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



%{

mid = [11422, 19706, 4689]; % glomB
mid = [14815, 18606, 2581] %glomG
mid = [13496, 19858, 3743] %glomA
mid = [12492, 21454, 5875] %massive input on 30001 225

mid = round([mid(1)/8 mid(2)/8 mid(3)/4]);
rad = [300 300 300];
viewProps.viewWindow = [mid(2) - rad(1)  mid(1) - rad(2)  mid(3) - rad(3) ; ...
    mid(2) + rad(1)  mid(1) + rad(2)  mid(3) + rad(3) ];
%}


% viewProps.viewWindow = [2300  1400  600 ; ...
%     2400  1600  1200 ];

%Target 2350 1500 700, Y 18800, x 1200, z = 2800

viewProps.maxScaleFactor = .1;
viewProps.sumScaleFactor = 3;
viewProps.viewWindow = double([1 1 1; fsize]);
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;


I_topSum = showCellsAndMore(viewProps);
image(uint8(I_topSum*1))
%{
imwrite(uint8(I_topSum),[cellPicDir 'axes201_4.png'])
%}


%%Synapse image


dSamp = [8 8 4];

dim = viewProps.dim;
if dim == 1
    dims = [3 2];
elseif dim == 2
    dims = [3 1];
elseif dim == 3
    dims = [1 2];
end


%%Ball
ballRad = 8;
ball = ones(ballRad*2+1,ballRad*2+1,ballRad*2+1);
[y x z] = ind2sub(size(ball),find(ball));
dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+ (z-mean(z)).^2);
ball(:) = dists;
ball(ball>ballRad) = 0;
ballSum = sum(ball>0,3);
ballSum = ballSum * 200/max(ballSum(:));



anchors = obI.colStruc.anchors;
anchors(:,1) = anchors(:,1)/dSamp(1);
anchors(:,2) = anchors(:,2)/dSamp(2);
anchors(:,3) = anchors(:,3)/dSamp(3);
anchors = round(anchors);
anchors(anchors<1) = 1;
synapses;


showSyn  = I_topSum *0;
for p = 1:length(prePop)
    
    usePre = [];
    for i = 1:length(prePop{p})
        usePre = cat(1,usePre,find(synapses(:,2)== prePop{p}(i)));
    end
    
    usePost = [];
    for i = 1:length(postPop{p})
        usePost = cat(1,usePost,find(synapses(:,1)== postPop{p}(i)));
    end
    
    foundSyn = intersect(usePre,usePost);
    
    useSynID = synapses(foundSyn,3);
    useAnchors = anchors(useSynID,:);
    
    Isize = [viewProps.fsize(dims(1)) viewProps.fsize(dims(2))];
    
    anchorInd = sub2ind(Isize, useAnchors(:,dims(1)),useAnchors(:,dims(2)));
    uAnchInd = unique(anchorInd);
    if length(uAnchInd)>1
        valAnch = histc(anchorInd,uAnchInd);
    else
        valAnch = 1;
    end

% 
% for d = 1:length(anchorInd)
%    
%     
% end

synImage = zeros(Isize); 
synImage(uAnchInd) = valAnch;
synImage = convn(synImage,ballSum,'same');
showSyn(:,:,p) =  synImage;

image(uint8(showSyn*3));


end

maskSyn = repmat(sum(showSyn,3),[1 1 3]);

showSyn(~maskSyn) = I_topSum(~maskSyn)/2;
%image(uint8(showSyn))



%{
imwrite(uint8(showSyn),[cellPicDir 'cell_10synMyelin_view1.png'])
%}

cropSyn = showSyn;
[y x] = find(sum(cropSyn,3)>0);
cropSyn = cropSyn(min(y):max(y),min(x):max(x),:);
image(uint8(cropSyn * 4)),pause(.01)


%{

imwrite(uint8(cropSyn),sprintf('%ssynPos_%06.0f.png',synDir,showCell))
%}

