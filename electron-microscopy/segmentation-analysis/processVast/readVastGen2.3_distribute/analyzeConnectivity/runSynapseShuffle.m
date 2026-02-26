%%Test for redundant synapses and shared synapses using simple shuffle of
%%edges

MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\export_14+04+25_matOut\'
load([MPN 'obI.mat'])


%% select edges

targCell = 108;
postHits =  find(obI.nameProps.edges(:,1) == targCell);
preCellList = obI.nameProps.edges(postHits,2);

targCells = unique(preCellList(preCellList>0));
edges = [];
for i = 1: length(targCells)
    foundEdge =  find(obI.nameProps.edges(:,2) == targCells(i));
    edges = cat(1,edges,obI.nameProps.edges(foundEdge,1:2));
end

goodEdge = find(edges(:,1)>0);
edges = edges(goodEdge,:);


%% Get Real data

con = edge2con(edges);
numEdges = size(edges,1);
numCombos = size(con,1) * size(con,2);

sharedPreReal = minSharedMat(con);
sharedPostReal = minSharedMat(con');
redundantReal = countRedundant(con)/numEdges * 100;
occupiedReal = sum(con(:)>0)/numCombos * 100;
maxConPreReal = sumMaxCon(con);
maxConPostReal = sumMaxCon(con');

tic

%% run Monte
reps = 10000;

redundant = zeros(reps,1);

sharedPre = redundant;
sharedPost = redundant;
occupied = redundant;
maxConPre = redundant;
maxConPost = redundant;

for i = 1:reps
    
    newEdges = shuffleEdges(edges);
        
    con = edge2con(newEdges);
    sharedPre(i) = minSharedMat(con);
    sharedPost(i) = minSharedMat(con');
    redundant(i) = countRedundant(con)/numEdges * 100;
    occupied(i) = sum(con(:)>0)/numCombos * 100;
    maxConPre(i) = sumMaxCon(con);
    maxConPost(i) = sumMaxCon(con');
    
    %image(con*5),pause(.1);
    
end
toc

%%

subplot(6,1,1)
scatterMonte(redundantReal,redundant,'% redundant');
subplot(6,1,2)
scatterMonte(occupiedReal,occupied,'% occupied');
subplot(6,1,3)
scatterMonte(maxConPreReal,maxConPre,'% syn in Pre max con');
subplot(6,1,4)
scatterMonte(maxConPostReal,maxConPost,'% syn in Post max con');
subplot(6,1,5)
scatterMonte(sharedPostReal,sharedPost,'sharedPost');
subplot(6,1,6)
scatterMonte(sharedPreReal,sharedPre,'sharedPre');


