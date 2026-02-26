clear all
load('MPN.mat');
load([MPN 'obI.mat'])
useImage  = 1;
edges = obI.nameProps.edges;

%%
seedList = [108 201 903 907];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
preList = useList.preList;
postList = useList.postList;

%% Spectral analysis

%%Make symmetric
nodeList = [preList setdiff(postList,seedList)]
symGraph = zeros(length(nodeList),length(nodeList));
for i = 1:size(edges,1)
   targPre = find(nodeList == edges(i,2) );
   targPost = find(nodeList == edges(i,1));
   if ~isempty(targPre) & ~isempty(targPost);
      symGraph(targPre,targPost) = symGraph(targPre,targPost)+1; 
      symGraph(targPost,targPre) = symGraph(targPost,targPre)+1; 
   end
end

isedge = find(sum(symGraph)>0);
symGraph = symGraph(isedge,:);
symGraph = symGraph(:,isedge);
nodeList = nodeList(isedge);

sum(symGraph(:))
% image(symGraph*100)

%[nodeIDX,C] = jwSpec(symGraph,8,4)
    
U = sim2U2(symGraph,8);

%U = symGraph;
%[nodeIDX,C] = kmeans(U,3,'replicates',100,'emptyaction','drop');
nodeIDX = clusterdata(U,12);












specRes = useList
specRes.nodeIDX = nodeIDX;
specRes.nodeList = nodeList;

save('.\data\otherClustering\aglomRes.mat','specRes')


for i = 1:size(U,2)
   pU = U(:,i);
   pU = pU-min(pU);
   pU = pU * 200/max(pU);
   plotU(:,i) = pU;
end
result.nodeX = plotU(:,1);
result.nodeY = plotU(:,2);
result.nodeIDs = nodeList';
result.snapTime = 0;
result.cellGroups = nodeIDX
save('.\data\otherClustering\aglom_result.mat','result')



%% save('.\data\otherClustering\specResGood1.mat','specRes')


return
%% huristic
%{
corners = [108 201];

postSim1 = postGraph(find(listPost==corners(1)),:);
postSim2 = postGraph(find(listPost==corners(2)),:);
goodPost = postSim1 | postSim2;

postVal = (postSim1-postSim2)./(postSim1+postSim2+.01);
postType = (postSim1>0) + (postSim2>0)*2;

[newPostVal postIDX] = sort(postVal,'descend')

sortPost = listPost(postIDX);
goodPost = goodPost(postIDX);
postType = postType(postIDX);
sortPost = sortPost(goodPost);
postType = postType(goodPost);


preVal = listPre*0;
for i = 1:length(listPre)
    conPost = protoGraph(i,:);
    preVal(i) = sum(conPost.*postVal)/sum(conPost);
end

preSim1 = protoGraph(:,find(listPost==corners(1)));
preSim2 = protoGraph(:,find(listPost==corners(2)));
goodPre = preSim1 | preSim2;

[newPreVal preIDX] = sort(preVal,'descend');
sortPre = listPre(preIDX);
goodPre = goodPre(preIDX);
sortPre = sortPre(goodPre);


%%hGraph
hGraph = zeros(length(sortPre),length(sortPost));
for i = 1:size(edges,1)
   targPre = find(sortPre == edges(i,2) );
   targPost = find(sortPost == edges(i,1));
   if ~isempty(targPre) & ~isempty(targPost);
      hGraph(targPre,targPost) = hGraph(targPre,targPost)+1; 
   end
end
sum(hGraph(:))

image(hGraph * 10)
%}
% 
% 
% imwrite(uint8(divGraph),[GPN 'divGraph.png']);

% %%
% axType = sortPre*0+1;
% axType(sortPre>2000) = 2;
% axType(sortPre>2031) = 3;
% 
% axGraph = repmat(axType',[1 length(sortPost)]);
% 
% hGraph = hGraph + axGraph * 2;

%% Random swap

%%recalculate links

[preNum postNum] = size(protoGraph);

preGraph2 = zeros(preNum);
for y = 1:preNum
    for x = 1:preNum
        preGraph2(y,x) = sum(min(protoGraph(y,:)/sum(protoGraph(y,:)),protoGraph(x,:)/sum(protoGraph(x,:))));
        preGraph2(y,x) = sum(min(protoGraph(y,:),protoGraph(x,:)));

    end
end

postGraph2 = zeros(postNum);
for y = 1:postNum
    for x = 1:postNum
        postGraph2(y,x) = sum(min(protoGraph(:,y)/sum(protoGraph(:,y)),protoGraph(:,x)/sum(protoGraph(:,x))));
         postGraph2(y,x) = sum(min(protoGraph(:,y),protoGraph(:,x)));

    end
end


[postGraph3 postIDX] = reshuffleSimilarity(postGraph2);
sortPost2 = listPost(postIDX);


[preGraph3 preIDX] = reshuffleSimilarity(preGraph2);
sortPre2 = listPre(preIDX);



%% rGraph
rGraph = zeros(length(sortPre2),length(sortPost2));
for i = 1:size(edges,1)
   targPre = find(sortPre2 == edges(i,2) );
   targPost = find(sortPost2 == edges(i,1));
   if ~isempty(targPre) & ~isempty(targPost);
      rGraph(targPre,targPost) = rGraph(targPre,targPost)+1; 
   end
end
sum(rGraph(:))

image(rGraph * 10)


%% save




%% sGraph
for i = 1:size(rGraph,1)
   hits = find(rGraph(i,:));
   vals = rGraph(i,hits);
   matchVal(i) = sum(hits.*vals)/sum(vals);
end

[sortMatch idx] = sort(matchVal,'ascend');
sortPre3 = sortPre2(idx);

sGraph = zeros(length(sortPre3),length(sortPost2));
for i = 1:size(edges,1)
   targPre = find(sortPre3 == edges(i,2) );
   targPost = find(sortPost2 == edges(i,1));
   if ~isempty(targPre) & ~isempty(targPost);
      sGraph(targPre,targPost) = sGraph(targPre,targPost)+1; 
   end
end
sum(sGraph(:))

image(sGraph * 10)


%% use Graph
usePostOrder = sortPost2;
usePreOrder = sortPre2;

uGraph = zeros(length(usePreOrder),length(usePostOrder));
for i = 1:size(edges,1)
   targPre = find(usePreOrder == edges(i,2) );
   targPost = find(usePostOrder == edges(i,1));
   if ~isempty(targPre) & ~isempty(targPost);
      uGraph(targPre,targPost) = uGraph(targPre,targPost)+1; 
   end
end
sum(uGraph(:))

image(uGraph * 10)




%% springs
%{
pos = rand(preNum,1)*100;
pushVal = 3;
pushMax = 10;
pullVal = 1;
pullMax = 10;
showPos = zeros(preNum,1);
reps = 10000;

for rep = 1:reps

for i = 1:preNum
   relPos = pos-pos(i);
   direct = relPos./abs(relPos);
   direct(isnan(direct)) = 0;
   
   push = pushVal./relPos.^2;
   push(i) = 0;
   push(push>pushMax) = pushMax;
   totPush = sum(push.*direct);
   
   weight = preGraph(:,i);
   weightPos = relPos .* weight;
   pull = pullVal .* weight .* direct;
   totPull = sum(pull);
   
   shift = totPush + totPull;
   pos(i) = pos(i) + shift/preNum;
   pos(pos<1) = 1;
   pos(pos>100) = 100;
   
   rPos = round(pos);
   uPos = unique(rPos);
   hPos = histc(rPos,uPos);
   showPos = showPos * 0;
   showPos(uPos) = hPos*30;
   
end

image(showPos),pause(.01)

end

%}

%%
%preGraph3 postGraph3 rGraph
eGraph = rGraph;
%eGraph = 1:25
cMap = jet(256)
cMap(1,:) = 0;
scaleSynColor = 10;
colormap(cMap)
colGraph = zeros(size(eGraph,1), size(eGraph,2),3);
showGraph = eGraph * scaleSynColor;
showGraph(showGraph>(size(cMap,1)-1)) = (size(cMap,1)-1);
for c = 1:3
    colLook = cMap(:,c);
    colTemp = colLook(showGraph+1);
    colGraph(:,:,c) = colTemp * 256;
end
image(uint8(colGraph))


preType = sortPre3 *0 + 1;
preType(sortPre3>1500) = 2;
preType(sortPre3>2031) = 3;
graphPreType = repmat(preType',[1,size(sGraph,2)]);



%{
imwrite(uint8(colGraph),[GPN 'rGraph.png']);


imwrite(uint8(graphPreType*10),[GPN 'sGraphPreType.png']);

%}

%% Marker


markG = colGraph;
markTarg = 109;
markG(:,(sortPost2 ==markTarg),:) = markG(:,(sortPost2==markTarg),:) +30;
image(uint8(markG))







