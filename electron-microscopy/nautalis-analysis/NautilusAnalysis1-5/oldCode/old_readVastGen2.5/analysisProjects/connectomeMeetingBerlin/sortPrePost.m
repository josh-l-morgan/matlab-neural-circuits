

eGraph;
%{

[5 1] = 1
[1 1] = 1
[2 2] = 2

%}

%%
tcrList = unique(obI.nameProps.cellNum(obI.nameProps.tcr));
tcrList = tcrList(tcrList>0);
linList = unique(obI.nameProps.cellNum(obI.nameProps.lin));
linList = linList(linList>0);
linList = [ 256 125 138 142 166 168 204 208 209 214] % cells with cell bodies
rgcList = unique(obI.nameProps.cellNum(obI.nameProps.rgc));
rgcList = rgcList(rgcList>0);

%%
listPre = unique(edges(:,2));
listPre = listPre(listPre>0);
listPost = unique(edges(:,1));
listPost = listPost(listPost>0);


listPre = intersect(listPre,rgcList);
%listPost = intersect(listPost,tcrList);


preNum = length(listPre);
postNum = length(listPost);





%% protoGraph
protoGraph = zeros(length(listPre),length(listPost));
for i = 1:size(edges,1)
   targPre = find(listPre == edges(i,2) );
   targPost = find(listPost == edges(i,1));
   if ~isempty(targPre) & ~isempty(targPost);
      protoGraph(targPre,targPost) = protoGraph(targPre,targPost)+1; 
   end
end
sum(protoGraph(:))


%% calculate links
preGraph = zeros(preNum);
for y = 1:preNum
    for x = 1:preNum
        preGraph(y,x) = sum(min(protoGraph(y,:),protoGraph(x,:)));
    end
end
    

postGraph = zeros(postNum);
for y = 1:postNum
    for x = 1:postNum
        postGraph(y,x) = sum(min(protoGraph(:,y),protoGraph(:,x)));
    end
end


%% huristic

corners = [108 201];
preSim1 = protoGraph(:,find(listPost==corners(1)));
preSim2 = protoGraph(:,find(listPost==corners(2)));
preVal = (preSim1 - preSim2);
goodPre = preSim1 | preSim2;

postSim1 = postGraph(find(listPost==corners(1)),:);
postSim2 = postGraph(find(listPost==corners(2)),:);
goodPost = postSim1 | postSim2;

postVal = (postSim1-postSim2)./(postSim1+postSim2+.01);
postType = (postSim1>0) + (postSim2>0)*2;

[newPostVal postIDX] = sort(postVal,'descend');
[newPreVal preIDX] = sort(preVal,'descend');

sortPost = listPost(postIDX);
goodPost = goodPost(postIDX);
postType = postType(postIDX);
sortPost = sortPost(goodPost);
postType = postType(goodPost);

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

%% springs

pos = rand(preNum,1)*100;
pushVal = 1;
pushMax = 1;
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







