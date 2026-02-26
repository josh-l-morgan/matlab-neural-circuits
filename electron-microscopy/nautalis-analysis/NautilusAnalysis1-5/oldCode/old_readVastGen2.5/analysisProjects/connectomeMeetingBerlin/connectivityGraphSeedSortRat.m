%{
 clear all
MPN = GetMyDir;
load([MPN 'obI.mat'])
GPN = [MPN 'graphs2\']
if ~exist(GPN),mkdir(GPN),end
%}

%% Decide on pres and posts



edges = obI.nameProps.edges(:,1:2);

preList = unique(edges(:,2));
preList = preList(preList>0);
postList = unique(edges(:,1));
postList = postList(postList>0);

tcrList = unique(obI.nameProps.cellNum(obI.nameProps.tcr));
tcrList = tcrList(tcrList>0);
linList = unique(obI.nameProps.cellNum(obI.nameProps.lin));
linList = linList(linList>0);
linList = [ 256 125 138 142 166 168 204 208 209 214] % cells with cell bodies
rgcList = unique(obI.nameProps.cellNum(obI.nameProps.rgc));
rgcList = rgcList(rgcList>0);


seed1 = 108;
seed2 = 201;

preSeed1 = edges(edges(:,1) == seed1,2);
preSeed2 = edges(edges(:,1) == seed2,2);

rgcSeed1 = intersect(preSeed1,rgcList);
rgcSeed2 = intersect(preSeed2,rgcList);
allRGC = unique([rgcSeed1,rgcSeed2]);

countRGC1 = zeros(size(allRGC));
countRGC2 = zeros(size(allRGC));
for i = 1:length(allRGC)
   countRGC1(i) = length(find(preSeed1 == allRGC(i)));
   countRGC2(i) = length(find(preSeed2 == allRGC(i)));

end

difRGC = countRGC1-countRGC2;
[sortedDifs idx] = sort(difRGC,'descend');
sortRGC = allRGC(idx);

weightTCR = zeros(max(tcrList),1);
countTCR = weightTCR;
for i = 1:length(sortRGC)
    postRGC = edges(edges(:,2) == sortRGC(i),1);
    tcrPost = intersect(postRGC,tcrList);
    for t = 1:length(tcrPost)
        weightTCR(tcrPost(t)) = weightTCR(tcrPost(t)) + sortedDifs(i);
        countTCR(tcrPost(t)) = countTCR(tcrPost(t))+1;
    end
end

useTCR = find(countTCR);
scaleTCR = weightTCR(useTCR)./countTCR(useTCR);
[sortTCRweights  idx] = sort(weightTCR(useTCR),'descend');
sortTCR = useTCR(idx);







%% repeat for inhibitory

weightLIN = zeros(max(linList),1);
countLIN = weightLIN;
for i = 1:length(sortRGC)
    postRGC = edges(edges(:,2) == sortRGC(i),1);
    linPost = intersect(postRGC,linList);
    for t = 1:length(linPost)
        weightLIN(linPost(t)) = weightLIN(linPost(t)) + sortedDifs(i);
        countLIN(linPost(t)) = countLIN(linPost(t))+1;
    end
end

useLIN = linList;
scaleLIN = weightLIN(useLIN)./countLIN(useLIN);
[sortLINweights  idx] = sort(weightLIN(useLIN),'descend');
sortLIN = useLIN(idx);


%%
markDivision = 2209823535;
bufGraph = 231085237098;

sortPre = [bufGraph sortRGC bufGraph markDivision bufGraph sortLIN bufGraph];
markPre = [bufGraph sortRGC *0+1 bufGraph markDivision bufGraph sortLIN*0+2 bufGraph];

sortPost = [bufGraph sortTCR' bufGraph markDivision bufGraph sortLIN bufGraph];
markPost = [sortTCR'*0+1 markDivision sortLIN *0+2];
markGraph = repmat(markPre',[1,length(markPost)]) + ...
    repmat(markPost,[length(markPre),1]);

imwrite(uint8(markGraph * 50),[GPN 'markGraph.png']);

%% Make graph

eGraph = zeros(length(sortPre),length(sortPost));
eGraphCell = cell(length(sortPre),length(sortPost));
for i = 1:size(edges,1)
   targPre = find(sortPre == edges(i,2) );
   targPost = find(sortPost == edges(i,1));
   if ~isempty(targPre) & ~isempty(targPost);
       if isempty(eGraphCell{targPre,targPost})
           eGraphCell{targPre,targPost} =1 ;
       else
                      eGraphCell{targPre,targPost} =eGraphCell{targPre,targPost}+1 ;
       end
      eGraph(targPre,targPost) = eGraph(targPre,targPost)+1; 
   end
end
sum(eGraph(:))

eGraphChar = cell(length(sortPre),length(sortPost));
for x = 1:size(eGraph,2)
    for y = 1:size(eGraph,1)
        gotNum = eGraph(y,x);
        if gotNum
            eGraphChar{y,x} = num2str(gotNum);
        else
            eGraphChar{y,x} = ' ';
        end
        
    end
end


%%
cMap = jet(256)
cMap(1,:) = 0;
scaleSynColor = 15;
colormap(cMap)
image(eGraph*scaleSynColor)
colGraph = zeros(size(eGraph,1), size(eGraph,2),3);
showGraph = eGraph * scaleSynColor;
showGraph(showGraph>(size(cMap,1)-1)) = (size(cMap,1)-1);
for c = 1:3
    colLook = cMap(:,c);
    colTemp = colLook(showGraph+1);
    colGraph(:,:,c) = colTemp * 256;
end
image(uint8(colGraph))
imwrite(uint8(colGraph),[GPN 'showGraph.png']);

% xAxLab = [0:20];
% image(xAxLab*scaleSynColor)
% xlabel(xAxLab)

%% identify cell
labelGraph = eGraph;
   targPre = find(sortPre == 1033);
   targPost = find(sortPost == 217);
   
   
   labelGraph(targPre,targPost) = labelGraph(targPre,targPost) + 1112;
   
   image(labelGraph*scaleSynColor)
   
%% divide graph

divGraph = colGraph;
preTarg = find(sortPre == markDivision);
divGraph(preTarg,:,:) = 20;
postTarg = find(sortPost == markDivision);
divGraph(:,postTarg,:) = 20;
image(uint8(divGraph))
imwrite(uint8(divGraph),[GPN 'divGraph.png']);

   
   
   
   
   
      
