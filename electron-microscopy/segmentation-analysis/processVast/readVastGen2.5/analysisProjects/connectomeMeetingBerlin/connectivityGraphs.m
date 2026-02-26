
   %{ 
MPN = GetMyDir;
load([MPN 'obI.mat'])
GPN = [MPN 'graphs\']
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
rgcList = unique(obI.nameProps.cellNum(obI.nameProps.rgc));
rgcList = rgcList(rgcList>0);

sortPre = [sort(rgcList) sort(linList)];
markPre = [sort(rgcList)*0+1 sort(linList)*0+2];

sortPost = [sort(tcrList) sort(linList)];
markPost = [sort(tcrList)*0+1 sort(linList)*0+2];
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

 
table(eGraph)
table(eGraphCell)


hist(eGraph(eGraph(:)>0))

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

xAxLab = [0:20];
image(xAxLab*scaleSynColor)
xlabel(xAxLab)

%% identify cell
labelGraph = eGraph;
   targPre = find(sortPre == edges(i,2) );
   targPost = find(sortPost == 237);
   
   labelGraph(:,targPost) = labelGraph(:,targPost) + 2;
   
   image(labelGraph*scaleSynColor)
   
   
   
   
   
   
   
      
