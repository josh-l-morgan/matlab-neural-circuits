
load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])


edges = obI.nameProps.edges(:,1:2);


%% select pre and post cells
targCell = 108;
[con1 a] = find(edges == targCell);
allPreCells = sort(unique(edges(con1,2)));
allPostCellsFull = [];
for i = 1:length(allPreCells)
    [con1 a] = find(edges == allPreCells(i));
    allPostCellsFull = cat(1,allPostCellsFull,edges(con1,1));
end

allPostCells = sort(unique(allPostCellsFull));
histPost = hist(allPostCellsFull,allPostCells);
[sortedPost postOrder] = sort(histPost,'descend');
Big2LittlePost = allPostCells(postOrder);


%%

con = zeros(length(allPreCells),length(allPostCells));
for i = 1:size(edges,1)
    preTarg = find(allPreCells == edges(i,2));
    postTarg = find(allPostCells == edges(i,1));
    con(preTarg,postTarg) =  con(preTarg,postTarg) + 1;
end
subplot(1,1,1)

image(con*30)


%% make hists grah
preHist = con*0;
subplot(1,1,1)
for i = 1:length(allPreCells)
    preHist(i,:) = sort(con(i,:),'descend');
    subplot(length(allPreCells)/4,4,i)
    bar(preHist(i,1:5))
    ylim([0 20])
    text(4,15,num2str(allPreCells(i)))
end 

%% make post hist
postHist = con'*0;
subplot(1,1,1)
for i = 1:length(allPostCells)
    postHist(i,:) = sort(con(:,i),'descend');
    subplot(length(allPostCells)/4,4,i)
    bar(postHist(i,1:20))
    ylim([0 20])
    text(4,15,num2str(allPostCells(i)))
end 

%% shared post

for i = 1:length(allPreCells)
    for p = 1:length(allPreCells)
        sharedPost(i,p) = sum(min(con(i,:),con(p,:)))/sum(con(i,:)+con(p,:));
        sharedPost(i,p) = sum(min(con(i,:),con(p,:)));

    end
end
subplot(1,1,1)
image(sharedPost)


%% shared pre
for i = 1:length(allPostCells)
    for p = 1:length(allPostCells)
        sharedPre(i,p) = sum(min(con(:,i),con(:,p)));
    end
end
subplot(1,1,1)
image(sharedPre)




%% shared pre

%plot(preHist')
