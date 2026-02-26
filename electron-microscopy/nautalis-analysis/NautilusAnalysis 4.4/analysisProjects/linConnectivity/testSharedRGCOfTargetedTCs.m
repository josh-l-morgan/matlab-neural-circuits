%% Checks whether the TCs with the highest innervation by LIN 125 ...
%% also share RGC input by creating a connectivity matrix

clear all
load('MPN.mat');
load([MPN 'obI.mat']);

seedList = 125;
allEdges = obI.nameProps.edges;
Pre = postTo(allEdges,seedList);
Post = preTo(allEdges,seedList);
pre125 = Pre(:,1);


[synNum bx] = sort(Post(:,2),'descend');
mostTC = Post(bx,1);

topTC = mostTC(2:9);

allPre = [];
for i = 1:length(topTC)
    
   Pre = postTo(allEdges,topTC(i));
   allPre = cat(1,allPre,Pre(:,1));
    
end

conPre =  unique(allPre);
conPost = [125; topTC];

con = zeros(length(conPre),length(conPost));
for y = 1:length(conPre)
    for x = 1:length(conPost)
        
        con(y,x) = sum((allEdges(:,2) == conPre(y)) & ...
            (allEdges(:,1) == conPost(x)));
        
    end
end

image(con*10)