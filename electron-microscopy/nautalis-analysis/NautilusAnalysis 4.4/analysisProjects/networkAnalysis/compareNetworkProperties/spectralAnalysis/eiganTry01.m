%%Notes


%{
Aij = Aji
eig(B B_)
U,S,V = svd(B)
B = U S V t
U(i,2)
S(2)
U S^2U^T = B B^t

U(i,1) * S(1)
%}



%% Get data
loadData = 1;
if loadData
    %clear all
    MPN = GetMyDir;
    load([MPN 'obI.mat']);
    seedList = [ 108 201 109 907 903];
    %seedList = [ 108  201 109 ];
    
    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    %seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
end




%% Reshape directed to undirected with pre first

conRaw = useList.con;
[s1 s2] = size(conRaw);
con = zeros(s1+s2);
con(1:s1,s1+1:end) = conRaw;
con(s1+1:end,1:s1) = conRaw';
nodeIDs = [useList.preList useList.postList];
image(con*20)

use1 = con(find(nodeIDs == 108),:)>0;
use2 = con(find(nodeIDs == 301),:)>0;
use3 = (con(find(nodeIDs == 109),:)>0) | ...
    (con(find(nodeIDs == 903),:)>0) | ...
    (con(find(nodeIDs == 907),:)>0) ;




%% get eigan modes
B = con;
[U S V] = svd(B);
diagS = diag(S);


scatter(U(:,1),U(:,2))

for k = 1:5
scatter(U(:,k) * diagS(k),U(:,k+1)* diagS(k+1))
hold on
scatter(U(use1,k) * diagS(k),U(use1,k+1)* diagS(k+1),'r')
scatter(U(use2,k) * diagS(k),U(use2,k+1)* diagS(k+1),'g')
scatter(U(use3,k) * diagS(k),U(use3,k+1)* diagS(k+1),'c')
hold off
pause
end




