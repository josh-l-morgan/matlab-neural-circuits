figure();
dS=2;
curColMap=hsv(max(cellSub.brhID(:))+1);
colInd=cellSub.brhID+1;
colMat=curColMap(colInd,:);
scatter3(cellSub.subs(1:dS:end,1),cellSub.subs(1:dS:end,2),cellSub.subs(1:dS:end,3),3,colMat(1:dS:end,:));
hold on

if 0
for i=1:length(arbor.branches)
    curBrPos=arbor.branches(i).nodePos;
    scatter3(curBrPos(:,1),curBrPos(:,2),curBrPos(:,3),5,curColMap(i,:));
    
    
end
end